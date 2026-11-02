#!/bin/bash
# sasa_one.sh <condition> <system> <slot>
# Version-A buried interfacial SASA for one replica, from the complex trajectory alone:
#   buried = SASA(receptor with ligand removed) + SASA(ligand with receptor removed) - SASA(complex)
# The partner must be physically stripped: cpptraj's `surf <mask>` sums only the masked atoms' contributions
# but still sees the whole system as neighbours, so masks alone give buried = 0 by construction.
# Frames match MM-PBSA exactly (1..PROD_FRAMES stride 10 -> 400 frames), so reused 100 ns files contribute
# only their first 40 ns. Writes analysis/sasa/<condition>/<system>_<slot>.csv
. "$(dirname "$0")/lib.sh"
COND=$1; SYS=$2; S=$(printf %03d $3)
TRAJ=$CAMPAIGN/trajectories/$COND/${SYS}_complex/${SYS}_complex_${S}.nc
PRM=$CAMPAIGN/topologies/$COND/${SYS}_complex.prmtop
OUT=$CAMPAIGN/analysis/sasa/$COND; mkdir -p "$OUT"
CSV=$OUT/${SYS}_${S}.csv
[ -f "$CSV" ] && exit 0
[ -f "$TRAJ" ] && [ -f "$PRM" ] || { echo "missing $TRAJ or $PRM"; exit 1; }
read -r RECRES _ _ <<<"$($PYTHON "$CAMPAIGN/scripts/check_topologies.py" --ranges "$PRM")"
NRES=$($PYTHON "$CAMPAIGN/scripts/solute_residues.py" "$PRM" | cut -d- -f2)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
for job in "cpx:" "rec:strip :$((RECRES+1))-$NRES" "lig:strip :1-$RECRES"; do
  nm=${job%%:*}; act=${job#*:}
  $CPPTRAJ > "$T/$nm.log" 2>&1 <<CP
parm $PRM
trajin $TRAJ 1 $PROD_FRAMES 10
$act
surf $nm out $T/$nm.dat
go
CP
done
$PYTHON - "$T" "$CSV" "$COND" "$SYS" "$S" <<'PY'
import sys, numpy as np
T, csvout, cond, sysn, slot = sys.argv[1:6]
c, r, l = (np.loadtxt(f"{T}/{n}.dat", skiprows=1)[:, 1] for n in ("cpx", "rec", "lig"))
b = r + l - c
with open(csvout, "w") as f:
    f.write("frame,complex,receptor_alone,ligand_alone,buried\n")
    for i in range(len(c)):
        f.write(f"{i+1},{c[i]:.3f},{r[i]:.3f},{l[i]:.3f},{b[i]:.3f}\n")
print(f"{cond}/{sysn}_{slot}: {len(c)} frames, buried {b.mean():.0f} +/- {b.std(ddof=1):.0f} A^2"
      + ("   *** some frames near zero - check imaging ***" if b.min() < 200 else ""))
PY
