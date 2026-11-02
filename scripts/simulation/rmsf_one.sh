#!/bin/bash
# rmsf_one.sh <system> <state: bound|unbound> <slot>   -- OPC reference condition only
. "$(dirname "$0")/lib.sh"
SYS=$1; STATE=$2; S=$(printf %03d $3)
OUT=$CAMPAIGN/analysis/rmsf/$STATE; mkdir -p "$OUT"
CSV=$OUT/${SYS}_${S}.csv; [ -s "$CSV" ] && exit 0
if [ "$STATE" = unbound ]; then
  TRJ=$CAMPAIGN/trajectories/opc-salt/${SYS}_kinase/${SYS}_kinase_${S}.nc
  PRM=$CAMPAIGN/topologies/opc-salt/${SYS}_kinase.prmtop; STRIP=""
else
  TRJ=$CAMPAIGN/trajectories/opc-salt/${SYS}_complex/${SYS}_complex_${S}.nc
  PRM=$CAMPAIGN/topologies/opc-salt/${SYS}_complex.prmtop; STRIP="strip :1-181"
fi
[ -f "$TRJ" ] || exit 1
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
$CPPTRAJ -p "$PRM" >/dev/null 2>&1 <<CP
trajin $TRJ 1 $PROD_FRAMES 10
$STRIP
rms first :1-59@CA
atomicfluct out $T/rmsf.dat :1-59@CA byres
secstruct :1-59 sumout $T/ss.sum
go
CP
$PYTHON - "$T" "$CSV" <<'PY'
import sys,numpy as np
T,out=sys.argv[1:3]
r=np.loadtxt(f'{T}/rmsf.dat',skiprows=1)
rows=[l.split() for l in open(f'{T}/ss.sum') if l.strip() and not l.startswith('#')]
hdr=[l for l in open(f'{T}/ss.sum') if l.startswith('#')][0].split()[2:]
d=np.array([[float(x) for x in q[1:]] for q in rows])
hel=d[:,[i for i,n in enumerate(hdr) if n in ('3-10','Alpha','Pi')]].sum(1)
with open(out,'w') as f:
    f.write('resid,rmsf,helix\n')
    for i in range(len(r)):
        f.write(f'{int(r[i,0])},{r[i,1]:.4f},{hel[i] if i<len(hel) else float("nan"):.4f}\n')
PY
