#!/bin/bash
# make_topologies.sh <condition>
# Builds the MM-PBSA topologies for a condition. The receptor and ligand topologies are cut FROM the stripped
# complex topology, never taken from the unbound builds: unbound topologies were built in separate tleap runs
# and carry different ZAFF-derived Lennard-Jones types, which produced an 18.6 kcal/mol artefact in the
# original analysis. Cutting from the complex guarantees matched parameters.
. "$(dirname "$0")/lib.sh"
COND=$1; TOP=$CAMPAIGN/topologies/$COND; mkdir -p "$TOP"
for SYS in lats1 ndr2; do
  CPX=$TOP/${SYS}_complex.prmtop
  [ -f "$CPX" ] || { echo "missing $CPX (produced by the first strip of ${SYS}_complex)"; continue; }
  NRES=$($PYTHON "$CAMPAIGN/scripts/solute_residues.py" "$CPX" | cut -d- -f2)
  # derive the receptor/ligand boundary from the molecule table instead of hard-coding it: solute molecule 1
  # is MOB1A with the Zn bonded in, molecule 2 is the kinase
  read -r RECRES RECAT LIGAT <<<"$($PYTHON "$CAMPAIGN/scripts/check_topologies.py" --ranges "$CPX")"
  REC=1-$RECRES; LIG=$((RECRES+1))-$NRES
  [ "$NRES" -gt "$RECRES" ] || { echo "$SYS: unexpected residue count $NRES vs receptor $RECRES"; continue; }
  $CPPTRAJ <<CP > "$CAMPAIGN/logs/topologies_${COND}_${SYS}.log" 2>&1
parm $CPX
parmstrip !(:$REC)
parmwrite out $TOP/${SYS}_receptor.prmtop
CP
  $CPPTRAJ <<CP >> "$CAMPAIGN/logs/topologies_${COND}_${SYS}.log" 2>&1
parm $CPX
parmstrip !(:$LIG)
parmwrite out $TOP/${SYS}_ligand.prmtop
CP
  if $PYTHON "$CAMPAIGN/scripts/check_topologies.py" "$CPX" "$TOP/${SYS}_receptor.prmtop" "$TOP/${SYS}_ligand.prmtop"; then
    echo "$COND $SYS: complex 1-$NRES, receptor :$REC ($RECAT atoms), ligand :$LIG ($LIGAT atoms) -> verified"
  else
    echo "$COND $SYS: TOPOLOGY CHECK FAILED - do not run MM-PBSA with these"; fi
done
