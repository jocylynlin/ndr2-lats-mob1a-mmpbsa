#!/bin/bash
# check_frames.sh [condition ...]
# Reports the frame count of every trajectory in the campaign. Reused files are 10,000 frames (100 ns) and
# new ones 4,000 (40 ns); analyses must restrict reused files to their first 4,000 frames.
. "$(dirname "$0")/lib.sh"
printf "%-12s%-15s%8s%8s%10s%10s\n" condition tag files "4000" "10000" other
for COND in ${*:-$CONDITIONS}; do for TAG in $TAGS; do
  d=$(traj_dir "$COND" "$TAG"); [ -d "$d" ] || continue
  n4=0; n10=0; no=0
  for f in "$d"/*.nc; do
    [ -e "$f" ] || continue
    nf=$($PYTHON -c "from scipy.io import netcdf_file;f=netcdf_file('$f','r',mmap=True);print(f.variables['coordinates'].shape[0])" 2>/dev/null)
    case "$nf" in 4000) n4=$((n4+1));; 10000) n10=$((n10+1));; *) no=$((no+1));; esac
  done
  tot=$((n4+n10+no)); [ "$tot" -gt 0 ] && printf "%-12s%-15s%8s%8s%10s%10s\n" "$COND" "$TAG" "$tot" "$n4" "$n10" "$no"
done; done
echo "reused = 10,000 frames (100 ns, use first 4,000); new = 4,000 frames (40 ns)"
