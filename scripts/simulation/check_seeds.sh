#!/bin/bash
# check_seeds.sh [--dir <replica tree>] [condition ...]
# With ig=-1 AMBER draws the heating seed from the wall clock, so replicas launched in the same second could
# in principle collide and produce duplicate trajectories. Each heat.out records the seed it used; this
# collects them and reports any duplicate. Heating is the stage that matters: it builds the initial velocity
# distribution from that seed, while eq and production restart from the previous stage's velocities.
# Writes logs/seeds_<condition>.csv. Exit status 1 if any duplicate is found.
. "$(dirname "$0")/lib.sh"
seed_of() { grep -ah "based on wallclock" "$1" 2>/dev/null | sed 's/.*seed to *\([0-9]\+\).*/\1/' | head -1; }
rc=0
if [ "${1:-}" = "--dir" ]; then
  DIR=$2; echo "scanning $DIR"
  tmp=$(mktemp)
  for d in "$DIR"/replica*/; do s=$(seed_of "$d/heat.out"); [ -n "$s" ] && echo "$(basename $d),$s" >> "$tmp"; done
  n=$(wc -l < "$tmp"); u=$(cut -d, -f2 "$tmp" | sort -u | wc -l)
  dup=$(cut -d, -f2 "$tmp" | sort | uniq -d | tr '\n' ' ')
  printf "  %s heat.out with a seed, %s unique%s\n" "$n" "$u" "$([ -n "$dup" ] && echo "  DUPLICATES: $dup" || echo "")"
  if [ -n "$dup" ]; then rc=1
    for s in $dup; do echo "     seed $s used by: $(awk -F, -v s="$s" '$2==s {printf "%s ", $1}' "$tmp")"; done
  fi
  rm -f "$tmp"; exit $rc
fi
CONDS=${*:-$CONDITIONS}
for COND in $CONDS; do
  out=$CAMPAIGN/logs/seeds_${COND}.csv; echo "tag,replica,seed" > "$out"
  total=0; miss=0
  for TAG in $TAGS; do
    RUN=$(runs_dir "$COND" "$TAG"); [ -d "$RUN" ] || continue
    for d in "$RUN"/replica*/; do
      [ -d "$d" ] || continue
      s=$(seed_of "$d/heat.out")
      if [ -n "$s" ]; then echo "$TAG,$(basename $d),$s" >> "$out"; total=$((total+1)); else miss=$((miss+1)); fi
    done
  done
  dup=$(tail -n +2 "$out" | cut -d, -f3 | sort | uniq -d)
  printf "%-12s %4s seeds recorded, %4s replicas not yet heated, %4s unique" \
         "$COND" "$total" "$miss" "$(tail -n +2 "$out" | cut -d, -f3 | sort -u | wc -l)"
  if [ -n "$dup" ]; then
    rc=1; echo "   *** DUPLICATE SEEDS ***"
    for s in $dup; do echo "     seed $s used by: $(awk -F, -v s="$s" '$3==s {printf "%s/%s ", $1, $2}' "$out")"; done
  else echo "   no duplicates"; fi
  echo "             -> $out"
done
exit $rc
