# shared helpers
set -u
CONF="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/campaign.conf"; . "$CONF"
log() { echo "$(date '+%F %T')  $*" | tee -a "${LOG:-$CAMPAIGN/logs/campaign.log}"; }
inputs_dir() { echo "$CAMPAIGN/inputs/$1/$2"; }
runs_dir()   { echo "$CAMPAIGN/runs/$1/$2"; }
traj_dir()   { echo "$CAMPAIGN/trajectories/$1/$2"; }
prmtop_of()  { ls "$(inputs_dir "$1" "$2")"/*.prmtop 2>/dev/null | head -1; }
# how many slots are already filled by reused trajectories
reused_count() { awk -F, -v c="$1" -v t="$2" 'NR>1 && $1==c && $4==t {print $13}' "$CAMPAIGN/MANIFEST.csv"; }
free_gb() { df -B1G --output=avail /home3 | tail -1 | tr -d ' '; }
