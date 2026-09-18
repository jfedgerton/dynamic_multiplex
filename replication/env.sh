# source replication/env.sh   (on Roar, from the repo root)
# Sets DM_ROOT to this checkout and loads R. Every script and sbatch file
# reads DM_ROOT; nothing else is configured.
export DM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if command -v module >/dev/null 2>&1; then module load r/4.5.0; fi
echo "DM_ROOT=$DM_ROOT"
