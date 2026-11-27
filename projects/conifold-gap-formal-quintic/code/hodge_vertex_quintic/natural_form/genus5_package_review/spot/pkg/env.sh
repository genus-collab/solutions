#!/bin/bash
# env.sh -- sourced (never run) by every job script: exports PY, a Python >= 3.8 interpreter. The package is pure Python (no packages).
# Order: the module named in G5_PYTHON_MODULE (if set), then python3 on the PATH, then the newest Python/3.x module that Lmod offers.
if ! command -v module >/dev/null 2>&1; then
  for f in /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
    if [ -f "$f" ]; then source "$f"; break; fi
  done
fi
g5_pyok() { "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' >/dev/null 2>&1; }
PY=""
if [ -n "$G5_PYTHON_MODULE" ]; then module load "$G5_PYTHON_MODULE"; fi
if command -v python3 >/dev/null 2>&1 && g5_pyok python3; then PY=python3; fi
if [ -z "$PY" ] && command -v module >/dev/null 2>&1; then
  for m in $(module -t avail Python 2>&1 | grep -E '^Python/3\.(8|9|1[0-9])\.' | sort -V -r); do
    if module load "$m" >/dev/null 2>&1 && g5_pyok python3; then PY=python3; echo "env.sh: loaded module $m" >&2; break; fi
    module unload "$m" >/dev/null 2>&1
  done
fi
if [ -z "$PY" ]; then
  echo "env.sh: no Python >= 3.8 found. Run 'module spider Python', pick a 3.x module and export G5_PYTHON_MODULE=<name> before sbatch." >&2
  return 1 2>/dev/null || exit 1
fi
export PY
echo "env.sh: using $($PY -c 'import sys; print(sys.executable, sys.version.split()[0])')" >&2
