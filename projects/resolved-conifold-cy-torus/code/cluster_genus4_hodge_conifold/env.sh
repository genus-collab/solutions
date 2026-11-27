#!/bin/bash
# env.sh -- Mathematica environment for every job script and for run_workers.sh (sourced, never run).
# Initialises the module system if the shell has not done so (non-interactive batch shells), then loads Mathematica.
if ! command -v module >/dev/null 2>&1; then
  for f in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /usr/share/lmod/lmod/init/bash; do
    if [ -f "$f" ]; then source "$f"; break; fi
  done
fi
module load Mathematica/13.2.1
if ! command -v wolframscript >/dev/null 2>&1; then echo "env.sh: wolframscript not found after 'module load Mathematica/13.2.1'" >&2; fi
