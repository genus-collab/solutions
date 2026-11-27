#!/bin/bash
# run_workers.sh -- launch one single-kernel wolframscript worker per allocated core inside a SLURM array task.
#   usage (from a job script):  bash run_workers.sh <script.m> <prefix args...> [-- <suffix args...>]
#   runs  wolframscript -file <script.m> <prefix args> <worker> <nWorkers> <suffix args>
#   with worker = (SLURM_ARRAY_TASK_ID - 1) * PER + k, k = 1..PER, PER = SLURM_CPUS_PER_TASK (6), nWorkers = NTASKS * PER (120),
#   NTASKS = number of array tasks (export NTASKS in the job script, must equal the --array size).
# All stage scripts are restartable: a resubmitted array simply skips the items already in the output files.
source "$(dirname "$0")/env.sh"
script=$1; shift
prefix=(); suffix=()
seen=0
for a in "$@"; do
  if [ "$a" == "--" ]; then seen=1; continue; fi
  if [ $seen -eq 0 ]; then prefix+=("$a"); else suffix+=("$a"); fi
done
PER=${SLURM_CPUS_PER_TASK:-6}
NT=${NTASKS:-20}
TID=${SLURM_ARRAY_TASK_ID:-1}
NW=$((NT * PER))
for k in $(seq 1 $PER); do
  w=$(( (TID - 1) * PER + k ))
  wolframscript -file "$script" "${prefix[@]}" $w $NW "${suffix[@]}" > "worker_${script%.m}_w${w}.out" 2>&1 &
done
wait
echo "task $TID: workers $(( (TID - 1) * PER + 1 ))..$(( TID * PER )) of $NW finished"
