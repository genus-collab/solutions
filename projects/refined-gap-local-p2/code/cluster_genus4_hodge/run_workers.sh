#!/bin/bash
# run_workers.sh -- launch one single-kernel wolframscript worker per allocated core inside a SLURM (array) task.
#   usage (from a job script):  bash run_workers.sh <script.m> <prefix args...> [-- <suffix args...>]
#   runs  wolframscript -file <script.m> <prefix args> <worker> <nWorkers> <suffix args>
#   with worker = (SLURM_ARRAY_TASK_ID - 1) * PER + k, k = 1..PER, PER = SLURM_CPUS_PER_TASK, nWorkers = NTASKS * PER
#   (export NTASKS in the job script = the --array size, 1 for a single job; default 20 x 6 = 120 workers).
# All stage scripts are restartable (a worker skips the items already in its output files); this launcher relies on that:
#   - the PER kernels of a task are started STAGGER seconds apart (on 2026-10-29, 120 kernels launched within 5 s on Stanage
#     and 9 of them never connected -- LinkConnect::linkc -- while their wolframscript processes hung for hours and kept the
#     array tasks alive);
#   - a worker that has printed nothing within START_TIMEOUT seconds (kernel never came up) is killed and relaunched;
#   - a worker that exits without its final line ("job finished" for hv_sample.m, "finished: ..." for hv_recon.m) is
#     relaunched; at most MAX_TRIES attempts per worker, after which the task exits with status 1 so that the afterok
#     dependency holds the rest of the chain.
# Output: worker_<script>_w<w>.out per worker (an existing one is moved to old_worker_out/ first), with [launcher] lines.
source "$(dirname "$0")/env.sh"
script=$1; shift
prefix=(); suffix=(); seen=0
for a in "$@"; do
  if [ "$a" == "--" ]; then seen=1; continue; fi
  if [ $seen -eq 0 ]; then prefix+=("$a"); else suffix+=("$a"); fi
done
PER=${SLURM_CPUS_PER_TASK:-6}
NT=${NTASKS:-20}
TID=${SLURM_ARRAY_TASK_ID:-1}
NW=$((NT * PER))
STAGGER=${STAGGER:-5}; START_TIMEOUT=${START_TIMEOUT:-600}; MAX_TRIES=${MAX_TRIES:-4}; POLL=${POLL:-10}
base=${script%.m}

has_output() { grep -v -E 'Failed to open configu|^[[:space:]]*$|LinkConnect' "$1" 2>/dev/null | grep -q .; }
is_complete() { grep -q -E 'job finished|^finished: ' "$1" 2>/dev/null; }
kill_tree() { # kill_tree <pid> <signal>: the process and all its descendants (wolframscript and its WolframKernel)
  local p=$1 sig=$2 c
  if command -v pgrep >/dev/null 2>&1; then for c in $(pgrep -P "$p" 2>/dev/null); do kill_tree "$c" "$sig"; done; fi
  kill -"$sig" "$p" 2>/dev/null
}

run_worker() { # run_worker <w>: at most MAX_TRIES attempts; returns 0 once the worker's output ends with its final line
  local w=$1 out="worker_${base}_w${w}.out" try pid waited tmp
  if [ -s "$out" ]; then mkdir -p old_worker_out; mv "$out" "old_worker_out/worker_${base}_w${w}.$(date +%Y%m%d_%H%M%S).out"; fi
  for try in $(seq 1 "$MAX_TRIES"); do
    tmp="${out}.try${try}"
    echo "[launcher] worker $w of $NW, attempt $try, $(date '+%F %T')" >> "$out"
    wolframscript -file "$script" "${prefix[@]}" "$w" "$NW" "${suffix[@]}" > "$tmp" 2>&1 &
    pid=$!
    waited=0
    while kill -0 "$pid" 2>/dev/null && ! has_output "$tmp"; do
      if [ "$waited" -ge "$START_TIMEOUT" ]; then
        echo "[launcher] worker $w: no kernel output after ${START_TIMEOUT} s, killing attempt $try ($(date '+%F %T'))" >> "$out"
        kill_tree "$pid" TERM; sleep 5; kill_tree "$pid" KILL
        break
      fi
      sleep "$POLL"; waited=$((waited + POLL))
    done
    wait "$pid" 2>/dev/null
    cat "$tmp" >> "$out" 2>/dev/null
    if is_complete "$tmp"; then
      rm -f "$tmp"
      echo "[launcher] worker $w finished on attempt $try, $(date '+%F %T')" >> "$out"
      return 0
    fi
    rm -f "$tmp"
    echo "[launcher] worker $w: attempt $try ended without the final line ($(date '+%F %T'))" >> "$out"
    sleep $((15 * try))
  done
  echo "[launcher] worker $w FAILED after $MAX_TRIES attempts, $(date '+%F %T')" >> "$out"
  return 1
}

lanes=(); fails=0
for k in $(seq 1 "$PER"); do
  w=$(( (TID - 1) * PER + k ))
  run_worker "$w" &
  lanes+=("$!")
  [ "$k" -lt "$PER" ] && sleep "$STAGGER"
done
for p in "${lanes[@]}"; do wait "$p" || fails=$((fails + 1)); done
echo "task $TID: workers $(( (TID - 1) * PER + 1 ))..$(( TID * PER )) of $NW finished, $fails failed after $MAX_TRIES attempts"
[ "$fails" -eq 0 ]
