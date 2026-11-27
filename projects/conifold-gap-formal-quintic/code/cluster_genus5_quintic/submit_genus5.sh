#!/bin/bash
# submit_genus5.sh -- submit genus-5 runs, each as one SLURM array plus a collect job that starts when the array has ended (afterany).
#   bash submit_genus5.sh main [NTASKS]       natural M at 2^89-1 and at 2^61-1 (default NTASKS = 100 each)
#   bash submit_genus5.sh controls [NTASKS]   bare symbol M_0 and tau_7 shifted (2^89-1, outer-edge sets of size <= 2: 20 tasks each),
#                                             and the natural M in the second conifold frame (1/3, -2/7) at 2^89-1 (NTASKS tasks)
#   bash submit_genus5.sh all [NTASKS]        both
# Run V_validate.sh first and check logs/validate.log.
set -e
cd "$(dirname "$0")"
mkdir -p logs/slurm
WHICH="${1:-main}"; NT="${2:-100}"
submit() {
  local nt="$1"; shift
  local jid cid
  jid=$(sbatch --parsable --array=0-$((nt - 1)) S_array.sh "$nt" "$@")
  cid=$(sbatch --parsable --dependency=afterany:"${jid%%;*}" C_collect.sh "$@")
  echo "$*: array job ${jid%%;*} ($nt tasks), collect job ${cid%%;*}"
}
case "$WHICH" in
  main|all)
    submit "$NT" --prime 89 --variant nat
    submit "$NT" --prime 61 --variant nat ;;
esac
case "$WHICH" in
  controls|all)
    submit 20 --prime 89 --variant M0 --maxout 2
    submit 20 --prime 89 --variant nat_t7 --maxout 2
    submit "$NT" --prime 89 --variant nat --ab 1/3,-2/7 ;;
esac
case "$WHICH" in main|controls|all) ;; *) echo "usage: bash submit_genus5.sh main|controls|all [NTASKS]" >&2; exit 1 ;; esac
