#!/bin/bash
# Runs every job of jobs.txt through ref_point.pl, P jobs in parallel; one output file per job in runs/.
cd "$(dirname "$0")"
mkdir -p runs
P=${1:-22}
grep -v '^\s*$' jobs.txt | xargs -P "$P" -I{} sh -c 'f=$(echo "{}" | tr " " "_"); perl ref_point.pl {} > "runs/pt_$f.out" 2>&1'
echo "batch finished"
