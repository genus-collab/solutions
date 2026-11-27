cluster_genus5_quintic -- genus 5 of the formal quintic with the natural factorisation M = M_0(tau(u) zhat), on Stanage (SLURM)
=================================================================================================================================

WHAT IS COMPUTED
  The Givental-Teleman graph sum of the formal quintic at genus 5 near the conifold point u = 1 - 5^5 q = 0, reorganised as the outer
  graph sum R'.(M.omega^) with M = M_0(tau(u) zhat), 1/tau = 2F1(1, 4/5; 9/5; u) (hodge_vertex_quintic/NATURAL_FORM.md). Method: edge
  splitting (hodge_vertex_quintic/fq4_split.py, NATURAL_FORM.md s.5) with one vertex elimination for all outer-edge subsets (g5_engine.py).
  A genus-5 term is u^-8 G(u); the run gives [u^0..u^7] G for each of the 4555 outer graphs. Conjecture N at genus 5 = every outer graph
  with edges is O(u^8) and the vertex graph carries the whole polar part c_5 t_c^-8.
  Exact arithmetic modulo 2^89-1 ("89") and 2^61-1 ("61"); pure Python >= 3.8, no packages; single core per array task.

ALREADY DONE ON THE WORKSTATION (shipped in parts/ and results/)
  Point-frame total (variant id) at both primes: total = (1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1 (the gap at genus 5;
  1/660000 = |B_10|/(10 . 8 . 5^4)). These runs are also the reference for the ADDITIVITY check of the collect step (for every inner graph
  the pieces of the natural run must add up to its point-frame term). Do not delete parts/g5_89_id/ and parts/g5_61_id/.

FILES
  genus5_run.py      driver: --validate | --task I --ntasks N | --run | --collect      (see its docstring)
  g5_engine.py       Split with the fast graph list, graph_dp (one elimination for all subsets), graph_id
  g5_graphs.py       fast canonical form; graph lists in data/graphs_<g>_0_<Emax>.json (genus 2..5; do NOT regenerate: indices matter)
  g5_checks.py       tags, task assignment (longest first, from data/cost_table_g5.json), part files, merge, checks
  g5_validate.py     validation + calibration
  fq*.py, gsum_q.py  unchanged copies of the project modules (hodge_vertex_quintic, 2026-11-11)
  make_refdata.py    how data/ref_g23.json, data/ref_g4.json were made (needs the project folder; NOT for the cluster)
  env.sh             finds a Python >= 3.8 (python3 on the PATH, else the newest Python/3.x module; override: export G5_PYTHON_MODULE=...)
  V_validate.sh      validation job (16 cores, about 10 minutes)
  S_array.sh         one array task;  C_collect.sh  merge + checks;  submit_genus5.sh  submits arrays + collect jobs

STEPS ON STANAGE
  1. Validate:   sbatch V_validate.sh        then   tail -3 logs/validate.log   -> "G5_VALIDATE ALL CHECKS PASSED"
     It re-checks genus 2, 3 against the old vertex-tensor code, eight genus-4 runs against the tables of record, the genus-5 set-up,
     and times three genus-5 graphs: the line "machine factor f" scales all times below (workstation: f = 1 by definition).
  2. Main runs:  bash submit_genus5.sh main          natural M at both primes, 100 one-core tasks each + a collect job each.
     Estimate on the workstation: 39-46 CPU hours per prime; with 100 tasks the heaviest task is about 0.6 h (x f), most about 0.4 h.
  3. Results:    tail -3 logs/g5_89_nat/collect.log  and  logs/g5_61_nat/collect.log   -> "GENUS5_RUN COLLECT (...) ALL CHECKS PASSED"
     The collect log lists: coverage of the 4555 inner graphs, E_out = O(u), the vertex graph in t_c (leading coefficient compared with
     1/660000), REGULARITY of all outer graphs with edges (any non-regular graph is printed with its coefficients), the GAP of the total,
     ADDITIVITY against the shipped id run. results/<tag>.json holds the per-outer-graph table.
  4. Controls (optional; bash submit_genus5.sh controls): bare symbol M_0 and tau_7 shifted by 1, both with --maxout 2 (exact for the
     vertex graph and the 1- and 2-edge outer graphs; these MUST FAIL the regularity check: that shows the test discriminates), and the
     natural M in the second conifold frame (1/3, -2/7) at 2^89-1 (full size, must pass). In the second frame the additivity check is
     skipped (no id run in that frame).

FAILURES AND RESUBMISSION
  A task that hits its time limit or dies keeps every inner graph it finished (parts/<tag>/task<I>of<N>.json, rewritten after each
  graph). Resubmit only the missing indices with the SAME N, e.g.
      sbatch --array=17,42 S_array.sh 100 --prime 89 --variant nat
  then rerun the collect job:  sbatch C_collect.sh --prime 89 --variant nat
  The collect log names the missing inner graphs; the task of inner graph t is found in logs/<tag>/task<I>.log ("assigned").
  Never mix different N for the same tag unless all tasks of the old N are complete (the merge tolerates duplicates with identical content).

COPY BACK
  The results needed on the workstation are small: results/, logs/ (without logs/slurm if large), and parts/ (per-inner-graph pieces,
  tens of MB per run).
