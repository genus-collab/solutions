# -*- coding: utf-8 -*-
"""g5_sample.py -- per-inner-graph cost estimates by timed random sampling (no full graph is computed), or full timing of chosen graphs.

For an inner graph Gt the work of Split.graph is  ncombo * nsub  colour sums, ncombo = prod_v (#ordered half-edge decorations of the
vertex v), nsub = #{outer-edge subsets S with |S| < K}.  (At K = 2g - 2 every E_out(a, b) has u-valuation exactly 1 -- see g5_setup.log --
so the valuation pruning of Split.graph never removes anything beyond |S| < K.)  SplitTools.sample times random (decoration, S) pairs
and returns  est = t_prepare + ncombo * nsub * mean  (an unbiased estimator of the colour-sum time).
Modes:
  sample  [--slice i/n]   all inner graphs of the slice, a process pool, results appended to g5_sample_<tag>.jsonl
  resample --list t1,t2,..   the same for chosen graphs with a larger budget -> g5_resample_<tag>.jsonl
  full --list t1,t2,..    Split.graph on chosen graphs (one process each), measured seconds -> g5_full_<tag>.jsonl
Global wall-clock deadline (--deadline seconds, default 1380): the pool is terminated at the deadline; whatever is done is on disk.
Usage:  py -I g5_sample.py sample --genus 5 --K 8 --nproc 5 --budget 0.4
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

_SP = {}


def get(spec):
    sp = _SP.get(spec)
    if sp is None:
        from g5_lib import make_split, SplitTools
        from fq_fp import P1, P2
        prime, g, K, variant = spec
        s = make_split(P2 if prime == "89" else P1, g, (0, 0), K, variant); sp = (s, SplitTools(s)); _SP[spec] = sp
    return sp


def job_sample(args):
    spec, t, budget, mins, maxs, seed = args; sp, tools = get(spec); t0 = time.time()
    r = tools.sample(t, budget=budget, min_samples=mins, max_samples=maxs, seed=seed); r['wall'] = time.time() - t0
    return r


def job_full(args):
    spec, t = args; sp, tools = get(spec); G = sp.inner_graphs[t]; t0 = time.time(); res = sp.graph(t); dt = time.time() - t0
    return dict(t=t, E=G['E'], V=G['V'], gv=list(G['gv']), val=[len(h) for h in G['hsAt']], full=dt, nouter=len(res))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("--genus", type=int, default=5); ap.add_argument("--prime", default="89"); ap.add_argument("--K", type=int, default=None)
    ap.add_argument("--variant", default="nat"); ap.add_argument("--nproc", type=int, default=5); ap.add_argument("--budget", type=float, default=0.4)
    ap.add_argument("--min", type=int, default=3); ap.add_argument("--max", type=int, default=200); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--slice", default="0/1"); ap.add_argument("--list", default=""); ap.add_argument("--deadline", type=float, default=1380.0)
    ap.add_argument("--tag", default="")
    a = ap.parse_args(); assert a.nproc <= 6, "at most 6 processes"
    K = a.K if a.K is not None else 2 * a.genus - 2; spec = (a.prime, a.genus, K, a.variant)
    tag = "g%d_%s_%s_K%d%s" % (a.genus, a.prime, a.variant, K, a.tag)
    import multiprocessing as mp
    from g5_lib import fast_graphs
    n = len(fast_graphs(a.genus, 0, 3 * a.genus - 3)); t0 = time.time()
    if a.mode == "sample":
        i, m = [int(x) for x in a.slice.split("/")]; todo = [t for t in range(n) if t % m == i]
        path = os.path.join(HERE, "g5_sample_%s.jsonl" % tag); jobs = [(spec, t, a.budget, a.min, a.max, a.seed) for t in todo]; fn = job_sample
    elif a.mode == "resample":
        todo = [int(x) for x in a.list.split(",")]
        path = os.path.join(HERE, "g5_resample_%s.jsonl" % tag); jobs = [(spec, t, a.budget, a.min, a.max, a.seed) for t in todo]; fn = job_sample
    else:
        todo = [int(x) for x in a.list.split(",")]
        path = os.path.join(HERE, "g5_full_%s.jsonl" % tag); jobs = [(spec, t) for t in todo]; fn = job_full
    done = set()
    if os.path.exists(path) and a.mode == "sample":
        for line in open(path):
            if line.strip(): done.add(json.loads(line)['t'])
        jobs = [j for j in jobs if j[1] not in done]
    print("%s: %d graphs to do (%d already on disk), %d processes, deadline %.0f s -> %s" % (a.mode, len(jobs), len(done), a.nproc, a.deadline, path), flush=True)
    if not jobs: return
    out = open(path, "a"); pool = mp.Pool(processes=min(a.nproc, len(jobs))); ndone = 0; finished = False
    try:
        it = pool.imap_unordered(fn, jobs, chunksize=1)
        while True:
            left = a.deadline - (time.time() - t0)
            if left <= 0: break
            try: r = it.next(timeout=left)
            except StopIteration: finished = True; break
            except mp.TimeoutError: break
            out.write(json.dumps(r) + "\n"); out.flush(); ndone += 1
            if ndone % max(1, len(jobs) // 20) == 0: print("   %d/%d done (%.0f s)" % (ndone, len(jobs), time.time() - t0), flush=True)
    finally:
        pool.terminate(); pool.join(); out.close()
    print("%s %s: %d of %d graphs, wall %.0f s" % (a.mode, "COMPLETE" if finished else "STOPPED AT THE DEADLINE", ndone, len(jobs), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
