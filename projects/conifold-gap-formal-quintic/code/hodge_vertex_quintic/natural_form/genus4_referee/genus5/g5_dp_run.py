# -*- coding: utf-8 -*-
"""g5_dp_run.py -- command-line driver for the subset-DP prototype g5_dp.py (modes check4, check5, sample, full; see g5_dp.py).

Global wall-clock deadline (--deadline, default 1380 s): the pool is terminated at the deadline; finished graphs are on disk.
Usage:  py -I g5_dp_run.py check4 --nproc 5
        py -I g5_dp_run.py check5 --list 0,6,12 --nproc 3
        py -I g5_dp_run.py sample --nproc 5 --budget 0.5 [--slice i/n]
        py -I g5_dp_run.py full --list t1,t2 --nproc 5
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g5_dp


def run_pool(fn, jobs, nproc, deadline, t0, on_result):
    import multiprocessing as mp
    pool = mp.Pool(processes=min(nproc, len(jobs))); ndone = 0; finished = False
    try:
        it = pool.imap_unordered(fn, jobs, chunksize=1)
        while True:
            left = deadline - (time.time() - t0)
            if left <= 0: break
            try: r = it.next(timeout=left)
            except StopIteration: finished = True; break
            except mp.TimeoutError: break
            on_result(r); ndone += 1
            if ndone % max(1, len(jobs) // 10) == 0: print("   %d/%d done (%.0f s)" % (ndone, len(jobs), time.time() - t0), flush=True)
    finally:
        pool.terminate(); pool.join()
    return ndone, finished


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("--genus", type=int, default=5); ap.add_argument("--prime", default="89"); ap.add_argument("--K", type=int, default=None)
    ap.add_argument("--nproc", type=int, default=5); ap.add_argument("--budget", type=float, default=0.5); ap.add_argument("--min", type=int, default=3)
    ap.add_argument("--max", type=int, default=200); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--slice", default="0/1")
    ap.add_argument("--list", default=""); ap.add_argument("--deadline", type=float, default=1380.0); ap.add_argument("--tag", default="")
    a = ap.parse_args(); assert a.nproc <= 6, "at most 6 processes"; t0 = time.time()
    LOG = open(os.path.join(HERE, "g5_dp.log"), "a", encoding="utf-8")
    def log(s):
        print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    from g5_lib import fast_graphs, ROOT
    if a.mode == "check4":
        from fq3_graphs import legged_graphs, canon
        from fq_fp import P1, P2
        p = P2 if a.prime == "89" else P1
        old = legged_graphs(4, 0, 9); old_pos = {G['key']: i for i, G in enumerate(old)}
        new = fast_graphs(4, 0, 9); to_old = [old_pos[canon(*G['key'])] for G in new]
        stored = json.load(open(os.path.join(ROOT, "parts", "g4_%s_nat_run.json" % a.prime)))
        spec = (a.prime, 4, 6, 'nat'); bad = []; times = {}
        def on(r):
            t, res, dt = r; times[t] = dt
            mine = {to_old[int(k)]: v for k, v in res.items()}
            ref = {int(k): [x % p for x in v] for k, v in stored['res'][str(to_old[t])].items() if any(x % p for x in v)}
            if mine != ref: bad.append(t)
        nd, fin = run_pool(g5_dp.job_check4, [(spec, t, None) for t in range(379)], a.nproc, a.deadline, t0, on)
        cpu_old = sum(stored['times'].values())
        log("=" * 100)
        log("%s genus 4, 2^%s-1, natural M: subset-DP graph_dp = stored per-inner-graph tables of genus4_split for %d of 379 inner graphs (%d mismatches %s); "
            "CPU %.0f s (stored run: %.0f s under 13-fold load), slowest graph %.1f s, wall %.0f s"
            % ("PASS" if (not bad and nd == 379) else "FAIL", a.prime, nd, len(bad), bad[:10], sum(times.values()), cpu_old, max(times.values()), time.time() - t0))
        return
    K = a.K if a.K is not None else 2 * a.genus - 2; spec = (a.prime, a.genus, K, 'nat'); tag = "g%d_%s_nat_K%d%s" % (a.genus, a.prime, K, a.tag)
    n = len(fast_graphs(a.genus, 0, 3 * a.genus - 3))
    if a.mode == "check5":
        todo = [int(x) for x in a.list.split(",")]; rows = []
        nd, fin = run_pool(g5_dp.job_check5, [(spec, t) for t in todo], a.nproc, a.deadline, t0, rows.append)
        ok = all(r['equal'] for r in rows) and nd == len(todo)
        log("=" * 100)
        log("%s genus %d, 2^%s-1, natural M, K = %d: graph_dp == Split.graph on %d of %d inner graphs; sum of times: Split.graph %.0f s, graph_dp %.0f s (ratio %.1f)"
            % ("PASS" if ok else "FAIL", a.genus, a.prime, K, nd, len(todo), sum(r['std'] for r in rows), sum(r['dp'] for r in rows),
               sum(r['std'] for r in rows) / max(1e-9, sum(r['dp'] for r in rows))))
        for r in sorted(rows, key=lambda r: r['t']):
            log("   t=%4d E=%2d V=%d val=%s: Split.graph %.1f s, graph_dp %.2f s, ratio %.1f, %d outer graphs, equal %s" % (r['t'], r['E'], r['V'], r['val'], r['std'], r['dp'], r['std'] / max(r['dp'], 1e-9), r['nouter'], r['equal']))
        return
    if a.mode == "sample":
        i, m = [int(x) for x in a.slice.split("/")]; todo = [t for t in range(n) if t % m == i]
        path = os.path.join(HERE, "g5_dpsample_%s.jsonl" % tag); done = set()
        if os.path.exists(path):
            for line in open(path):
                if line.strip(): done.add(json.loads(line)['t'])
        jobs = [(spec, t, a.budget, a.min, a.max, a.seed) for t in todo if t not in done]; fn = g5_dp.job_sample
    elif a.mode == "resample":
        path = os.path.join(HERE, "g5_dpresample_%s.jsonl" % tag)
        jobs = [(spec, int(x), a.budget, a.min, a.max, a.seed) for x in a.list.split(",")]; fn = g5_dp.job_sample
    else:
        path = os.path.join(HERE, "g5_dpfull_%s.jsonl" % tag); jobs = [(spec, int(x)) for x in a.list.split(",")]; fn = g5_dp.job_full
    print("%s: %d graphs, %d processes -> %s" % (a.mode, len(jobs), a.nproc, path), flush=True)
    if not jobs: return
    out = open(path, "a")
    def on(r): out.write(json.dumps(r) + "\n"); out.flush()
    nd, fin = run_pool(fn, jobs, a.nproc, a.deadline, t0, on); out.close()
    print("%s %s: %d of %d graphs, wall %.0f s" % (a.mode, "COMPLETE" if fin else "STOPPED AT THE DEADLINE", nd, len(jobs), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
