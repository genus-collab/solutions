# -*- coding: utf-8 -*-
"""genus5_run.py -- driver of the genus-5 package (formal quintic, natural factorisation, edge splitting). Python >= 3.8, no packages.

  --validate [--nproc N]        genus 2, 3, 4 against the stored results of the project, genus-5 set-up, timing calibration (V_validate.sh)
  --task I --ntasks N           SLURM array task I (0-based) of N: the inner graphs of slice I (longest-processing-time assignment from
                                data/cost_table_g5.json), resumable graph by graph -> parts/<tag>/task<I>of<N>.json, log logs/<tag>/task<I>of<N>.log
  --run [--nproc N] [--list ..] local pool over all (or the listed) inner graphs -> parts/<tag>/local.json
  --collect                     merge parts/<tag>/*, run the checks -> results/<tag>.json, logs/<tag>/collect.log
Run options: --genus (5)  --prime 89|61  --variant nat|id|M0|nat_t<k>  --ab alpha,beta (0,0)  --maxout m (only outer-edge sets of size
<= m: exact for the vertex graph and the outer graphs with <= m edges; much cheaper; for controls).
tag = g<genus>_<prime>_<variant>[_ab..][_out<m>].  A genus-g term is u^(-(2g-2)) G(u); results are [u^0..u^(2g-3)] per outer graph.
"""
import os, sys, json, argparse, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import P1, P2
import g5_checks as C

_CACHE = {}


def get_split(spec):
    sp = _CACHE.get(spec)
    if sp is None:
        from g5_engine import make_split
        prime, g, ab, K, variant = spec
        sp = make_split(P2 if prime == "89" else P1, g, tuple(ab.split(",")), K, variant); _CACHE[spec] = sp
    return sp


def compute(spec, t, maxout):
    from g5_engine import graph_dp, graph_id
    sp = get_split(spec); t0 = time.time()
    if sp.variant == 'id': res = {str(t): graph_id(sp, t)}
    else: res = {str(q): v for q, v in graph_dp(sp, t, maxout).items() if any(v)}
    return res, time.time() - t0


def job(args):
    spec, t, maxout = args; res, dt = compute(spec, t, maxout)
    return t, res, dt


def pool_map(jobs, nproc, log=None, label=""):
    import multiprocessing as mp
    out = []; t0 = time.time()
    if nproc <= 1:
        for j in jobs: out.append(job(j))
        return out
    with mp.Pool(processes=min(nproc, len(jobs))) as pool:
        for k, r in enumerate(pool.imap_unordered(job, jobs)):
            out.append(r)
            if log and ((k + 1) % max(1, len(jobs) // 5) == 0 or k + 1 == len(jobs)): log("      %s %d/%d (%.0f s)" % (label, k + 1, len(jobs), time.time() - t0))
    return out


def opener(path):
    os.makedirs(os.path.dirname(path), exist_ok=True); f = open(path, "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); f.write(s + "\n"); f.flush()
    return log, f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", type=int, default=5); ap.add_argument("--prime", default="89"); ap.add_argument("--variant", default="nat")
    ap.add_argument("--ab", default="0,0"); ap.add_argument("--maxout", type=int, default=None); ap.add_argument("--nproc", type=int, default=1)
    ap.add_argument("--validate", action="store_true"); ap.add_argument("--run", action="store_true"); ap.add_argument("--collect", action="store_true")
    ap.add_argument("--task", type=int, default=None); ap.add_argument("--ntasks", type=int, default=None); ap.add_argument("--list", default="")
    a = ap.parse_args()
    if sys.version_info < (3, 8): sys.exit("genus5_run.py needs Python >= 3.8 (found %s)" % sys.version.split()[0])
    if a.validate:
        import g5_validate
        sys.exit(0 if g5_validate.main(a.nproc) else 1)
    if a.variant == 'id' and a.maxout is not None: sys.exit("--maxout is meaningless for the variant id")
    if (a.task is None) != (a.ntasks is None): sys.exit("--task and --ntasks go together")
    tag = C.tag_of(a.genus, a.prime, a.variant, a.ab, a.maxout); K = 2 * a.genus - 2
    spec = (a.prime, a.genus, a.ab, K, a.variant); pdir = os.path.join(HERE, "parts", tag); os.makedirs(pdir, exist_ok=True)
    if a.task is not None:
        log, f = opener(os.path.join(HERE, "logs", tag, "task%dof%d.log" % (a.task, a.ntasks)))
        log("=" * 90); log("genus5_run %s: task %d of %d, python %s, host %s" % (tag, a.task, a.ntasks, sys.version.split()[0], __import__("socket").gethostname()))
        t0 = time.time(); sp = get_split(spec); log("   set-up %.0f s: %d inner graphs, K = %d" % (time.time() - t0, len(sp.inner_graphs), sp.K))
        slices, load = C.lpt_assign(C.costs_for(a.genus, len(sp.inner_graphs), a.variant), a.ntasks); mine = slices[a.task]
        log("   %d inner graphs assigned, estimated %.2f h on the reference machine (heaviest task %.2f h)" % (len(mine), load[a.task] / 3600, max(load) / 3600))
        path = os.path.join(pdir, "task%dof%d.json" % (a.task, a.ntasks)); pieces = {}; times = {}
        if os.path.exists(path):
            prev = json.load(open(path)); assert prev['tag'] == tag; pieces = prev['pieces']; times = prev['times']; log("   resuming: %d done" % len(pieces))
        for t in mine:
            if str(t) in pieces: continue
            res, dt = compute(spec, t, a.maxout); pieces[str(t)] = res; times[str(t)] = dt; C.write_part(path, tag, pieces, times)
            G = sp.inner_graphs[t]; log("   inner graph %4d (E=%2d, V=%d): %7.1f s, %d outer graphs" % (t, G['E'], G['V'], dt, len(res)))
        log("TASK %d OF %d (%s) COMPLETE: %d inner graphs, %.0f s" % (a.task, a.ntasks, tag, len(mine), time.time() - t0)); f.close(); return
    if a.run:
        log, f = opener(os.path.join(HERE, "logs", tag, "local.log")); log("=" * 90); log("genus5_run %s: local run, %d processes" % (tag, a.nproc))
        t0 = time.time(); sp = get_split(spec); nin = len(sp.inner_graphs)
        path = os.path.join(pdir, "local.json"); pieces = {}; times = {}
        if os.path.exists(path): prev = json.load(open(path)); pieces = prev['pieces']; times = prev['times']
        todo = [int(x) for x in a.list.split(",")] if a.list else sorted(range(nin), key=lambda t: -C.costs_for(a.genus, nin, a.variant)[t])
        todo = [t for t in todo if str(t) not in pieces]
        for t, res, dt in pool_map([(spec, t, a.maxout) for t in todo], a.nproc, log, "inner graphs"):
            pieces[str(t)] = res; times[str(t)] = dt
        C.write_part(path, tag, pieces, times)
        log("   local run: %d inner graphs, wall %.0f s, CPU %.0f s" % (len(todo), time.time() - t0, sum(times[str(t)] for t in todo))); f.close()
        if a.list: return
    if a.collect or a.run:
        from fq3_rec import Rec
        log, f = opener(os.path.join(HERE, "logs", tag, "collect.log")); log("=" * 90); log("genus5_run %s: collect" % tag)
        sp = get_split(spec); names, pieces, times = C.read_parts(tag)
        log("   part files: %s; %d inner graphs; CPU %.1f h" % (names, len(pieces), sum(times.values()) / 3600))
        tab = C.merge(pieces, sp.K, sp.p)
        idtag = C.tag_of(a.genus, a.prime, 'id', a.ab, None); id_pieces = None
        if a.variant != 'id' and a.maxout is None:
            _, idp, _ = C.read_parts(idtag)
            if len(idp) == len(sp.inner_graphs): id_pieces = idp
            else: log("   (no complete %s run: the additivity check is skipped)" % idtag)
        FAILED, info = C.checks(sp, tab, set(int(t) for t in pieces), a.maxout, log, Rec(sp.R), id_pieces, pieces)
        os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
        json.dump(dict(tag=tag, table={str(i): tab[i] for i in sorted(tab)}, info={k: v for k, v in info.items()}, failed=FAILED,
                       outer=[dict(E=G['E'], gv=list(G['gv']), ed=[list(e) for e in G['ed']], aut=G['aut']) for G in sp.inner_graphs]),
                  open(os.path.join(HERE, "results", "%s.json" % tag), "w"))
        log("GENUS5_RUN COLLECT (%s) %s" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); f.close()


if __name__ == "__main__":
    main()
