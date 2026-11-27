# -*- coding: utf-8 -*-
"""g5l_fullrun.py -- local FULL-maxout (maxout = None) run of the package engine on the cheapest inner graphs that fit a budget.

Why: with --maxout m the pieces of an inner graph are only a partial sum over the outer-edge subsets, so the per-inner-graph
additivity of the method (sum over all subsets = point-frame term) cannot be tested.  This driver recomputes complete piece
families for as many inner graphs as fit, cheapest first by data/cost_table_g5.json, so that additivity can be tested exactly
on them (g5_additivity.py, check (A)) and so that the maxout-m pieces can be shown to be their restriction (check (C)).

It calls the package's own genus5_run.compute / pool_map (hence g5_engine.graph_dp with maxout = None) and writes the package's
own part format to parts/g5_89_nat/local.json (resumable, merged with whatever is already there).  It never runs a collect step,
so results/ is not touched.

  py -I -B g5l_fullrun.py [--prime 89] [--variant nat] [--budget 17000] [--nproc 13]
--budget is in seconds of the REFERENCE machine (the cost table's unit); the log prints the measured local time.
"""
import os, sys, json, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g5l_lib as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", default="89")
    ap.add_argument("--variant", default="nat")
    ap.add_argument("--budget", type=float, default=17000.0)
    ap.add_argument("--nproc", type=int, default=13)
    a = ap.parse_args()
    log, fh = L.log_opener(os.path.join(OUT, "fullrun_%s_%s.log" % (a.prime, a.variant)))
    import genus5_run as GR
    import g5_checks as C

    tag = C.tag_of(5, a.prime, a.variant, "0,0", None)
    spec = (a.prime, 5, "0,0", L.K, a.variant)
    pdir = os.path.join(L.PKG, "parts", tag)
    os.makedirs(pdir, exist_ok=True)
    path = os.path.join(pdir, "local.json")
    pieces, times = {}, {}
    if os.path.exists(path):
        prev = json.load(open(path))
        assert prev["tag"] == tag
        pieces, times = prev["pieces"], prev["times"]
        log("resuming from %s: %d inner graphs already stored" % (path, len(pieces)))

    costs = C.costs_for(5, L.NIN, a.variant)
    order = sorted(range(L.NIN), key=lambda t: (costs[t], t))
    todo, cum = [], 0.0
    for t in order:
        if str(t) in pieces:
            continue
        if cum + costs[t] > a.budget:
            break
        todo.append(t)
        cum += costs[t]
    log("tag %s (maxout = None), %d inner graphs to do, %.0f reference-seconds budgeted, %d processes" % (
        tag, len(todo), cum, a.nproc))
    if not todo:
        log("nothing to do")
        fh.close()
        return 0
    t0 = time.time()
    done = GR.pool_map([(spec, t, None) for t in todo], a.nproc, log, "inner graphs")
    for t, res, dt in done:
        pieces[str(t)] = res
        times[str(t)] = dt
    C.write_part(path, tag, pieces, times)
    cpu = sum(times[str(t)] for t in todo)
    log("local full-maxout run: %d inner graphs, wall %.0f s, CPU %.0f s (%.2f h); part file now holds %d inner graphs" % (
        len(todo), time.time() - t0, cpu, cpu / 3600, len(pieces)))
    log("   local/reference time ratio on this slice: %.3f" % (cpu / cum))
    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
