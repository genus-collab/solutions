# -*- coding: utf-8 -*-
"""ig_plan.py -- choose the sample of inner graphs for the independent genus-5 verification, under a CPU budget.

Rules (deliberate, not random):
  * COMPLETE STRATA at the low end: every inner graph with E <= 5 edges, so that nothing with few edges is left out.
    These are also the graphs whose pieces reach the low-E outer graphs, i.e. the ones that are regular only by cancellation.
  * the loop-heaviest graph 149 (5 loops on one genus-0 vertex, |Aut| = 3840) is included but split into decoration shards:
    it is one of the three graphs (149, 245, 2668) that no independent code had ever covered.
  * COMPLETE STRATUM at the top: every inner graph with E = 12 (the trivalent frontier), which is cheap for this algorithm.
  * the PACKAGE-HEAVIEST graphs: the n_heavy inner graphs with the largest cluster wall time in parts/g5_89_nat.
  * a COST-SPREAD sample in every stratum E = 6..11: the stratum is sorted by predicted independent cost and sampled at
    evenly spaced quantiles, so cheap, median and expensive graphs all occur, until the stratum budget is spent.
Predicted costs come from the proxy of ig_cost.py (rate 5e-6 s per proxy unit, calibrated on the 17 spot graphs).
  py -I -B ig_plan.py [--budget 56000] [--p61budget 20000] [--n-heavy 20]
Writes sample_nat89.txt, sample_nat61.txt, sample_id89.txt, sample_id61.txt and prints the accounting.
"""
import os, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RATE = 5e-6
SHARD_GRAPHS = {149: 10}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=float, default=56000.0)
    ap.add_argument("--p61budget", type=float, default=20000.0)
    ap.add_argument("--n-heavy", type=int, default=20)
    a = ap.parse_args()
    ix = json.load(open(os.path.join(OUT, "index.json")))
    ci = json.load(open(os.path.join(OUT, "cost_ind.json")))
    G = {int(t): v for t, v in ix["graphs"].items()}
    pred = {t: max(1.0, RATE * ci[str(t)][3]) for t in G}
    pid = {t: max(2.0, RATE * ci[str(t)][0] * ci[str(t)][2] * 2) for t in G}
    byE = {}
    for t in G:
        byE.setdefault(G[t]["E"], []).append(t)

    chosen = set()
    why = {}

    def take(t, tag):
        if t not in chosen:
            chosen.add(t)
            why[t] = tag

    for e in range(0, 6):
        for t in byE[e]:
            take(t, "E<=5 complete")
    for t in byE[12]:
        take(t, "E=12 complete")
    heavy = sorted(G, key=lambda t: -G[t]["t89"])[:a.n_heavy]
    for t in heavy:
        take(t, "package-heaviest")
    base = sum(pred[t] for t in chosen)
    print("complete strata E <= 5 (%d graphs) + E = 12 (%d) + %d package-heaviest: %d graphs, predicted %.0f s"
          % (sum(len(byE[e]) for e in range(6)), len(byE[12]), a.n_heavy, len(chosen), base))
    print("   (of which graph 149 alone: %.0f s, run in %d shards)" % (pred[149], SHARD_GRAPHS[149]))

    left = a.budget - base
    strata = [6, 7, 8, 9, 10, 11]
    share = {6: 0.16, 7: 0.17, 8: 0.19, 9: 0.19, 10: 0.17, 11: 0.12}
    for e in strata:
        bud = left * share[e]
        pool = sorted((t for t in byE[e] if t not in chosen), key=lambda t: pred[t])
        n = len(pool)
        spent = 0.0
        k = 0
        # evenly spaced quantiles of the cost order, widening until the budget is spent
        step = max(1, n // 200)
        order = []
        for off in range(step):
            order += pool[off::step]
        for t in order:
            if spent + pred[t] > bud:
                continue
            take(t, "E=%d cost-spread" % e)
            spent += pred[t]
            k += 1
        print("   E=%2d: %4d graphs in the stratum, %3d taken, predicted %.0f s of the %.0f s stratum budget"
              % (e, len(byE[e]), k, spent, bud))
    tot = sum(pred[t] for t in chosen)
    print("SAMPLE p89 nat: %d inner graphs (%.1f%% of 4555), predicted %.0f CPU-s = %.1f CPU-h"
          % (len(chosen), 100.0 * len(chosen) / len(G), tot, tot / 3600))
    npieces = sum(G[t]["npieces"] for t in chosen)
    print("   pieces to be compared: %d (of 168475 in the full run); outer graphs reached: %d"
          % (npieces, len(set(q for t in chosen for q in G[t]["qs"]))))
    byEc = {}
    for t in chosen:
        byEc.setdefault(G[t]["E"], []).append(t)
    print("   by E: " + ", ".join("E=%d:%d/%d" % (e, len(byEc.get(e, [])), len(byE[e])) for e in sorted(byE)))

    # second prime: the cheap complete strata plus a spread, under its own budget
    s61 = set(t for t in chosen if G[t]["E"] <= 4)
    s61 |= set(heavy[:8])
    spent = sum(pred[t] for t in s61)
    for e in range(5, 13):
        pool = sorted((t for t in byEc.get(e, []) if t not in s61), key=lambda t: pred[t])
        step = max(1, len(pool) // 40)
        for t in pool[::step]:
            if spent + pred[t] <= a.p61budget:
                s61.add(t)
                spent += pred[t]
    print("SAMPLE p61 nat: %d inner graphs, predicted %.0f CPU-s = %.1f CPU-h; pieces %d"
          % (len(s61), spent, spent / 3600, sum(G[t]["npieces"] for t in s61)))

    def write(name, items, note):
        with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
            f.write("# %s\n" % note)
            for t in sorted(items):
                f.write("%d   # E=%d V=%d aut=%d pred=%.0fs %s\n" % (t, G[t]["E"], G[t]["V"], G[t]["aut"], pred[t], why.get(t, "")))
        print("wrote %s (%d graphs)" % (name, len(items)))

    nat89 = sorted(t for t in chosen if t not in SHARD_GRAPHS)
    write("sample_nat89.txt", nat89, "p89 natural-M sample, unsharded (graph 149 is run separately with --shards)")
    write("sample_nat61.txt", sorted(s61 - set(SHARD_GRAPHS)), "p61 natural-M sample")
    id89 = sorted(chosen)
    write("sample_id89.txt", id89, "p89 point-frame (variant id) sample = the whole p89 sample")
    write("sample_id61.txt", sorted(s61), "p61 point-frame sample")
    print("predicted id cost: p89 %.0f s, p61 %.0f s" % (sum(pid[t] for t in id89), sum(pid[t] for t in s61)))
    print("predicted grand total: %.1f CPU-h" % ((tot + spent + sum(pid[t] for t in id89) + sum(pid[t] for t in s61)) / 3600))
    json.dump(dict(nat89=nat89, nat61=sorted(s61 - set(SHARD_GRAPHS)), id89=id89, id61=sorted(s61),
                   shards={str(k): v for k, v in SHARD_GRAPHS.items()}, why={str(t): why[t] for t in chosen},
                   pred={str(t): pred[t] for t in chosen}),
              open(os.path.join(OUT, "sample.json"), "w"))


if __name__ == "__main__":
    main()
