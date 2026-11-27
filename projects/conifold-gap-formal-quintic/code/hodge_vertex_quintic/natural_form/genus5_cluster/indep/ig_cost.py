# -*- coding: utf-8 -*-
"""ig_cost.py -- static cost proxy for the INDEPENDENT code (genus4_referee/indep) at genus 5, calibrated on the
17 inner graphs that natural_form/genus5_package_review/spot measured.

The inner loop of ind_sum.Summer.graph runs over decorations (one ordered exponent tuple per half-edge, from
ind_tables.vertex_table) and, inside each decoration, over the subsets S of outer edges and the 5^(E-V+1) momentum
configurations.  Proxy:
      cost(t) = Ndec(t) * ( Nsub(t) * 5^(E-V+1) + c0 )
with Ndec = prod_v #{permutations of the exponent tuples of vertex_table(g_v, valence_v)} and Nsub = #subsets with
|S| <= K-1 = 7.  Only the vertex tables are needed, so this is cheap (no graph is summed).
  py -I -B ig_cost.py [--prime 89]
Writes out/cost_ind.json: {t: [Ndec, Nsub, pow5, proxy]} plus the calibration print-out.
"""
import os, sys, json, time, itertools, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
IND = os.path.normpath(os.path.join(HERE, "..", "..", "genus4_referee", "indep"))
PKG = os.environ.get("G5PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic"))
sys.path.insert(0, IND)

# the 17 graphs of spot/SPOT_G5.md s.3 with the independent-code seconds from out/spot3_89.log
SPOT = {1: 0.8, 3: 4.5, 5: 0.0, 25: 0.8, 28: 1.7, 118: 12.7, 119: 348.4, 383: 0.9, 444: 3.8, 1234: 19.1,
        2318: 13.8, 2539: 18.1, 3638: 186.3, 3847: 348.9, 4356: 72.0, 4484: 31.7, 4552: 45.6}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", default="89")
    a = ap.parse_args()
    from ind_tables import Tables
    t0 = time.time()
    T = Tables(a.prime, 5, 8, "nat")
    print("Tables(%s, g=5, K=8, nat) built in %.1f s; all assertions passed" % (a.prime, time.time() - t0))

    graphs = json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))
    cost = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))

    # number of ordered exponent tuples per (genus, valence)
    perms = {}

    def nperm(g, hs):
        key = (g, hs)
        if key not in perms:
            tab = T.vertex_table(g, hs)
            perms[key] = sum(len(set(itertools.permutations(hp))) for hp in tab)
        return perms[key]

    res = {}
    for t, (gk, aut) in enumerate(graphs):
        gv, ed, _ = gk
        V = len(gv)
        E = len(ed)
        val = [0] * V
        for u, w in ed:
            val[u] += 1
            val[w] += 1
        nd = 1
        for v in range(V):
            nd *= nperm(gv[v], val[v])
        nsub = sum(1 for S in range(1 << E) if bin(S).count("1") <= 7)
        p5 = 5 ** (E - V + 1)
        res[t] = [nd, nsub, p5, nd * (nsub * p5 + 1)]
    json.dump({str(t): res[t] for t in res}, open(os.path.join(OUT, "cost_ind.json"), "w"))

    print("\ncalibration on the 17 spot graphs (seconds per unit of proxy):")
    print("%6s %3s %3s %12s %7s %8s %12s %10s %8s" % ("t", "E", "V", "Ndec", "Nsub", "5^(E-V+1)", "proxy", "ind s", "s/Mproxy"))
    rates = []
    for t in sorted(SPOT):
        nd, nsub, p5, pr = res[t]
        s = SPOT[t]
        r = s / (pr / 1e6) if pr else 0
        rates.append((r, t))
        print("%6d %3d %3d %12d %7d %8d %12d %10.1f %8.3f" % (t, len(graphs[t][0][1]), len(graphs[t][0][0]), nd, nsub, p5, pr, s, r))
    rates.sort()
    print("\ns per 1e6 proxy: min %.3f (t=%d), median %.3f, max %.3f (t=%d)" %
          (rates[0][0], rates[0][1], rates[len(rates) // 2][0], rates[-1][0], rates[-1][1]))
    hi = max(r for r, _ in rates if r)
    print("conservative rate used for planning: %.3f s per 1e6 proxy" % hi)
    print("\nhow the proxy ranks the whole list (top 12 by proxy):")
    top = sorted(res, key=lambda t: -res[t][3])[:12]
    for t in top:
        print("   t=%4d E=%2d V=%d proxy=%.3e  est_dp=%.0f s  predicted ind %.0f s"
              % (t, len(graphs[t][0][1]), len(graphs[t][0][0]), res[t][3], cost[str(t)]["est_dp"], hi * res[t][3] / 1e6))


if __name__ == "__main__":
    main()
