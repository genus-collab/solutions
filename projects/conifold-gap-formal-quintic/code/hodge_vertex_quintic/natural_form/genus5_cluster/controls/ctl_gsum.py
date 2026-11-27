# -*- coding: utf-8 -*-
"""ctl_gsum.py -- full outer graph sums at a small genus for every tau control, to test the prediction
   "nat_tk makes exactly the outer graphs with k + E - 1 <= 2g-3 fail, from the order u^(k+E-1) upwards".

   py -I -B ctl_gsum.py <genus> [prime] [nproc]      e.g.  py -I -B ctl_gsum.py 3 89 1

   Uses the shipped package cluster_genus5_quintic unchanged (genus5_run.compute, g5_checks).  Nothing is written.
"""
import os, sys, collections, time
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)


def run(g, prime, variant, nproc=1):
    import genus5_run as D, g5_checks as C
    K = 2 * g - 2
    spec = (prime, g, "0,0", K, variant)
    sp = D.get_split(spec)
    n = len(sp.inner_graphs)
    pieces = {}
    for t, res, dt in D.pool_map([(spec, t, None) for t in range(n)], nproc):
        pieces[str(t)] = res
    tab = C.merge(pieces, K, sp.p)
    return sp, tab


def polar(sp, row):
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3
    tc = tc_series(sp.C, sp.R, sp.K + 4)
    return tc_expansion3(row, tc, sp.R, sp.g)[: 2 * sp.g - 2]


def main():
    g = int(sys.argv[1]); prime = sys.argv[2] if len(sys.argv) > 2 else "89"
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    K = 2 * g - 2
    print("=" * 104)
    print("genus %d, 2^%s-1: full outer graph sums, regularity = row is O(u^%d)" % (g, prime, K))
    print("=" * 104)
    variants = ["nat", "M0"] + ["nat_t%d" % k for k in range(1, K)]
    first = True
    for variant in variants:
        t0 = time.time()
        sp, tab = run(g, prime, variant, nproc)
        E = [G['E'] for G in sp.inner_graphs]
        if first:
            print("graphs by E: %s (total %d)" % (sorted(collections.Counter(E).items()), len(E)))
            first = False
        bad = {i: [k for k in range(K) if tab.get(i, [0] * K)[k]] for i in range(1, len(E)) if any(tab.get(i, [0] * K))}
        V = polar(sp, tab.get(0, [0] * K))
        kk = 1 if variant == "M0" else (int(variant[5:]) if variant.startswith("nat_t") else None)
        pred = sorted(i for i in range(1, len(E)) if kk is not None and kk + E[i] - 1 <= K - 1)
        print()
        print("-" * 104)
        print("%-8s (%.0f s)  non-regular outer graphs with edges: %d" % (variant, time.time() - t0, len(bad)))
        print("   vertex graph polar part in t_c (t_c^-%d .. t_c^-1): %s" % (K, ["0" if x == 0 else "NONZERO" for x in V]))
        for i in sorted(bad):
            G = sp.inner_graphs[i]
            print("      %3d E=%d gv=%s ed=%s |Aut|=%d  failing u-orders %s  (predicted first order %s)"
                  % (i, E[i], list(G['gv']), [list(e) for e in G['ed']], G['aut'], bad[i],
                     (kk + E[i] - 1) if kk is not None else "-"))
        if kk is not None:
            agree = sorted(bad) == pred and all(min(bad[i]) == kk + E[i] - 1 and bad[i] == list(range(kk + E[i] - 1, K)) for i in bad)
            print("   PREDICTION (fail set = {E <= %d}, orders u^(k+E-1)..u^%d with k = %d): %s"
                  % (K - kk, K - 1, kk, "MATCHES" if agree else "DIFFERS, predicted fail set %s" % pred))


if __name__ == "__main__":
    main()
