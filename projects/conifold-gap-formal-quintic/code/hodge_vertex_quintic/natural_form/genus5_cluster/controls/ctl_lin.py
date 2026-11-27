# -*- coding: utf-8 -*-
"""ctl_lin.py -- is the response of the test to a shift of tau_k the predicted LINEAR one?

   The prediction of the report: with tau' = tau + c u^k one has, mod u^(2k),
       E_M' - E_M = (a+b+1) tau^(a+b) c u^k (Eh_0(a, b)),      E_out' - E_out = -(E_M' - E_M),
   so every visible perturbation of an outer graph row is LINEAR in c as long as 2k >= 2g-2 (no second-order term fits
   inside the window u^0..u^(2g-3)) and the row is perturbed in exactly one edge factor.
   Checked at genus 3 (window u^0..u^3) for k = 1, 2, 3 and c = 1, 2, 3: expect exact linearity for k = 2, 3 and
   a visible quadratic term for k = 1 (2k = 2 < 4).

   py -I -B ctl_lin.py [prime]            (no files are written; nothing in the package is modified)
"""
import os, sys
from fractions import Fraction as Fr
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)

G = 3
K = 2 * G - 2


def table(prime, k, c):
    """full genus-3 outer table with tau_k -> tau_k + c (c = 0: the natural tau), by patching tau_pred in memory"""
    import fq4_split, genus5_run as D, g5_checks as C
    orig = fq4_split.tau_pred
    if c:
        def patched(n, _o=orig, _k=k, _c=c):
            t = list(_o(n))
            if _k < n: t[_k] = t[_k] + Fr(_c)
            return t
        fq4_split.tau_pred = patched
    D._CACHE.clear()
    try:
        spec = (prime, G, "0,0", K, "nat")
        sp = D.get_split(spec)
        pieces = {str(t): D.compute(spec, t, None)[0] for t in range(len(sp.inner_graphs))}
        tab = C.merge(pieces, K, sp.p)
    finally:
        fq4_split.tau_pred = orig
        D._CACHE.clear()
    return sp, tab


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"
    print("=" * 100)
    print("genus 3, 2^%s-1: response of the outer graph rows to tau_k -> tau_k + c" % prime)
    print("=" * 100)
    sp0, base = table(prime, 0, 0)
    p = sp0.p
    E = [g['E'] for g in sp0.inner_graphs]
    nz = [i for i in range(1, len(E)) if any(base.get(i, [0] * K))]
    print("natural tau: non-regular outer graphs with edges = %s (expected none)" % (nz or "none"))
    for k in (1, 2, 3):
        print()
        print("-" * 100)
        print("k = %d  (prediction: visible perturbation linear in c iff 2k >= %d)" % (k, K))
        ref = None
        for c in (1, 2, 3):
            sp, tab = table(prime, k, c)
            d = {i: [(tab.get(i, [0] * K)[j] - base.get(i, [0] * K)[j]) % p for j in range(K)] for i in range(len(E))}
            bad = {i: [j for j in range(K) if tab.get(i, [0] * K)[j]] for i in range(1, len(E)) if any(tab.get(i, [0] * K))}
            if ref is None:
                ref = d
                print("   c = 1: non-regular graphs %s" % sorted(bad))
            lin = all(d[i][j] == ref[i][j] * c % p for i in range(len(E)) for j in range(K))
            firstbad = {i: min(bad[i]) for i in bad}
            print("   c = %d: non-regular graphs %s; all deltas = c x (delta at c=1): %s"
                  % (c, sorted(bad), "YES (exactly linear)" if lin else "NO (non-linear terms present)"))
            if c == 3:
                for i in sorted(bad):
                    print("      graph %2d E=%d gv=%s: delta/c at the first failing order u^%d linear: %s"
                          % (i, E[i], list(sp.inner_graphs[i]['gv']), firstbad[i],
                             d[i][firstbad[i]] == ref[i][firstbad[i]] * c % p))


if __name__ == "__main__":
    main()
