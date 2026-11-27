# -*- coding: utf-8 -*-
"""ref_run2b.py -- REFEREE: the weight decomposition of G and log G.

Weight of a monomial u^n zhat^k in the flat-bivector component K^a (x) K^b (a = -m-1, b = -m'-1 for the Fourier component (m, m')):
      j = 5n - 4k + 2 - w_a - w_b,   w = (0, 4, 3, 2, 1)     (j is a multiple of 5 by the Z/5 selection rule; rule W: j >= 0).
In the un-hatted variables the monomial is z^k s^j K^a (x) K^b with the frame K = (1, u H, H^2, H^3, H^4) (hatted normalisations
understood), so the weight-j part is the coefficient of s^j of a formal power series in z.
Prints, for X = G and X = log G: for every zhat-order k the components present at weight 0 (the 'leading part' Gamma_0(z) = X at
s = 0, z fixed) with their rational values (prime 2^521-1), and the table n_min for log G.
Run: py -I ref_run2b.py [--NZ 16] [--K 12] [--prime 521]
"""
import os, sys, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, ratrec, blog
from ref_hat import HatRows, tau_pred, M0_of, G_of, fourier

PR = {'89': P89, '61': P61, '521': P521}
W = (0, 4, 3, 2, 1)


def weights(X, z, p):
    """{(j, k, a, b): (n, value)}"""
    out = {}
    for (k, n), A in X.items():
        if (k, n) == (0, 0): continue
        for (m, mp), v in fourier(A, z, p).items():
            aa = (-m - 1) % 5; bb = (-mp - 1) % 5
            out[(5 * n - 4 * k + 2 - W[aa] - W[bb], k, aa, bb)] = (n, v)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="521"); ap.add_argument("--NZ", type=int, default=16)
    ap.add_argument("--K", type=int, default=12); a = ap.parse_args(); p = PR[a.prime]; NZ = a.NZ; K = a.K
    H = HatRows(p, NZ, K); z = H.zeta; M0, M0i = M0_of(H.Rt, NZ, p)
    tau = [of(x, p) for x in tau_pred(K)]; G = G_of(H.Rt, M0i, tau, p, NZ, K); Lg = blog(G, p, NZ, K, K)
    print("=" * 100); print("ref_run2b: p = 2^%s-1, NZ = %d, K = %d; checks %s" % (a.prime, NZ, K, H.checks))
    for name, X in (("G", G), ("log G", Lg)):
        wt = weights(X, z, p); js = sorted(set(key[0] for key in wt))
        print("%s: weights j present: %s" % (name, js))
        for j in (0, 5):
            print("  weight j = %d part of %s (z^k s^%d K^a (x) K^b), per zhat-order k: component (a,b) [u-order n]: value" % (j, name, j))
            for k in range(1, NZ + 1):
                items = []; complete = True
                for aa in range(5):
                    for bb in range(5):
                        if (aa + bb + 2 + k) % 5: continue
                        n5 = j + 4 * k - 2 + W[aa] + W[bb]
                        n = n5 // 5
                        if n >= K: complete = False; continue
                        if n < 1: continue
                        e = wt.get((j, k, aa, bb))
                        items.append("(%d,%d)[u^%d]: %s" % (aa, bb, n, "0" if e is None else str(ratrec(e[1], p))))
                print("    k = %2d%s: %s" % (k, "" if complete else " (incomplete: some components need u-order >= K)", ";  ".join(items)))


if __name__ == "__main__":
    main()
