# -*- coding: utf-8 -*-
"""ref_run2.py -- REFEREE item 2 (structure of G): support of G in the flat-bivector (Fourier) components.

Dictionary (checked in ref_run4.py): the Fourier component Xhat(m, m') of a colour matrix is the coefficient of
Hh^a (x) Hh^b with a = -m-1, b = -m'-1 (mod 5), Hh^a = s^(1-c_a) H^a the hatted flat classes, c = (0,-1,3,2,1).
A monomial u^n zhat^k Hh^a (x) Hh^b equals z^k s^(5n - 4k + 2 - c_a - c_b) H^a (x) H^b in the un-hatted variables (z, s = u^(1/5)).
For every component (a, b) and every zhat-order k the script prints the minimal u-order n_min(k; a, b) that occurs in G and
compares with the candidate rules
   rule C:  5n - 4k + 2 - c_a - c_b >= 0     (regularity in s of the un-hatted flat bivector, weights c)
   rule W:  5n - 4k + 2 - w_a - w_b >= 0,  w = (0, 4, 3, 2, 1)   (c with c_1 = -1 replaced by +4)
Run: py -I ref_run2.py --prime 89 [--NZ 16] [--K 12]
"""
import os, sys, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, fmt
from ref_hat import HatRows, tau_pred, M0_of, G_of, degrees, fourier

PR = {'89': P89, '61': P61, '521': P521}
C = (0, -1, 3, 2, 1); W = (0, 4, 3, 2, 1)


def ceil_div(a, b): return -((-a) // b)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--NZ", type=int, default=16)
    ap.add_argument("--K", type=int, default=12); ap.add_argument("--alpha", default="0"); ap.add_argument("--beta", default="0")
    a = ap.parse_args(); p = PR[a.prime]; NZ = a.NZ; K = a.K
    H = HatRows(p, NZ, K, Fr(a.alpha), Fr(a.beta)); z = H.zeta; Rt = H.Rt; M0, M0i = M0_of(Rt, NZ, p)
    tau = [of(x, p) for x in tau_pred(K)]; G = G_of(Rt, M0i, tau, p, NZ, K); deg = degrees(G, K)
    print("=" * 100); print("ref_run2: p = 2^%s-1, NZ = %d, K = %d, frame (%s, %s); checks %s" % (a.prime, NZ, K, a.alpha, a.beta, H.checks))
    print("degrees of G_n:", deg); print("floor((5n+2)/4):", [(5 * n + 2) // 4 for n in range(K)])
    nmin = {}
    for (k, n), A in G.items():
        if (k, n) == (0, 0): continue
        for (m, mp), v in fourier(A, z, p).items():
            assert (m + mp - k) % 5 == 0
            ab = ((-m - 1) % 5, (-mp - 1) % 5)
            key = (ab, k)
            if key not in nmin or n < nmin[key]: nmin[key] = n
    viol = {'C': 0, 'W': 0}; slack = {'C': 0, 'W': 0}; total = 0
    print("component (a,b) = H^a (x) H^b;  per zhat-order k: n_min observed | rule C bound | rule W bound   ('-' = absent for all n < K)")
    for aa in range(5):
        for bb in range(5):
            line = []
            for k in range(1, NZ + 1):
                if (aa + bb + 2 + k) % 5: continue                    # Z/5 selection rule: m + m' = k, i.e. a + b = -k - 2 (mod 5)
                bC = max(1, ceil_div(4 * k - 2 + C[aa] + C[bb], 5)); bW = max(1, ceil_div(4 * k - 2 + W[aa] + W[bb], 5))
                nm = nmin.get(((aa, bb), k))
                line.append("k=%d: %s|%d|%d" % (k, "-" if nm is None else nm, bC, bW))
                for name, b in (('C', bC), ('W', bW)):
                    if b < K:
                        if nm is None or nm > b: slack[name] += 1
                        elif nm < b: viol[name] += 1
                total += 1
            print("   (%d,%d): %s" % (aa, bb, "  ".join(line)))
    print("rule C: violations (n_min below the bound) %d, non-sharp cases (n_min above the bound, bound < K) %d, of %d (component, k) pairs" % (viol['C'], slack['C'], total))
    print("rule W: violations %d, non-sharp cases %d, of %d" % (viol['W'], slack['W'], total))


if __name__ == "__main__":
    main()
