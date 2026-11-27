# -*- coding: utf-8 -*-
"""ref_cross.py -- REFEREE: the factorisation at THREE primes (2^521-1, 2^89-1, 2^61-1) and the cross-prime consistency of G.

For each prime: G = Rt . M_0(tau zhat)^(-1) with the predicted tau, its degrees, the tau fitted from the top zhat-order alone, and all
Fourier components Xhat(m, m') of all coefficients [u^n zhat^k] G.  The components are reconstructed as rationals at 2^521-1
(bound sqrt(p/2) ~ 10^78) and their reductions are compared with the residues at the two smaller primes (every prime uses its own
zeta: agreement shows that the Fourier components are rational numbers, the same at all primes).
Run: py -I ref_cross.py [--NZ 16] [--K 12] [--alpha 0] [--beta 0]
"""
import os, sys, time, argparse, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, ratrec
from ref_hat import HatRows, tau_pred, M0_of, G_of, degrees, fourier, fit_tau


def one(p, NZ, K, al, be):
    H = HatRows(p, NZ, K, al, be); M0, M0i = M0_of(H.Rt, NZ, p); tau = [of(x, p) for x in tau_pred(K)]
    G = G_of(H.Rt, M0i, tau, p, NZ, K); tfit, rep = fit_tau(H.Rt, M0, M0i, p, NZ, K)
    F = {}
    for (k, n), A in G.items():
        for key, v in fourier(A, H.zeta, p).items(): F[(n, k) + key] = v
    return dict(checks=H.checks, deg=degrees(G, K), fit_ok=[tfit[m] == tau[m] for m in range(K)], rep=rep, F=F,
                idrow=[H.Q[0][0][n][0] for n in range(NZ + 1)])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--NZ", type=int, default=16); ap.add_argument("--K", type=int, default=12)
    ap.add_argument("--alpha", default="0"); ap.add_argument("--beta", default="0"); ap.add_argument("--out", default="")
    a = ap.parse_args(); NZ = a.NZ; K = a.K; al = Fr(a.alpha); be = Fr(a.beta); res = {}
    print("=" * 100); print("ref_cross: NZ = %d, K = %d, frame (alpha, beta) = (%s, %s)" % (NZ, K, a.alpha, a.beta))
    for name, p in (('521', P521), ('89', P89), ('61', P61)):
        t = time.time(); r = one(p, NZ, K, al, be); res[name] = r
        print("p = 2^%s-1 (%.0f s): checks %s" % (name, time.time() - t, r['checks']))
        print("   degrees of G_n (n = 0..%d): %s" % (K - 1, r['deg']))
        print("   tau fitted from the top zhat-order == predicted: %s" % r['fit_ok'])
        print("   vanishing matrix coefficients above the degree, per u-order: %s (total %d)" %
              ([x['vanishing_orders'] for x in r['rep']], sum(x['vanishing_orders'] for x in r['rep'])))
    big = res['521']['F']; rat = {}; bad = 0
    for key, v in big.items():
        q = ratrec(v, P521)
        if q is None: bad += 1
        else: rat[key] = q
    print("Fourier components of G: %d nonzero at 2^521-1; rational reconstruction failures: %d" % (len(big), bad))
    hmax = max(max(abs(q.numerator), q.denominator) for q in rat.values())
    print("   largest numerator/denominator: %d decimal digits (bound of the reconstruction: 78 digits)" % len(str(hmax)))
    for name, p in (('89', P89), ('61', P61)):
        F = res[name]['F']; same_support = set(F) == set(big)
        agree = all(of(q, p) == F.get(key, 0) for key, q in rat.items())
        print("   reductions agree with the residues at 2^%s-1: %s; same support: %s (%d components)" % (name, agree, same_support, len(F)))
    idr = [ratrec(x, P521) for x in res['521']['idrow']]
    print("identity row f(zhat) of the symbol (rational, from 2^521-1):"); print("   " + ", ".join(str(x) for x in idr[:9]))
    print("   reductions agree at 2^89-1, 2^61-1: %s" % [all(of(q, p) == x for q, x in zip(idr, res[nm]['idrow'])) for nm, p in (('89', P89), ('61', P61))])
    if a.out:
        json.dump(dict(NZ=NZ, K=K, deg=res['521']['deg'], G_fourier={"u%d z%d (%d,%d)" % key: str(q) for key, q in sorted(rat.items())},
                       identity_row=[str(x) for x in idr]), open(os.path.join(HERE, a.out), "w"), indent=0)
        print("written", a.out)


if __name__ == "__main__":
    main()
