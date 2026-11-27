# -*- coding: utf-8 -*-
"""ref_run1.py -- REFEREE item 1/2: independent recomputation of the Teleman form and the test of
      Rt(zhat; u) = G(u, zhat) M_0(tau(u) zhat),   tau = 1 / Sum_n 4 u^n/(5n+4).
Run: py -I ref_run1.py --prime 89|61|521 [--NZ 6] [--K 6] [--alpha 0] [--beta 0]
Prints: internal checks, identity row of the symbol, fitted tau (from the top zhat-order only) vs prediction, degrees of G_n,
number of vanishing matrix coefficients, Fourier components of the low orders.
"""
import os, sys, time, argparse, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, ratrec, fmt, mid, mzero, miszero
from ref_hat import HatRows, tau_pred, M0_of, G_of, degrees, fourier, fit_tau

PR = {'89': P89, '61': P61, '521': P521}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--NZ", type=int, default=6)
    ap.add_argument("--K", type=int, default=6); ap.add_argument("--alpha", default="0"); ap.add_argument("--beta", default="0")
    ap.add_argument("--dump", default=""); a = ap.parse_args(); p = PR[a.prime]; NZ = a.NZ; K = a.K
    t0 = time.time(); H = HatRows(p, NZ, K, Fr(a.alpha), Fr(a.beta)); z = H.zeta
    print("=" * 100); print("ref_run1: p = 2^%s-1, NZ = %d, K = %d, frame (alpha, beta) = (%s, %s)   [%.1f s]" % (a.prime, NZ, K, a.alpha, a.beta, time.time() - t0))
    print("internal checks:", H.checks)
    print("WKB: degrees of ph_n(u) (colour 0):", H.wkb_info[0]['deg'], " free coefficients (n, power of u):", H.wkb_info[0]['free'])
    F = lambda x: fmt(x, p, z)
    Rt = H.Rt; M0, M0i = M0_of(Rt, NZ, p)
    print("identity row of the symbol  f = R^_0^(0) = Qh[0][0][n](u = 0):", [F(H.Q[0][0][n][0]) for n in range(min(NZ, 6) + 1)])
    print("gamma (g = u + gamma u^2 + ..):", F(H.g[2]), "  f0:", [F(x) for x in H.f[0][:5]], "  t_c:", [F(x) for x in H.tc[:4]])
    # ---- prediction
    tp = tau_pred(K); tau = [of(x, p) for x in tp]
    G = G_of(Rt, M0i, tau, p, NZ, K); deg = degrees(G, K)
    print("PREDICTED tau:", [str(x) for x in tp])
    print("zhat-degree of [u^n] G, n = 0..%d (truncation NZ = %d): %s" % (K - 1, NZ, deg))
    # ---- fit from the top order
    tfit, rep = fit_tau(Rt, M0, M0i, p, NZ, K)
    print("FITTED tau (top zhat-order only):", [F(x) for x in tfit])
    print("fitted == predicted:", [tfit[m] == tau[m] for m in range(K)])
    for r in rep: print("   ", r)
    # ---- low-order Fourier components
    for n in range(1, min(K, 4)):
        for k in range(1, (deg[n] or 0) + 1):
            if (k, n) in G: print("   [u^%d zhat^%d] G Fourier (m,m'): %s" % (n, k, {key: F(v) for key, v in sorted(fourier(G[(k, n)], z, p).items())}))
    if a.dump:
        json.dump(dict(prime=a.prime, NZ=NZ, K=K, checks=H.checks, deg=deg, tau_fit=[str(x) for x in tfit],
                       G={"%d,%d" % key: A for key, A in G.items()}), open(os.path.join(HERE, a.dump), "w"))
    print("time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
