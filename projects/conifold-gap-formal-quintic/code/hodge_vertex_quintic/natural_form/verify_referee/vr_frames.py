# -*- coding: utf-8 -*-
"""vr_frames.py -- VERIFIER: frame dependence of G (conifold frame varpi_0 = f_0 + alpha g + beta f_2), on the PROJECT rows with the
verifier's own arithmetic.  For each frame: degrees of G_n, rule W (violations = component present below the bound; non-sharp =
component absent at the bound), chain relation at weight 0, the Fourier components of G_1, the H (x) H = (3,3) component of the
zhat-linear part at u^2, u^3, u^4, and the fit  (3,3)[u^2 zhat] = -127/15750 - (2/5) beta,  (2,4)[u^1 zhat] = (alpha - 8/45)/5.
Run: py -I -B vr_frames.py --prime 89 [--Nz 12 --K 9 --prec 50 --Nu 20]
"""
import os, sys, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from vr_proj import mmul, madd, mzero, mid, nz, bmul, tau_series, fourier, build, W, PSI

FRAMES = [("0", "0"), ("1/3", "-2/7"), ("0", "5"), ("8/45", "0"), ("0", "-127/6300"), ("8/45", "-127/6300"), ("-3/11", "2/13")]


def cdiv(a, b): return -((-a) // b)


def one(p, NZ, K, prec, Nu, al, be):
    R, H, Rt = build(p, NZ, K, prec, Nu, Fr(al), Fr(be)); z = R.zeta
    of = lambda x: Fr(x).numerator % p * pow(Fr(x).denominator, p - 2, p) % p
    M0 = [Rt.get((n, 0), mzero()) for n in range(NZ + 1)]; Mi = [mid()]
    for n in range(1, NZ + 1):
        acc = mzero()
        for m in range(1, n + 1): acc = madd(acc, mmul(M0[m], Mi[n - m], p), p)
        Mi.append([[(-x) % p for x in r] for r in acc])
    tq, tau = tau_series(K, p); tp = [[1] + [0] * (K - 1)]
    for n in range(1, NZ + 1):
        prev = tp[-1]; tp.append([sum(prev[j] * tau[k - j] for j in range(k + 1)) % p for k in range(K)])
    Mit = {(n, k): [[x * tp[n][k] % p for x in r] for r in Mi[n]] for n in range(NZ + 1) for k in range(K) if tp[n][k] and nz(Mi[n])}
    G = bmul(Rt, Mit, p, NZ, K)
    deg = [max([n for (n, k) in G if k == j] or [-1]) for j in range(K)]
    comp = {}
    for (k, n), A in G.items():
        for (m, mp), v in fourier(A, z, p).items(): comp[(n, k, m, mp)] = v
    viol = []; nonsharp = []; test = 0
    for aa in range(5):
        for bb in range(5):
            for k in range(1, NZ + 1):
                if (aa + bb + 2 + k) % 5: continue
                b = max(1, cdiv(4 * k - 2 + W[aa] + W[bb], 5))
                ns = [n for n in range(1, K) if (n, k, (-aa - 1) % 5, (-bb - 1) % 5) in comp]
                if ns and min(ns) < b: viol.append((aa, bb, k, min(ns), b))
                if b < K:
                    test += 1
                    if not ns or min(ns) > b: nonsharp.append((aa, bb, k, min(ns) if ns else None, b))
    X = {(n, k, (-m - 1) % 5, (-mp - 1) % 5): PSI[(4 + m) % 5] * PSI[(4 + mp) % 5] * v % p for (n, k, m, mp), v in comp.items()}
    cL = [0, 0]; cR = [0, 0]
    for n in range(1, K):
        for k in range(1, NZ):
            for b in range(5):
                if 5 * n - 4 * k + 2 - W[1] - W[b] != 0 or (1 + b + 2 + k) % 5: continue
                cL[0] += 1; cL[1] += X.get((n, k + 1, 0, b), 0) == X.get((n, k, 1, b), 0)
                cR[0] += 1; cR[1] += X.get((n, k + 1, b, 0), 0) == (-X.get((n, k, b, 1), 0)) % p
    g1 = {key[2:]: R.fmt(v) for key, v in sorted(comp.items()) if key[:2] == (1, 1)}
    hh = [R.fmt(comp.get((n, 1, 3, 3), 0)) for n in range(1, min(K, 5))]
    pred33 = comp.get((2, 1, 3, 3), 0) == of(Fr(-127, 15750) - Fr(2, 5) * Fr(be))
    pred24 = comp.get((1, 1, 2, 4), 0) == of((Fr(al) - Fr(8, 45)) / 5)
    return deg, viol, nonsharp, test, cL, cR, g1, hh, pred33, pred24, M0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--Nz", type=int, default=12)
    ap.add_argument("--K", type=int, default=9); ap.add_argument("--prec", type=int, default=50); ap.add_argument("--Nu", type=int, default=20)
    a = ap.parse_args(); p = {'89': (1 << 89) - 1, '61': (1 << 61) - 1}[a.prime]
    print("=" * 100); print("vr_frames: p = 2^%s-1, Nz = %d, K = %d; claimed degrees floor((5n+2)/4) = %s" % (a.prime, a.Nz, a.K, [(5 * n + 2) // 4 for n in range(a.K)]))
    M0ref = None
    for al, be in FRAMES:
        t0 = time.time(); deg, viol, nonsharp, test, cL, cR, g1, hh, p33, p24, M0 = one(p, a.Nz, a.K, a.prec, a.Nu, al, be)
        if M0ref is None: M0ref = M0
        print("frame (alpha, beta) = (%s, %s)  [%.0f s]" % (al, be, time.time() - t0))
        print("   degrees of G_n: %s;  symbol M_0 identical to that of the frame (0,0): %s" % (deg, M0 == M0ref))
        print("   rule W: violations %d, non-sharp %d %s of %d testable pairs;  chain L %d/%d, R %d/%d" % (
            len(viol), len(nonsharp), nonsharp[:6], test, cL[1], cL[0], cR[1], cR[0]))
        print("   G_1 Fourier: %s" % g1)
        print("   (3,3) of the zhat-linear part at u^1..u^4: %s;  u^2 value = -127/15750 - (2/5) beta: %s;  (2,4)[u zhat] = (alpha - 8/45)/5: %s" % (hh, p33, p24))


if __name__ == "__main__":
    main()
