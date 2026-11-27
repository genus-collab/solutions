# -*- coding: utf-8 -*-
"""vr_proj_main.py -- driver of vr_proj.py (see its docstring).  Run: py -I -B vr_proj_main.py --prime 89 --Nz 10 --K 8 --prec 46 --Nu 14"""
import os, sys, time, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from vr_proj import mmul, madd, mzero, mid, nz, bmul, tau_series, fourier, build, W, PSI


def cdiv(a, b): return -((-a) // b)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--Nz", type=int, default=10)
    ap.add_argument("--K", type=int, default=8); ap.add_argument("--prec", type=int, default=46); ap.add_argument("--Nu", type=int, default=14)
    a = ap.parse_args(); p = {'89': (1 << 89) - 1, '61': (1 << 61) - 1}[a.prime]; NZ = a.Nz; K = a.K
    t0 = time.time(); R, H, Rt = build(p, NZ, K, a.prec, a.Nu); z = R.zeta
    print("=" * 100); print("vr_proj: p = 2^%s-1, Nz = %d, K = %d, prec_abs = %d, Nu = %d; project rows built in %.0f s; min u-precision %d" % (
        a.prime, NZ, K, a.prec, a.Nu, time.time() - t0, H.min_prec))
    of = lambda x: Fr(x).numerator % p * pow(Fr(x).denominator, p - 2, p) % p
    M0 = [Rt.get((n, 0), mzero()) for n in range(NZ + 1)]; assert M0[0] == mid()
    Mi = [mid()]
    for n in range(1, NZ + 1):
        acc = mzero()
        for m in range(1, n + 1): acc = madd(acc, mmul(M0[m], Mi[n - m], p), p)
        Mi.append([[(-x) % p for x in r] for r in acc])
    tq, tau = tau_series(K, p)
    tp = [[1] + [0] * (K - 1)]
    for n in range(1, NZ + 1):
        prev = tp[-1]; tp.append([sum(prev[j] * tau[k - j] for j in range(k + 1)) % p for k in range(K)])
    Minv_tau = {(n, k): [[x * tp[n][k] % p for x in r] for r in Mi[n]] for n in range(NZ + 1) for k in range(K) if tp[n][k] and nz(Mi[n])}
    G = bmul(Rt, Minv_tau, p, NZ, K)
    deg = [max([n for (n, k) in G if k == j] or [-1]) for j in range(K)]
    claim = [(5 * n + 2) // 4 for n in range(K)]
    print("(1) degrees of G_n: %s" % deg); print("    floor((5n+2)/4):  %s" % claim)
    print("    vanishing matrix coefficients above the claimed degree (within zhat <= %d): %d; nonzero ones: %d" % (
        NZ, sum(max(0, NZ - claim[n]) for n in range(1, K)), sum(1 for (n, k) in G if n > claim[k])))
    # (2) comparison with the referee's table
    d = json.load(open(os.path.join(os.path.dirname(HERE), "referee", "G_fourier_NZ16_K12.json")))
    ref = {}
    for key, v in d['G_fourier'].items():
        u_, z_, mm = key.split(); n = int(u_[1:]); k = int(z_[1:]); m, mp = (int(x) for x in mm.strip("()").split(","))
        if n < K and k <= NZ: ref[(n, k, m, mp)] = of(Fr(v))
    mine = {}
    for (k, n), A in G.items():
        for (m, mp), v in fourier(A, z, p).items(): mine[(n, k, m, mp)] = v
    print("(2) Fourier components of G (u < %d, zhat <= %d): verifier %d nonzero, referee table %d; identical as residues: %s" % (
        K, NZ, len(mine), len(ref), mine == ref))
    # (3) rule W as an iff
    miss = []; extra = []; tot = 0
    for n in range(1, K):
        for k in range(1, NZ + 1):
            for aa in range(5):
                for bb in range(5):
                    if (aa + bb + 2 + k) % 5: continue
                    j = 5 * n - 4 * k + 2 - W[aa] - W[bb]; key = (n, k, (-aa - 1) % 5, (-bb - 1) % 5); tot += 1
                    if j >= 0 and key not in mine: miss.append((n, k, aa, bb, j))
                    if j < 0 and key in mine: extra.append((n, k, aa, bb, j))
    print("(3) rule W as an IFF over %d Z/5-allowed (n, k, a, b): weight >= 0 but ZERO: %d %s; weight < 0 but NONZERO: %d" % (
        tot, len(miss), miss[:12], len(extra)))
    # (4) chain relation at weight 0 (flat components X^{ab} = psi_{3-a} psi_{3-b} Xhat(-a-1,-b-1))
    X = {(n, k, (-m - 1) % 5, (-mp - 1) % 5): PSI[(3 + m + 1) % 5] * PSI[(3 + mp + 1) % 5] * v % p for (n, k, m, mp), v in mine.items()}
    cL = [0, 0]; cR = [0, 0]
    for n in range(1, K):
        for k in range(1, NZ):
            for b in range(5):
                if 5 * n - 4 * k + 2 - W[1] - W[b] != 0 or (1 + b + 2 + k) % 5: continue
                cL[0] += 1; cL[1] += X.get((n, k + 1, 0, b), 0) == X.get((n, k, 1, b), 0) != 0
                cR[0] += 1; cR[1] += X.get((n, k + 1, b, 0), 0) == (-X.get((n, k, b, 1), 0)) % p != 0
    print("(4) chain relation, weight 0, k >= 1 (both sides nonzero required): L %d of %d, R %d of %d" % (cL[1], cL[0], cR[1], cR[0]))
    # (5) M_0 on the left at order u^1 with a free c_1:  [u^1] M_0(c zhat)^-1 Rt = A - c_1 B,  A = M_0^-1 Rt_1, B = M_0^-1 zhat M_0'
    A = [mzero() for _ in range(NZ + 1)]; B = [mzero() for _ in range(NZ + 1)]
    for i in range(NZ + 1):
        for j in range(NZ + 1 - i):
            A[i + j] = madd(A[i + j], mmul(Mi[i], Rt.get((j, 1), mzero()), p), p)
            if j: B[i + j] = madd(B[i + j], mmul(Mi[i], M0[j], p), p, j)
    for k in (NZ, NZ - 1, NZ - 2):
        fa = fourier(A[k], z, p); fb = fourier(B[k], z, p)
        cands = {key: fa.get(key, 0) * pow(fb[key], p - 2, p) % p for key in fb}
        bad = [key for key in fa if key not in fb]
        print("(5) left placement, u^1 zhat^%d: c_1 required per Fourier component: %d distinct values among %d components (tau_1 = -4/9 among them: %s); components with B = 0, A != 0: %d" % (
            k, len(set(cands.values())), len(cands), of(Fr(-4, 9)) in cands.values(), len(bad)))
    # right placement for comparison: [u^1] Rt M_0(c zhat)^-1 = A' - c_1 B',  A' = Rt_1 M_0^-1, B' = zhat M_0' M_0^-1
    A2 = [mzero() for _ in range(NZ + 1)]; B2 = [mzero() for _ in range(NZ + 1)]
    for i in range(NZ + 1):
        for j in range(NZ + 1 - i):
            A2[i + j] = madd(A2[i + j], mmul(Rt.get((i, 1), mzero()), Mi[j], p), p)
            if i: B2[i + j] = madd(B2[i + j], mmul(M0[i], Mi[j], p), p, i)
    for k in (NZ, NZ - 1):
        fa = fourier(A2[k], z, p); fb = fourier(B2[k], z, p)
        cands = set(fa.get(key, 0) * pow(fb[key], p - 2, p) % p for key in fb)
        print("    right placement, u^1 zhat^%d: distinct c_1: %d (= -4/9: %s)" % (k, len(cands), cands == {of(Fr(-4, 9))}))
    # (6) A = M_0 Z M_0^-1
    Zm = [[pow(z, i, p) if i == j else 0 for j in range(5)] for i in range(5)]
    Aw = [mzero() for _ in range(NZ + 1)]
    for i in range(NZ + 1):
        for j in range(NZ + 1 - i): Aw[i + j] = madd(Aw[i + j], mmul(mmul(M0[i], Zm, p), Mi[j], p), p)
    print("(6) A = M_0 Z M_0^-1: nonzero coefficients at w^k for k = %s (not a polynomial within the truncation)" % [k for k in range(NZ + 1) if nz(Aw[k])])
    print("    identity row f:", [R.fmt(H.Q[0][0][n].coef(0)) for n in range(min(NZ, 5) + 1)])
    print("time %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
