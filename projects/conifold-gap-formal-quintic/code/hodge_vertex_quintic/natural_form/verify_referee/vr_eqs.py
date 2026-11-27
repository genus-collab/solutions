# -*- coding: utf-8 -*-
"""vr_eqs.py -- VERIFIER: the three differential equations of the referee (R-eq, M-eq, G-eq), the A-equation used in the uniqueness
proof, the spectrum of B_0, log G, and frame independence -- all recomputed on the PROJECT rows with the verifier's own arithmetic
(nothing imported from referee/).

  nabla := zhat (1-u) ((4/5) theta_zhat - theta_u),  Z = diag(zeta^i),  Psi_{a,i}(u) = Qhat_{a,0}^{(i)}(u),
  B(u) := Psi^-1 (E/5 + theta_u) Psi = (1/5) Sum_a Psi_{3-a,j} (e_a/5 + theta_u) Psi_{a,i},  E = diag(-1,-2,2,1,0)
  (R-eq) nabla Rt = [Z, Rt] + zhat (1-u) B(u) Rt
  (M-eq) (4/5) w^2 M_0' = [Z, M_0] + w B_0 M_0
  (G-eq) nabla G = [Z, G] + zhat ((1-u) B(u) G - tau G B_0)
  (A-eq) (4/5) w^2 A' = [Z + w B_0, A],  A = M_0 Z M_0^-1
Also: eigenvalues of B(u) on the Fourier vectors through u^1 (hand computation of the verifier: B f_m = t_m(u) f_m with
t = (0, 1/5, 2/5 + (2/5)u, -2/5 - (2/5)u, -1/5) + O(u^2)), and the hand-derived G_1 (see VERIFY_REFEREE notes in the final report):
x_0 = x_1 = -1/75, x_2 = x_4 = -8/225, x_3 = 0 for x_m = Ghat_1(m, 1-m).
Run: py -I -B vr_eqs.py --prime 61 --Nz 12 --K 9 --prec 50 --Nu 20 [--alpha 1/3 --beta -2/7]
"""
import os, sys, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from vr_proj import mmul, madd, mzero, mid, nz, bmul, tau_series, fourier, build, W

E = (-1, -2, 2, 1, 0)


def bsub(X, Y, p):
    out = dict(X)
    for key, B in Y.items(): out[key] = madd(out.get(key, mzero()), B, p, p - 1)
    return {k: v for k, v in out.items() if nz(v)}


def nabla(X, p, NZ, K, c45):
    out = {}
    for (k, n), A in X.items():
        c = (c45 * k - n) % p
        if not c or k + 1 > NZ: continue
        for nn, sg in ((n, 1), (n + 1, p - 1)):
            if nn < K: out[(k + 1, nn)] = madd(out.get((k + 1, nn), mzero()), A, p, c * sg % p)
    return {k: v for k, v in out.items() if nz(v)}


def commZ(X, z, p):
    return {key: [[(pow(z, j, p) - pow(z, i, p)) * A[j][i] % p for i in range(5)] for j in range(5)] for key, A in X.items()}


def blog(G, p, NZ, K):
    Y = {k: v for k, v in G.items() if k != (0, 0)}; out = {}; pw = {(0, 0): mid()}
    for m in range(1, K + NZ + 1):
        pw = bmul(pw, Y, p, NZ, K)
        if not pw: break
        c = pow(m, p - 2, p) * (1 if m % 2 else p - 1) % p
        for key, A in pw.items(): out[key] = madd(out.get(key, mzero()), A, p, c)
    return {k: v for k, v in out.items() if nz(v)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--Nz", type=int, default=12)
    ap.add_argument("--K", type=int, default=9); ap.add_argument("--prec", type=int, default=50); ap.add_argument("--Nu", type=int, default=20)
    ap.add_argument("--alpha", default="0"); ap.add_argument("--beta", default="0")
    a = ap.parse_args(); p = {'89': (1 << 89) - 1, '61': (1 << 61) - 1}[a.prime]; NZ = a.Nz; K = a.K
    t0 = time.time(); R, H, Rt = build(p, NZ, K, a.prec, a.Nu, Fr(a.alpha), Fr(a.beta)); z = R.zeta; i5 = pow(5, p - 2, p)
    of = lambda x: Fr(x).numerator % p * pow(Fr(x).denominator, p - 2, p) % p; c45 = of(Fr(4, 5))
    print("=" * 100); print("vr_eqs: p = 2^%s-1, Nz = %d, K = %d, frame (%s, %s); rows %.0f s" % (a.prime, NZ, K, a.alpha, a.beta, time.time() - t0))
    Psi = [[[H.Q[aa][i][0].coef(k) for k in range(K)] for i in range(5)] for aa in range(5)]          # Psi[a][i][k]
    Bu = [mzero() for _ in range(K)]
    for j in range(5):
        for i in range(5):
            for aa in range(5):
                x = Psi[(3 - aa) % 5][j]; y = [(E[aa] * i5 + k) % p * Psi[aa][i][k] % p for k in range(K)]
                for k1 in range(K):
                    for k2 in range(K - k1): Bu[k1 + k2][j][i] = (Bu[k1 + k2][j][i] + x[k1] * y[k2] * i5) % p
    antis = all(Bu[k][j][i] == (-Bu[k][i][j]) % p for k in range(K) for i in range(5) for j in range(5))
    # eigenvalues on Fourier vectors: B f_m = t_m f_m  <=>  Bhat(m', m) = t_m delta_{m+m'=0}/5
    eig = []
    for k in range(min(K, 3)):
        fo = fourier(Bu[k], z, p)
        assert all((m + mp) % 5 == 0 for (m, mp) in fo), "B(u) not Fourier-diagonal"
        eig.append([R.fmt(fo.get(((-m) % 5, m), 0) * 5 % p) for m in range(5)])
    print("(1) B(u) antisymmetric: %s; eigenvalues t_m on f_m (m = 0..4): u^0 %s, u^1 %s, u^2 %s" % (antis, eig[0], eig[1], eig[2] if len(eig) > 2 else None))
    zB = {(1, k): madd(Bu[k], Bu[k - 1] if k else mzero(), p, p - 1) for k in range(K)}; zB = {k: v for k, v in zB.items() if nz(v)}
    d = bsub(bsub(nabla(Rt, p, NZ, K, c45), commZ(Rt, z, p), p), bmul(zB, Rt, p, NZ, K), p)
    print("    (R-eq) residual zero through zhat^%d, u^%d: %s" % (NZ, K - 1, not d))
    M0 = [Rt.get((n, 0), mzero()) for n in range(NZ + 1)]; Mi = [mid()]
    for n in range(1, NZ + 1):
        acc = mzero()
        for m in range(1, n + 1): acc = madd(acc, mmul(M0[m], Mi[n - m], p), p)
        Mi.append([[(-x) % p for x in r] for r in acc])
    B0 = Bu[0]; okM = True
    for k in range(NZ):
        lhs = [[x * c45 * k % p for x in r] for r in M0[k]]
        rhs = madd([[(pow(z, j, p) - pow(z, i, p)) * M0[k + 1][j][i] % p for i in range(5)] for j in range(5)], mmul(B0, M0[k], p), p)
        okM &= lhs == rhs
    print("    (M-eq) through w^%d: %s" % (NZ, okM))
    tq, tau = tau_series(K, p)
    tp = [[1] + [0] * (K - 1)]
    for n in range(1, NZ + 1):
        prev = tp[-1]; tp.append([sum(prev[j] * tau[k - j] for j in range(k + 1)) % p for k in range(K)])
    Mit = {(n, k): [[x * tp[n][k] % p for x in r] for r in Mi[n]] for n in range(NZ + 1) for k in range(K) if tp[n][k] and nz(Mi[n])}
    G = bmul(Rt, Mit, p, NZ, K)
    print("    degrees of G_n: %s" % [max([n for (n, k) in G if k == j] or [-1]) for j in range(K)])
    tB0 = {(1, k): [[x * tau[k] % p for x in r] for r in B0] for k in range(K)}
    d = bsub(bsub(nabla(G, p, NZ, K, c45), commZ(G, z, p), p), bmul(zB, G, p, NZ, K), p)
    d = bsub(d, {k: [[(-x) % p for x in r] for r in v] for k, v in bmul(G, tB0, p, NZ, K).items()}, p)
    print("    (G-eq) residual zero: %s" % (not d))
    Zm = [[pow(z, i, p) if i == j else 0 for j in range(5)] for i in range(5)]; Aw = [mzero() for _ in range(NZ + 1)]
    for i in range(NZ + 1):
        for j in range(NZ + 1 - i): Aw[i + j] = madd(Aw[i + j], mmul(mmul(M0[i], Zm, p), Mi[j], p), p)
    okA = True
    for k in range(NZ):                                              # w^(k+1): (4/5) k A_k = [Z, A_(k+1)] + [B_0, A_k]
        lhs = [[x * c45 * k % p for x in r] for r in Aw[k]]
        rhs = madd(madd(mmul(Zm, Aw[k + 1], p), mmul(Aw[k + 1], Zm, p), p, p - 1), madd(mmul(B0, Aw[k], p), mmul(Aw[k], B0, p), p, p - 1), p)
        okA &= lhs == rhs
    print("    (A-eq) through w^%d: %s; A_k nonzero for k = %s; A_1 - B_0 diagonal: %s" % (
        NZ, okA, [k for k in range(NZ + 1) if nz(Aw[k])], all(i == j or Aw[1][i][j] == B0[i][j] for i in range(5) for j in range(5))))
    g1 = fourier(G[(1, 1)], z, p)
    print("(2) G_1 Fourier components: %s   (hand derivation: (0,1),(1,0): -1/75; (2,4),(4,2): -8/225; (3,3): 0)" % {k: R.fmt(v) for k, v in sorted(g1.items())})
    Lg = blog(G, p, NZ, K)
    print("(3) log G: degrees %s" % [max([n for (n, k) in Lg if k == j] or [-1]) for j in range(K)])
    for (n, k) in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (3, 4)):
        print("    [u^%d zhat^%d] log G: %s" % (n, k, {key: R.fmt(v) for key, v in sorted(fourier(Lg.get((k, n), mzero()), z, p).items())}))
    print("    zhat-linear parts of G and log G equal at every u-order: %s" % all(G.get((1, n)) == Lg.get((1, n)) for n in range(K)))
    print("    (3,3) of the zhat-linear part, u^1..: %s" % [R.fmt(fourier(G.get((1, n), mzero()), z, p).get((3, 3), 0)) for n in range(1, K)])
    print("time %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
