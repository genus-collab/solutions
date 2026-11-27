# -*- coding: utf-8 -*-
"""ref_run5.py -- REFEREE: (1) the differential equations behind the factorisation, (2) the 'chain' relation of the leading part.

(1) With nabla := zhat (1-u) ((4/5) theta_zhat - theta_u)  (= s z D at fixed z), Z = diag(zeta^i), Psi_{ai}(u) = Qh_{a,0}^{(i)}(u),
    B(u) := Psi^(-1) (E/5 + theta_u) Psi,  E = diag(-1, -2, 2, 1, 0),  the referee derives (ref_rows.py conventions)
        (R-eq)  nabla Rt = [Z, Rt] + zhat (1-u) B(u) Rt                       (hatted quantum differential equation)
        (M-eq)  (4/5) w^2 dM_0/dw = [Z, M_0] + w B(0) M_0                      (its u = 0 limit: the model)
        (G-eq)  nabla G = [Z, G] + zhat ( (1-u) B(u) G - tau(u) G B(0) )       (G = Rt M_0(tau zhat)^(-1))
    using  nabla f(tau(u) zhat) = (4/5) w^2 f'(w) at w = tau zhat, which holds iff (4/5) F + theta_u F = (4/5)/(1-u), F = 1/tau, i.e. iff
    F = Sum_n 4 u^n/(5n+4): THE PREDICTED tau IS THE UNIQUE RESCALING FOR WHICH M_0(tau zhat) SOLVES AN EQUATION OF THE SAME FORM.
    All three equations are checked here coefficient by coefficient (orders zhat^1..zhat^NZ, u^0..u^(K-1)).
(2) Weight-0 part (see ref_run2b.py): with X^{ab} the flat components (Hh^a (x) Hh^b) of the coefficient [u^n zhat^k] of G,
        chain L:  X^{0,b}[u^n zhat^(k+1)] = + X^{1,b}[u^n zhat^k],      chain R:  X^{a,0}[u^n zhat^(k+1)] = - X^{a,1}[u^n zhat^k]
    whenever the right-hand component has weight 0 (5n - 4k - 2 - w_b = 0): the leading part contains Hh^1 and Hh^0 only through
    Hh^1 + zhat Hh^0 = s^(-3) (u H + z) on the left leg and Hh^1 - zhat Hh^0 = s^(-3)(u H - z) on the right leg.
    Tested for G (k >= 1) and for G_flat = Psi(0)^(-1) Psi(u) G (k >= 0).
Run: py -I -B ref_run5.py --prime 89 [--NZ 16] [--K 12]
"""
import os, sys, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, inv, fmt, mid, mzero, mmul, madd, mscale, miszero, bmul, bclean, bsub
from ref_hat import HatRows, tau_pred, M0_of, G_of, fourier

PR = {'89': P89, '61': P61, '521': P521}
E = (-1, -2, 2, 1, 0); W = (0, 4, 3, 2, 1)


def nabla(X, p, NZ, K):
    """zhat (1-u) ((4/5) theta_zhat - theta_u) X, exact for zhat-orders <= NZ"""
    out = {}; c45 = of(Fr(4, 5), p)
    for (k, n), A in X.items():
        c = (c45 * k - n) % p
        if not c or k + 1 > NZ: continue
        for nn, sgn in ((n, 1), (n + 1, p - 1)):
            if nn >= K: continue
            B = out.get((k + 1, nn), mzero()); out[(k + 1, nn)] = madd(B, A, p, c * sgn % p)
    return bclean(out)


def commZ(X, zeta, p):
    return bclean({key: [[(pow(zeta, j, p) - pow(zeta, i, p)) * A[j][i] % p for i in range(5)] for j in range(5)] for key, A in X.items()})


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--NZ", type=int, default=16)
    ap.add_argument("--K", type=int, default=12); a = ap.parse_args(); p = PR[a.prime]; NZ = a.NZ; K = a.K
    H = HatRows(p, NZ, K); z = H.zeta; Rt = H.Rt; M0, M0i = M0_of(Rt, NZ, p); F = lambda x: fmt(x, p, z)
    tau = [of(x, p) for x in tau_pred(K)]; G = G_of(Rt, M0i, tau, p, NZ, K); i5 = inv(5, p)
    print("=" * 100); print("ref_run5: p = 2^%s-1, NZ = %d, K = %d; checks %s" % (a.prime, NZ, K, H.checks))
    # B(u)_{ji} = Sum_a (1/5) Psi_{3-a,j} (e_a/5 + theta_u) Psi_{a,i}
    Bu = [mzero() for _ in range(K)]
    for j in range(5):
        for i in range(5):
            acc = [0] * K
            for aa in range(5):
                x = H.Q[(3 - aa) % 5][j][0]; y = H.Q[aa][i][0]
                y2 = [(E[aa] * i5 + k) % p * y[k] % p for k in range(K)]
                for k1 in range(K):
                    if x[k1]:
                        for k2 in range(K - k1): acc[k1 + k2] = (acc[k1 + k2] + x[k1] * y2[k2]) % p
            for k in range(K): Bu[k][j][i] = acc[k] * i5 % p
    B0 = Bu[0]
    print("(1) B(u) antisymmetric and circulant at every u-order: %s; B(0) Fourier components (diagonal in the flat basis): %s" % (
        all(Bu[k][j][i] == (-Bu[k][i][j]) % p and Bu[k][(j + 1) % 5][(i + 1) % 5] == Bu[k][j][i] for k in range(K) for i in range(5) for j in range(5)),
        {key: F(v * 25 % p * i5 % p) for key, v in sorted(fourier(B0, z, p).items())}))
    zB = bclean({(1, k): madd(Bu[k], Bu[k - 1] if k else mzero(), p, p - 1) for k in range(K)})          # zhat (1-u) B(u)
    d = bsub(nabla(Rt, p, NZ, K), commZ(Rt, z, p), p); d = bsub(d, bmul(zB, Rt, p, NZ, K), p)
    d = {key: A for key, A in d.items() if key[0] >= 1}
    print("    (R-eq) residual vanishes at all orders zhat^1..zhat^%d, u^0..u^%d: %s" % (NZ, K - 1, not d))
    ztB0 = bclean({(1, k): mscale(B0, tau[k], p) for k in range(K)})                                   # tau(u) zhat B(0)
    d = bsub(nabla(G, p, NZ, K), commZ(G, z, p), p); d = bsub(d, bmul(zB, G, p, NZ, K), p)
    d = bsub(d, {key: mscale(A, p - 1, p) for key, A in bmul(G, ztB0, p, NZ, K).items()}, p)        # + tau zhat G B(0)
    d = {key: A for key, A in d.items() if key[0] >= 1}
    print("    (G-eq) residual vanishes at all orders: %s" % (not d))
    ok = True
    for k in range(NZ):                                              # (M-eq): (4/5) k M_k at w^(k+1) = [Z, M_(k+1)] + B0 M_k
        lhs = mscale(M0[k], of(Fr(4, 5), p) * k % p, p)
        rhs = madd([[(pow(z, j, p) - pow(z, i, p)) * M0[k + 1][j][i] % p for i in range(5)] for j in range(5)], mmul(B0, M0[k], p), p)
        ok &= (lhs == rhs)
    print("    (M-eq) holds for the symbol through w^%d: %s" % (NZ, ok))
    # ---------------------------------------------------------------- (2)
    psi = [H.Q[b][0][0][0] for b in range(5)]
    Pl = bclean({(0, k): [[sum(H.Q[(3 - aa) % 5][j][0][0] * H.Q[aa][i][0][k] for aa in range(5)) * i5 % p for i in range(5)] for j in range(5)] for k in range(K)})
    Gflat = bmul(Pl, G, p, NZ, K)
    for name, X in (("G (Teleman normalisation)", G), ("G_flat = Psi(0)^-1 Psi(u) G", Gflat)):
        comp = {}
        for (k, n), A in X.items():
            for (m, mp), v in fourier(A, z, p).items():
                aa = (-m - 1) % 5; bb = (-mp - 1) % 5
                comp[(n, k, aa, bb)] = psi[(3 - aa) % 5] * psi[(3 - bb) % 5] * v % p
        for j in (0, 5):
            res = {'L': [0, 0, 0], 'R': [0, 0, 0]}                    # tested, holds, of which k = 0
            for n in range(1, K):
                for k in range(0, NZ):
                    for b in range(5):
                        if 5 * n - 4 * k + 2 - W[1] - W[b] != j or (1 + b + 2 + k) % 5: continue
                        x1 = comp.get((n, k, 1, b), 0); x0 = comp.get((n, k + 1, 0, b), 0)
                        y1 = comp.get((n, k, b, 1), 0); y0 = comp.get((n, k + 1, b, 0), 0)
                        if k == 0 and (n, 0, 1, b) not in comp and name.startswith("G ("): continue      # G(zhat = 0) = 1: nothing to test
                        res['L'][0] += 1; res['L'][1] += (x0 == x1); res['L'][2] += (k == 0)
                        res['R'][0] += 1; res['R'][1] += (y0 == (-y1) % p); res['R'][2] += (k == 0)
            print("(2) %-28s weight %d: chain L holds in %d of %d cases, chain R in %d of %d (cases with k = 0: %d)" % (
                name, j, res['L'][1], res['L'][0], res['R'][1], res['R'][0], res['L'][2]))


if __name__ == "__main__":
    main()
