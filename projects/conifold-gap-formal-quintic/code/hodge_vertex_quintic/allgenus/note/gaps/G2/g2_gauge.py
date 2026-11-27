# -*- coding: utf-8 -*-
"""g2_gauge.py -- the graded regular gauge of Theorem U built along the proof written in GAPS_G2 (projection formula,
lifting lemma, scalar lemma), exact over Q, for the programme frames; checks of G8 (lem:gau-diag step 4), of the graded
normalisation (P_11(0,z) in Q[[z^5]]) and of the uniqueness corollary G9 through (G), each with a failing control.
Usage: py -I -B g2_gauge.py [K] [N]   (z-orders 0..K, u-precision N).  Stdlib only; imports g2_core only.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2_core import Ser, systems, prod, frobenius, log_obstruction

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def gauge(b, bt, K):
    """graded P_k (entries p[k][b] at (b+k, b)) with z P' = A P - P At, following steps (i)-(v) of the proof"""
    d = [Ser.const(1, bt[0].N)]
    for a in range(4):
        d.append(d[a] * bt[a] * b[a].inv())
    close = d[4] * bt[4] * b[4].inv() - d[0]
    assert close.iszero(), "Theta_U not consistent (products differ)"
    Dl = [d[a].d() * d[a].inv() for a in range(5)]
    trDl = Dl[0] + Dl[1] + Dl[2] + Dl[3] + Dl[4]
    B = prod(bt)
    dlogB = B.d() * B.inv()

    def W(a, m):  # calB_a^(m) = bt_a .. bt_(a+m-1)
        out = Ser.const(1, bt[0].N)
        for i in range(m):
            out = out * bt[(a + i) % 5]
        return out

    def scalar(k, Pi):
        """Lemma 2 (lem:gau-scal), graded: c_k = gam y^j, j = k mod 5, with pr(c' + Dl c) = -pr(Pi' + Dl Pi)"""
        j = k % 5
        Z = [Pi[bb].d() + Dl[(bb + k) % 5] * Pi[bb] for bb in range(5)]
        phi5 = Z[0] * W(0, j)
        for bb in range(1, 5):
            phi5 = phi5 + Z[bb] * W(bb, j)        # 5 phi_j = Sum_b Z_(b+k,b) calB_b^(j)
        assert phi5.v >= -1, "pr(Z) has a pole of order > 1"
        kap5 = dlogB.scale(-j) + trDl             # 5 kappa_j = j/u + holomorphic
        assert kap5.v >= -1 and kap5.resid() == j, "residue of kappa_j"
        n_max = min(phi5.N, kap5.N)
        g = []
        for n in range(0, n_max + 1):           # [u^(n-1)] of 5 c' + kap5 c + 5 phi = 0
            rhs = -phi5.coef(n - 1) - sum(kap5.coef(m) * g[n - 1 - m] for m in range(0, n))
            piv = 5 * n + j
            if piv == 0:
                assert rhs == 0
                g.append(Fr(1) if k == 0 else Fr(0))   # free constant: 1 at k = 0, 0 at k = 5q > 0
            else:
                g.append(rhs / piv)
        return Ser(0, g, n_max + 1)

    def lift(k, X):
        """Lemma 1 (lem:gau-lift): Pi_(k+1) holomorphic with [At, Pi_(k+1)] = X (entries of X at (b+k, b))"""
        dd = (-(k + 1)) % 5
        rho = {1: Ser(bt[0].N, [], bt[0].N)}
        for a in (1, 2, 3, 4):
            rho[(a + 1) % 5] = rho[a] + X[(a - k) % 5] * W(a, dd + 1).inv()
        closing = rho[0] + X[(0 - k) % 5] * W(0, dd + 1).inv()
        assert closing.iszero(), "pr(X) != 0"
        Pi = [rho[(bb + k + 1) % 5] * W((bb + k + 1) % 5, dd) for bb in range(5)]
        assert all(x.v >= 0 for x in Pi), "lifting lemma: pole in Pi"
        return Pi

    Pi = [Ser(bt[0].N, [], bt[0].N) for _ in range(5)]
    P = []
    for k in range(K + 1):
        j = k % 5
        gam = scalar(k, Pi)
        Ph = [Pi[bb] + gam * W(bb, j).inv() for bb in range(5)]
        assert all(x.v >= 0 for x in Ph), "Ph_k not holomorphic"
        P.append([d[(bb + k) % 5] * Ph[bb] for bb in range(5)])
        X = [Ph[bb].d() + Dl[(bb + k) % 5] * Ph[bb] for bb in range(5)]
        if k < K:
            Pi = lift(k, X)
    return P, d


def times_c(P, c):
    """P(z) c(z) for a scalar series c = {power: Fr}; graded storage requires c in Q[[z^5]]"""
    K = len(P) - 1
    out = []
    for k in range(K + 1):
        row = None
        for e, ce in c.items():
            if e <= k and ce != 0:
                assert e % 5 == 0
                t = [x.scale(ce) for x in P[k - e]]
                row = t if row is None else [row[i] + t[i] for i in range(5)]
        out.append(row if row is not None else [x.scale(0) for x in P[k]])
    return out


def inverse_series(c, K):
    """1/c(z) for c with c[0] != 0, truncated at z^K"""
    a0 = c[0]; out = {0: 1 / a0}
    for k in range(1, K + 1):
        out[k] = -sum(c.get(j, 0) * out[k - j] for j in range(1, k + 1)) / a0
    return out


def gauge_eq_ok(P, b, bt):
    """P_(k-1)' = A P_k - P_k At at every z-order k = 0..K, all entries, every known u-coefficient"""
    bad = 0
    for k in range(len(P)):
        for bb in range(5):
            # graded: A P_k - P_k At has entries (b+k-1, b) = b_(b+k-1) p_(k,b) - p_(k,b-1) bt_(b-1); P_(k-1)' has p_(k-1,b)'
            lhs = P[k - 1][bb].d() if k >= 1 else None
            rhs = b[(bb + k - 1) % 5] * P[k][bb] - P[k][(bb - 1) % 5] * bt[(bb - 1) % 5]
            diff = rhs if lhs is None else rhs - lhs
            if not diff.iszero():
                bad += 1
    return bad


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    print("g2_gauge.py: K = %d, N = %d (exact over Q)" % (K, N))
    fb = frobenius(10)
    obs = log_obstruction(fb)
    print("  log obstruction matrix (rows u^-3, u^-2, u^-1; columns f_0, f_1, f_2):", [[str(x) for x in r] for r in obs])
    check("obstruction matrix = [[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]], kernel (0,1,7/10)",
          obs == [[-2, 0, 0], [Fr(17, 5), 0, 0], [Fr(-7, 5), Fr(-7, 5), 2]]
          and all(sum(obs[i][j] * [0, 1, Fr(7, 10)][j] for j in range(3)) == 0 for i in range(3)))
    for (al, be) in [(0, 0), (Fr(1, 3), Fr(-2, 7))]:
        S = systems(N, al, be)
        b, bt = S["b"], S["bt"]
        print(" frame (alpha, beta) = (%s, %s); t_c = u + %s u^2 + ..., v = u + %s u^2 + ..." % (al, be, S["tc"].coef(2), S["v"].coef(2)))
        check("equal products prod b = prod bt = -1/(u(1-u)^5)", (prod(b) - prod(bt)).iszero()
              and (prod(b) * Ser(1, [-1], N) * Ser(0, [1, -1], N).powq(5) - 1).iszero())
        check("values at 0: b_a(0) = bt_a(0) (a != 1), (u b_1)(0) = (u bt_1)(0) = -1",
              all(b[a].coef(0) == bt[a].coef(0) for a in (0, 2, 3, 4)) and b[1].resid() == bt[1].resid() == -1)
        P, d = gauge(b, bt, K)
        print("  constructed graded P_0..P_%d; precisions %s" % (K, [min(x.N for x in Pk) for Pk in P]))
        check("gauge equation z P' = A P - P At, every entry and known coefficient", gauge_eq_ok(P, b, bt) == 0)
        check("P_k holomorphic (power series)", all(x.v >= 0 for Pk in P for x in Pk))
        check("P_0 diagonal with unit entries", all(P[0][a].coef(0) != 0 for a in range(5)))
        c11 = {5 * q: P[5 * q][1].coef(0) for q in range(K // 5 + 1)}
        print("  P_11(0,z) coefficients at z^0, z^5, ..:", [str(c11[e]) for e in sorted(c11)])
        Pc = times_c(P, inverse_series(c11, K))
        check("P_can: P_11(0,z) = 1", all(Pc[5 * q][1].coef(0) == (1 if q == 0 else 0) for q in range(K // 5 + 1)))
        check("P_can satisfies the gauge equation", gauge_eq_ok(Pc, b, bt) == 0)
        p0 = [Pc[0][a] for a in range(5)]
        check("G8: p_(0,a)(0) = 1 for every a (P_can)", all(x.coef(0) == 1 for x in p0), [str(x.coef(0)) for x in p0])
        pr0 = prod(p0)
        check("G8: prod_a p_(0,a) = 1 identically (to u^%d)" % (pr0.N - 1), (pr0 - 1).iszero())
        tau = S["tau"]; var0 = S["var0"]
        check("lem:gau-diag: p_(0,0) = tau^(-1/4)/varpi_0", (p0[0] - tau.powq(Fr(-1, 4)) * var0.inv()).iszero())
        c2 = times_c(Pc, {0: Fr(2)})
        check("control: P_can * 2 has prod_a p_(0,a) = 32 (so the normalisation matters)",
              (prod([c2[0][a] for a in range(5)]) - 32).iszero())
        btw = [bt[0].scale(2), bt[1].scale(Fr(1, 2)), bt[2], bt[3], bt[4]]
        Pw, _ = gauge(b, btw, 0)
        cw = Pw[0][1].coef(0)
        prw = prod([Pw[0][a].scale(1 / cw) for a in range(5)])
        check("control: with bt_0 -> 2 bt_0, bt_1 -> bt_1/2 (bt(0) != b(0)) the canonical product is 1/16, not 1",
              (prw - Fr(1, 16)).iszero(), "const = %s" % prw.coef(0))
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
