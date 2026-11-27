# -*- coding: utf-8 -*-
"""v2_gauge.py -- adversarial re-implementation of the regular gauge of UNIFORM_REDUCTION.md s.3 (Theorem U), own code and own
series class (Fractions), following the PROOF (reduction by D_0, commutant O[y], Lemma 1, Lemma 2), not the algorithm of
natural_form/theory/ur_gauge.py (particular solution + kernel + ODE for rho).  Holomorphy of every step is asserted, so a run of
this module is also a test of Lemma 1 and Lemma 2 order by order.

Systems: th Y_a = b_a Y_{a+1} (quintic, conifold frame (alpha, beta)), bt_a = v' bmod_a(v), v = u tau^(-5/4), bmod = (1,-1/v,1,-1,-1).
Output: p[k][b] = entry (b + k, b) of P_k, P = Sum z^k P_k, z P' = A P - P At.  Normalisation: p_{5,0}(0) = 0 (as in the note).
Only input from the project: the period series varpi_0, t_c of fq_conifold.Conifold (as rational coefficient lists).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)


class S(object):
    """Laurent series: coefficients c[i] of u^(v+i), absolute precision prec = v + len(c)"""
    def __init__(self, v, c):
        c = list(c); prec = v + len(c)
        while c and c[0] == 0: c.pop(0); v += 1
        self.v = v if c else prec; self.c = c; self.prec = prec

    def coef(self, k):
        assert k < self.prec, "precision"
        return self.c[k - self.v] if k >= self.v else Fr(0)

    def __add__(self, o):
        prec = min(self.prec, o.prec); v = min(self.v, o.v)
        return S(v, [self.coef(k) + o.coef(k) for k in range(v, prec)]) if v < prec else S(prec, [])

    def scale(self, a): return S(self.v, [a * x for x in self.c])
    def __neg__(self): return self.scale(-1)
    def __sub__(self, o): return self + (-o)

    def __mul__(self, o):
        prec = min(self.prec + o.v, o.prec + self.v); v = self.v + o.v
        if not self.c or not o.c or v >= prec: return S(prec, [])
        out = [Fr(0)] * (prec - v)
        for i, x in enumerate(self.c):
            if x == 0: continue
            for j, y in enumerate(o.c):
                if i + j >= len(out): break
                out[i + j] += x * y
        return S(v, out)

    def inv(self):
        assert self.c, "inverse of zero"
        n = len(self.c); a0 = self.c[0]; out = [1 / a0]
        for k in range(1, n):
            out.append(-sum(self.c[j] * out[k - j] for j in range(1, k + 1)) / a0)
        return S(-self.v, out)

    def d(self):
        return S(self.v - 1, [(self.v + i) * x for i, x in enumerate(self.c)])

    def shift(self, k): return S(self.v + k, self.c)

    def log1(self):
        """log of a series with v = 0, c[0] = 1"""
        assert self.v == 0 and self.c[0] == 1
        f = self.d() * self.inv()
        return S(0, [Fr(0)] + [f.coef(k) / (k + 1) for k in range(f.prec)])

    def exp0(self):
        """exp of a series with v >= 1 (or zero)"""
        assert self.v >= 1
        n = self.prec; out = [Fr(1)]
        for k in range(1, n):
            out.append(sum(j * self.coef(j) * out[k - j] for j in range(1, k + 1)) / k)
        return S(0, out)

    def powq(self, q):
        return self.log1().scale(q).exp0()


def const(a, N): return S(0, [Fr(a)] + [Fr(0)] * (N - 1))


def periods(N, alpha=0, beta=0):
    from fq_fp import FpRing, P2
    from fq_conifold import Conifold
    C = Conifold(FpRing(P2), Nu=N, alpha=Fr(alpha), beta=Fr(beta))
    return S(0, [Fr(C.varpi0.coef(n)) for n in range(N)]), S(0, [Fr(C.tc.coef(n)) for n in range(N)])


def systems(N, alpha=0, beta=0):
    var0, tc = periods(N, alpha, beta)
    one_m = S(0, [Fr(1), Fr(-1)] + [Fr(0)] * (N - 2))
    tcp = tc.d()
    u = S(1, [Fr(1)] + [Fr(0)] * (N - 1))
    b1 = -(u * var0 * var0 * one_m * one_m * one_m * tcp * tcp).inv()
    b3 = -(var0 * one_m.inv())
    bq = [tcp, b1, tcp, b3, b3]
    tau = S(0, [Fr(4, 5 * n + 4) for n in range(N)]).inv()
    v = tau.powq(Fr(-5, 4)).shift(1)
    vp = v.d()
    bt = [vp, -(vp * v.inv()), vp, -vp, -vp]
    return bq, bt, tau, v, var0, tc


def gauge(K, N, alpha=0, beta=0, verbose=None):
    bq, bt, tau, v, var0, tc = systems(N, alpha, beta)
    prod = lambda lst: lst[0] * lst[1] * lst[2] * lst[3] * lst[4]
    dd = prod(bq) - prod(bt)
    assert not dd.c, "prod b != prod bt"
    d = [const(1, N)]
    for a in range(4): d.append(d[a] * bt[a] * bq[a].inv())
    chk = d[4] * bt[4] * bq[4].inv() - d[0]
    assert not chk.c, "D_0 inconsistent"
    Dl = [d[a].d() * d[a].inv() for a in range(5)]
    trDl = Dl[0] + Dl[1] + Dl[2] + Dl[3] + Dl[4]
    B = prod(bt); dlogB = B.d() * B.inv()

    def W(a, j):
        out = None
        for i in range(j):
            out = bt[(a + i) % 5] if out is None else out * bt[(a + i) % 5]
        return out if out is not None else const(1, N)

    def lemma2(k, pi):
        """gamma with pr(Ph' + Dl Ph) = 0 for Ph = pi + gamma y^j, j = k mod 5"""
        j = k % 5
        Sx = None
        for b in range(5):
            t = W(b, j) * (pi[b].d() + Dl[(b + k) % 5] * pi[b])
            Sx = t if Sx is None else Sx + t
        assert Sx.v >= -1, "pr(Z) has a pole of order > 1"
        kap = dlogB.scale(-j) + trDl
        assert kap.v >= -1 and (kap.coef(-1) if kap.v == -1 else 0) == j, "residue of kappa"
        nmax = min(Sx.prec, kap.prec) - 1
        g = []
        for n in range(0, nmax + 1):
            rhs = -Sx.coef(n - 1) - sum(kap.coef(m) * g[n - 1 - m] for m in range(0, n))
            piv = 5 * n + j
            if piv == 0:
                assert rhs == 0, "resonance with non-zero right-hand side"
                g.append(Fr(1) if k == 0 else Fr(0))
            else:
                g.append(rhs / piv)
        return S(0, g)

    def lemma1(k, X):
        """holomorphic Pi_{k+1} (entries at E_{b+k+1,b}) with [At, Pi] = X (entries at E_{b+k,b})"""
        dg = (-(k + 1)) % 5
        rho = {1: S(N, [])}
        for a in (1, 2, 3, 4):
            rho[(a + 1) % 5] = rho[a] + X[(a - k) % 5] * W(a, dg + 1).inv()
        close = rho[0] + X[(0 - k) % 5] * W(0, dg + 1).inv()
        assert not close.c, "pr(X) != 0 (cycle does not close)"
        pi = [rho[(b + k + 1) % 5] * W((b + k + 1) % 5, dg) for b in range(5)]
        assert all(x.v >= 0 for x in pi), "Lemma 1: a pole in Pi_%d: %s" % (k + 1, [x.v for x in pi])
        return pi

    pi = [S(N, []) for _ in range(5)]
    P = []
    for k in range(K + 1):
        j = k % 5
        gam = lemma2(k, pi)
        ph = [pi[b] + gam * W(b, j).inv() for b in range(5)]
        assert all(x.v >= 0 for x in ph), "Ph_%d not holomorphic" % k
        P.append([d[(b + k) % 5] * ph[b] for b in range(5)])
        X = [ph[b].d() + Dl[(b + k) % 5] * ph[b] for b in range(5)]
        if verbose: verbose("   k = %d: valuations of p_{k,b}: %s, precision O(u^%d)" % (k, [x.v for x in P[-1]], min(x.prec for x in P[-1])))
        if k < K: pi = lemma1(k, X)
    # normalisation p_{0,0}(0) = 1 is automatic; impose p_{5m,0}(0) = 0 for m >= 1 by P -> P (1 - c z^(5m))
    for m in range(1, K // 5 + 1):
        c = P[5 * m][0].coef(0) / P[0][0].coef(0)
        for k in range(5 * m, K + 1):
            P[k] = [P[k][b] - P[k - 5 * m][b].scale(c) for b in range(5)]
    return dict(p=P, tau=tau, v=v, bq=bq, bt=bt, var0=var0, tc=tc)
