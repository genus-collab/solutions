# -*- coding: utf-8 -*-
"""tr_core.py -- truncated series over a field (tr_field) and the local data of genus-0 spectral curves with y = p.

Curves (global coordinate p on P^1, B = dp1 dp2/(p1-p2)^2, y = p):
  ModelCurve(F, N, kappa):  x(p) = kappa W(p),  W = p^4/4 + 1/p  (the W-form of the LG model at u = 1, s = 1);
                  kappa = 1: x = W;  kappa = 4: x = p^4 + 4/p (lambda = p^4 + v/p at v = 4);  kappa = -1: x = -W.
                  x' = kappa (p^5 - 1)/p^2; branch points c_j = zeta^(-j) (colour j), x(c_j) = (5/4) kappa zeta^j,
                  x''(c_j) = 5 kappa c_j^2;  Z/5 symmetric: x(zeta^-1 p) = zeta x(p), y(zeta^-1 p) = zeta^-1 y(p).
  AiryCurve(F):   x(p) = -p^2/2, one branch point 0 (regression: Witten-Kontsevich numbers).
Series: list c[0..N] = coefficients of t^0..t^N (field elements); Laurent series: (v, c) = Sum_k c[k] t^(v+k).
Local data at a branch point a (t = p - a):  l[k] = [t^k] x(a + t); tau(t) = t sqrt(1 + Sum_{k>=1} (l_{k+2}/l_2) t^k), so
x - x(a) = l_2 tau^2; trev = t(tau); sigma(t) = trev(-tau(t)) (the deck transformation); dx[k] = [t^k] x'(a + t).
New file of the allgenus/tr task (replaces the unreviewed version of the first attempt).
"""
from fractions import Fraction as Fr


def smul(a, b, N, F):
    out = [F.zero] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x.iszero(): continue
        for j, y in enumerate(b[:N + 1 - i]):
            out[i + j] = out[i + j] + x * y
    return out


def sinv(a, N, F):
    a0i = F.one / a[0]; out = [a0i] + [F.zero] * N
    for n in range(1, N + 1):
        s = F.zero
        for k in range(1, min(n, len(a) - 1) + 1):
            s = s + a[k] * out[n - k]
        out[n] = -(s * a0i)
    return out


def ssqrt1(a, N, F):
    assert (a[0] - F.one).iszero()
    b = [F.one] + [F.zero] * N; half = F.of(Fr(1, 2))
    for n in range(1, N + 1):
        s = a[n] if n < len(a) else F.zero
        for k in range(1, n): s = s - b[k] * b[n - k]
        b[n] = s * half
    return b


def scomp(f, g, N, F):
    """f(g(t)), g[0] = 0 (Horner)"""
    assert g[0].iszero()
    out = [F.zero] * (N + 1)
    for c in reversed(f[:N + 1]):
        out = smul(out, g, N, F); out[0] = out[0] + c
    return out


def srevert(f, N, F):
    """g with f(g(t)) = t;  f[0] = 0, f[1] != 0"""
    assert f[0].iszero() and not f[1].iszero()
    f1i = F.one / f[1]; g = [F.zero, f1i] + [F.zero] * (N - 1)
    for n in range(2, N + 1):
        c = scomp(f, g, n, F)[n]
        g[n] = -(c * f1i)
    return g


def sder(a, F):
    return [a[k] * k for k in range(1, len(a))] or [F.zero]


def lmul(A, B, T, F):
    """Laurent product truncated at absolute degree <= T"""
    va, ca = A; vb, cb = B; v = va + vb; n = T - v + 1
    if n <= 0: return (v, [])
    out = [F.zero] * n
    for i, x in enumerate(ca[:n]):
        if x.iszero(): continue
        for j, y in enumerate(cb[:n - i]):
            out[i + j] = out[i + j] + x * y
    return (v, out)


def lcoef(A, n, F):
    v, c = A; k = n - v
    return c[k] if 0 <= k < len(c) else F.zero


def binom(e, m):
    """binomial(e, m) for integer e (possibly negative), m >= 0, as a Fraction"""
    r = Fr(1)
    for i in range(m): r = r * (e - i) / (i + 1)
    return r


class Local(object):
    """local data of x at a branch point a (field element), series order N"""
    def __init__(self, F, lcoefs, a, N):
        self.F = F; self.a = a; self.N = N
        l = lcoefs                                                     # l[k] = [t^k] x(a+t), k = 0..N+2
        assert l[1].iszero() and not l[2].iszero(), "not a simple branch point"
        self.l = l; l2i = F.one / l[2]
        h = [l[k + 2] * l2i for k in range(N + 1)]
        r = ssqrt1(h, N, F)
        self.tau = [F.zero] + r[:N]
        self.trev = srevert(self.tau, N, F)
        self.sigma = scomp(self.trev, [-c for c in self.tau], N, F)
        assert self.sigma[0].iszero() and (self.sigma[1] + F.one).iszero()
        self.dx = [l[k + 1] * (k + 1) for k in range(N + 1)]
        self.dsig = sder(self.sigma, F)[:N + 1]

    def check(self, M=None):
        """x(a + sigma(t)) = x(a + t) and sigma(sigma(t)) = t through t^M"""
        F = self.F; M = (self.N - 2) if M is None else M
        f = [F.zero] + self.l[1:M + 1]
        lhs = scomp(f, self.sigma, M, F)
        ok1 = all((lhs[k] - f[k]).iszero() for k in range(M + 1))
        ss = scomp(self.sigma, self.sigma, M, F)
        ok2 = all((ss[k] - (F.one if k == 1 else F.zero)).iszero() for k in range(M + 1))
        return ok1 and ok2


class ModelCurve(object):
    """x = kappa (p^4/4 + 1/p), y = p, B standard; branch points c_j = zeta^-j"""
    nbp = 5
    symmetric = True

    def __init__(self, F, N, kappa=1):
        self.F = F; self.N = N; self.kappa = Fr(kappa); self.kap = F.of(self.kappa)
        self.name = "x = %s (p^4/4 + 1/p), y = p" % self.kappa
        self.c = [F.z(-j) for j in range(5)]
        self._loc = {}

    def x_taylor(self, a, N):
        """[t^k] x(a + t) = kappa ( binom(4,k) a^(4-k)/4 + (-1)^k a^(-k-1) )"""
        F = self.F; ai = F.one / a; out = []
        for k in range(N + 1):
            v = F.zero
            if k <= 4: v = v + F.of(binom(4, k) / 4) * self._pw(a, 4 - k)
            v = v + F.of((-1) ** k) * self._pw(ai, k + 1)
            out.append(v * self.kap)
        return out

    def _pw(self, a, e):
        r = self.F.one
        for _ in range(e): r = r * a
        return r

    def local(self, j):
        L = self._loc.get(j)
        if L is None:
            L = Local(self.F, self.x_taylor(self.c[j], self.N + 2), self.c[j], self.N); self._loc[j] = L
        return L

    def ydx_taylor(self, j, N):
        """[t^k] of y x' at c_j + t (y x' = kappa (p^4 - 1/p)), k = 0..N"""
        F = self.F; a = self.c[j]; ai = F.one / a; out = []
        for k in range(N + 1):
            v = F.zero
            if k <= 4: v = v + F.of(binom(4, k)) * self._pw(a, 4 - k)
            v = v - F.of((-1) ** k) * self._pw(ai, k + 1)
            out.append(v * self.kap)
        return out


class AiryCurve(object):
    """x = -p^2/2, y = p: TR gives omega_{g,n} = Sum <tau_d1..tau_dn>_g prod (2d_j+1)!! dp_j/p_j^(2d_j+2)"""
    name = "Airy x = -p^2/2, y = p"
    nbp = 1
    symmetric = False

    def __init__(self, F, N):
        self.F = F; self.N = N; self.c = [F.zero]; self._loc = {}

    def x_taylor(self, a, N):
        F = self.F
        return [F.zero, F.zero, F.of(Fr(-1, 2))] + [F.zero] * (N - 2)

    def local(self, j):
        L = self._loc.get(j)
        if L is None:
            L = Local(self.F, self.x_taylor(self.c[j], self.N + 2), self.c[j], self.N); self._loc[j] = L
        return L

    def ydx_taylor(self, j, N):
        F = self.F
        return [F.zero, F.zero, F.of(-1)] + [F.zero] * (N - 2)  # y x' = -p^2


if __name__ == "__main__":
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from tr_field import QzField, make_fp, P89
    for F in (QzField(), make_fp(P89)):
        for kap in (1, -1, 4):
            C = ModelCurve(F, 16, kap)
            for j in range(5):
                L = C.local(j); assert L.check(12)
                # x''(c_j) = 5 kappa c_j^2, x(c_j) = (5/4) kappa zeta^j
                assert (L.l[2] * 2 - C.c[j] * C.c[j] * 5 * kap).iszero() and (L.l[0] - F.z(j) * F.of(Fr(5, 4) * kap)).iszero()
        A = AiryCurve(F, 10); assert A.local(0).check(8)
    L0 = ModelCurve(QzField(), 10).local(0)
    print("tr_core ok; sigma at c_0 = 1:", [L0.sigma[k].rational() for k in range(1, 8)])

