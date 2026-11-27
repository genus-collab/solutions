"""Referee's INDEPENDENT topological-recursion engine (no import of worker / tr modules).

Curve x = p^4/4 + v/p, y = p at v = 1, B = dp1 dp2/(p1-p2)^2, mod a prime P = 1 mod 5.
Method (different from the worker's e-basis engine): local coordinate zeta_b at each branch point c_b with
x - x(c_b) = alpha_b zeta^2 (alpha_b = x''(c_b)/2, no square roots), sigma: zeta -> -zeta; omega_{g,n} stored as a
SYMMETRIC tensor in the odd basis dxi_{b,k}(p) = Res_{q->c_b} B(p,q) zeta_b(q)^-(2k+1)/(2k+1) = zeta^-(2k+2) dzeta + reg.
Kernel K(p0,q) = (1/2) Int_{sigma q}^{q} B(p0,.)/((y(q)-y(sigma q)) dx(q)) = Sum_k dxi_{b,k}(p0) kappa_k(zeta)/dzeta,
kappa_k = zeta^(2k)/(2 alpha (y(zeta)-y(-zeta))).  (Orientation q-bar -> q: the OPPOSITE of arXiv:math-ph/0702045 (4-3).)
omega(sigma q) is the pull-back sigma^*: f(zeta) dzeta -> -f(-zeta) dzeta.
B(q, p_l) in a J-slot is projected on its odd part Sum_k (2k+1) zeta^(2k) dzeta dxi_{b,k}(p_l) (legitimate because the
result is odd in p_l; the symmetry of every computed tensor is CHECKED, rf_tr.check_sym).
Usage: build(P, KM, targets) -> dict (g,n) -> {sorted index tuple: value}; index = b*(KM+1)+k.
"""
import sys

def inv(a, P):
    return pow(a % P, P - 2, P)

def zeta5(P):
    assert (P - 1) % 5 == 0
    for g in range(2, 1000):
        z = pow(g, (P - 1) // 5, P)
        if z != 1:
            return z
    raise ValueError

# ---------- truncated power series (lists, length N) ----------
def smul(a, b, N, P):
    r = [0] * N
    for i, ai in enumerate(a[:N]):
        if ai:
            for j in range(min(len(b), N - i)):
                r[i + j] = (r[i + j] + ai * b[j]) % P
    return r

def sinv(a, N, P):  # a[0] != 0
    r = [0] * N
    i0 = inv(a[0], P)
    r[0] = i0
    for n in range(1, N):
        s = 0
        for k in range(1, min(n, len(a) - 1) + 1):
            s += a[k] * r[n - k]
        r[n] = (-s * i0) % P
    return r

def ssqrt1(a, N, P):  # a[0] == 1
    r = [0] * N
    r[0] = 1
    h = inv(2, P)
    for n in range(1, N):
        s = a[n] if n < len(a) else 0
        for k in range(1, n):
            s -= r[k] * r[n - k]
        r[n] = (s * h) % P
    return r

def spow(a, m, N, P):
    r = [1] + [0] * (N - 1)
    for _ in range(m):
        r = smul(r, a, N, P)
    return r

def sderiv(a, P):
    return [(i * a[i]) % P for i in range(1, len(a))] + [0]

def ssub_neg(a, P):  # a(-z)
    return [(-c) % P if i % 2 else c for i, c in enumerate(a)]

class Local:
    """Data at branch point c (v=1): t(zeta) with p = c + t, alpha."""
    def __init__(self, c, N, P):
        self.c, self.N, self.P = c, N, P
        M = N + 4
        from math import comb
        xs = [0] * (M + 3)
        for j in range(M + 3):
            a = comb(4, j) * pow(c, 4 - j, P) * inv(4, P) if j <= 4 else 0
            b = (-1) ** j * inv(pow(c, j + 1, P), P)
            xs[j] = (a + b) % P
        assert xs[1] == 0, "c must be a branch point"
        self.alpha = xs[2]
        ia = inv(xs[2], P)
        g = [(xs[j + 2] * ia) % P for j in range(M + 1)]   # (x-x_b)/(alpha t^2)
        u = ssqrt1(g, M, P)                                  # zeta = t u(t)
        # invert zeta = t u(t): t = zeta / u(t), fixed point iteration
        t = [0, 1] + [0] * (M - 2)
        for _ in range(M + 1):
            ut = self.compose(u, t, M)
            q = sinv(ut, M, P)
            t = [0] + q[:M - 1]
        self.t = t[:N + 2]
        # check: x(c+t(zeta)) - x(c) == alpha zeta^2
        tt = self.t
        chk = [0] * (N + 2)
        tp = [1] + [0] * (N + 1)
        for j in range(1, N + 2):
            tp = smul(tp, tt, N + 2, P)
            for i in range(N + 2):
                chk[i] = (chk[i] + xs[j] * tp[i]) % P
        assert chk[2] == self.alpha and all(chk[i] == 0 for i in range(N + 2) if i != 2), "local coordinate"

    def compose(self, f, t, M):  # f(t(zeta)), t[0]=0
        P = self.P
        r = [0] * M
        tp = [1] + [0] * (M - 1)
        for j in range(M):
            if j > 0:
                tp = smul(tp, t, M, P)
            if j < len(f) and f[j]:
                for i in range(M):
                    r[i] = (r[i] + f[j] * tp[i]) % P
        return r
