# -*- coding: utf-8 -*-
"""xy_local.py -- local data of the spectral curve x = p^4 + 4/p, y = p at the branch point a_0 = 1, and Laurent-list
arithmetic over a coefficient field (QZ or FP of xy_field.py).

Curve: x = p^4 + v/p with v = 4 (homogeneity restores v: p -> c p, v -> c^5 v), branch points p^5 = v/4 = 1, a_i = zeta^i.
Local coordinate t = p - 1 at a_0.  Everything below is over Q (Fractions); the field enters only through 1/(1 - zeta^i).
A Laurent series is a pair (val, list) meaning Sum_j list[j] t^(val + j); lists are truncated at a stated absolute order.
"""
from fractions import Fraction as Fr
from math import comb


# ---------------- rational power series (lists, index = power) ----------------
def ps_mul(a, b, n):
    out = [Fr(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai:
            for k, bk in enumerate(b[: n + 1 - i]):
                if bk: out[i + k] += ai * bk
    return out


def ps_inv(a, n):
    out = [Fr(0)] * (n + 1); out[0] = 1 / a[0]
    for k in range(1, n + 1):
        s = Fr(0)
        for i in range(1, min(k, len(a) - 1) + 1):
            if a[i]: s += a[i] * out[k - i]
        out[k] = -s / a[0]
    return out


def ps_sqrt1(a, n):
    """b with b^2 = a, a[0] = b[0] = 1"""
    assert a[0] == 1
    b = [Fr(1)] + [Fr(0)] * n
    for k in range(1, n + 1):
        s = a[k] if k < len(a) else Fr(0)
        for i in range(1, k): s -= b[i] * b[k - i]
        b[k] = s / 2
    return b


def ps_compose(f, g, n):
    """f(g(t)) mod t^(n+1), g[0] = 0 (Horner)"""
    assert g[0] == 0
    out = [Fr(0)] * (n + 1)
    for c in reversed(f[: n + 1]):
        out = ps_mul(out, g, n); out[0] += c
    return out


def ps_revert(f, n):
    """g with f(g(t)) = t, f[0] = 0, f[1] != 0 (coefficient by coefficient)"""
    g = [Fr(0), 1 / f[1]] + [Fr(0)] * (n - 1)
    for k in range(2, n + 1):
        h = ps_compose(f, g, k)
        g[k] = -h[k] / f[1]
    return g


def ps_deriv(a):
    return [k * a[k] for k in range(1, len(a))]


# ---------------- the curve at a_0 = 1 (rational) ----------------
class LocalQ:
    """c[k] = [t^k](x(1+t) - x(1));  dx[j] = [t^j] x'(1+t);  sig = sigma(t) (deck transformation), dsig = sigma'(t);
    kap[m] = (val, list): kappa_m(t) = (t^m - sigma^m)/((t - sigma) 2 x'(1+t)), the kernel K_0 = Sum_m e_{0,m}(p0) kappa_m(t)/dt;
    Bss = (val, list): B(p, sigma p)/dt^2 = sigma'/(t - sigma)^2;  sp[n] = sigma^n sigma' (power series)."""

    def __init__(self, N, mmax):
        self.N = N
        c = [Fr(0)] * (N + 4)
        for k in range(2, N + 4): c[k] = Fr(comb(4, k)) + 4 * (-1) ** k
        self.c = c
        self.dx = [Fr(k) * c[k] for k in range(1, N + 4)]           # x'(1+t) = Sum dx[j] t^j, dx[0] = 0
        u = [Fr(1)] + [c[k + 2] / c[2] for k in range(1, N + 1)]   # f/(c2 t^2)
        sq = ps_sqrt1(u, N)
        tau = [Fr(0)] + sq[:N]                                       # f = c2 tau^2
        inv = ps_revert(tau, N)
        sig = ps_compose(inv, [-a for a in tau], N)
        assert sig[0] == 0 and sig[1] == -1
        assert ps_compose(sig, sig, N)[:N + 1] == [Fr(0), Fr(1)] + [Fr(0)] * (N - 1), "sigma is not an involution"
        self.sig = sig
        self.dsig = ps_deriv(sig) + [Fr(0)]
        self.su = sig[1:] + [Fr(0)]                                 # sigma = t * su
        # 1/(2 x'(1+t)) = t^-1 * inv2dx
        dxu = self.dx[1:] + [Fr(0)]
        inv2dx = [a / 2 for a in ps_inv(dxu, N)]
        self.kap = {}
        for m in range(1, mmax + 1):
            q = [Fr(0)] * (N + 1)                                      # Sum_{j<m} t^(m-1-j) sigma^j
            sj = [Fr(1)] + [Fr(0)] * N
            for j in range(m):
                sh = m - 1 - j
                for i in range(N + 1 - sh): q[i + sh] += sj[i]
                sj = ps_mul(sj, sig, N)
            self.kap[m] = (-1, ps_mul(q, inv2dx, N))
        tms = [Fr(0)] + [(1 if k == 1 else 0) - sig[k] for k in range(1, N + 1)]   # t - sigma = t * (2 + ..)
        w = ps_inv(tms[1:] + [Fr(0)], N)
        self.Bss = (-2, ps_mul(ps_mul(w, w, N), self.dsig, N))
        self.sp = [self.dsig[:N + 1]]
        for n in range(1, N + 1): self.sp.append(ps_mul(self.sp[-1], sig, N))
        self.sup = {}                                                # su^(-k) cached

    def su_neg(self, k):
        if k not in self.sup:
            self.sup[k] = ps_inv(self.su, self.N) if k == 1 else ps_mul(self.su_neg(k - 1), self.su_neg(1), self.N)
        return self.sup[k]


# ---------------- Laurent lists over a field F ----------------
def lconv(F, ser, top):
    """rational Laurent (val, list) -> F-valued, truncated at absolute order top"""
    v, c = ser
    return (v, [F.frac(x) for x in c[: max(top - v + 1, 0)]])


def lmul(A, B, top, zero):
    va, ca = A; vb, cb = B
    v = va + vb; n = top - v + 1
    if n <= 0: return (v, [])
    out = [zero] * n
    for i, a in enumerate(ca[:n]):
        if a.is_zero(): continue
        for j, b in enumerate(cb[: n - i]):
            out[i + j] = out[i + j] + a * b
    return (v, out)


def ladd_into(acc, A, coef=None):
    """acc: dict exponent -> F (sparse accumulator)"""
    v, c = A
    for j, a in enumerate(c):
        if coef is not None: a = a * coef
        e = v + j
        acc[e] = acc[e] + a if e in acc else a


def acc_to_l(acc, zero, top):
    if not acc: return (0, [])
    v = min(acc)
    return (v, [acc.get(e, zero) for e in range(v, top + 1)])
