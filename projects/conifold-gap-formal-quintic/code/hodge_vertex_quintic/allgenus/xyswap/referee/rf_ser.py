# -*- coding: utf-8 -*-
"""rf_ser.py -- referee's own exact series helpers (Fractions only; no worker module imported).

Power series: lists of Fractions, index = exponent.  Laurent series: pairs (val, list) = Sum list[j] t^(val+j).
Laurent polynomials in one variable: dicts {exponent: Fraction}.
"""
from fractions import Fraction as Fr


def gbinom(m, k):
    """binomial(m, k) for any integer m and k >= 0"""
    r = Fr(1)
    for i in range(k):
        r = r * (m - i) / (i + 1)
    return r


def pmul(a, b, n):
    """product of power series a, b mod t^(n+1)"""
    out = [Fr(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai:
            for j, bj in enumerate(b[: n + 1 - i]):
                if bj:
                    out[i + j] += ai * bj
    return out


def pinv(a, n):
    """1/a mod t^(n+1), a[0] != 0"""
    out = [Fr(0)] * (n + 1)
    out[0] = 1 / a[0]
    for k in range(1, n + 1):
        s = Fr(0)
        for i in range(1, min(k, len(a) - 1) + 1):
            if a[i]:
                s += a[i] * out[k - i]
        out[k] = -s / a[0]
    return out


def lmul(A, B, top):
    """product of Laurent series, keeping exponents <= top"""
    va, ca = A
    vb, cb = B
    v = va + vb
    n = top - v + 1
    if n <= 0:
        return (v, [])
    out = [Fr(0)] * n
    for i, a in enumerate(ca[:n]):
        if a:
            for j, b in enumerate(cb[: n - i]):
                if b:
                    out[i + j] += a * b
    return (v, out)


def coeff(A, e):
    v, c = A
    j = e - v
    return c[j] if 0 <= j < len(c) else Fr(0)


def acc_add(acc, A, cf=Fr(1)):
    """acc: dict exponent -> Fraction"""
    v, c = A
    for j, a in enumerate(c):
        if a:
            acc[v + j] = acc.get(v + j, Fr(0)) + a * cf


def acc_to_laurent(acc, top):
    ks = [e for e, c in acc.items() if c and e <= top]
    if not ks:
        return (0, [])
    v = min(ks)
    return (v, [acc.get(e, Fr(0)) for e in range(v, top + 1)])


def taylor_laurentpoly(L, a, N):
    """[t^k] L(a + t), k = 0..N, for a Laurent polynomial L = {e: c} at a rational point a != 0"""
    a = Fr(a)
    out = [Fr(0)] * (N + 1)
    for e, c in L.items():
        for k in range(N + 1):
            out[k] += c * gbinom(e, k) * a ** (e - k)
    return out


def taylor_log(a, N):
    """[t^k] log(a + t) for k >= 1 (k = 0 term, log a, dropped), a rational != 0"""
    a = Fr(a)
    return [Fr(0)] + [Fr((-1) ** (k + 1), k) / a ** k for k in range(1, N + 1)]


def bernoulli(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    for m in range(1, n + 1):
        s = Fr(0)
        for j in range(m):
            s += gbinom(m + 1, j) * B[j]
        B[m] = -s / (m + 1)
    return B


def lp_mul(a, b):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            out[e1 + e2] = out.get(e1 + e2, Fr(0)) + c1 * c2
    return {e: c for e, c in out.items() if c}


def lp_deriv(a, k=1):
    out = dict(a)
    for _ in range(k):
        out = {e - 1: c * e for e, c in out.items() if e != 0 and c}
    return out


def partitions(n, maxpart=None):
    """partitions of n as lists of parts (non-increasing)"""
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest
