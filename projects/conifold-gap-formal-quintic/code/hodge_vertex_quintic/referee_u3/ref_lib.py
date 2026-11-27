# -*- coding: utf-8 -*-
"""ref_lib.py -- referee ELIMINATION: independent univariate / multivariate polynomial tools over F_p (pure Python).
Univariate: coefficient lists low -> high.  Multivariate: dict {exponent tuple: coefficient}.
None of the project helpers (interp_simple, sylvester_res, pgcd, det_fp, roots_fp) is used here."""
import random


# ------------------------------------------------------------------ univariate
def utrim(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0: a.pop()
    return a


def uadd(a, b, p):
    n = max(len(a), len(b)); return utrim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)], p)


def usub(a, b, p):
    n = max(len(a), len(b)); return utrim([(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(n)], p)


def uscale(a, c, p): return utrim([x * c for x in a], p)


def umul(a, b, p):
    if not a or not b: return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] += x * y
    return utrim(out, p)


def udivmod(a, b, p):
    a = utrim(a, p); b = utrim(b, p)
    if not b: raise ZeroDivisionError
    if len(a) < len(b): return [], a
    inv = pow(b[-1], p - 2, p); q = [0] * (len(a) - len(b) + 1); a = list(a); db = len(b) - 1
    for s in range(len(a) - len(b), -1, -1):
        f = a[s + db] * inv % p
        if f:
            q[s] = f
            for i in range(db + 1): a[s + i] = (a[s + i] - f * b[i]) % p
    return utrim(q, p), utrim(a[:db], p)


def umod(a, b, p): return udivmod(a, b, p)[1]


def udiv_exact(a, b, p):
    q, r = udivmod(a, b, p)
    if r: raise ValueError("inexact division")
    return q


def umonic(a, p):
    a = utrim(a, p)
    if not a: return a
    inv = pow(a[-1], p - 2, p); return [x * inv % p for x in a]


def ugcd(a, b, p):
    a = utrim(a, p); b = utrim(b, p)
    while b: a, b = b, umod(a, b, p)
    return umonic(a, p)


def ueval(a, x, p):
    r = 0
    for c in reversed(a): r = (r * x + c) % p
    return r


def uderiv(a, p): return utrim([k * a[k] for k in range(1, len(a))], p)


def upowmod(base, e, mod, p):
    res = [1]; b = umod(base, mod, p)
    while e:
        if e & 1: res = umod(umul(res, b, p), mod, p)
        e >>= 1
        if e: b = umod(umul(b, b, p), mod, p)
    return res


def usqfree(a, p):
    """squarefree part (characteristic is huge compared with the degrees)"""
    a = utrim(a, p)
    if len(a) <= 1: return a
    return umonic(udiv_exact(a, ugcd(a, uderiv(a, p), p), p), p)


def umult(a, r, p):
    """multiplicity of the root r in a"""
    m = 0; a = utrim(a, p)
    while a and ueval(a, r, p) == 0: a = udiv_exact(a, [(-r) % p, 1], p); m += 1
    return m


def uremove(a, f, p):
    """remove from a every factor it shares with f (repeatedly); returns the monic remainder"""
    a = umonic(a, p)
    while True:
        g = ugcd(a, f, p)
        if len(g) <= 1: return a
        a = umonic(udiv_exact(a, g, p), p)


def newton_interp(xs, ys, p):
    """Newton divided differences (different algorithm from the project's Lagrange/master polynomial)"""
    n = len(xs); c = [y % p for y in ys]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow(xs[i] - xs[i - j], p - 2, p) % p
    out = [c[n - 1]]
    for i in range(n - 2, -1, -1):
        new = [0] * (len(out) + 1)
        for k, v in enumerate(out):
            new[k + 1] = (new[k + 1] + v) % p; new[k] = (new[k] - xs[i] * v) % p
        new[0] = (new[0] + c[i]) % p; out = new
    return utrim(out, p)


def det_modp(M, p):
    M = [[x % p for x in r] for r in M]; n = len(M); d = 1
    for c in range(n):
        pr = None
        for r in range(c, n):
            if M[r][c]: pr = r; break
        if pr is None: return 0
        if pr != c: M[c], M[pr] = M[pr], M[c]; d = -d
        pv = M[c][c]; d = d * pv % p; inv = pow(pv, p - 2, p); rc = M[c]
        for r in range(c + 1, n):
            f = M[r][c]
            if f:
                f = f * inv % p; rr = M[r]; M[r] = [(rr[k] - f * rc[k]) % p for k in range(n)]
    return d % p


def ures_true(f, g, p):
    """resultant of f, g at their TRUE degrees by the Euclidean algorithm (Sylvester sign convention: f rows first)"""
    f = utrim(f, p); g = utrim(g, p)
    if not f or not g: return 0
    res = 1
    while True:
        m = len(f) - 1; n = len(g) - 1
        if n == 0: return res * pow(g[0], m, p) % p
        if m == 0: return res * pow(f[0], n, p) % p
        if m < n:
            f, g = g, f
            if (m * n) & 1: res = (-res) % p
            continue
        r = umod(f, g, p)
        if not r: return 0
        k = len(r) - 1                                     # Res(f, g) = (-1)^(m n) lc(g)^(m - k) Res(g, r) (-1)^(n k) ... handled by the swap
        res = res * pow(g[-1], m - k, p) % p
        if (m * n) & 1: res = (-res) % p
        f, g = g, r


def ures(f, g, df, dg, p):
    """resultant with FORMAL degrees df, dg (the determinant of the (df + dg) Sylvester matrix, f rows first)"""
    f = utrim(f, p); g = utrim(g, p); m = len(f) - 1; n = len(g) - 1
    if m > df or n > dg: raise ValueError("degree exceeds the formal degree")
    if m < df and n < dg: return 0
    if not f or not g: return 0
    if m == df: return pow(f[-1], dg - n, p) * ures_true(f, g, p) % p
    s = -1 if (dg * (df - m)) & 1 else 1
    return s * pow(g[-1], df - m, p) * ures_true(f, g, p) % p


def ddf(f, p, dmax=None):
    """distinct-degree factorisation of a squarefree monic f: list of (d, product of the irreducible factors of degree d)"""
    f = umonic(f, p); out = []; xq = [0, 1]; d = 0
    while len(f) - 1 >= 2 * (d + 1):
        d += 1; xq = upowmod(xq, p, f, p); g = ugcd(f, usub(xq, [0, 1], p), p)
        if len(g) > 1: out.append((d, g)); f = umonic(udiv_exact(f, g, p), p); xq = umod(xq, f, p)
        if dmax and d >= dmax: break
    if len(f) > 1: out.append((len(f) - 1, f))
    return out


def edf(g, d, p, rng):
    """split a product g of irreducible factors of degree d into the factors"""
    g = umonic(g, p)
    if len(g) - 1 == d: return [g]
    while True:
        r = [rng.randrange(p) for _ in range(len(g) - 1)]
        t = usub(upowmod(r, (p ** d - 1) // 2, g, p), [1], p); h = ugcd(g, t, p)
        if 1 < len(h) < len(g): return edf(h, d, p, rng) + edf(udiv_exact(g, h, p), d, p, rng)


def factor(f, p, rng):
    """irreducible factors (monic) of the squarefree part of f, sorted by degree"""
    out = []
    for d, g in ddf(usqfree(f, p), p): out += edf(g, d, p, rng)
    return sorted(out, key=lambda h: (len(h), h))


def uroots(f, p, rng):
    f = usqfree(f, p)
    if len(f) <= 1: return []
    g = ugcd(f, usub(upowmod([0, 1], p, f, p), [0, 1], p), p)
    return sorted((-h[0]) % p for h in (edf(g, 1, p, rng) if len(g) > 1 else []))


# ------------------------------------------------------------------ multivariate (dict)
def madd(a, b, p):
    r = dict(a)
    for k, v in b.items():
        x = (r.get(k, 0) + v) % p
        if x: r[k] = x
        else: r.pop(k, None)
    return r


def mscale(a, c, p):
    c %= p
    return {k: v * c % p for k, v in a.items()} if c else {}


def msub(a, b, p): return madd(a, mscale(b, p - 1, p), p)


def mmul(a, b, p):
    r = {}
    if len(a) > len(b): a, b = b, a
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb)); r[k] = (r.get(k, 0) + va * vb) % p
    return {k: v for k, v in r.items() if v}


def mconst(c, nv, p): return {(0,) * nv: c % p} if c % p else {}


def mvar(i, nv):
    e = [0] * nv; e[i] = 1; return {tuple(e): 1}


def mdeg(a, i=None):
    if not a: return -1
    return max(sum(k) for k in a) if i is None else max(k[i] for k in a)


def msupport(a): return sorted(set(i for k in a for i, e in enumerate(k) if e))


def mcoeffs(a, i):
    """coefficients with respect to variable i: list (low -> high) of polynomials (variable i set to exponent 0)"""
    d = mdeg(a, i); out = [dict() for _ in range(d + 1)]
    for k, v in a.items(): out[k[i]][k[:i] + (0,) + k[i + 1:]] = v
    return out


def msubs(a, vals, p):
    """substitute numbers for some variables: vals = {index: value}; the exponents of those variables become 0"""
    r = {}
    for k, v in a.items():
        kk = list(k)
        for i, x in vals.items():
            if k[i]: v = v * pow(x, k[i], p) % p; kk[i] = 0
        kk = tuple(kk); r[kk] = (r.get(kk, 0) + v) % p
    return {k: v for k, v in r.items() if v}


def muni(a, i, p):
    """a polynomial involving only variable i -> univariate list"""
    out = [0] * (mdeg(a, i) + 1)
    for k, v in a.items():
        assert sum(k) == k[i], "not univariate"
        out[k[i]] = v
    return utrim(out, p)


def mdet(M, p, nv):
    n = len(M)
    if n == 1: return M[0][0]
    tot = {}
    for c in range(n):
        if not M[0][c]: continue
        sub = mdet([r[:c] + r[c + 1:] for r in M[1:]], p, nv)
        if not sub: continue
        t = mmul(M[0][c], sub, p); tot = madd(tot, t, p) if c % 2 == 0 else msub(tot, t, p)
    return tot
