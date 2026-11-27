"""Critic helper: exact Laurent-log polynomials in p with a symbolic parameter v.

A term is keyed by (j, l, m): coefficient * p^j * (log p)^l * v^m, coefficient a Fraction.
Used by cr_dictionary.py (TR legs of the model curve) and cr_conv.py.  Pure stdlib, exact over Q.
"""
from fractions import Fraction as Fr


def clean(A):
    return {k: c for k, c in A.items() if c != 0}


def add(A, B, s=1):
    C = dict(A)
    for k, c in B.items():
        C[k] = C.get(k, 0) + s * c
    return clean(C)


def scal(A, s):
    return clean({k: s * c for k, c in A.items()})


def mul(A, B):
    C = {}
    for (j1, l1, m1), c1 in A.items():
        for (j2, l2, m2), c2 in B.items():
            k = (j1 + j2, l1 + l2, m1 + m2)
            C[k] = C.get(k, 0) + c1 * c2
    return clean(C)


def _int_term(j, l):
    """Antiderivative in p of p^j (log p)^l, as a dict {(j', l'): coeff} (no constant)."""
    if j == -1:
        return {(0, l + 1): Fr(1, l + 1)}
    out = {}
    # int p^j L^l = p^(j+1) L^l/(j+1) - l/(j+1) int p^j L^(l-1)
    sgn = Fr(1)
    for t in range(l, -1, -1):
        out[(j + 1, t)] = out.get((j + 1, t), 0) + sgn / (j + 1)
        if t == 0:
            break
        sgn = sgn * (-Fr(t, j + 1))
    return out


def integ(A):
    """Antiderivative in p (constants of integration dropped: a constant is a function of x)."""
    C = {}
    for (j, l, m), c in A.items():
        for (j2, l2), d in _int_term(j, l).items():
            if j2 == 0 and l2 == 0:
                continue
            k = (j2, l2, m)
            C[k] = C.get(k, 0) + c * d
    return clean(C)


def deriv(A):
    C = {}
    for (j, l, m), c in A.items():
        if j != 0:
            k = (j - 1, l, m)
            C[k] = C.get(k, 0) + c * j
        if l != 0:
            k = (j - 1, l - 1, m)
            C[k] = C.get(k, 0) + c * l
    return clean(C)


def split_log(A):
    """{l: Laurent dict {(j, m): c}}"""
    out = {}
    for (j, l, m), c in A.items():
        out.setdefault(l, {})[(j, m)] = c
    return out


def lpow(X, n):
    R = {(0, 0, 0): Fr(1)}
    for _ in range(n):
        R = mul(R, X)
    return R


def as_poly_in_x(F, X, top_coeff):
    """Write a log-free Laurent dict F = {(j, m): c} as Sum b_{n,m'} v^m' x^n exactly, x = X (Laurent dict
    with keys (j, 0, m)) whose top term is top_coeff * p^4.  Returns {(n, m'): b} or None if impossible."""
    G = {(j, 0, m): c for (j, m), c in F.items()}
    res = {}
    powers = {}
    guard = 0
    while G:
        guard += 1
        if guard > 10000:
            return None
        jt = max(j for (j, l, m) in G)
        if jt < 0 or jt % 4 != 0:
            return None
        n = jt // 4
        for (j, l, m), c in list(G.items()):
            if j == jt:
                b = c / top_coeff ** n
                res[(n, m)] = res.get((n, m), 0) + b
                if n not in powers:
                    powers[n] = lpow(X, n)
                G = add(G, mul({(0, 0, m): b}, powers[n]), -1)
                break
    return clean(res)


def l_reduce(F, X, w):
    """L-membership test for a log-free Laurent dict F = {(j, m): c} of weight w (p weight 1, v weight 5):
    remove negative powers of p by subtracting v^m'' x^n (lowest term v^n p^-n of x^n has coefficient 1);
    m'' may be negative (functions of x pair to zero).  Returns (ok, remainder): ok iff the remainder lies
    in Q[v][p] (all v-exponents >= 0)."""
    G = {(j, 0, m): c for (j, m), c in F.items()}
    powers = {}
    guard = 0
    while True:
        guard += 1
        if guard > 10000:
            return False, G
        neg = [j for (j, l, m) in G if j < 0]
        if not neg:
            break
        jb = min(neg)
        n = -jb
        for (j, l, m), c in list(G.items()):
            if j == jb:
                if n not in powers:
                    powers[n] = lpow(X, n)
                G = add(G, mul({(0, 0, m - n): c}, powers[n]), -1)
                break
    for (j, l, m) in G:
        assert j + 5 * m == w, ("inhomogeneous", j, m, w)
    ok = all(m >= 0 for (j, l, m) in G)
    return ok, G
