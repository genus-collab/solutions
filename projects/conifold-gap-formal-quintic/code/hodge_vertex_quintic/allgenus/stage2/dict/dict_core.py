"""Stage 2, D(iii): exact core for the dictionary K_0 <-> L/(functions of x) on the model curve
x = W = p^4/4 + v/p, y = p.  Pure stdlib, exact over Q.  Written from scratch: imports nothing from critic/ or
any other worker folder (independent of critic/cr_laurent.py and critic/cr_dictionary.py).

Objects.
  Laurent-log element: dict {(j, m, l): Fraction} = sum c p^j v^m (log p)^l   (route A, generic calculus).
  Weight: p -> 1, v -> 5, x -> 4, q = v/p -> 4, log p -> 0.
  Class M = Q[p, q] (+) log p . Q[v, x]   (q = v/p; Q[p, q] = {Laurent polys with every term v^m p^j, j >= -m}).
  Normal form NF (route B, independent of route A's L-test): Q[p, q] = Q[p, x] is free over Q[v, x] with basis
  1, p, .., p^4 because p^5 = 4 x p - 4 v; NF(F) = the coordinates (a_1, .., a_4), a_i in Q[v, x], i.e. F mod Q[v, x].
Legs (DNOPS-type, as in critic/cr_dictionary.py docstring, up to one constant per b):
  Phi_0: b=0 -> p, b=1 -> log p, b=2 -> v/p, b=3 -> p^3, b=4 -> p^2;  Phi_(k+1) = I(Phi_k) := -Int Phi_k dx (p-constants dropped).
"""
from fractions import Fraction as Fr

ZERO = Fr(0)


def clean(A):
    return {k: c for k, c in A.items() if c != 0}


def add(A, B, s=1):
    C = dict(A)
    for k, c in B.items():
        C[k] = C.get(k, ZERO) + s * c
    return clean(C)


def mul(A, B):
    C = {}
    for (j1, m1, l1), c1 in A.items():
        for (j2, m2, l2), c2 in B.items():
            k = (j1 + j2, m1 + m2, l1 + l2)
            C[k] = C.get(k, ZERO) + c1 * c2
    return clean(C)


def weight_set(A):
    return {j + 5 * m for (j, m, l) in A}


# ---------------- route A: generic antiderivative of p^j (log p)^l in p ----------------
def antider(j, l):
    """Int p^j L^l dp, L = log p, as {(j', l'): c}; recursion by parts, no constant."""
    if j == -1:
        return {(0, l + 1): Fr(1, l + 1)}
    out = {(j + 1, l): Fr(1, j + 1)}
    if l > 0:
        for (jj, ll), c in antider(j, l - 1).items():
            out[(jj, ll)] = out.get((jj, ll), ZERO) - Fr(l, j + 1) * c
    return out


def integ_p(A):
    C = {}
    for (j, m, l), c in A.items():
        for (j2, l2), d in antider(j, l).items():
            if (j2, l2) == (0, 0):
                continue                      # drop p-constants (a constant is a function of x)
            k = (j2, m, l2)
            C[k] = C.get(k, ZERO) + c * d
    return clean(C)


DX = {(3, 0, 0): Fr(1), (-2, 1, 0): Fr(-1)}    # dx/dp = p^3 - v p^-2
XL = {(4, 0, 0): Fr(1, 4), (-1, 1, 0): Fr(1)}  # x as a Laurent polynomial


def I_A(A):
    """Descendant step -Int A dx, route A."""
    return {k: -c for k, c in integ_p(mul(A, DX)).items()}


PHI0 = {0: {(1, 0, 0): Fr(1)}, 1: {(0, 0, 1): Fr(1)}, 2: {(-1, 1, 0): Fr(1)}, 3: {(3, 0, 0): Fr(1)},
        4: {(2, 0, 0): Fr(1)}}
WT0 = {0: 1, 1: 0, 2: 4, 3: 3, 4: 2}


def legs_A(kmax, phi0=None):
    phi0 = phi0 or PHI0
    T = {}
    for b in range(5):
        cur = phi0[b]
        for k in range(kmax + 1):
            T[(b, k)] = cur
            cur = I_A(cur)
    return T


# ---------------- membership tests ----------------
def split(A):
    out = {}
    for (j, m, l), c in A.items():
        out.setdefault(l, {})[(j, m)] = c
    return out


def in_Qpq(F):
    """F = {(j, m): c} log-free Laurent: in Q[p, q] iff every term has m >= 0 and j >= -m."""
    return all(m >= 0 and j >= -m for (j, m) in F)


_XP = {}


def xpow(n):
    if n not in _XP:
        R = {(0, 0, 0): Fr(1)}
        for _ in range(n):
            R = mul(R, XL)
        _XP[n] = R
    return _XP[n]


def as_poly_x(F):
    """Log coefficient F = {(j, m): c} as sum b v^m x^n (top-down in p).  Returns {(n, m): b} or None."""
    G = {(j, m, 0): c for (j, m), c in F.items()}
    res = {}
    while G:
        jt = max(j for (j, m, l) in G)
        if jt < 0 or jt % 4:
            return None
        n = jt // 4
        (j, m, l), c = next((k, c) for k, c in G.items() if k[0] == jt)
        b = c * 4 ** n
        res[(n, m)] = res.get((n, m), ZERO) + b
        G = add(G, mul({(0, m, 0): b}, xpow(n)), -1)
    return clean(res)


def in_M(A):
    """(ok, reason): A in Q[p, q] (+) log p Q[v, x]?"""
    sp = split(A)
    if any(l not in (0, 1) for l in sp):
        return False, "log power > 1"
    if not in_Qpq(sp.get(0, {})):
        return False, "non-log part not in Q[p,q]"
    if 1 in sp:
        f = as_poly_x(sp[1])
        if f is None:
            return False, "log coefficient not a polynomial in x"
        if any(m < 0 for (n, m) in f):
            return False, "log coefficient has negative v-power"
    return True, ""
