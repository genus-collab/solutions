"""Stage 2, D(iii): normal forms and an L-membership test by linear algebra (independent of dict_core.in_Qpq).

NF (route B): Q[p, q] = Q[p, x] = (+)_{i=0..4} Q[v, x] p^i via p^5 = 4 x p - 4 v; F mod Q[v, x] <-> (a_1, .., a_4).
L_test: is a log-free Laurent F of weight w in Q[v][p] + Q(v)[x]?  Solve F - sum_n a_n v^((w-4n)/5) x^n in Q[v][p]
for rational a_n by exact Gaussian elimination (v-exponents (w-4n)/5 may be negative).
"""
from fractions import Fraction as Fr
from math import comb

from dict_core import ZERO, xpow, split, as_poly_x


def to_px(F):
    """Log-free F = {(j, m): c} in Q[p, q] -> {(i, m, n): c} = sum c p^i v^m x^n, then reduce p-degree to <= 4."""
    B = {}

    def put(i, m, n, c):
        if c:
            B[(i, m, n)] = B.get((i, m, n), ZERO) + c

    for (j, m), c in F.items():
        if j >= 0:
            put(j, m, 0, c)
        else:
            s = -j                            # v^m p^-s = v^(m-s) q^s,  q = x - p^4/4
            assert m - s >= 0, "not in Q[p,q]"
            for t in range(s + 1):
                put(4 * t, m - s, s - t, c * comb(s, t) * Fr(-1, 4) ** t)
    while True:
        top = [k for k, c in B.items() if c != 0 and k[0] >= 5]
        if not top:
            break
        imax = max(k[0] for k in top)
        for (i, m, n) in [k for k in top if k[0] == imax]:
            c = B.pop((i, m, n))
            put(i - 4, m, n + 1, 4 * c)       # p^i = 4 x p^(i-4) - 4 v p^(i-5)
            put(i - 5, m + 1, n, -4 * c)
    return {k: c for k, c in B.items() if c != 0}


def NF(A):
    """Normal form of A in M modulo Q[v, x]: dict with keys ('p', i, m, n) (i = 1..4) and ('L', m, n)."""
    sp = split(A)
    out = {}
    for (i, m, n), c in to_px(sp.get(0, {})).items():
        if i >= 1:
            out[('p', i, m, n)] = c
    if 1 in sp:
        f = as_poly_x(sp[1])
        assert f is not None
        for (n, m), c in f.items():
            out[('L', m, n)] = c
    return out


def rank_Q(rows):
    """Rank over Q of a list of dicts (column key -> Fraction)."""
    rows = [dict(r) for r in rows if r]
    rk = 0
    while rows:
        r = rows.pop()
        r = {k: c for k, c in r.items() if c != 0}
        if not r:
            continue
        col = min(r)
        piv = r[col]
        rk += 1
        new = []
        for s in rows:
            if s.get(col, 0) != 0:
                f = s[col] / piv
                s = dict(s)
                for k, c in r.items():
                    s[k] = s.get(k, ZERO) - f * c
                s = {k: c for k, c in s.items() if c != 0}
            if s:
                new.append(s)
        rows = new
    return rk


def kernel_Q(rows):
    """Left kernel {coeff vectors y : sum y_i rows_i = 0} as a list of dicts index -> Fraction."""
    n = len(rows)
    M = [({k: c for k, c in r.items() if c != 0}, {i: Fr(1)}) for i, r in enumerate(rows)]
    done = []
    kern = []
    while M:
        r, e = M.pop()
        if not r:
            kern.append(e)
            continue
        col = min(r)
        piv = r[col]
        done.append((col, r, e))
        new = []
        for s, f in M:
            if s.get(col, 0) != 0:
                t = s[col] / piv
                s = dict(s)
                f = dict(f)
                for k, c in r.items():
                    s[k] = s.get(k, ZERO) - t * c
                for k, c in e.items():
                    f[k] = f.get(k, ZERO) - t * c
                s = {k: c for k, c in s.items() if c != 0}
            new.append((s, f))
        M = new
    assert len(kern) + len(done) == n
    return kern


def L_test(F, w):
    """F = {(j, m): c} log-free Laurent, weight w.  True iff F in Q[v][p] + Q(v)[x]."""
    for (j, m) in F:
        assert j + 5 * m == w
    N = max([0] + [-j for (j, m) in F])
    ns = [n for n in range(N + 1) if (w - 4 * n) % 5 == 0]
    # columns = p-exponents j outside [0, w] (these must vanish); unknown a_n multiplies v^e x^n
    cols = {}
    rows = []
    for n in ns:
        e = (w - 4 * n) // 5
        r = {}
        for (j, m, l), c in xpow(n).items():
            jj = j
            if jj < 0 or jj > w:              # the term is c v^(m+e) p^j; j > w <=> v-exponent < 0
                r[jj] = r.get(jj, ZERO) + c
        rows.append(r)
    rhs = {j: c for (j, m), c in F.items() if j < 0 or j > w}
    # solve sum a_n rows_n = rhs: rhs in the row span  <=>  rank(rows) == rank(rows + [rhs])
    return rank_Q(rows + [rhs]) == rank_Q(rows)
