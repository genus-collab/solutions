# -*- coding: utf-8 -*-
"""cyc_ref_jet.py -- Stage-5 referee 'cyc_formal': jets, the symbolic-n calculus, PF4, small linear algebra (stdlib only).

  Jet  bilinear differential polynomial Sum q_ab f_a g_b in two differential indeterminates f, g (f_a = D^a f), q_ab exact
       polynomials in u (Ser with N = EX); D(q f_a g_b) = (Dq) f_a g_b + q f_(a+1) g_b + q f_a g_(b+1).
  PN   polynomial in a symbol n (list of Fractions, low degree first); a shifted monomial is a dict s -> PN meaning
       Sum_s p_s(n) u^(n+s), and D(p(n) u^(n+s)) = p(n) (n+s) (u^(n+s) - u^(n+s-1)) EXACTLY in n.
  PF4  A_k of CYC s.0 / fq_conifold (u-form), built here from the large-radius form, not typed in (see cyc_ref_alg.py).
Imports only cyc_ref_core and cyc_ref_ops (this referee's own modules).
"""
from fractions import Fraction as Fr
from math import comb
from itertools import permutations
from cyc_ref_core import Ser, ONE, U, EX
from cyc_ref_ops import Op, REl


def poly(d):
    return Ser.poly({e: Fr(x) for e, x in d.items()})


# ---------------------------------------------------------------- PF4 from D^4 - 5 q prod_(k=1..4) (5D + k), q = (1-u)/3125
def pf4_from_large_radius(order='left'):
    q = poly({0: Fr(1, 3125), 1: Fr(-1, 3125)})
    prod = Op.fn(ONE)
    for k in range(1, 5):
        prod = prod * Op({1: poly({0: 5}), 0: poly({0: k})})
    qprod = Op.fn(q.sc(5)) * prod if order == 'left' else prod * Op.fn(q.sc(5))
    return Op.Dk(4) - qprod


def ops_equal_exact(a, b):
    c = a - b
    return all(not s.c for s in c.A.values())


# ---------------------------------------------------------------- jets
def jet_add(J, K):
    d = dict(J)
    for k, s in K.items():
        d[k] = d[k] + s if k in d else s
    return d


def jet_neg(J):
    return {k: -s for k, s in J.items()}


def jet_D(J):
    out = {}
    for (a, b), q in J.items():
        for key, s in (((a, b), q.D()), ((a + 1, b), q), ((a, b + 1), q)):
            out[key] = out[key] + s if key in out else s
    return out


def jet_iszero(J):
    return all(not s.c for s in J.values())


def jet_swap(J):
    return {(b, a): s for (a, b), s in J.items()}


def concomitant_jet(A, signs=True, skip=()):
    """P(f, g) = Sum_(k=1..4) Sum_(j<k) (-1)^j D^j(A_k g) D^(k-1-j) f, expanded by Leibniz; key (f-order, g-order)."""
    out = {}
    for k, a in A.items():
        if k == 0 or k in skip:
            continue
        for j in range(k):
            dk = [a]
            for _ in range(j):
                dk.append(dk[-1].D())
            for i in range(j + 1):          # D^j(A_k g) = Sum_i C(j,i) D^(j-i)(A_k) g_i
                s = dk[j - i].sc(Fr((-1) ** j if signs else 1) * comb(j, i))
                key = (k - 1 - j, i)
                out[key] = out[key] + s if key in out else s
    return out


def lagrange_lhs(A, Aad):
    """g L f - f L* g as a jet: Sum_k A_k f_k g_0 - Sum_k (L*)_k f_0 g_k."""
    out = {}
    for k, a in A.items():
        out[(k, 0)] = out[(k, 0)] + a if (k, 0) in out else a
    for k, a in Aad.items():
        out[(0, k)] = out[(0, k)] - a if (0, k) in out else -a
    return out


def jet_eval(J, fj, gj):
    """Sum_ab D^a f . q_ab . D^b g on REl jets fj[a] = D^a f, gj[b] = D^b g"""
    out = None
    for (a, b), q in J.items():
        t = (fj[a] * gj[b]).mulS(q)
        out = t if out is None else out + t
    return out


def jets_of(y, n=4):
    out = [y]
    for _ in range(n - 1):
        out.append(out[-1].D())
    return out


# ---------------------------------------------------------------- polynomials in n and shifted monomials
def pn_trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def pn_add(p, q):
    n = max(len(p), len(q))
    return pn_trim([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])


def pn_mul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1) if p and q else []
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return pn_trim(out)


def pn_eval(p, x):
    return sum(Fr(c) * Fr(x) ** i for i, c in enumerate(p))


def sm_D(X):
    out = {}
    for s, p in X.items():
        t = pn_mul(p, [Fr(s), Fr(1)])            # p(n) (n + s)
        out[s] = pn_add(out.get(s, []), t)
        out[s - 1] = pn_add(out.get(s - 1, []), [-c for c in t])
    return {s: p for s, p in out.items() if p}


def sm_apply(Aop):
    """PF u^n = Sum_s c_s(n) u^(n+s) for an operator with exact polynomial coefficients (symbolic n)."""
    out, cur = {}, {0: [Fr(1)]}
    for k in range(max(Aop.A) + 1):
        if k in Aop.A:
            a = Aop.A[k]
            for i, x in enumerate(a.c):
                e = a.v + i
                for s, p in cur.items():
                    out[s + e] = pn_add(out.get(s + e, []), [x * c for c in p])
        cur = sm_D(cur)
    return {s: p for s, p in out.items() if p}


# ---------------------------------------------------------------- linear algebra over Q, and det over R
def det_inv(M):
    n = len(M); a = [[Fr(x) for x in r] + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]; det = Fr(1)
    for col in range(n):
        piv = next((r for r in range(col, n) if a[r][col] != 0), None)
        if piv is None:
            return Fr(0), None
        if piv != col:
            a[col], a[piv] = a[piv], a[col]; det = -det
        det *= a[col][col]; pv = a[col][col]; a[col] = [x / pv for x in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                fac = a[r][col]; a[r] = [x - fac * y for x, y in zip(a[r], a[col])]
    return det, [row[n:] for row in a]


def pf4(B):
    return B[0][1] * B[2][3] - B[0][2] * B[1][3] + B[0][3] * B[1][2]


def perm_sign(p):
    s, p = 1, list(p)
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


def det_R(W):
    """determinant of a 4x4 matrix of REl (Leibniz formula)"""
    out = None
    for p in permutations(range(4)):
        t = W[0][p[0]] * W[1][p[1]] * W[2][p[2]] * W[3][p[3]]
        t = t if perm_sign(p) > 0 else -t
        out = t if out is None else out + t
    return out


def matT_B_mat(T, B):
    """T^T B T over Q"""
    n = len(B)
    return [[sum(T[a][i] * B[a][b] * T[b][j] for a in range(n) for b in range(n)) for j in range(n)] for i in range(n)]
