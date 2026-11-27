# -*- coding: utf-8 -*-
"""fq3_quad.py -- the O(u^2) quadrics of the two-edge graphs at genus 3 (genus3_quad programme).

For M = M_0 exp(u N_1) (right placement) every primed edge factor is E' = u E'_1 + O(u^2) with E'_1 AFFINE in N_1 (R' = R^(u) exp(-u N_1)
M_0^{-1} = Psi + u R'_1(N_1) + O(u^2), R'_1(N_1) = R^_1 M_0^{-1} - R^_0 N_1 M_0^{-1}); leaves are O(u); the vertex tensors are Omega^{M_0} + O(u)
and xhat = 1 + O(u^3).  Hence for a two-edge graph Gamma (no leaves can contribute)
    [u^2] G_Gamma[M_0 exp(u N_1)] = Q_Gamma(N_1) := Sum_dec (1/|Aut|) Omega^{M_0}-vertices (x) E'_1 (x) E'_1,
a quadratic polynomial in the 75 sp coordinates of N_1, and Q_Gamma(N_1) = q_Gamma(N_1 - N_1^taut) is HOMOGENEOUS in Delta = N_1 - N_1^taut
because E'_1 vanishes at the tautological correction N_1^taut = M_0^{-1} R~_1 (R' = Psi + O(u^2) there).
Computation: the existing GraphSummer6 with the hatted rows at u-precision K = 3, the u-dependent M = MU(R, M_0, [N_1], 3) and the
CONSTANT vertex provider MVertex6(R, M_0) (u^0 tensors only; exact for the u^2 coefficient of two-edge graphs, wrong beyond): one
evaluation of Q_Gamma for the five graphs 3, 6, 7, 8, 17 costs a few seconds.  Polarisation on the 2926 points N_1 = 0, e_i, 2 e_i,
e_i + e_j gives the exact coefficients  Q = c_0 + Sum_i b_i n_i + Sum_i A_ii n_i^2 + Sum_{i<j} A_ij n_i n_j.
API: job_quad((frame3, n_vector_or_None, idxs)) -> (tag, {idx: [u^0, u^1, u^2]}); points(); coefficients(vals, p) -> QuadForm;
     QuadForm(c0, b, A) with .eval(n), .restrict(x0, B) (affine change n = x0 + B t), .hom_part_matrix(); fp rank helpers.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq3_m import NSP, sp_mats

IDX_2EDGE = [3, 6, 7, 8, 17]


def job_quad(args):
    """args = (frame (with K = 3), par (75 coordinates of N_1) or None for N_1 = 0, idxs, tag); returns (tag, {idx: [u^0, u^1, u^2]})"""
    frame, par, idxs, tag = args
    from fq3_par import setup
    from fq3_m import MU
    from fq3_mvertex import MVertex6
    from fq3_msum import GraphSummer6
    st = setup(frame); R = st['R']; p = R.p; M0 = st['M0']
    N = sp_mats(par if par is not None else [0] * NSP, p)
    M = MU(R, M0, [N], st['K'], 'right')
    mv = st.get('mv0')
    if mv is None: mv = MVertex6(R, M0); st['mv0'] = mv
    GS = GraphSummer6(st['H'], M, mv); out = {}
    for i in idxs:
        s = GS.graph(st['graphs'][i]); out[i] = [s.coef(k) if k < s.prec else None for k in range(3)]
    return tag, out


def points():
    """the 2926 evaluation points as (tag, coordinate vector): tag = ('0',), ('e', i), ('2e', i), ('ee', i, j)"""
    pts = [(('0',), [0] * NSP)]
    for i in range(NSP):
        v = [0] * NSP; v[i] = 1; pts.append((('e', i), v))
        w = [0] * NSP; w[i] = 2; pts.append((('2e', i), w))
    for i in range(NSP):
        for j in range(i + 1, NSP):
            v = [0] * NSP; v[i] = 1; v[j] = 1; pts.append((('ee', i, j), v))
    return pts


class QuadForm(object):
    """Q(n) = c0 + b . n + n^T A n with A symmetric (A_ii on the diagonal, A_ij/2 off the diagonal), over F_p"""
    def __init__(self, p, c0, b, A):
        self.p = p; self.c0 = c0 % p; self.b = [x % p for x in b]; self.A = [[x % p for x in r] for r in A]; self.n = len(b)

    def eval(self, n):
        p = self.p; v = self.c0
        for i in range(self.n):
            if n[i]: v = (v + self.b[i] * n[i]) % p
        for i in range(self.n):
            if not n[i]: continue
            row = self.A[i]; s = 0
            for j in range(self.n):
                if n[j]: s = (s + row[j] * n[j]) % p
            v = (v + n[i] * s) % p
        return v

    def gradient(self, n):
        """b + 2 A n"""
        p = self.p
        return [(self.b[i] + 2 * sum(self.A[i][j] * n[j] for j in range(self.n))) % p for i in range(self.n)]

    def restrict(self, x0, B):
        """n = x0 + B t (B: list of columns, each a vector of length self.n): QuadForm in t"""
        p = self.p; m = len(B)
        Ax0 = [sum(self.A[i][j] * x0[j] for j in range(self.n)) % p for i in range(self.n)]
        c0 = (self.c0 + sum(self.b[i] * x0[i] for i in range(self.n)) + sum(x0[i] * Ax0[i] for i in range(self.n))) % p
        AB = [[sum(self.A[i][j] * col[j] for j in range(self.n)) % p for i in range(self.n)] for col in B]      # A B_k
        b = [(sum(self.b[i] * B[k][i] for i in range(self.n)) + 2 * sum(x0[i] * AB[k][i] for i in range(self.n))) % p for k in range(m)]
        A = [[sum(B[k][i] * AB[l][i] for i in range(self.n)) % p for l in range(m)] for k in range(m)]
        return QuadForm(p, c0, b, A)


def coefficients(vals, p):
    """vals: {tag: value}; returns QuadForm on F_p^75"""
    inv2 = pow(2, p - 2, p); c0 = vals[('0',)]
    b = [0] * NSP; A = [[0] * NSP for _ in range(NSP)]
    for i in range(NSP):
        q1 = vals[('e', i)]; q2 = vals[('2e', i)]
        A[i][i] = (q2 - 2 * q1 + c0) * inv2 % p
        b[i] = (4 * q1 - q2 - 3 * c0) * inv2 % p
    for i in range(NSP):
        for j in range(i + 1, NSP):
            aij = (vals[('ee', i, j)] - vals[('e', i)] - vals[('e', j)] + c0) % p
            A[i][j] = A[j][i] = aij * inv2 % p
    return QuadForm(p, c0, b, A)


# ------------------------------------------------------------------ F_p linear algebra helpers
def rank_fp(M, p):
    from fq_mfit import fp_rref
    if not M: return 0
    return len(fp_rref([list(r) for r in M], len(M[0]), p)[1])


def nullspace_fp(rows, ncols, p):
    from fq_mfit import fp_nullspace
    return fp_nullspace([list(r) for r in rows], ncols, p)


def sym_rank(A, p): return rank_fp(A, p)


def hom_eval(A, v, p):
    """v^T A v"""
    n = len(v); s = 0
    for i in range(n):
        if not v[i]: continue
        row = A[i]; t = 0
        for j in range(n):
            if v[j]: t = (t + row[j] * v[j]) % p
        s = (s + v[i] * t) % p
    return s


def sqrt_fp(a, p):
    """square root mod p for p = 3 mod 4 (both primes used), or None"""
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    r = pow(a, (p + 1) // 4, p)
    return r if r * r % p == a else None
