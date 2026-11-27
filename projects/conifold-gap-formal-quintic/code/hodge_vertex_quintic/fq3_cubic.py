# -*- coding: utf-8 -*-
"""fq3_cubic.py -- the O(u^3) cubics of the three-edge graphs at genus 3 (genus3_cubic programme).

For M = M_0 exp(u N_1 + ...) (right placement) every primed edge factor is E' = u E'_1 + O(u^2) with E'_1 affine in N_1 and vanishing at
N_1 = N_1^taut; leaves are O(u), the vertex tensors are Omega^{M_0} + O(u), xhat = 1 + O(u^3).  Hence for a THREE-edge graph Gamma
    [u^3] G_Gamma[M_0 exp(u N_1 + u^2 N_2 + ...)] = c_Gamma(Delta),   Delta = N_1 - N_1^taut,
a homogeneous cubic form in Delta, independent of N_2, N_3 (three edge factors E'_1, vertex tensors at M_0, no leaves).
A three-edge graph has Sum_v dim = 6 - E = 3, so an edge of zhat-order k (psi-degree k - 1) can enter only with Sum (k_e - 1) <= 3: the cubics
see the orders 1..4 of Delta only, and order 4 only linearly (monomials of type (1, 1, 4)).
Computation as in fq3_quad: GraphSummer6 with hatted rows at u-precision K = 4, M = MU(R, M_0, [N_1], 4), CONSTANT vertex provider
MVertex6(R, M_0) (exact for the u^3 coefficient of three-edge graphs and for the u^2 coefficient of two-edge graphs; wrong beyond).
API: job_cubic((frame4, par75 or None, idxs, tag)) -> (tag, {idx: [u^0..u^3]});
     CubicForm(p, nvar, coef) homogeneous cubic on F_p^nvar (dict {(i, j, k) i <= j <= k: coefficient}); fit_cubic(evalf, nvar, p).
"""
import os, sys, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq3_m import NSP, sp_mats

IDX_3EDGE = [4, 9, 10, 11, 12, 18, 19, 20, 32]


def job_cubic(args):
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
        s = GS.graph(st['graphs'][i]); out[i] = [s.coef(k) if k < s.prec else None for k in range(4)]
    return tag, out


class CubicForm(object):
    """homogeneous cubic c(t) = Sum_{i<=j<=k} a_ijk t_i t_j t_k over F_p"""
    def __init__(self, p, nvar, coef):
        self.p = p; self.n = nvar; self.coef = {k: v % p for k, v in coef.items() if v % p}

    def eval(self, t):
        p = self.p; s = 0
        for (i, j, k), a in self.coef.items():
            if t[i] and t[j] and t[k]: s = (s + a * t[i] * t[j] * t[k]) % p
        return s

    def gradient(self, t):
        p = self.p; g = [0] * self.n
        for (i, j, k), a in self.coef.items():
            g[i] = (g[i] + a * t[j] * t[k]) % p; g[j] = (g[j] + a * t[i] * t[k]) % p; g[k] = (g[k] + a * t[i] * t[j]) % p
        return g

    def restrict(self, B):
        """t = Sum_m s_m B[m] (B: list of vectors of length self.n): CubicForm in s"""
        p = self.p; m = len(B); out = {}
        for (i, j, k), a in self.coef.items():
            for x in range(m):
                bi = B[x][i]
                if not bi: continue
                for y in range(m):
                    bj = B[y][j]
                    if not bj: continue
                    ab = a * bi * bj % p
                    for z in range(m):
                        bk = B[z][k]
                        if not bk: continue
                        key = tuple(sorted((x, y, z))); out[key] = (out.get(key, 0) + ab * bk) % p
        return CubicForm(p, m, out)

    def support_vars(self):
        return sorted(set(v for key in self.coef for v in key))


def cubic_points(nvar):
    """polarisation points for a homogeneous cubic: e_i; e_i + e_j, e_i + 2 e_j (i < j); e_i + e_j + e_k (i < j < k)"""
    pts = []
    for i in range(nvar): pts.append((('e', i), {i: 1}))
    for i in range(nvar):
        for j in range(i + 1, nvar):
            pts.append((('ee', i, j), {i: 1, j: 1})); pts.append((('e2e', i, j), {i: 1, j: 2}))
    for i, j, k in itertools.combinations(range(nvar), 3): pts.append((('eee', i, j, k), {i: 1, j: 1, k: 1}))
    return pts


def cubic_from_values(vals, nvar, p):
    """vals: {tag: value of the homogeneous cubic at the polarisation point}; returns CubicForm"""
    inv2 = pow(2, p - 2, p); a = {}
    for i in range(nvar): a[(i, i, i)] = vals[('e', i)] % p
    for i in range(nvar):
        for j in range(i + 1, nvar):
            # c(e_i + e_j) = a_iii + a_jjj + a_iij + a_ijj ;  c(e_i + 2 e_j) = a_iii + 8 a_jjj + 2 a_iij + 4 a_ijj
            s1 = (vals[('ee', i, j)] - a[(i, i, i)] - a[(j, j, j)]) % p
            s2 = (vals[('e2e', i, j)] - a[(i, i, i)] - 8 * a[(j, j, j)]) % p
            aijj = (s2 - 2 * s1) * inv2 % p; aiij = (s1 - aijj) % p
            a[(i, i, j)] = aiij; a[(i, j, j)] = aijj
    for i, j, k in itertools.combinations(range(nvar), 3):
        v = vals[('eee', i, j, k)]
        for x, y in ((i, j), (i, k), (j, k)):
            v -= a[(x, x, y)] + a[(x, y, y)]
        v -= a[(i, i, i)] + a[(j, j, j)] + a[(k, k, k)]
        a[(i, j, k)] = v % p
    return CubicForm(p, nvar, a)
