# -*- coding: utf-8 -*-
"""fq3_circsys.py -- the complete system of N_1-conditions of genus 2 and 3 on the circulant class, in the 15 circulant coordinates of
Delta = N_1 - N_1^taut:  c1_0 c1_1 c1_2 | c2_1 c2_2 | c3_0 c3_1 c3_2 | c4_1 c4_2 | c5_0 c5_1 c5_2 | c6_1 c6_2   (indices 0..14).

CircSys(prime_tag) loads: the two genus-2 O(u) rows L2 (orders <= 3), the two genus-3 rows J3, the five two-edge quadrics (15 x 15 matrices,
homogeneous), the nine three-edge cubics (on the coordinates 0..9), the linear form m_1 (c_32 = kappa m_1^3).
  .low(c1, tau, w)        the 8 low coordinates with c2_2 from m_1 = 0 and (c3_0, c3_1) from the genus-2 rows (linear in (c1, tau, w))
  .cub(i, v), .quad(i, v), .quad_grad(i, v), .residuals(v)   evaluations at a 15-vector
  .fibre(low8)            for a given low part: numeric data of the linear systems for c4 (cubics 9, 10, 12, 11, 19) and for (c5, c6)
Univariate helpers over F_p (coefficient lists, low to high): interp_simple, sylvester_res, pgcd, ptrim, peval.
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_rref, fp_solve
from fq3_quad import QuadForm, IDX_2EDGE
from fq3_cubic import IDX_3EDGE, CubicForm
from fq3_m import NSP, class_columns
import fq3_quadlib as L


def ptrim(a, p):
    a = [x % p for x in a]
    while a and a[-1] == 0: a.pop()
    return a


def peval(a, x, p):
    r = 0
    for c in reversed(a): r = (r * x + c) % p
    return r


def pgcd(a, b, p):
    a = ptrim(a, p); b = ptrim(b, p)
    while b:
        inv = pow(b[-1], p - 2, p); a = list(a)
        while len(a) >= len(b):
            f = a[-1] * inv % p; s = len(a) - len(b)
            for i, x in enumerate(b): a[s + i] = (a[s + i] - f * x) % p
            a = ptrim(a, p)
            if not a: break
        a, b = b, a
    if a:
        inv = pow(a[-1], p - 2, p); a = [x * inv % p for x in a]
    return a


def interp_simple(xs, ys, p):
    """Lagrange interpolation, O(n^2), via the master polynomial"""
    n = len(xs); master = [1]
    for x in xs:
        new = [0] * (len(master) + 1)
        for j, c in enumerate(master):
            new[j] = (new[j] - x * c) % p; new[j + 1] = (new[j + 1] + c) % p
        master = new
    out = [0] * n
    for i in range(n):
        # q = master / (x - xs[i]) by synthetic division
        q = [0] * n; carry = 0
        for j in range(n, 0, -1):
            carry = (master[j] + carry * xs[i]) % p; q[j - 1] = carry
        den = peval(q, xs[i], p); w = ys[i] * pow(den, p - 2, p) % p
        if w:
            for j in range(n): out[j] = (out[j] + w * q[j]) % p
    return out


def sylvester_res(f, g, df, dg, p):
    """resultant of f, g regarded as polynomials of formal degrees df, dg (coefficient lists low to high, padded)"""
    f = (list(f) + [0] * (df + 1))[:df + 1]; g = (list(g) + [0] * (dg + 1))[:dg + 1]; n = df + dg
    M = []
    for i in range(dg):
        row = [0] * n
        for j in range(df + 1): row[i + j] = f[df - j]
        M.append(row)
    for i in range(df):
        row = [0] * n
        for j in range(dg + 1): row[i + j] = g[dg - j]
        M.append(row)
    return L.det_fp(M, p)


class CircSys(object):
    def __init__(self, tag):
        self.tag = tag; p = self.p = P2 if tag == "89" else P1; R = self.R = FpRing(p)
        qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % tag))); g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % tag)))
        cs = json.load(open(os.path.join(HERE, "g3_cubic_%s.json" % tag))); self.nt = qs['n_taut']
        self.cols, self.names = class_columns(R, 'circ'); assert cs['vis_names'] == self.names[:10]
        self.C = {i: CubicForm(p, 10, {tuple(int(t) for t in k.split(",")): v for k, v in cs['C'][str(i)].items()}) for i in IDX_3EDGE}
        self.Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']).restrict(self.nt, self.cols).A for i in IDX_2EDGE}
        dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
        self.L2 = [[dot(r, self.cols[j]) for j in range(15)] for r in (g3['J2']['2'], g3['J2']['4'])]
        self.J3 = [[dot(r, self.cols[j]) for j in range(15)] for r in (g3['J']['2'], g3['J']['5'])]
        self.m1 = qs['m1_circ_coords']; assert all(self.m1[j] == 0 for j in (0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) and self.m1[4]
        rref, piv = fp_rref([[self.L2[r][j] for j in (5, 6, 7)] + [self.L2[r][j] for j in range(5)] for r in range(2)], 3, p)
        assert list(piv) == [0, 1], piv; self.rref = rref

    def low(self, c1, tau, w):
        p = self.p; m = self.m1
        c22 = (-(m[1] * c1[1] + m[2] * c1[2] + m[3] * tau)) * pow(m[4], p - 2, p) % p
        x = [c1[0] % p, c1[1] % p, c1[2] % p, tau % p, c22]
        c3 = [(-(self.rref[r][2] * w + sum(self.rref[r][3 + j] * x[j] for j in range(5)))) % p for r in range(2)]
        return x + c3 + [w % p]

    def cub(self, i, v): return self.C[i].eval(v[:10])

    def quad(self, i, v):
        p = self.p; A = self.Q[i]; s = 0
        for a in range(15):
            if v[a]: s = (s + v[a] * sum(A[a][b] * v[b] for b in range(15) if v[b])) % p
        return s

    def quad_grad(self, i, v):
        p = self.p; A = self.Q[i]
        return [2 * sum(A[a][b] * v[b] for b in range(15) if v[b]) % p for a in range(15)]

    def residuals(self, v):
        """all conditions at a 15-vector: dict name -> value"""
        p = self.p; out = {}
        for r in range(2): out['L2_%d' % r] = sum(self.L2[r][j] * v[j] for j in range(15)) % p
        for r in range(2): out['J3_%d' % r] = sum(self.J3[r][j] * v[j] for j in range(15)) % p
        for i in IDX_2EDGE: out['q_%d' % i] = self.quad(i, v)
        for i in IDX_3EDGE: out['c_%d' % i] = self.cub(i, v)
        return out

    def to_sp(self, v):
        p = self.p
        return [(self.nt[k] + sum(self.cols[j][k] * v[j] for j in range(15) if v[j])) % p for k in range(NSP)]

    def high_system(self, v10, rows_q=(3, 7, 8)):
        """linear system for u = (c5_0, c5_1, c5_2, c6_1, c6_2) given the coordinates 0..9: rows J3 (2) and the quadrics rows_q; returns (A, rhs)"""
        p = self.p; v = list(v10) + [0] * 5; A = []; rhs = []
        for r in range(2): A.append([self.J3[r][j] for j in range(10, 15)]); rhs.append((-sum(self.J3[r][j] * v[j] for j in range(10))) % p)
        for i in rows_q:
            g = self.quad_grad(i, v); A.append([g[j] for j in range(10, 15)]); rhs.append((-self.quad(i, v)) % p)
        return A, rhs
