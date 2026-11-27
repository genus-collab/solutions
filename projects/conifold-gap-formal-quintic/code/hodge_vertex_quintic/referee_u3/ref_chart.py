# -*- coding: utf-8 -*-
"""ref_chart.py -- the chart c1_0 = 1 (c1 = (1, a, b), tau = c2_1, w = c3_2): symbolic matrices in the 3 variables (a, b, t).
  E-matrix (the script's): six c4-rows after w = -B4/a4, unknown a4 c4:  cubics 9, 10, 11, 12, 19: [A1, A2, a4 P - B4 S];
                           q_3: [a4 g1, a4 g2, a4^2 q0 - a4 B4 q1 + B4^2 q2]                                   (6 x 3, rank <= 2 needed)
  D-matrix (mine, no scaling): the six cubics 4, 9, 10, 11, 12, 19 are affine in (c4_1, c4_2, w): rows [A1, A2, S, P]   (6 x 4, rank <= 3
                           needed; valid also on a4 = 0)."""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *

CUB6 = ['c_4', 'c_9', 'c_10', 'c_11', 'c_12', 'c_19']


def chart_data(P, p):
    """returns dict with 3-variable polynomials (a, b, t)"""
    C1 = {k: msubs(f, {0: 1}, p) for k, f in P.items()}
    def split(f):
        """f affine in c4_1, c4_2 (indices 5, 6) and polynomial in w (index 4): returns (A1, A2, [w^0, w^1, ...]) as 3-variable polys"""
        assert all(e[5] + e[6] <= 1 and not any(e[7:]) for e in f)
        A1 = {e: c for e, c in f.items() if e[5]}; A2 = {e: c for e, c in f.items() if e[6]}; rest = {e: c for e, c in f.items() if not e[5] and not e[6]}
        assert all(e[4] == 0 for e in A1) and all(e[4] == 0 for e in A2), "c4-coefficient depends on w"
        t3 = lambda g: {(e[1], e[2], e[3]): c for e, c in g.items()}
        wc = mcoeffs(rest, 4) if rest else [{}]
        return t3(A1), t3(A2), [t3(c) for c in wc]
    D = {}; rows = {}
    for k in CUB6 + ['q_3']:
        A1, A2, wc = split(C1[k]); rows[k] = (A1, A2, wc)
    D['rows'] = rows; D['C1'] = C1
    a4 = rows['c_4'][2][1]; B4 = rows['c_4'][2][0]; assert not rows['c_4'][0] and not rows['c_4'][1] and len(rows['c_4'][2]) == 2
    D['a4'] = a4; D['B4'] = B4
    mm = lambda x, y: mmul(x, y, p)
    Em = []
    for k in ('c_9', 'c_10', 'c_11', 'c_12', 'c_19'):
        A1, A2, wc = rows[k]; assert len(wc) <= 2; Sg = wc[1] if len(wc) > 1 else {}
        Em.append([A1, A2, msub(mm(a4, wc[0]), mm(B4, Sg), p)])
    A1, A2, wc = rows['q_3']; wc = wc + [{}] * (3 - len(wc))
    Em.append([mm(a4, A1), mm(a4, A2), madd(msub(mm(mm(a4, a4), wc[0]), mm(mm(a4, B4), wc[1]), p), mm(mm(B4, B4), wc[2]), p)])
    D['E'] = Em
    Dm = []
    for k in CUB6:
        A1, A2, wc = rows[k]; Dm.append([A1, A2, wc[1] if len(wc) > 1 else {}, wc[0]])
    D['D'] = Dm
    A1, A2, wc = rows['q_3']; assert len(wc) <= 2, "q_3 not affine in w"
    D['D7'] = Dm + [[A1, A2, wc[1] if len(wc) > 1 else {}, wc[0]]]            # q_3 is affine in (c4, w) as well on {m_1 = 0, L2}
    return D


def minors(M, size, p, nv=3):
    """all maximal minors of a matrix with `size` columns: dict {row tuple: poly}"""
    out = {}
    for T in itertools.combinations(range(len(M)), size):
        out[T] = mdet([M[r] for r in T], p, nv)
    return out


def combos(ms, lam, p):
    """random combinations: lam = list of coefficient lists (one per combination) over the sorted minors"""
    keys = sorted(ms); out = []
    for l in lam:
        s = {}
        for c, k in zip(l, keys):
            if ms[k]: s = madd(s, mscale(ms[k], c, p), p)
        out.append(s)
    return out


def tau_coeffs(f):
    """coefficients in t (index 2) of a 3-variable poly: list of 2-variable-supported dicts (still 3-tuples, t exponent 0)"""
    return mcoeffs(f, 2) if f else [{}]


def eval_ab(cf, a, b, p):
    """evaluate a list of (a, b)-polynomials (3-tuples with t exponent 0) at numbers -> list of numbers"""
    out = []
    for c in cf:
        s = 0
        for (i, j, _), v in c.items(): s += v * pow(a, i, p) * pow(b, j, p)
        out.append(s % p)
    return out


class Fast2(object):
    """fast repeated evaluation of a list of bivariate polynomials at many b for fixed a: specialise a first"""
    def __init__(self, cf, p):
        self.p = p; self.polys = []
        for c in cf:
            byj = {}
            for (i, j, _), v in c.items(): byj.setdefault(j, []).append((i, v))
            self.polys.append(byj)

    def at_a(self, a):
        p = self.p; out = []; pw = [1]
        for _ in range(40): pw.append(pw[-1] * a % p)
        for byj in self.polys:
            if not byj: out.append([]); continue
            u = [0] * (max(byj) + 1)
            for j, lst in byj.items(): u[j] = sum(v * pw[i] for i, v in lst) % p
            out.append(u)
        return out
