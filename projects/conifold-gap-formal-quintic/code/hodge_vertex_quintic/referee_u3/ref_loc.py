# -*- coding: utf-8 -*-
"""ref_loc.py -- shared pieces for the referee scripts: the stored symbolic system (ref_polys_<prime>.json, written by ref1_structure.py),
ring evaluation of dict polynomials, the linear system for the high unknowns, Jacobian ranks."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sys import *

NV = 12; LOWV = [0, 1, 2, 3, 4]; HIGHV = [5, 6, 7, 8, 9, 10, 11]
CONDS = ['J3_0', 'J3_1', 'q_3', 'q_6', 'q_7', 'q_8', 'q_17', 'c_4', 'c_9', 'c_10', 'c_11', 'c_12', 'c_18', 'c_19', 'c_20', 'c_32']


def load_polys(tag):
    d = json.load(open(os.path.join(HERE, "ref_polys_%s.json" % tag)))
    return {k: {tuple(e): c for e, c in f} for k, f in d['P'].items()}


def mevalR(f, vals, R):
    """evaluate a dict polynomial at ring elements vals (list indexed by variable)"""
    s = R.zero; pw = {}
    for e, c in f.items():
        t = None
        for i, k in enumerate(e):
            for _ in range(k): t = vals[i] if t is None else R.mul(t, vals[i])
        s = R.add(s, R.const(c) if t is None else R.scal(c, t))
    return s


def mdiff(f, i, p):
    out = {}
    for e, c in f.items():
        if e[i]:
            k = e[:i] + (e[i] - 1,) + e[i + 1:]; out[k] = c * e[i] % p
    return {k: v for k, v in out.items() if v}


def high_system(P, lowvals, R, p):
    """every condition is affine in the 7 high variables: rows (coefficients over R) and right-hand sides at the given low values
    lowvals = [z, x, y, t, w] (ring elements)"""
    vals = list(lowvals) + [R.zero] * 7; rows = []; rhs = []; names = []
    for k in CONDS:
        f = P[k]
        if not f: continue
        assert all(sum(e[i] for i in HIGHV) <= 1 for e in f), "not affine in the high variables: " + k
        row = []
        for i in HIGHV:
            cs = mcoeffs(f, i); row.append(mevalR(cs[1], vals, R) if len(cs) > 1 else R.zero)
        const = {e: c for e, c in f.items() if not any(e[i] for i in HIGHV)}
        rows.append(row); rhs.append(R.sub(R.zero, mevalR(const, vals, R))); names.append(k)
    return rows, rhs, names


def rank_field(R, M):
    M = [list(r) for r in M]; r = 0; n = len(M[0])
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if not R.iszero(M[i][c])), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; inv = R.inv(M[r][c])
        for i in range(r + 1, len(M)):
            if not R.iszero(M[i][c]):
                f = R.mul(M[i][c], inv); M[i] = [R.sub(x, R.mul(f, y)) for x, y in zip(M[i], M[r])]
        r += 1
        if r == len(M): break
    return r


def jac_rank(P, vals12, R, p):
    """rank of the Jacobian of the 16 conditions (on {m_1 = 0, genus-2 rows}) with respect to the 12 variables at a point (ring elements)"""
    J = [[mevalR(mdiff(P[k], i, p), vals12, R) for i in range(NV)] for k in CONDS if P[k]]
    return rank_field(R, J)


def full_lift(Y, P, lowvals, R, p):
    """solve for the high part and verify EVERY condition with the direct evaluator Sys.residuals (not through the symbolic polynomials).
    returns (15-vector or None, rank, number of free high directions, list of non-zero residual names)"""
    rows, rhs, names = high_system(P, lowvals, R, p); sol, rk, free = solve_field(R, rows, rhs)
    if sol is None: return None, rk, len(free), ['inconsistent']
    v = Y.low(R, *lowvals) + sol; res = Y.residuals(R, v)
    return v, rk, len(free), sorted(k for k, x in res.items() if not R.iszero(x))


def to3(f, idx=(1, 2, 3)):
    """12-variable dict with support in idx -> 3-variable dict"""
    out = {}
    for e, c in f.items():
        assert all(e[i] == 0 for i in range(len(e)) if i not in idx), "unexpected variable"
        out[tuple(e[i] for i in idx)] = c
    return out


def locus0_forms(P, p):
    """on z = 0: a4, B4, a10, B10, Hw, H0 (H = Hw w + H0 = det [c4-coefficients | value] of c_9, q_3, q_8) as forms in (x, y, t);
    G0 = a4 B10 - a10 B4 (degree 5), HH = a4 H0 - B4 Hw (degree 7)"""
    Z = {k: msubs(f, {0: 0}, p) for k, f in P.items()}
    c4 = mcoeffs(Z['c_4'], 4); c10 = mcoeffs(Z['c_10'], 4); assert len(c4) == 2 and len(c10) == 2
    def coef(f, i):
        cs = mcoeffs(f, i); return cs[1] if len(cs) == 2 else {}
    A = [[coef(Z[k], 5), coef(Z[k], 6), msubs(Z[k], {5: 0, 6: 0}, p)] for k in ('c_9', 'q_3', 'q_8')]
    H = mdet(A, p, NV); h = mcoeffs(H, 4); assert len(h) == 2
    F = dict(a4=to3(c4[1]), B4=to3(c4[0]), a10=to3(c10[1]), B10=to3(c10[0]), Hw=to3(h[1]), H0=to3(h[0]))
    F['G0'] = msub(mmul(F['a4'], F['B10'], p), mmul(F['a10'], F['B4'], p), p)
    F['HH'] = msub(mmul(F['a4'], F['H0'], p), mmul(F['B4'], F['Hw'], p), p)
    F['A'] = [[to3(r[0]), to3(r[1])] for r in A]
    return F, Z
