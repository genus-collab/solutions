# -*- coding: utf-8 -*-
"""fq_leadline.py -- the common zero set of the three pure-c_1 equations (graphs 3, 6, 7) of the CIRCULANT class in the 3-dimensional
c_1-space (genus2_symbol step 3 (ii)), when the pairwise resultants eliminate to zero (positive-dimensional locus).

  solve_on_plane(eqs, lin, p): the three t-polys restricted to an affine plane c_1 = a + s b + t c (variables 0 = s, 1 = t);
     Res_s of the pairs, gcd in t, F_p roots, back substitution -> points of the locus on that plane, plus the gcd degree
     (= the degree of the locus if it is a curve meeting a general plane in gcd-degree points over the closure).
  line_through(P, Q): parametrisation; on_line_check(eqs, P, d, p): the univariate restrictions along P + t d are identically zero.
  residual(eqs, P, d, e1, e2, p): coordinates c_1 = P + u d + v e1 + w e2 (line = {v = w = 0}); write each L = v A + w B and
     solve the strict-transform systems of the two blow-up charts (w = lam v: A(u,v,lam v) + lam B(u,v,lam v); v = lam w) with
     fq_leadmodes.solve_cubics3 (general degrees): solutions with v != 0 (resp. w != 0) are points of the locus OFF the line;
     an identically vanishing eliminant signals a further positive-dimensional component.  Returns the two chart results.
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_leadclass import tp_vars, tp_coeff_in, tp_to_uni, tp_eval, uni_gcd_all, uni_roots, dd_pattern, tp_add, tp_scale
from fq_leadmodes import tp_resultant, solve_cubics3, tp_eval_partial


def affine_substitute(tp, lin, p):
    """lin: list over the t-poly's variables of {new monomial tuple: coeff} (affine forms: () = constant); returns t-poly"""
    out = {}
    for m, c in tp.items():
        terms = {(): c}
        for j in m:
            new = {}
            for mm, cc in terms.items():
                for k, coef in lin[j].items():
                    m2 = tuple(sorted(mm + k)); new[m2] = (new.get(m2, 0) + cc * coef) % p
            terms = new
        for mm, cc in terms.items(): out[mm] = (out.get(mm, 0) + cc) % p
    return {m: c for m, c in out.items() if c}


def solve_on_plane(eqs, a, b, c, p):
    lin = [{(): a[j], (0,): b[j], (1,): c[j]} for j in range(3)]
    e = [affine_substitute(E, lin, p) for E in eqs]
    rs = [tp_resultant(e[i], e[j], 0, p) for (i, j) in ((0, 1), (0, 2), (1, 2))]
    us = [tp_to_uni(r, 1) for r in rs if r]
    g = uni_gcd_all(us, p) if us else []
    pts = []
    for t in (uni_roots(g, p) if len(g) > 1 else []):
        cu = [tp_to_uni(tp_eval_partial(E, {1: t}, p), 0) for E in e]
        g0 = uni_gcd_all([x for x in cu if x], p)
        for s in (uni_roots(g0, p) if len(g0) > 1 else []):
            if all(tp_eval(E, {0: s, 1: t}, p) == 0 for E in e):
                pts.append([(a[j] + s * b[j] + t * c[j]) % p for j in range(3)])
    return dict(points=pts, gcd_degree=(len(g) - 1) if g else None, dd=dd_pattern(g, p) if len(g) > 1 else {}, zero_resultants=sum(1 for r in rs if not r))


def on_line_check(eqs, P, d, p):
    """each equation restricted to P + t d is identically zero?"""
    lin = [{(): P[j], (0,): d[j]} for j in range(3)]
    return [not affine_substitute(E, lin, p) for E in eqs]


def residual(eqs, P, d, e1, e2, p, charts=(2,)):
    """blow-up charts around the line {P + u d} in the coordinates c_1 = P + u d + v e1 + w e2.  Chart 2 (v = lam w) finds every
    point of the locus off the line with w != 0; an isolated point has w = 0 only for special (e1, e2), so random coordinates
    (repeated by the caller) see all of them.  Chart 1 (w = lam v) eliminates to the lam-axis and is degenerate whenever a point of
    the line is a singular point of all three surfaces (its whole exceptional fibre lies in the strict transform), so it is off by default.
    Returns dict(chart -> dict(result, off_line, n_solutions), orders = vanishing orders of the equations along the line)."""
    lin = [{(): P[j], (0,): d[j], (1,): e1[j], (2,): e2[j]} for j in range(3)]
    e = [affine_substitute(E, lin, p) for E in eqs]          # variables u = 0, v = 1, w = 2
    out = {}
    orders = [min(m.count(1) + m.count(2) for m in E) for E in e]
    for chart in charts:
        strict = []
        for E in e:
            k = min(m.count(1) + m.count(2) for m in E)          # vanishing order along the line
            assert k >= 1, "equation does not vanish on the line"
            S = {}
            for m, c in E.items():
                nv = m.count(1); nw = m.count(2)
                rest = tuple(x for x in m if x == 0)
                if chart == 1:   # w = lam v (lam = variable 2): u^a v^nv w^nw -> u^a v^{nv+nw-k} lam^nw
                    m2 = tuple(sorted(rest + (1,) * (nv + nw - k) + (2,) * nw))
                else:            # v = lam w (lam = variable 1): u^a v^nv w^nw -> u^a w^{nv+nw-k} lam^nv
                    m2 = tuple(sorted(rest + (2,) * (nv + nw - k) + (1,) * nv))
                S[m2] = (S.get(m2, 0) + c) % p
            strict.append({m: c for m, c in S.items() if c})
        res = solve_cubics3(strict, p)
        key = 1 if chart == 1 else 2
        off = [s for s in res['solutions'] if s[key] != 0]
        pts = []
        for s in off:
            u, v, lam = s[0], s[1], s[2]
            vv, ww = (v, lam * v % p) if chart == 1 else (s[2] * s[1] % p, s[2])
            if chart == 2: vv, ww = (s[1] * s[2]) % p, s[2]
            pts.append([(P[j] + u * d[j] + vv * e1[j] + ww * e2[j]) % p for j in range(3)])
        out[chart] = dict(result=res, off_line=pts, n_solutions=len(res['solutions']))
    out['orders'] = orders
    return out
