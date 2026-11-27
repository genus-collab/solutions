# -*- coding: utf-8 -*-
"""fq_leadcol.py -- the per-COLOURED-graph version of the leading-order system for the diagonal classes (genus2_symbol step 3 (ii)).

For diagonal M_0 = diag(exp(x_i zhat + y_i zhat^3)) (mu_2 diagonal antisymmetric = 0) colours are preserved and the coloured
graph terms are defined (fq_msum colmask on the FROZEN theory, fq_leadsolve.FrozenRows): the leading coefficient of the coloured
graph (Gamma, c) is a polynomial of weight <= 3 in (x_0..x_4, y_0..y_4) (weights 1, 3), fitted on the 61 monomials of weight <= 3
from random diagonal frozen evaluations (fq_leadclass.fit_monos).  Structure expected from DIAG_SCAN/SCALAR_SCAN s.5: monochrome
graphs at colour c: T0 + T1 x_c + T2 x_c^2 + T3 x_c^3 + T4 y_c with (T0, T1, T2, T3, T4) = (T0_-10, T1_-6, T2_-2, T3_2, T4_2) of the
SCALAR_SCAN table (the leading s-coefficients, since the frozen theory keeps exactly the leading coefficients); mixed colourings of
the two-vertex graphs 4, 5, 7: 1, x_c, x_c', x_c x_c'; mixed theta colourings: constants.
Systems: per colour c the 7 monochrome equations in (x_c, y_c) (eliminate y_c with graph 1, gcd of the 6 cubics in x_c; also the edge
graphs only: eliminate with graph 2, gcd of 5); the 20 mixed theta colourings (constant equations, inconsistent iff non-zero);
the covariant restriction x_i = zeta^{-i} a_1, y_i = zeta^{-3i} a_3 of all 115 coloured polynomials (SCALAR_SCAN s.5: 23 orbit
representatives, gcd 1 of the 22 eliminated polynomials; reproduced from the orbit representatives, all 110 edge equations).
CAVEAT (fq_msum, not a bug for its M = 1 use): colmask cannot restrict the leafless genus-2 vertex tensor, so the coloured vertex graph
at colour c carries the other colours' leafless terms zeta^{3j}(x_j^3/10 + y_j/48) for diagonal M != 1; coloured_polys removes them
(and checks that the removed part is exactly that) -- see the comment in the code.
API: coloured_polys(fr, graphs, R, rng, nsamp=75) -> ({(idx, col): t-poly}, rank, nbad, drop_ok); monochrome_table(cp, R); solve_percolour(cp, R);
     theta_constants(cp); covariant_percoloured(cp, R).
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_sser import SympM
from fq_msum import GraphSummer
from fq_leadclass import (fit_monos, monos_weighted, tp_vars, tp_coeff_in, tp_to_uni, eliminate_linear, uni_gcd_all, uni_roots,
                          tp_eval, tp_add, tp_mul, tp_scale)

NV = 10; WEIGHTS = [1] * 5 + [3] * 5
MONOS = monos_weighted(NV, WEIGHTS, 3)      # 61
assert len(MONOS) == 61
SCALAR_TABLE = {1: (Fr(1807, 15360000), Fr(37, 51200), Fr(-119, 7680), Fr(-131, 3456), Fr(23, 1152)),
                2: (Fr(-811, 1280000), Fr(-343, 64000), Fr(89, 1920), Fr(-1, 288), Fr(-1, 48)),
                3: (Fr(-57, 320000), Fr(87, 16000), Fr(-7, 160), Fr(1, 8), Fr(0)),
                4: (Fr(-31, 15360000), Fr(-73, 256000), Fr(71, 7680), Fr(-397, 3456), Fr(1, 1152)),
                5: (Fr(-17, 1280000), Fr(67, 64000), Fr(-53, 1920), Fr(23, 96), Fr(0)),
                6: (Fr(1, 96000), Fr(-1, 1600), Fr(1, 80), Fr(-1, 12), Fr(0)),
                7: (Fr(1, 64000), Fr(-3, 3200), Fr(3, 160), Fr(-1, 8), Fr(0))}


def colourings(G):
    return list(itertools.product(range(5), repeat=G['V']))


def coloured_polys(fr, graphs, R, rng, nsamp=75):
    p = R.p; samples = []
    for t in range(nsamp):
        x = [rng.randrange(p) for _ in range(5)]; y = [rng.randrange(p) for _ in range(5)]
        M = SympM.diagonal(R, x, y); GS = GraphSummer(fr, M); vals = {}
        for G in graphs:
            for col in colourings(G):
                s = GS.graph(G, colmask=col); vals[(G['idx'], col)] = s.coef(0) if s.prec > 0 else 0
        samples.append((x + y, vals))
    polys, rank, nbad = fit_monos(samples, MONOS, p)
    # fq_msum's colmask cannot restrict the LEAFLESS genus-2 vertex tensor (n = 0: nothing to mask), so the 'coloured' vertex
    # graph at colour c also contains the other colours' leafless terms Sum_{j != c} zeta^{3j} (x_j^3/10 + y_j/48)
    # (Delta^_j^{(2g-2)/2} I(2; empty), I(2; empty) = a_1^3/10 + a_2/48, SCALAR_SCAN s.1).  At M = 1 these vanish (phase sum), which is
    # why genus2_symplectic's per-colour check passed.  Here they are removed, after checking that they are exactly those terms.
    inv10 = R.inv(10); inv48 = R.inv(48); drop_ok = True
    for c in range(5):
        tp = polys[(1, (c,))]; own = {m: v for m, v in tp.items() if set(m) <= {c, 5 + c}}
        foreign = {m: v for m, v in tp.items() if not set(m) <= {c, 5 + c}}
        exp = {}
        for j in range(5):
            if j == c: continue
            exp[(j, j, j)] = R.z(3 * j) * inv10 % p; exp[(5 + j,)] = R.z(3 * j) * inv48 % p
        if foreign != exp: drop_ok = False
        polys[(1, (c,))] = own
    return polys, rank, nbad, drop_ok


def monochrome_table(cp, R):
    """{(idx, c): (T0..T4)} from the monochrome polynomials; also whether each involves only (x_c, y_c) and matches SCALAR_TABLE"""
    p = R.p; out = {}; ok_vars = True; ok_table = True
    for (idx, col), tp in cp.items():
        if len(set(col)) != 1: continue
        c = col[0]
        if not set(tp_vars(tp)) <= {c, 5 + c}: ok_vars = False
        T = (tp.get((), 0), tp.get((c,), 0), tp.get((c, c), 0), tp.get((c, c, c), 0), tp.get((5 + c,), 0))
        out[(idx, c)] = T
        if c == 0 and tuple(R.of(x) for x in SCALAR_TABLE[idx]) != T: ok_table = False
    return out, ok_vars, ok_table


def solve_percolour(cp, R):
    """per colour c: eliminate y_c (pivot graph 1, and separately pivot graph 2 for the edge-only system); gcds of the univariate
    polynomials in x_c; returns {c: dict(gcd_all, gcd_edge, roots_all, roots_edge)}"""
    p = R.p; out = {}
    for c in range(5):
        eqs = [cp[(idx, (c,) * (1 if idx in (1, 2, 3) else 2))] for idx in range(1, 8)]
        red, cpiv, epiv = eliminate_linear(eqs, 5 + c, 0, p)
        unis = [tp_to_uni(e, c) for e in red if e]
        g_all = uni_gcd_all(unis, p) if unis else []
        red_e, _, _ = eliminate_linear(eqs[1:], 5 + c, 0, p)
        unis_e = [tp_to_uni(e, c) for e in red_e if e]
        g_edge = uni_gcd_all(unis_e, p) if unis_e else []
        out[c] = dict(gcd_all=g_all, gcd_edge=g_edge, roots_all=uni_roots(g_all, p) if len(g_all) > 1 else [],
                      roots_edge=uni_roots(g_edge, p) if len(g_edge) > 1 else [], n_eqs=(len(unis), len(unis_e)))
    return out


def theta_constants(cp, R):
    """the mixed theta colourings: {col: constant} and whether every one is a non-zero constant (no variable dependence)"""
    out = {}; ok = True
    for (idx, col), tp in cp.items():
        if idx != 6 or len(set(col)) == 1: continue
        if tp_vars(tp): ok = False
        out[col] = tp.get((), 0)
        if not tp.get((), 0): ok = False
    return out, ok


def mixed_structure(cp, R):
    """monomials occurring in the mixed colourings of graphs 4, 5, 7 (expected subsets of {1, x_c, x_c', x_c x_c'})"""
    bad = []
    for (idx, col), tp in cp.items():
        if idx in (4, 5, 7) and len(set(col)) == 2:
            c, d = col
            allowed = {(), (c,), (d,), tuple(sorted((c, d)))}
            if not set(tp) <= allowed: bad.append((idx, col, sorted(set(tp) - allowed)))
    return bad


def covariant_percoloured(cp, R, edge_only=True):
    """restrict all coloured polynomials to x_i = zeta^{-i} a_1, y_i = zeta^{-3i} a_3 (t = (a_1, a_3)); eliminate a_3 with the
    vertex graph at colour 0 (or, edge_only, with graph 2 at colour 0); gcd of the univariate polynomials in a_1"""
    p = R.p; sub = [{0: R.z(-i)} for i in range(5)] + [{1: R.z(-3 * i)} for i in range(5)]
    def restrict(tp):
        out = {}
        for m, c in tp.items():
            terms = {(): c}
            for j in m:
                new = {}
                for mm, cc in terms.items():
                    for k, coef in sub[j].items():
                        m2 = tuple(sorted(mm + (k,))); new[m2] = (new.get(m2, 0) + cc * coef) % p
                terms = new
            for mm, cc in terms.items(): out[mm] = (out.get(mm, 0) + cc) % p
        return {m: c for m, c in out.items() if c}
    keys = sorted(cp); eqs = [restrict(cp[k]) for k in keys]
    piv = keys.index((2, (0,)) if edge_only else (1, (0,)))
    sel = [i for i, k in enumerate(keys) if (k[0] >= 2 or not edge_only)]
    eqs_sel = [eqs[i] for i in sel]; piv_sel = sel.index(piv)
    red, cpiv, epiv = eliminate_linear(eqs_sel, 1, piv_sel, p)
    unis = [tp_to_uni(e, 0) for e in red if e]
    g = uni_gcd_all(unis, p) if unis else []
    # orbit representatives count: distinct polynomials
    distinct = len(set(tuple(sorted(e.items())) for e in eqs_sel))
    return dict(gcd=g, roots=uni_roots(g, p) if len(g) > 1 else [], n_eqs=len(unis), n_distinct=distinct, pivot=keys[piv])
