# -*- coding: utf-8 -*-
"""fq_leadmodes.py -- exact solvers for the structured classes of the leading-order system (genus2_symbol step 3 (ii)-(iii)).

Resultants (Groebner-free elimination):
  sylvester(f, g, u, p): Sylvester matrix of two t-polys in the variable t_u (entries: t-polys in the other variables);
  tp_resultant(f, g, u, p): its determinant (Laplace, fq_leadclass.pm-style on t-polys) -- a t-poly without t_u.
  uni_resultant_eval(F, G, v, w, p): resultant in t_v of two t-polys in (t_v, t_w) by evaluation at w-points and interpolation
     (numeric F_p determinants; degree bound deg_v F deg_w G + deg_v G deg_w F).
  solve_cubics3(eqs, p, R): all F_p solutions of three polynomials in three variables (0, 1, 2) of degree <= 3, via
     Res_0 of the pairs, Res_1 of the resulting bivariate polynomials, gcd of the three eliminants (eliminant of t_2),
     F_p roots and back-substitution with verification; also the distinct-degree pattern of the eliminant's squarefree part
     (an upper bound on the number of solutions over the algebraic closure, spurious factors from vanishing leading coefficients
     included).  Returns dict(solutions, elim (gcd of eliminants), dd, identically_zero).
Diagonal class (fq_leadsolve.modes_diagonal; verified structure: every graph is  c + alpha X_1 + beta Y_3 + delta X_0X_2 +
  eps X_1^2 + phi X_3X_4 + gamma C(X),  C(X) = X_0^2X_3 + 2X_0X_1X_2 + X_0X_4^2 + X_1^3/3 + 2X_1X_3X_4 + X_2^2X_4 + X_2X_3^2):
  diagonal_solve(tps, R): the 6 edge equations are affine in w = (Y_3, P_02, P_34, C) with a CONSTANT coefficient matrix M and a
  right-hand side r(X_1) of degree 2; solutions exist iff u . r(X_1) = 0 for every left null vector u of M has a common root X_1
  (gcd of the quadratics).  Realisability of any (P_02, P_34, C) by (X_0, X_2, X_3, X_4) over the algebraic closure is shown in
  the report (the map (X_0,X_2,X_3,X_4) -> (X_0X_2, X_3X_4, Q) is dominant with 1-dimensional fibres).
Circulant class (fq_leadsolve.class_circulant): graphs 3, 6, 7 involve only c_1 = (c1_0, c1_1, c1_2) -> solve_cubics3; for each
  root the equations of graphs 2, 4, 5 are affine in (c2_1, c2_2, c3_0, c3_1, c3_2): fibre = affine space of dim 5 - rank.
Random points of W = {mu_1 : L_3 = L_6 = L_7 = 0} in A^15 (full system): restrict to a random affine 3-plane and solve_cubics3.
"""
import itertools, os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_leadclass import (tp_add, tp_sub, tp_mul, tp_scale, tp_vars, tp_degree, tp_coeff_in, tp_to_uni, tp_eval, tp_const,
                          uni_gcd_all, uni_roots, dd_pattern, fit_monos)
from fq_scalar_lib import p_trim, p_add, p_sub, p_mul, p_scale, p_gcd, p_monic, p_eval, p_roots
from fq_mfit import fp_rref, fp_nullspace, fp_solve


# ------------------------------------------------------------------ resultants on t-polys
def tp_det(M, p):
    n = len(M)
    if n == 1: return dict(M[0][0])
    out = {}
    for j in range(n):
        if not M[0][j]: continue
        sub = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        t = tp_mul(M[0][j], tp_det(sub, p), p)
        out = tp_add(out, t, p) if j % 2 == 0 else tp_sub(out, t, p)
    return out


def sylvester(f, g, u, p):
    cf = tp_coeff_in(f, u); cg = tp_coeff_in(g, u)
    df = max(cf); dg = max(cg); n = df + dg
    rows = []
    for i in range(dg):
        rows.append([cf.get(df - (c - i), {}) if 0 <= c - i <= df else {} for c in range(n)])
    for i in range(df):
        rows.append([cg.get(dg - (c - i), {}) if 0 <= c - i <= dg else {} for c in range(n)])
    return rows, df, dg


def tp_resultant(f, g, u, p):
    if not f or not g: return {}
    cf = tp_coeff_in(f, u); cg = tp_coeff_in(g, u)
    if max(cf) == 0 or max(cg) == 0: return {}
    rows, df, dg = sylvester(f, g, u, p)
    return tp_det(rows, p)


def det_fp(M, p):
    n = len(M); A = [list(r) for r in M]; d = 1
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] % p), None)
        if pr is None: return 0
        if pr != c: A[c], A[pr] = A[pr], A[c]; d = -d
        d = d * A[c][c] % p; inv = pow(A[c][c], p - 2, p)
        for i in range(c + 1, n):
            f = A[i][c] * inv % p
            if f: A[i] = [(x - f * y) % p for x, y in zip(A[i], A[c])]
    return d % p


def interpolate(xs, ys, p):
    """Newton interpolation -> coefficient list"""
    n = len(xs); coef = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) * pow((xs[i] - xs[i - j]) % p, p - 2, p) % p
    poly = []
    for i in range(n - 1, -1, -1):
        poly = p_add(p_mul(poly, [(-xs[i]) % p, 1], p), [coef[i]], p)
    return p_trim(poly)


def uni_resultant_eval(F, G, v, w, p):
    """Res_{t_v}(F, G) for t-polys in (t_v, t_w): coefficient list in t_w"""
    cF = tp_coeff_in(F, v); cG = tp_coeff_in(G, v)
    dv_F = max(cF); dv_G = max(cG)
    if dv_F == 0 or dv_G == 0: return []
    dw_F = max((tp_degree(c) for c in cF.values()), default=0); dw_G = max((tp_degree(c) for c in cG.values()), default=0)
    bound = dv_F * dw_G + dv_G * dw_F; N = bound + 1
    rows, df, dg = sylvester(F, G, v, p)
    xs = list(range(1, N + 1)); ys = []
    for x in xs:
        M = [[tp_eval(e, {w: x}, p) if e else 0 for e in r] for r in rows]
        ys.append(det_fp(M, p))
    return interpolate(xs, ys, p)


def solve_cubics3(eqs, p, R=None):
    """eqs: three t-polys in the variables 0, 1, 2.  Returns dict(solutions=[(t0,t1,t2)], elim, dd, identically_zero)"""
    f, g, h = eqs
    R01 = tp_resultant(f, g, 0, p); R02 = tp_resultant(f, h, 0, p); R12 = tp_resultant(g, h, 0, p)   # in t_1, t_2
    elims = []
    for (A, B) in ((R01, R02), (R01, R12), (R02, R12)):
        if not A or not B: elims.append([]); continue
        elims.append(uni_resultant_eval(A, B, 1, 2, p))
    nonzero = [e for e in elims if e]
    out = dict(identically_zero=(len(nonzero) < len(elims)), solutions=[], elim=[], dd={})
    if not nonzero: return out
    E = uni_gcd_all(nonzero, p); out['elim'] = E
    if len(E) <= 1: return out
    out['dd'] = dd_pattern(E, p)
    for t2 in uni_roots(E, p):
        # t_1 from the gcd of R01, R02, R12 at t2
        cands = [tp_to_uni(tp_eval_partial(Rr, {2: t2}, p), 1) for Rr in (R01, R02, R12) if Rr]
        cands = [c for c in cands if c]
        g1 = uni_gcd_all(cands, p) if cands else []
        for t1 in (uni_roots(g1, p) if len(g1) > 1 else []):
            cu = [tp_to_uni(tp_eval_partial(e, {1: t1, 2: t2}, p), 0) for e in eqs]
            g0 = uni_gcd_all([c for c in cu if c], p)
            for t0 in (uni_roots(g0, p) if len(g0) > 1 else []):
                if all(tp_eval(e, {0: t0, 1: t1, 2: t2}, p) == 0 for e in eqs): out['solutions'].append((t0, t1, t2))
    return out


def tp_eval_partial(tp, vals, p):
    """substitute some variables (dict var -> value); returns a t-poly in the remaining ones"""
    out = {}
    for m, c in tp.items():
        x = c; rest = []
        for j in m:
            if j in vals: x = x * vals[j] % p
            else: rest.append(j)
        if x:
            rest = tuple(sorted(rest)); out[rest] = (out.get(rest, 0) + x) % p
    return {m: c for m, c in out.items() if c}


# ------------------------------------------------------------------ diagonal class
CUBIC_C = {(0, 0, 3): 1, (0, 1, 2): 2, (0, 4, 4): 1, (1, 1, 1): None, (1, 3, 4): 2, (2, 2, 4): 1, (2, 3, 3): 1}   # X_1^3 coefficient 1/3


def diagonal_structure(tp, R):
    """decompose a diagonal-mode t-poly as c + alpha X_1 + beta Y_3 + delta X_0X_2 + eps X_1^2 + phi X_3X_4 + gamma C(X);
    returns (coeffs dict, ok) with ok = False if other monomials occur or the cubic part is not proportional to C"""
    p = R.p; inv3 = R.inv(3); co = dict(c=tp.get((), 0), alpha=tp.get((1,), 0), beta=tp.get((8,), 0), delta=tp.get((0, 2), 0),
                                        eps=tp.get((1, 1), 0), phi=tp.get((3, 4), 0))
    known = {(), (1,), (8,), (0, 2), (1, 1), (3, 4)} | set(CUBIC_C)
    ok = all(m in known for m in tp)
    gamma = tp.get((0, 0, 3), 0); co['gamma'] = gamma
    for m, k in CUBIC_C.items():
        want = gamma * (inv3 if k is None else k) % p
        if tp.get(m, 0) != want: ok = False
    return co, ok


def diagonal_solve(tps, R, edge=range(2, 8)):
    """the 6 edge equations as affine in w = (Y_3, P_02, P_34, C) with constant matrix M and rhs r(X_1) (degree 2):
    returns dict(M, rank, left_null, conditions (univariate polys in X_1), gcd, roots, structure_ok)"""
    p = R.p; M = []; rvec = []; ok = True
    for G in edge:
        co, good = diagonal_structure(tps[G], R); ok = ok and good
        M.append([co['beta'], co['delta'], co['phi'], co['gamma']])
        rvec.append(p_trim([co['c'], co['alpha'], co['eps']]))
    rank = len(fp_rref(M, 4, p)[0])
    MT = [[M[i][j] for i in range(len(M))] for j in range(4)]
    U = fp_nullspace(MT, len(M), p)
    conds = []
    for u in U:
        c = []
        for g, r in enumerate(rvec):
            if u[g]: c = p_add(c, p_scale(r, u[g], p), p)
        conds.append(c)
    gcd = uni_gcd_all([c for c in conds if c], p) if any(conds) else []
    roots = uni_roots(gcd, p) if len(gcd) > 1 else []
    return dict(M=M, rank=rank, left_null=U, conditions=conds, gcd=gcd, roots=roots, structure_ok=ok, rhs=rvec)


# ------------------------------------------------------------------ circulant class
def circulant_solve(tps, R, rng):
    """tps: {G: t-poly in c1_0, c1_1, c1_2, c2_1, c2_2, c3_0, c3_1, c3_2}.  Graphs 3, 6, 7 are pure c_1 (checked);
    returns dict(c1_only_ok, cubics result, fibres: [(c1 root, rank, consistent, dim)])"""
    p = R.p
    pure = [3, 6, 7]; ok = all(set(tp_vars(tps[G])) <= {0, 1, 2} for G in pure)
    res = solve_cubics3([tps[G] for G in pure], p, R)
    fibres = []
    for sol in res['solutions']:
        A = []; b = []
        for G in (2, 4, 5):
            tp = tp_eval_partial(tps[G], {0: sol[0], 1: sol[1], 2: sol[2]}, p)
            row = [tp.get((j,), 0) for j in range(3, 8)]; const = tp.get((), 0)
            assert all(len(m) <= 1 for m in tp), "not affine in (c2, c3) at fixed c1"
            A.append(row); b.append((-const) % p)
        x, N, rank = fp_solve(A, b, p)
        fibres.append((sol, rank, x is not None, (5 - rank) if x is not None else None))
    return dict(c1_only_ok=ok, cubics=res, fibres=fibres)


# ------------------------------------------------------------------ random points of W (full system)
def random_W_points(L367, p, rng, ntries=6):
    """restrict the three mu_1-only equations (t-polys in the 15 mu_1 variables, indices 0..14) to random affine 3-planes
    mu_1 = a + Sum_k t_k b_k and solve; returns list of mu_1 points (15-vectors) found"""
    from fq_leadsolve import restrict_poly
    pts = []
    for _ in range(ntries):
        a = [rng.randrange(p) for _ in range(15)]; bs = [[rng.randrange(p) for _ in range(15)] for _ in range(3)]
        A = [dict() for _ in range(40)]
        for v in range(15): A[v] = {k: bs[k][v] for k in range(3)}; A[v][3] = a[v]
        # affine: use a fourth variable t_3 = 1 -> substitute afterwards
        eqs = [tp_eval_partial(restrict_poly(L, A, p), {3: 1}, p) for L in L367]
        res = solve_cubics3(eqs, p)
        for (t0, t1, t2) in res['solutions']:
            pts.append([(a[v] + t0 * bs[0][v] + t1 * bs[1][v] + t2 * bs[2][v]) % p for v in range(15)])
        if pts: break
    return pts
