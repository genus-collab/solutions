# -*- coding: utf-8 -*-
"""fq_leadclass.py -- exact (Groebner-free) solving tools for the structured classes of the leading-order system
(genus2_symbol step 3 (ii)-(iii)), on t-polynomials {sorted exponent tuple: coeff} over F_p (fq_leadsolve.restrict_poly output).

  t-polynomial arithmetic: tp_add, tp_sub, tp_mul, tp_scale, tp_vars, tp_coeff_in(tp, j) -> {degree: t-poly without t_j}.
  Univariate: tp_to_uni(tp, j) (a t-poly in the single variable t_j -> coefficient list), uni_gcd_all(list) (fq_scalar_lib.p_gcd),
     uni_roots (p_roots), uni_fmt.
  eliminate_linear(eqs, j, pivot): for equations LINEAR in t_j, E_G c_piv - E_piv c_G (c = the t_j-coefficient t-poly) -- the
     resulting equations do not contain t_j; solutions of the original system are solutions of these together with E_piv = 0;
     conversely a common zero is genuine iff t_j can be solved back (checked by back_substitute).
  Polynomial matrices (entries = univariate coefficient lists): pm_det (Laplace, n <= 6), pm_minor; used for the Cramer
     solution of the mode systems over F_p[X_1] (fq_leadmodes) and the circulant class.
  Generic small fitter fit_monos(samples, monos, p): exact interpolation on an arbitrary monomial list (RREF), for the
     per-coloured-graph diagonal polynomials (61 monomials of weight <= 3 in x_0..x_4, y_0..y_4).
  Distinct-degree pattern dd_pattern(f, p): degrees of the irreducible factors of a squarefree univariate polynomial over F_p
     (gcd with x^{p^k} - x), to count roots over the algebraic closure when they are not in F_p.
"""
import itertools, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_scalar_lib import p_trim, p_add, p_sub, p_scale, p_mul, p_divmod, p_gcd, p_monic, p_eval, p_powmod, p_roots, p_fmt
from fq_mfit import fp_rref


# ------------------------------------------------------------------ t-polynomials
def tp_add(a, b, p):
    out = dict(a)
    for m, c in b.items():
        x = (out.get(m, 0) + c) % p
        if x: out[m] = x
        elif m in out: del out[m]
    return out
def tp_scale(a, c, p):
    c %= p
    return {m: v * c % p for m, v in a.items() if v * c % p}
def tp_sub(a, b, p): return tp_add(a, tp_scale(b, p - 1, p), p)
def tp_mul(a, b, p):
    out = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            m = tuple(sorted(m1 + m2)); x = (out.get(m, 0) + c1 * c2) % p
            if x: out[m] = x
            elif m in out: del out[m]
    return out
def tp_const(c, p): return {(): c % p} if c % p else {}
def tp_vars(tp): return sorted(set(v for m in tp for v in m))
def tp_degree(tp): return max((len(m) for m in tp), default=-1)
def tp_coeff_in(tp, j):
    """{e: t-poly} with tp = Sum_e t_j^e (coefficient poly without t_j)"""
    out = {}
    for m, c in tp.items():
        e = m.count(j); rest = tuple(v for v in m if v != j)
        out.setdefault(e, {})[rest] = (out.setdefault(e, {}).get(rest, 0) + c)
    return out
def tp_subs_zero(tp, vars_zero):
    return {m: c for m, c in tp.items() if not any(v in vars_zero for v in m)}
def tp_eval(tp, t, p):
    r = 0
    for m, c in tp.items():
        x = c
        for j in m: x = x * t[j] % p
        r = (r + x) % p
    return r


def tp_to_uni(tp, j):
    """t-poly in the single variable t_j -> coefficient list (low -> high)"""
    assert all(set(m) <= {j} for m in tp), "not univariate in t_%d: %s" % (j, tp_vars(tp))
    d = tp_degree(tp)
    out = [0] * (d + 1)
    for m, c in tp.items(): out[len(m)] = c
    return p_trim(out)


def eliminate_linear(eqs, j, pivot, p):
    """eqs: list of t-polys linear in t_j; returns (new equations without t_j for every non-pivot equation, pivot coefficient)"""
    cp = tp_coeff_in(eqs[pivot], j); c_piv = cp.get(1, {}); e_piv = cp.get(0, {})
    assert max(cp) <= 1 and c_piv, "pivot not linear in t_%d or coefficient zero" % j
    out = []
    for g, E in enumerate(eqs):
        if g == pivot: continue
        cg = tp_coeff_in(E, j); assert max(cg) <= 1, "equation %d not linear in t_%d" % (g, j)
        out.append(tp_sub(tp_mul(cg.get(0, {}), c_piv, p), tp_mul(e_piv, cg.get(1, {}), p), p))
    return out, c_piv, e_piv


def uni_gcd_all(polys, p):
    g = []
    for f in polys:
        f = p_trim(list(f))
        if not f: continue
        g = p_gcd(g, f, p) if g else p_monic(f, p)
        if len(g) == 1: return g
    return g


def uni_roots(f, p): return p_roots(list(f), p)
def uni_fmt(f, R, var="x"): return p_fmt(f, R, var)


def dd_pattern(f, p):
    """degrees of the irreducible factors of the squarefree part of f (distinct-degree factorisation): {degree: count}"""
    f = p_monic(list(f), p)
    df = p_trim([(k * c) % p for k, c in enumerate(f)][1:])
    f = p_divmod(f, p_gcd(f, df, p), p)[0] if df else f          # squarefree part
    out = {}; k = 1; xp = [0, 1]
    while len(f) > 1:
        if 2 * k > len(f) - 1: out[len(f) - 1] = out.get(len(f) - 1, 0) + 1; break
        xp = p_powmod(xp, p, f, p)                                  # x^{p^k} mod f
        g = p_gcd(p_sub(xp, [0, 1], p), f, p)
        if len(g) > 1:
            out[k] = (len(g) - 1) // k; f = p_divmod(f, g, p)[0]; xp = p_divmod(xp, f, p)[1] if len(f) > 1 else xp
        k += 1
    return out


# ------------------------------------------------------------------ polynomial matrices (univariate entries)
def pm_det(M, p):
    n = len(M)
    if n == 1: return p_trim(list(M[0][0]))
    if n == 2: return p_sub(p_mul(M[0][0], M[1][1], p), p_mul(M[0][1], M[1][0], p), p)
    out = []
    for j in range(n):
        if not p_trim(list(M[0][j])): continue
        sub = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        t = p_mul(M[0][j], pm_det(sub, p), p)
        out = p_add(out, t, p) if j % 2 == 0 else p_sub(out, t, p)
    return out


def pm_minor(M, rows, cols, p):
    return pm_det([[M[i][j] for j in cols] for i in rows], p)


def tp_matrix_to_uni(rows, j):
    """rows of t-polys in the single variable t_j -> rows of coefficient lists"""
    return [[tp_to_uni(e, j) for e in r] for r in rows]


# ------------------------------------------------------------------ generic small fitter
def fit_monos(samples, monos, p):
    """samples: [(t vector, {key: value})]; monos: list of exponent tuples; returns ({key: t-poly}, rank, n_inconsistent)"""
    keys = sorted(set(k for (_, v) in samples for k in v)); nm = len(monos)
    rows = []
    for t, vals in samples:
        row = []
        for m in monos:
            x = 1
            for j in m: x = x * t[j] % p
            row.append(x)
        rows.append(row + [vals[k] % p for k in keys])
    Rr, piv = fp_rref(rows, nm + len(keys), p)
    rank = sum(1 for c in piv if c < nm); nbad = sum(1 for c in piv if c >= nm)
    out = {}
    if rank == nm and nbad == 0:
        for kk, key in enumerate(keys):
            tp = {}
            for r, c in enumerate(piv):
                v = Rr[r][nm + kk]
                if v: tp[monos[c]] = v
            out[key] = tp
    return out, rank, nbad


def monos_weighted(nvars, weights, wmax):
    """all exponent tuples (sorted, with repetition) of total weight <= wmax"""
    out = [()]
    for d in range(1, wmax + 1):
        for m in itertools.combinations_with_replacement(range(nvars), d):
            if sum(weights[j] for j in m) <= wmax: out.append(m)
    return out
