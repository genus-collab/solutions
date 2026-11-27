# -*- coding: utf-8 -*-
"""fq_leadfit.py -- interpolation of the colour-summed genus-2 graph s-series as polynomials in the constant symplectic M over
ALL monomials of weight <= 3 (genus2_symbol programme, step 2).  Extends fq_mfit (weight <= 2 plus m3) by the weight-3
monomials m1 m2 (150) and m1^3 (680), which matter for the leading-order (q-dependent, M = M_0(zhat)) system.

Variables (fq_mfit numbering): 0..14 m1_{ij} (i <= j), 15..24 m2_{ij} (i < j), 25..39 m3_{ij} (i <= j); weights 1, 2, 3.
MONOS3 (991): 1; 40 variables; m1 m1 (120); m1 m2 (150); m1 m1 m1 (680).  MONO3_W = weight.
fit3(samples, keys, p): exact fit by forward Gaussian elimination + back substitution on the Vandermonde-type matrix
   [mono_mu(par_t)]_{t, mu} | [val_t(key)]_{t, key}; rank must be 991 and the surplus rows consistent (zero LHS -> zero RHS).
   Returns ({key: {mono idx: coeff}}, rank, n_inconsistent).
Galois (sigma_m: zeta -> zeta^m = colour relabelling i -> m i, with the sign of m2 when the pair flips): conj_coeffs3(poly, m, p);
   with fq_scalar_lib.qzeta_coords this reconstructs every coefficient in Q(zeta) (poly_qzeta3).
Colour shift i -> i+1 (z -> zeta^{-1} z): c_{shift mu} = zeta^w c_mu for a weight-w monomial (shift_check).
Weight bookkeeping (SYMPLECTIC_SCAN s.2): the weight-w part of every colour-summed graph lies in s^{-10+4w} Q(zeta)[[u]], so the
   coefficient of s^k can only contain monomials of weight w with -w = k mod 5: s^-10, -5, 0: weight 0; s^-6, -1: 1; s^-2: 2;
   s^+2: 3; s^-9..-7, -4, -3, 1: nothing.  expected_weights(k) encodes this; weights_of(poly) reads the fitted content.
"""
import itertools, time
from fractions import Fraction as Fr
from fq_mfit import PAIRS1, PAIRS2, NV, VNAME, VWEIGHT, params_of, var_perm
from fq_scalar_lib import qzeta_coords

M1V = list(range(15)); M2V = list(range(15, 25)); M3V = list(range(25, 40))
MONOS3 = ([()] + [(v,) for v in range(NV)]
          + [(v, w) for v in M1V for w in M1V if v <= w]
          + [(v, w) for v in M1V for w in M2V]
          + [(u, v, w) for u in M1V for v in M1V for w in M1V if u <= v <= w])
MONO3_W = [sum(VWEIGHT[v] for v in mu) for mu in MONOS3]
MONO3_INDEX = {mu: k for k, mu in enumerate(MONOS3)}
assert len(MONOS3) == 991 and len(MONO3_INDEX) == 991


def mono_name3(mu):
    return "*".join(VNAME[v] for v in mu) if mu else "1"


def mono_eval3(par, mu, p):
    r = 1
    for v in mu: r = r * par[v] % p
    return r


def expected_weights(k):
    """the weights w <= 3 with -w = k mod 5 (the only ones that can occur at s^k)"""
    return [w for w in range(4) if (k + w) % 5 == 0]


def weights_of(poly):
    return sorted(set(MONO3_W[k] for k in poly))


def weight_part(poly, w):
    return {k: c for k, c in poly.items() if MONO3_W[k] == w}


def poly_eval3(poly, par, p):
    return sum(c * mono_eval3(par, MONOS3[k], p) for k, c in poly.items()) % p


# ------------------------------------------------------------------ the fit
def fit3(samples, keys, p):
    """samples: [(par (40 F_p entries), {key: value})]; forward elimination to row echelon form, then back substitution.
    Returns ({key: {mono idx: coeff}}, rank, n_inconsistent_rows)."""
    nm = len(MONOS3); nk = len(keys)
    rows = [[mono_eval3(par, mu, p) for mu in MONOS3] + [vals[k] % p for k in keys] for (par, vals) in samples]
    nr = len(rows); piv_cols = []; r = 0
    for c in range(nm):
        pr = next((i for i in range(r, nr) if rows[i][c]), None)
        if pr is None: continue
        rows[r], rows[pr] = rows[pr], rows[r]
        inv = pow(rows[r][c], p - 2, p); rows[r] = [x * inv % p for x in rows[r]]
        pivot = rows[r]
        for i in range(r + 1, nr):
            f = rows[i][c]
            if f: rows[i] = [(x - f * y) % p for x, y in zip(rows[i], pivot)]
        piv_cols.append(c); r += 1
        if r == nr: break
    rank = len(piv_cols)
    # surplus rows (r..nr-1) have zero LHS now; they must have zero RHS
    nbad = sum(1 for i in range(rank, nr) if any(rows[i][nm:]))
    out = {}
    if rank == nm and nbad == 0:
        # back substitution: rows[t] has pivot 1 at column piv_cols[t] = t (full rank => piv_cols = 0..nm-1)
        sol = [[0] * nk for _ in range(nm)]
        for t in range(nm - 1, -1, -1):
            row = rows[t]; rhs = row[nm:]
            acc = list(rhs)
            for c in range(t + 1, nm):
                f = row[c]
                if f:
                    sc = sol[c]
                    acc = [(a - f * b) % p for a, b in zip(acc, sc)]
            sol[t] = acc
        for kk, key in enumerate(keys):
            poly = {}
            for mu_i in range(nm):
                v = sol[mu_i][kk]
                if v: poly[mu_i] = v
            out[key] = poly
    return out, rank, nbad


# ------------------------------------------------------------------ Galois and colour-shift structure
def conj_coeffs3(poly, m, p):
    """coefficients of sigma_m(poly): sigma_m(c_mu) = sgn * c_{pi_m(mu)} (pi_m: relabelling i -> m i on the variables)"""
    vp = var_perm(m); out = {}
    for k, mu in enumerate(MONOS3):
        img = [vp[v] for v in mu]
        mu2 = tuple(sorted(x[0] for x in img)); sgn = 1
        for x in img: sgn *= x[1]
        k2 = MONO3_INDEX.get(mu2)
        if k2 is None: continue
        c = poly.get(k2, 0)
        if c: out[k] = c * sgn % p
    return out


def poly_qzeta3(R, poly, var_filter=None):
    """{mono idx: Q(zeta) coords (c0..c3)} of the fitted polynomial via its Galois conjugates"""
    p = R.p
    conjs = [conj_coeffs3(poly, m, p) for m in range(1, 5)]
    out = {}
    for k in poly:
        if var_filter is not None and not var_filter(MONOS3[k]): continue
        out[k] = qzeta_coords(R, [conjs[m - 1].get(k, 0) for m in range(1, 5)])
    return out


def shift_var(v):
    """colour shift i -> i+1 on a variable index (sign for m2 when the pair wraps)"""
    if v < 15:
        i, j = PAIRS1[v]; return PAIRS1.index(tuple(sorted(((i + 1) % 5, (j + 1) % 5)))), 1
    if v < 25:
        i, j = PAIRS2[v - 15]; a, b = (i + 1) % 5, (j + 1) % 5
        return (15 + PAIRS2.index((a, b)), 1) if a < b else (15 + PAIRS2.index((b, a)), -1)
    i, j = PAIRS1[v - 25]; return 25 + PAIRS1.index(tuple(sorted(((i + 1) % 5, (j + 1) % 5)))), 1


def shift_check(poly, R):
    """c_{shift mu} = zeta^w sgn c_mu for every monomial: returns the number of violations"""
    p = R.p; bad = 0
    for k, c in poly.items():
        mu = MONOS3[k]; img = [shift_var(v) for v in mu]; sgn = 1
        for x in img: sgn *= x[1]
        k2 = MONO3_INDEX[tuple(sorted(x[0] for x in img))]
        if poly.get(k2, 0) != c * R.z(MONO3_W[k]) * sgn % p: bad += 1
    return bad


def collect_samples(S, graphs, n, rng, orders, SympM, GraphSummer, log=None):
    """n random symplectic M -> [(par, {(idx, k): coef})]"""
    out = []; t0 = time.time()
    for t in range(n):
        Mt = SympM.random(S.R, rng); r = GraphSummer(S, Mt).all_graphs(graphs)
        out.append((params_of(Mt), {(idx, k): r[idx].coef(k) for idx in r for k in orders}))
        if log and (t + 1) % 200 == 0: log("      %d evaluations, %.1f s" % (t + 1, time.time() - t0))
    return out
