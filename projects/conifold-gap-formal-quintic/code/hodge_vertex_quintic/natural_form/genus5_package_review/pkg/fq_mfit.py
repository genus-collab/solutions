# -*- coding: utf-8 -*-
"""fq_mfit.py -- interpolation and linear algebra over F_p for genus2_symplectic.py (steps 4-5).

Variables (40): v = 0..14  m1_{ij} (i <= j, row-major), 15..24  m2_{ij} (i < j), 25..39  m3_{ij} (i <= j).
Monomials (161): 1; the 40 variables; m1 m1 (120, v <= w).  Weight: m1 -> 1, m2 -> 2, m3 -> 3.
  fit(samples, keys): samples = list of (param vector (40 F_p entries), {key: value}); least-squares-free exact fit by
  Gaussian elimination on the Vandermonde-type matrix (rank must be 161, extra rows must be consistent).
Galois action (sigma_m: zeta -> zeta^m, colour i -> m i): sigma_m(c_mu) = c_{pi_m(mu)}, pi_m relabelling the variables
  i -> m i (m2 picks up a sign when the order of the pair flips): conj_coeffs(poly, m).  With fq_scalar_lib.qzeta_coords this reconstructs coefficients in Q(zeta).
Linear algebra: fp_rref(rows, ncols, p) -> (R, pivots); fp_nullspace(rows, ncols, p) -> RREF basis of the kernel;
  fp_solve(A, b, p) -> (particular solution or None, nullspace basis); rank.
"""
import itertools, random
from fractions import Fraction as Fr

PAIRS1 = [(i, j) for i in range(5) for j in range(i, 5)]        # 15
PAIRS2 = [(i, j) for i in range(5) for j in range(i + 1, 5)]    # 10
NV = 40
VNAME = (["m1_%d%d" % ij for ij in PAIRS1] + ["m2_%d%d" % ij for ij in PAIRS2] + ["m3_%d%d" % ij for ij in PAIRS1])
VWEIGHT = [1] * 15 + [2] * 10 + [3] * 15
MONOS = [()] + [(v,) for v in range(NV)] + [(v, w) for v in range(15) for w in range(v, 15)]
MONO_W = [sum(VWEIGHT[v] for v in mu) for mu in MONOS]


def mono_name(mu):
    if not mu: return "1"
    return "*".join(VNAME[v] for v in mu)


def params_of(M):
    v1, v2, v3 = M.params()
    return list(v1) + list(v2) + list(v3)


def mono_eval(par, mu, p):
    r = 1
    for v in mu: r = r * par[v] % p
    return r


def var_perm(m):
    """(index map, sign) of the relabelling colour i -> m i on the 40 variables"""
    out = []
    for v in range(NV):
        if v < 15:
            i, j = PAIRS1[v]; a, b = sorted(((m * i) % 5, (m * j) % 5)); out.append((PAIRS1.index((a, b)), 1))
        elif v < 25:
            i, j = PAIRS2[v - 15]; a, b = (m * i) % 5, (m * j) % 5
            if a < b: out.append((15 + PAIRS2.index((a, b)), 1))
            else: out.append((15 + PAIRS2.index((b, a)), -1))
        else:
            i, j = PAIRS1[v - 25]; a, b = sorted(((m * i) % 5, (m * j) % 5)); out.append((25 + PAIRS1.index((a, b)), 1))
    return out


def conj_coeffs(poly, m, p):
    """poly: {mono index: coeff}; returns the coefficients of sigma_m(poly): sigma_m(c_mu) = sgn * c_{pi_m(mu)}
    (sigma_m(F(M)) = F(P_m M P_m^T) with (P_m M P_m^T)_{mi, mj} = M_{ij}; e.g. c_{m1_ii} = zeta^i c_{m1_00} for a covariant
    family, and sigma_m(zeta^i) = zeta^{mi} = c_{m1_{mi,mi}}/c_{00})."""
    vp = var_perm(m)
    out = {}
    for k, mu in enumerate(MONOS):
        img = [vp[v] for v in mu]
        mu2 = tuple(sorted(x[0] for x in img)); sgn = 1
        for x in img: sgn *= x[1]
        k2 = MONOS.index(mu2)
        c = poly.get(k2, 0)
        if c: out[k] = c * sgn % p
    return out


# ------------------------------------------------------------------ F_p linear algebra
def fp_rref(rows, ncols, p):
    A = [list(r) for r in rows]; piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(A)) if A[i][c] % p), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c], p - 2, p); A[r] = [x * inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] % p:
                f = A[i][c]; A[i] = [(x - f * y) % p for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == len(A): break
    return A[:r], piv


def fp_nullspace(rows, ncols, p):
    """RREF basis of {x: rows . x = 0}"""
    Rr, piv = fp_rref(rows, ncols, p)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        x = [0] * ncols; x[f] = 1
        for r, c in enumerate(piv): x[c] = (-Rr[r][f]) % p
        basis.append(x)
    return basis


def fp_solve(A, b, p):
    """A x = b (A list of rows, b list); returns (x particular or None if inconsistent, nullspace basis, rank)"""
    n = len(A[0]) if A else 0
    aug = [list(r) + [bb] for r, bb in zip(A, b)]
    Rr, piv = fp_rref(aug, n + 1, p)
    if n in piv: return None, fp_nullspace(A, n, p), len(piv) - 1
    x = [0] * n
    for r, c in enumerate(piv): x[c] = Rr[r][n]
    return x, fp_nullspace(A, n, p), len(piv)


def fit(samples, keys, p):
    """samples: [(par, {key: val})]; returns ({key: {mono idx: coeff}}, rank, n_inconsistent_rows)"""
    nm = len(MONOS)
    rows = []
    for par, vals in samples:
        rows.append([mono_eval(par, mu, p) for mu in MONOS] + [vals[k] for k in keys])
    Rr, piv = fp_rref(rows, nm + len(keys), p)
    rank = sum(1 for c in piv if c < nm)
    bad = [c for c in piv if c >= nm]           # a pivot in a value column = inconsistent system
    out = {}
    if rank == nm and not bad:
        for kk, key in enumerate(keys):
            poly = {}
            for r, c in enumerate(piv):
                v = Rr[r][nm + kk]
                if v: poly[c] = v
            out[key] = poly
    return out, rank, len(bad)


def poly_eval(poly, par, p):
    return sum(c * mono_eval(par, MONOS[k], p) for k, c in poly.items()) % p


def poly_weights(poly):
    return sorted(set(MONO_W[k] for k in poly))


def linear_part(poly, var_range):
    """coefficients of the linear monomials in the variables var_range (list), as a vector over var_range"""
    out = []
    for v in var_range:
        k = MONOS.index((v,)); out.append(poly.get(k, 0))
    return out


def quad_part(poly, p):
    """15x15 symmetric matrix Qm with  Sum_{v<=w} c_vw m1_v m1_w = m1^T Qm m1"""
    Q = [[0] * 15 for _ in range(15)]
    inv2 = pow(2, p - 2, p)
    for k, c in poly.items():
        mu = MONOS[k]
        if len(mu) == 2:
            v, w = mu
            if v == w: Q[v][v] = (Q[v][v] + c) % p
            else: Q[v][w] = (Q[v][w] + c * inv2) % p; Q[w][v] = (Q[w][v] + c * inv2) % p
    return Q


def quad_restrict(Q, basis, p):
    """the quadratic form on the span of basis (list of vectors): matrix B[a][b] = k_a^T Q k_b"""
    r = len(basis)
    out = [[0] * r for _ in range(r)]
    for a in range(r):
        Qk = [sum(Q[i][j] * basis[a][j] for j in range(15)) % p for i in range(15)]
        for b in range(r): out[a][b] = sum(basis[b][i] * Qk[i] for i in range(15)) % p
    return out


def vec_fmt(R, x, names=None):
    parts = []
    for i, c in enumerate(x):
        if c: parts.append("%s %s" % (R.fmt(c), names[i] if names else "e%d" % i))
    return " + ".join(parts) if parts else "0"
