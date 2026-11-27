# -*- coding: utf-8 -*-
"""fq_msolve.py -- the solve of genus2_symplectic.py (step 5) on the fitted polynomials of fq_mfit.

Equations: for the graphs in `graphs` and the orders in `orders`,  P_{Gamma,k}(m) = target_{Gamma,k}, target_{Gamma,-5} = -(7/5) c_Gamma
(c_Gamma = the s^-10 coefficient), all other targets 0.  Structure (fq_mfit weights): P_{Gamma,-10}, P_{Gamma,-5} constant;
P_{Gamma,-6}, P_{Gamma,-1} linear homogeneous in m1; P_{Gamma,-2} = m1^T Q m1 + b . m2; all other orders identically zero.
solve(polys, graphs, orders, c10, R) returns a dict with
  obstructed: [(Gamma, k, residual)] for the constant equations whose value differs from the target (each is independently
              unsatisfiable: the polynomial is constant),
  lin_rank, K (RREF basis of the m1 kernel of the linear equations), dimK,
  B (the m2 coefficient rows of the s^-2 equations), rankB, N2 (RREF basis of the m2 kernel), U (left null vectors of B),
  quad_ok (True iff u . (Q^Gamma restricted to K) = 0 for every left-null vector u: then every m1 in K extends to a solution),
  quad_obstruction (the restricted quadratic forms u . Q_K that are non-zero), dim_solution (= dimK + 10 - rankB when consistent).
Galois (Q(zeta)) reconstruction of an RREF basis: conj_basis(basis, m) = RREF of the relabelled vectors (sigma_m(K) =
{(y_{pi_m mu})_mu : y in K}), then fq_scalar_lib.qzeta_coords entrywise (basis_qzeta).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import (MONOS, MONO_W, PAIRS1, PAIRS2, VNAME, var_perm, fp_rref, fp_nullspace, fp_solve, linear_part, quad_part,
                     quad_restrict, poly_weights)
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt

M1 = list(range(15)); M2 = list(range(15, 25)); M3 = list(range(25, 40))


def const_part(poly):
    return poly.get(0, 0)


def solve(polys, graphs, orders, c10, R):
    p = R.p
    out = dict(graphs=list(graphs), orders=list(orders))
    obstructed = []; lin_rows = []; lin_rhs = []; quad = []
    for G in graphs:
        for k in orders:
            poly = polys[(G, k)]
            target = (-R.of(Fr(7, 5)) * c10[G]) % p if k == -5 else 0
            ws = poly_weights(poly)
            # weight-3 monomials must be absent (m3 never polar)
            assert all(w <= 2 for w in ws), (G, k, ws)
            if all(w == 0 for w in ws):
                r = (const_part(poly) - target) % p
                if r: obstructed.append((G, k, r))
                continue
            if max(ws) == 1:
                assert all(MONOS[m] and MONOS[m][0] < 15 for m in poly if m), (G, k)   # linear in m1 only
                lin_rows.append(linear_part(poly, M1)); lin_rhs.append((target - const_part(poly)) % p)
                continue
            # weight 2: m1 m1 and m2
            assert const_part(poly) == 0 and all(MONOS[m][0] < 25 for m in poly if m) and k == -2, (G, k)
            quad.append((G, quad_part(poly, p), linear_part(poly, M2), target))
    out['obstructed'] = obstructed
    # linear system in m1
    if lin_rows:
        xpart, K, rank = fp_solve(lin_rows, lin_rhs, p)
        out['lin_rank'] = rank; out['lin_consistent'] = xpart is not None
        out['lin_particular'] = xpart
    else:
        K = fp_nullspace([[0] * 15], 15, p); out['lin_rank'] = 0; out['lin_consistent'] = True; out['lin_particular'] = [0] * 15
    out['K'] = K; out['dimK'] = len(K)
    # s^-2 equations: B m2 = target - m1^T Q m1 on m1 = xpart + K t
    if quad:
        B = [q[2] for q in quad]
        Rr, piv = fp_rref(B, 10, p); rankB = len(piv)
        out['B'] = B; out['rankB'] = rankB; out['N2'] = fp_nullspace(B, 10, p)
        # left null vectors of B: u with u^T B = 0  <=>  B^T u = 0
        BT = [[B[i][j] for i in range(len(B))] for j in range(10)]
        U = fp_nullspace(BT, len(B), p); out['U'] = U
        # restricted quadratic forms (homogeneous part only is tested when the particular m1 solution is zero)
        QK = [quad_restrict(q[1], K, p) for q in quad] if K else [[] for q in quad]
        out['QK'] = QK
        obstr = []
        for u in U:
            tot = [[sum(u[g] * QK[g][a][b] for g in range(len(quad))) % p for b in range(len(K))] for a in range(len(K))]
            if any(x for row in tot for x in row): obstr.append((u, tot))
        # the inhomogeneous parts (targets - const) are zero for k = -2; the particular m1 solution is zero (homogeneous)
        out['quad_ok'] = not obstr and all(q[3] == 0 for q in quad) and not any(out['lin_particular'])
        out['quad_obstruction'] = obstr
        out['dim_solution'] = (len(K) + 10 - rankB) if out['quad_ok'] else None
    else:
        out['rankB'] = None; out['quad_ok'] = True; out['dim_solution'] = len(K) + 10
    return out


# ------------------------------------------------------------------ Galois reconstruction of RREF bases
def relabel_vec(y, m, var_range):
    """z_mu = y_{pi_m mu} on the variables var_range (indices into the 40-variable list)"""
    vp = var_perm(m)
    z = [0] * len(y)
    for a, v in enumerate(var_range):
        w, sgn = vp[v]
        z[a] = y[var_range.index(w)] * sgn
    return z


def conj_basis(basis, m, var_range, p):
    if not basis: return []
    rel = [[x % p for x in relabel_vec(y, m, var_range)] for y in basis]
    Rr, piv = fp_rref(rel, len(var_range), p)
    return Rr


def basis_qzeta(R, basis, var_range):
    """[(pivot structure ok, [Q(zeta) coordinate tuples per entry])] -- the RREF basis of sigma_m(K) for m = 1..4 must share
    the pivot columns; returns (ok, rows of qzeta coordinate tuples) or (False, None)"""
    p = R.p
    conj = [conj_basis(basis, m, var_range, p) for m in range(1, 5)]
    if any(len(c) != len(basis) for c in conj): return False, None
    def pivots(B): return [next((a for a, x in enumerate(row) if x), None) for row in B]
    if any(pivots(c) != pivots(conj[0]) for c in conj): return False, None
    rows = []
    for r in range(len(basis)):
        row = []
        for a in range(len(var_range)):
            xs = [conj[m - 1][r][a] for m in range(1, 5)]
            row.append(qzeta_coords(R, xs))
        rows.append(row)
    return True, rows


def fmt_basis_qzeta(R, rows, var_range, names=VNAME):
    lines = []
    for r, row in enumerate(rows):
        parts = []
        for a, cs in enumerate(row):
            if any(cs):
                parts.append("[%s] %s" % (qzeta_fmt(R, cs), names[var_range[a]]))
        lines.append("      k_%d = " % (r + 1) + (" + ".join(parts) if parts else "0"))
    return lines


def poly_qzeta(R, polys, key, var_filter=None):
    """{mono: Q(zeta) coords} of P_key via the Galois conjugates sigma_m(c_mu) = c_{pi_m mu}"""
    from fq_mfit import conj_coeffs
    p = R.p; poly = polys[key]
    conjs = [conj_coeffs(poly, m, p) for m in range(1, 5)]
    out = {}
    for k, c in poly.items():
        if var_filter is not None and not var_filter(MONOS[k]): continue
        out[k] = qzeta_coords(R, [conjs[m - 1].get(k, 0) for m in range(1, 5)])
    return out
