# -*- coding: utf-8 -*-
"""fq_mfourier.py -- Z/5 Fourier structure of the fitted polar polynomials (genus2_symplectic, step 5 analysis).

Colour shift i -> i+1 (lambda -> zeta lambda) with z -> zeta^{-1} z is a symmetry of the colour-summed graph terms, acting on
the constant matrices by m_k -> zeta^{-k} P m_k P^T.  Hence a weight-w coefficient polynomial satisfies
c_{shift^s(mu)} = zeta^{s w} c_mu, and in the Fourier modes
   L_{d,a}(m1) := Sum_i zeta^{a i} m1_{i,i+d}   (d = 0,1,2; indices mod 5, m1 symmetric),
   B_{d,a}(m2) := Sum_i zeta^{a i} m2_{i,i+d}   (d = 1,2; m2_{i,i+d} := -m2_{i+d,i} when i+d < i mod 5)
the weight-1 coefficients (s^-6, s^-1) are combinations of L_{d,1} only, the m2 part of the s^-2 coefficient of B_{d,2} only,
and its m1 m1 part involves only products L_{d,a} L_{d',a'} with a + a' = 2 mod 5.  (All three statements are CHECKED on the
fitted polynomials, not assumed.)
Consequences for the subleading system {s^-6 = s^-1 = s^-2 = 0 for all graphs}:
   K = {m1 : L_{0,1} = L_{1,1} = L_{2,1} = 0} (dim 12) is the solution set of the linear part (if the 14 forms span all three modes);
   on K the quadratic part is  Q^Gamma = x^T C^Gamma y + z^T E^Gamma w  with x_d = L_{d,0}, y_d = L_{d,2}, z_d = L_{d,3}, w_d = L_{d,4}
   (3x3 matrices C, E over F_p), and the s^-2 equations read  beta^Gamma . (B_{1,2}, B_{2,2})(m2) + Q^Gamma(m1) = 0.
API: modes(R) -> (F1, F2, MODES1, MODES2); lin_in_modes(coeffs, F); quad_in_modes(Q, F1); structure(polys, R) -> dict with the
mode decompositions, class checks, the matrices C, E, beta, kappa (if the m2 forms are proportional) and jacobian_dims (local
dimensions of the quadric locus at random points of the four 6-dim linear subspaces {x or y = 0} x {z or w = 0}).
"""
import os, sys, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import PAIRS1, PAIRS2, MONOS, linear_part, quad_part, fp_rref, fp_nullspace

M1 = list(range(15)); M2 = list(range(15, 25))
MODES1 = [(d, a) for d in range(3) for a in range(5)]
MODES2 = [(d, a) for d in (1, 2) for a in range(5)]


def modes(R):
    """F1 (15x15): m1_var = Sum_mode F1[var][mode] L_mode (inverse Fourier); F2 (10x10) likewise for m2"""
    p = R.p; inv5 = R.inv(5)
    F1 = [[0] * 15 for _ in range(15)]
    for mi, (d, a) in enumerate(MODES1):
        for i in range(5):
            v = PAIRS1.index(tuple(sorted((i, (i + d) % 5))))
            F1[v][mi] = (F1[v][mi] + inv5 * R.z(-a * i)) % p
    F2 = [[0] * 10 for _ in range(10)]
    for mi, (d, a) in enumerate(MODES2):
        for i in range(5):
            j = (i + d) % 5; sgn = 1 if i < j else -1
            v = PAIRS2.index((min(i, j), max(i, j)))
            F2[v][mi] = (F2[v][mi] + sgn * inv5 * R.z(-a * i)) % p
    return F1, F2


def lin_in_modes(coeffs, F, p):
    return [sum(coeffs[v] * F[v][mi] for v in range(len(coeffs))) % p for mi in range(len(F[0]))]


def quad_in_modes(Q, F1, p):
    n = len(F1[0])
    return [[sum(F1[v][mi] * Q[v][w] * F1[w][mj] for v in range(15) for w in range(15)) % p for mj in range(n)] for mi in range(n)]


def structure(polys, R, graphs=range(1, 8), rng=None):
    p = R.p; F1, F2 = modes(R)
    out = dict(lin={}, m2={}, quad={}, classes_ok=True)
    for G in graphs:
        for k in (-6, -1):
            c = lin_in_modes(linear_part(polys[(G, k)], M1), F1, p)
            out['lin'][(G, k)] = c
            if any(x for mi, x in enumerate(c) if MODES1[mi][1] != 1): out['classes_ok'] = False
        b = lin_in_modes(linear_part(polys[(G, -2)], M2), F2, p)
        out['m2'][G] = b
        if any(x for mi, x in enumerate(b) if MODES2[mi][1] != 2): out['classes_ok'] = False
        QF = quad_in_modes(quad_part(polys[(G, -2)], p), F1, p)
        out['quad'][G] = QF
        if any(QF[mi][mj] for mi in range(15) for mj in range(15) if (MODES1[mi][1] + MODES1[mj][1]) % 5 != 2): out['classes_ok'] = False
    # rank of the linear forms in the (d,1) modes
    rows = [[c[MODES1.index((d, 1))] for d in range(3)] for c in out['lin'].values()]
    out['lin_rank'] = len(fp_rref(rows, 3, p)[0])
    # beta^Gamma = (coefficient of B_{1,2}, B_{2,2}); kappa if rank 1
    beta = {G: (out['m2'][G][MODES2.index((1, 2))], out['m2'][G][MODES2.index((2, 2))]) for G in graphs}
    out['beta'] = beta
    out['beta_rank'] = len(fp_rref([list(beta[G]) for G in graphs], 2, p)[0])
    # C^Gamma, E^Gamma: coefficient of x_d y_d' (= QF[(d,0)][(d',2)] + QF[(d',2)][(d,0)]) and z_d w_d'
    C = {}; E = {}
    for G in graphs:
        QF = out['quad'][G]
        C[G] = [[(QF[MODES1.index((d, 0))][MODES1.index((dd, 2))] + QF[MODES1.index((dd, 2))][MODES1.index((d, 0))]) % p for dd in range(3)] for d in range(3)]
        E[G] = [[(QF[MODES1.index((d, 3))][MODES1.index((dd, 4))] + QF[MODES1.index((dd, 4))][MODES1.index((d, 3))]) % p for dd in range(3)] for d in range(3)]
    out['C'] = C; out['E'] = E
    # the reduced quadrics on K: eliminate B(m2) with a pivot graph g0 (beta != 0) if beta has rank 1
    gl = list(graphs); g0 = next((G for G in gl if any(beta[G])), None); out['pivot'] = g0
    if g0 is not None and out['beta_rank'] == 1:
        red = []
        for G in gl:
            if G == g0: continue
            b1 = beta[g0][0] if beta[g0][0] else beta[g0][1]; bG = beta[G][0] if beta[g0][0] else beta[G][1]
            Cr = [[(b1 * C[G][d][dd] - bG * C[g0][d][dd]) % p for dd in range(3)] for d in range(3)]
            Er = [[(b1 * E[G][d][dd] - bG * E[g0][d][dd]) % p for dd in range(3)] for d in range(3)]
            red.append((G, Cr, Er))
        out['reduced'] = red
        # independent quadrics: rank of the 18-vectors
        vecs = [[x for row in Cr for x in row] + [x for row in Er for x in row] for (G, Cr, Er) in red]
        out['n_quadrics'] = len(fp_rref(vecs, 18, p)[0])
        # local dimension of the locus at random points of the linear subspaces: Jacobian rank
        rng = rng or random.Random(7)
        dims = {}
        for (sx, sz) in (("x", "z"), ("x", "w"), ("y", "z"), ("y", "w")):
            ranks = []
            for _ in range(3):
                pt = {v: [rng.randrange(p) for _ in range(3)] for v in "xyzw"}
                pt[sx] = [0, 0, 0]; pt[sz] = [0, 0, 0]
                J = []
                for (G, Cr, Er) in red:
                    gx = [sum(Cr[d][dd] * pt['y'][dd] for dd in range(3)) % p for d in range(3)]
                    gy = [sum(Cr[d][dd] * pt['x'][d] for d in range(3)) % p for dd in range(3)]
                    gz = [sum(Er[d][dd] * pt['w'][dd] for dd in range(3)) % p for d in range(3)]
                    gw = [sum(Er[d][dd] * pt['z'][d] for d in range(3)) % p for dd in range(3)]
                    J.append(gx + gy + gz + gw)
                ranks.append(len(fp_rref(J, 12, p)[0]))
            dims[(sx, sz)] = [12 - r for r in ranks]
        out['jacobian_dims'] = dims
    return out


def fmt_mat(R, A):
    return "[" + "; ".join(" ".join(R.fmt(x) for x in row) for row in A) + "]"
