# -*- coding: utf-8 -*-
"""fq_ou.py -- the O(u) correction of the tautological factorisation at genus 2 (genus2_ou programme): M(zhat; u) = M_0 exp(u N(zhat)),
M_0 = R~ mod zhat^4 (fq_hat.tautological_M0), N = nu_1 zhat + nu_2 zhat^2 + nu_3 zhat^3 in sp (nu_1, nu_3 symmetric, nu_2 antisymmetric;
40 coordinates in the fq_leadsolve.mu_params / params_to_mats convention: 0..14 nu_1 (i <= j), 15..24 nu_2 (i < j), 25..39 nu_3).

Every genus-2 graph term is a polynomial in the entries of the coefficient matrices of M (total M-order <= 3), so with eps a scalar
   P_Gamma(eps; u) := G_Gamma[M_0 exp(eps N)](u)   (hatted sum, fq_hat)   is a polynomial in eps of degree <= 3 with u-series coefficients,
interpolated exactly from eps = 0..9 (eps^4..eps^9 must vanish: checked).  The genuine u-dependent factorisation is eps = u:
   G_Gamma[M_0 exp(u N)](u) = Sum_k u^k P_k(u),   [u^0] = L_Gamma(M_0) = (1/1200, 0, .., 0),   [u^1] = g1_Gamma + D_N L_Gamma,
with g1_Gamma := [u^1] G_Gamma[M_0] and D_N L_Gamma := d/d eps|_0 G^0_Gamma[M_0 exp(eps N)] (frozen theory, fq_leadsolve.FrozenRows), linear
in N: D_N L_Gamma = Sum_v J[Gamma][v] n_v.  J is computed from the 40 basis directions (jacobian_N).  All six edge graphs are regular iff
   J[2] . n = -g1_2,  J[4] . n = -g1_4        (J[3], J[5], J[6], J[7] vanish identically: checked).
Classes of N (RIGHT placement M_0 exp(uN); LEFT placement exp(uN') M_0, i.e. N = M_0^{-1} N' M_0 with N' in the class): covariant scalar
(fq_leadsolve.class_covariant_diagonal: nu_1 = b_1 diag(zeta^{-i}), nu_3 = b_3 diag(zeta^{-3i})), diagonal in Fourier modes
(modes_diagonal: x_i = (1/5) Sum_a zeta^{-ai} X_a; the restricted coefficients are then RATIONAL), circulant (class_circulant, 8 parameters),
full (identity, 40).  class_columns gives the 40-vectors of the class basis (left: conjugated by M_0), restricted gives the 2 x k matrix,
solve_restricted solves it over F_p (fq_mfit.fp_solve: particular solution, nullspace, rank).
Cross-check of J (optional): the leading polynomials L_Gamma(mu) of symbol_result.json (M_0 = exp(mu)) have mu-gradient J_mu at the
tautological point; exp(mu + eps Y) = exp(mu) exp(eps N) + O(eps^2) with N = dexp(Y) = Y - [mu,Y]/2 + [mu,[mu,Y]]/6 (mod zhat^4), so
J[G] . coords(dexp(Y_e)) must equal J_mu[G][e] for the 40 basis directions Y_e (jacobian_from_json, dexp).
verify_solution: the eps-polynomial from HATTED sums and its eps = u substitution for a given N and placement.
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_mfit import PAIRS1, PAIRS2, VNAME, fp_solve, fp_rref
from fq_leadsolve import (FrozenRows, frozen_sums, mu_params, params_to_mats, class_covariant_diagonal, class_diagonal, class_circulant,
                          modes_diagonal, jacobian)
from fq_leadfit import MONO3_INDEX
from fq_symbol import mat_series_mul, _mat_add
from fq_hat import (RawSympM, nu_series, series_in_sp, m0_exp, exp_m0, conj_series, hat_sums, interp_fp, interp_ser, vander_inv, ucoef)

EPS_NODES = list(range(10))
GRAPHS = list(range(1, 8))


def N_of_params(par, p):
    nu1, nu2, nu3 = params_to_mats(par, p); return nu_series(nu1, nu2, nu3)


def params_of_N(N, p):
    assert series_in_sp(N, p), "N not in sp"
    return mu_params(N[1], N[2], N[3], p)


def direction_N(v, p):
    par = [0] * 40; par[v] = 1; return N_of_params(par, p)


def placed(R, M0, M0inv, N, eps, placement):
    return m0_exp(R, M0, M0inv, N, eps, check=False) if placement == 'right' else exp_m0(R, M0, M0inv, N, eps, check=False)


def frozen_eps_poly(fr, M0, M0inv, N, graphs, R, Vinv, placement='right'):
    """{G: [c_0, .., c_9]} coefficients in eps of G^0_Gamma[M_0 exp(eps N)] (frozen theory)"""
    p = R.p; vals = {}
    for eps in EPS_NODES:
        res, _ = frozen_sums(fr, placed(R, M0, M0inv, N, eps, placement), graphs)
        for G, x in res.items(): vals.setdefault(G, []).append(x)
    return {G: interp_fp(vals[G], Vinv, p) for G in vals}


def jacobian_N(fr, M0, M0inv, graphs, R, Vinv, log=None):
    """J[G][v] = d/d eps|_0 G^0_Gamma[M_0 exp(eps N_v)]; also the maximal eps-degree seen and the constant terms (must be L_Gamma(M_0))"""
    p = R.p; J = {G: [0] * 40 for G in GRAPHS}; maxdeg = 0; consts = {}
    for v in range(40):
        polys = frozen_eps_poly(fr, M0, M0inv, direction_N(v, p), graphs, R, Vinv)
        for G, c in polys.items():
            deg = max([k for k, x in enumerate(c) if x] or [0]); maxdeg = max(maxdeg, deg)
            J[G][v] = c[1]; consts.setdefault(G, set()).add(c[0])
        if log and (v + 1) % 10 == 0: log("      direction %d/40 done" % (v + 1))
    return J, maxdeg, consts


def hat_eps_poly(H, M0, M0inv, N, graphs, R, Vinv, placement='right'):
    """{G: [P_0(u), .., P_9(u)]} (Ser in s) of P_Gamma(eps; u) = G_Gamma[M_0 exp(eps N)](u)"""
    sers = {}
    for eps in EPS_NODES:
        res = hat_sums(H, placed(R, M0, M0inv, N, eps, placement), graphs)
        for G, s in res.items(): sers.setdefault(G, []).append(s)
    return {G: interp_ser(sers[G], Vinv) for G in sers}


def eps_to_u(P):
    """Sum_k u^k P_k(u) as a Ser in s"""
    acc = None
    for k, Pk in enumerate(P):
        t = Pk.shift(5 * k); acc = t if acc is None else acc.add(t)
    return acc


def eps_degree(P):
    """highest k with P_k != 0 (Ser list) or coefficient list"""
    return max([k for k, c in enumerate(P) if (c.c if isinstance(c, Ser) else c)] or [0])


# ------------------------------------------------------------------ classes and placements
def class_columns(R, M0, M0inv, cls, placement):
    """list over the class parameters of 40-vectors (N-coordinates of the basis element), names; left: N = M_0^{-1} N' M_0"""
    p = R.p
    if cls == 'cov': A, names = class_covariant_diagonal(R)
    elif cls == 'diag': A, names = modes_diagonal(R)
    elif cls == 'circ': A, names = class_circulant(R)
    elif cls == 'full': A, names = [{v: 1} for v in range(40)], ["n%d" % v for v in range(40)]
    else: raise ValueError(cls)
    k = len(names); cols = []
    for j in range(k):
        col = [A[v].get(j, 0) % p for v in range(40)]
        if placement == 'left': col = params_of_N(conj_series(M0inv, N_of_params(col, p), M0, p), p)
        cols.append(col)
    return cols, names


def restricted(J, cols, p, graphs=(2, 4)):
    """rows J[G] . col_j for G in graphs: the restricted coefficient matrix"""
    return [[sum(J[G][v] * col[v] for v in range(40)) % p for col in cols] for G in graphs]


def solve_restricted(Jr, b, p):
    """(x particular or None, nullspace basis, rank) of Jr t = b"""
    return fp_solve(Jr, b, p)


def expand_solution(cols, t, p):
    """N-coordinates Sum_j t_j col_j"""
    return [sum(c[v] * t[j] for j, c in enumerate(cols)) % p for v in range(40)]


# ------------------------------------------------------------------ cross-check against the stored leading polynomials
def polys_from_json(path, R):
    """{G: {mono idx: F_p coefficient}} from symbol_result.json['lead_polys_qzeta'] (exact Q(zeta) coordinates c0..c3)"""
    p = R.p; d = json.load(open(path)); out = {}
    for Gs, mp in d['lead_polys_qzeta'].items():
        poly = {}
        for name, cs in mp.items():
            if not cs: raise ValueError("non-reconstructed coefficient in symbol_result.json")
            mono = () if name == "1" else tuple(sorted(VNAME.index(t) for t in name.split("*")))
            poly[MONO3_INDEX[mono]] = sum(R.of(Fr(cs[k])) * R.z(k) for k in range(4)) % p
        out[int(Gs)] = poly
    return out


def bracket(A, B, p):
    return [_mat_add(mat_series_mul(A, B, p, 3)[k], mat_series_mul(B, A, p, 3)[k], p, p - 1) for k in range(4)]


def dexp(mu, Y, R):
    """N with exp(mu + eps Y) = exp(mu) exp(eps N) + O(eps^2) mod zhat^4: N = Y - [mu,Y]/2 + [mu,[mu,Y]]/6 (ad_mu^3 Y = O(zhat^4))"""
    p = R.p; c1 = bracket(mu, Y, p); c2 = bracket(mu, c1, p)
    return [_mat_add(_mat_add(Y[k], c1[k], p, (-R.inv(2)) % p), c2[k], p, R.inv(6)) for k in range(4)]


def jacobian_from_json(path, R, mu):
    """(J_mu rows for graphs 1..7 at the tautological point, and the 40 x 40 matrix Phi with coords(dexp(Y_e)) as columns)"""
    p = R.p; polys = polys_from_json(path, R)
    tpar = mu_params(mu[1], mu[2], mu[3], p)
    Jmu = jacobian([polys[G] for G in GRAPHS], tpar, p)
    Phi = [params_of_N(dexp(mu, direction_N(e, p), R), p) for e in range(40)]       # Phi[e] = coords of dexp(Y_e)
    return {G: Jmu[i] for i, G in enumerate(GRAPHS)}, Phi, polys


def dexp_check(J, Jmu, Phi, p):
    """J[G] . Phi[e] == Jmu[G][e] for all G, e; returns the number of failures"""
    return sum(1 for G in GRAPHS for e in range(40) if sum(J[G][v] * Phi[e][v] for v in range(40)) % p != Jmu[G][e] % p)


def dexp_numeric_check(R, mu, Y, Vinv):
    """exp(mu + eps Y) (SympM at 4 nodes, Lagrange) has linear eps-coefficient exp(mu) . dexp(Y): returns ok"""
    from fq_sser import SympM
    p = R.p; nodes = [0, 1, 2, 3]; Vi = vander_inv(nodes, p)
    vals = [SympM(R, *[[[(mu[k][i][j] + e * Y[k][i][j]) % p for j in range(5)] for i in range(5)] for k in (1, 2, 3)], check=False).M for e in nodes]
    lin = [[[sum(Vi[1][t] * vals[t][k][i][j] for t in range(4)) % p for j in range(5)] for i in range(5)] for k in range(4)]
    M0 = SympM(R, mu[1], mu[2], mu[3], check=False).M
    target = mat_series_mul(M0, dexp(mu, Y, R), p, 3)
    return all((lin[k][i][j] - target[k][i][j]) % p == 0 for k in range(4) for i in range(5) for j in range(5))


# ------------------------------------------------------------------ verification of a solution by eps -> u
def verify_solution(H, M0, M0inv, N, graphs, R, Vinv, placement, kmax=4):
    """returns dict: P (eps-polynomial, Ser coefficients), Gu (eps = u substitution), maxdeg, edge_regular, vertex u-coefficients, total"""
    P = hat_eps_poly(H, M0, M0inv, N, graphs, R, Vinv, placement)
    Gu = {G: eps_to_u(P[G]) for G in P}
    maxdeg = max(eps_degree(P[G]) for G in P)
    uc = {G: [ucoef(Gu[G], k) for k in range(kmax + 1)] for G in Gu}
    edge_reg = all(uc[G][0] == 0 and uc[G][1] == 0 for G in Gu if G != 1)
    tot = [sum(uc[G][k] for G in uc) % R.p for k in range(kmax + 1)]
    return dict(P=P, Gu=Gu, maxdeg=maxdeg, uc=uc, edge_regular=edge_reg, total=tot)
