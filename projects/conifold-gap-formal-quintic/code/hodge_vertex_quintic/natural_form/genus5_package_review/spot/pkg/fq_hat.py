# -*- coding: utf-8 -*-
"""fq_hat.py -- the HATTED theory of the formal quintic at genus 2 (genus2_ou programme): the full u-dependence of the colour-summed
graph terms for a q-dependent symplectic factorisation R = R' M(zhat; u), M constant in zhat = z/s^4 (s = 1/L, u = s^5).

Rescaling (FORMAL_QUINTIC s.3, CONIFOLD_SYMBOL s.1): Hhat^a := s^{1-c_a} H^a (c = (0,-1,3,2,1); eta is preserved on the eta-pairs
since c_a + c_b = 2), z = s^4 zhat.  The normalised rows become u-series
    Qhat[a][i][n] := Q[a][i][n] s^{4n+1-c_a}    (Qhat|_{u=0} = the conifold symbol R^_{a,n}^{(i)}, fq_symbol.Symbol.Rh),
and the TFT Delta^hat_i^{1/2} = zeta^{4i} xhat, xhat := x/s = f_0(u) = varpi^cf_0(u) (xhat|_{u=0} = 1).  Valuation count: for a decoration of
a genus-2 graph with half-edge exponents a_h, leaves d_l, inner psi-deficits D_v, the s-exponent of the un-hatted product is
    (2g-2+m) - 4 Sum D_v - 4 Sum_e (a_h + a_h' + 1) - Sum_l (4(d_l-1) + 1) = -10(g-1) = -10,
so for ANY symplectic M(zhat) constant in u (a fq_sser.SympM or a RawSympM whose coefficient matrices are read as [zhat^k], i.e.
M(z) = Sum_k m_k s^{-4k} z^k) the fq_msum.GraphSummer run on HatRows (duck-typed SRows: Q = Qhat, x = xhat, xinv, R, Nz, eta_inv; the
same role as fq_leadsolve.FrozenRows, which is HatRows at u = 0) returns  G_Gamma(u) := s^{10} T_Gamma  as a Laurent series in s whose
only non-zero coefficients sit at s^{5k} (a u-series; this is CHECKED, not assumed: HatRows works in the variable s).  The genuine term
is T_Gamma = u^{-2} G_Gamma(u): REGULAR at the conifold point iff [u^0] G = [u^1] G = 0.  In t_c = u (1 + c_1 u + ...):
    T_Gamma = a_0 t_c^{-2} + (a_1 + 2 c_1 a_0) t_c^{-1} + ...,   a_k = [u^k] G_Gamma     (tc_expansion, exact via Ser.reverse).
RawSympM: a duck-typed SympM for a matrix series M(zhat) = Sum_k M_k zhat^k given by its coefficient matrices (.M, .Minv lists of 4,
.R, .p, .is_identity = False); symplecticity M(zhat)^T M(-zhat) = 1 and M M^{-1} = 1 mod zhat^4 are asserted.  rows_prime, EdgeLeaf
and MVertex use only these attributes.  m0_exp(R, M0, N, eps) builds M_0 exp(eps N) (RIGHT placement), exp_m0 builds exp(eps N) M_0
(LEFT), conj_series gives M_0^{-1} N' M_0 mod zhat^4; tautological_M0 returns M_0 = R~ mod zhat^4 as a SympM (from mu_k = [zhat^k] log R~)
together with the check that it reproduces the colour matrices C_1, C_2, C_3 of the symbol.
Interpolation in eps over F_p: vander_inv(nodes, p) (exact inverse of the Vandermonde matrix), interp_fp(vals, Vinv, p) -> coefficients,
interp_ser(sers, Vinv) -> coefficient Ser (linear, applied term by term).
API: HatRows(S, prec_s=None); hat_sums(H, M, graphs) -> {idx: Ser}; ucoef(ser, k) = [u^k]; ucoefs(ser, kmax); RawSympM; nu_series;
     m0_exp / exp_m0 / conj_series / tautological_M0; tc_series(C, R, N); tc_expansion(acoefs, tc, R) -> (A_{-2}, A_{-1}, A_0, ...);
     vander_inv / interp_fp / interp_ser; point_frame_check(S, H, graphs, kmax) -> (ok, detail).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_sser import SympM, mat_mul, mat_T, mat_id, mat_zero, mat_add, mat_eq
from fq_symbol import mat_series_mul, mat_series_exp, mat_series_log, CEXP
from fq_msum import GraphSummer, outer_graphs
from gsum_q import is_u_series


class HatRows(object):
    """duck-typed fq_sser.SRows for GraphSummer: Qhat[a][i][n] = Q[a][i][n] s^{4n+1-c_a}, xhat = x/s (Laurent series in s, u-series)"""
    def __init__(self, S, prec_s=None):
        R = S.R; self.R = R; self.Nz = S.Nz; self.eta_inv = S.eta_inv; self.S = S
        self.Q = [[[S.Q[a][i][n].shift(4 * n + 1 - CEXP[a]) for n in range(S.Nz + 1)] for i in range(5)] for a in range(5)]
        self.x = S.x.shift(-1)
        if prec_s is not None:
            self.Q = [[[q.truncate(prec_s) for q in r] for r in a] for a in self.Q]; self.x = self.x.truncate(prec_s)
        self.xinv = self.x.inv()
        allq = [self.Q[a][i][n] for a in range(5) for i in range(5) for n in range(S.Nz + 1)]
        self.useries_ok = all(is_u_series(q) for q in allq) and is_u_series(self.x)
        self.val_ok = all(q.c and q.v == 0 for q in allq) and self.x.v == 0 and self.x.c[0] == R.one
        self.min_prec = min(q.prec for q in allq + [self.x])

    def symbol_check(self, sym):
        """Qhat|_{u=0} = R^ (the conifold symbol of fq_symbol.Symbol)"""
        return all(self.Q[a][i][n].coef(0) == sym.Rh[a][i][n] for a in range(5) for i in range(5) for n in range(self.Nz + 1))


def hat_sums(H, M, graphs=None):
    """{idx: G_Gamma as a Laurent series in s} = s^{10} x (colour-summed term of the stable graph) for the factorisation R = R' M(zhat)"""
    return GraphSummer(H, M).all_graphs(graphs if graphs is not None else outer_graphs())


def ucoef(ser, k):
    """[u^k] of a u-series written in s"""
    return ser.coef(5 * k)


def ucoefs(ser, kmax):
    return [ucoef(ser, k) for k in range(kmax + 1)]


def uprec(ser):
    """number of known u-coefficients (u^0..u^{uprec-1})"""
    return (ser.prec + 4) // 5


# ------------------------------------------------------------------ RawSympM and the products
class RawSympM(object):
    """M(zhat) = Sum_k M[k] zhat^k mod zhat^4, symplectic (M(zhat)^T M(-zhat) = 1), with inverse Minv; duck-typed fq_sser.SympM"""
    def __init__(self, R, M, Minv, check=True):
        p = R.p; self.R = R; self.p = p; self.M = M; self.Minv = Minv; self.is_identity = False
        if check:
            for k in range(4):
                S1 = mat_zero(5); S2 = mat_zero(5)
                for a in range(k + 1):
                    S1 = mat_add(S1, mat_mul(mat_T(M[a]), M[k - a], p), p, (1 if (k - a) % 2 == 0 else p - 1))
                    S2 = mat_add(S2, mat_mul(M[a], Minv[k - a], p), p)
                target = mat_id(5, p) if k == 0 else mat_zero(5)
                assert mat_eq(S1, target, p), "RawSympM: not symplectic at zhat^%d" % k
                assert mat_eq(S2, target, p), "RawSympM: M M^-1 != 1 at zhat^%d" % k


def nu_series(nu1, nu2, nu3):
    """N(zhat) = nu1 zhat + nu2 zhat^2 + nu3 zhat^3 as a matrix series [N_0 = 0, N_1, N_2, N_3]"""
    return [mat_zero(5), nu1, nu2, nu3]


def series_scale(X, c, p):
    return [[[x * c % p for x in r] for r in Xk] for Xk in X]


def series_add(X, Y, p, c=1):
    return [mat_add(X[k], Y[k], p, c) for k in range(len(X))]


def series_in_sp(N, p):
    """N_k^T = (-1)^{k+1} N_k for k = 1..3 (N_0 = 0): the Lie algebra of the symplectic series"""
    return (not any(x for r in N[0] for x in r)) and all((N[k][i][j] - (-1) ** (k + 1) * N[k][j][i]) % p == 0 for k in (1, 2, 3) for i in range(5) for j in range(5))


def m0_exp(R, M0, M0inv, N, eps, check=True):
    """RIGHT placement: M_0 exp(eps N) mod zhat^4 and its inverse exp(-eps N) M_0^{-1}"""
    p = R.p; e = R.of(eps)
    E = mat_series_exp(series_scale(N, e, p), R, 3); Ei = mat_series_exp(series_scale(N, (-e) % p, p), R, 3)
    return RawSympM(R, mat_series_mul(M0, E, p, 3), mat_series_mul(Ei, M0inv, p, 3), check)


def exp_m0(R, M0, M0inv, N, eps, check=True):
    """LEFT placement: exp(eps N) M_0 mod zhat^4 and its inverse M_0^{-1} exp(-eps N)"""
    p = R.p; e = R.of(eps)
    E = mat_series_exp(series_scale(N, e, p), R, 3); Ei = mat_series_exp(series_scale(N, (-e) % p, p), R, 3)
    return RawSympM(R, mat_series_mul(E, M0, p, 3), mat_series_mul(M0inv, Ei, p, 3), check)


def conj_series(M0inv, N, M0, p):
    """M_0^{-1} N M_0 mod zhat^4 (the N of the RIGHT placement equivalent to N' = N in the LEFT placement)"""
    return mat_series_mul(mat_series_mul(M0inv, N, p, 3), M0, p, 3)


def tautological_M0(sym, R):
    """M_0 = R~ mod zhat^4 as a SympM (mu_k = [zhat^k] log R~, k <= 3); returns (SympM, mu, C, ok) with ok: SympM.M[k] == C[k], k <= 3"""
    C = sym.colour_matrix(); mu = mat_series_log(C, R, 3)
    Ms = SympM(R, mu[1], mu[2], mu[3])
    ok = all(mat_eq(Ms.M[k], C[k], R.p) for k in range(4))
    return Ms, mu, C, ok


# ------------------------------------------------------------------ t_c conversion
def tc_series(C, R, N=8):
    """t_c(u) = g/varpi^cf_0 of the fq_conifold.Conifold object C as a Ser over R, truncated at O(u^N)"""
    return Ser(R, C.tc.v, [R.of(c) for c in C.tc.c]).truncate(N)


def tc_expansion(acoefs, tc, R):
    """u^{-2} Sum_k a_k u^k rewritten in t_c: returns (coefficients of t_c^{-2}, t_c^{-1}, t_c^0, ...) as far as the data allow"""
    G = Ser(R, 0, [R.of(a) for a in acoefs])
    u_of_t = tc.reverse()
    T = G.compose(u_of_t).mul(u_of_t.inv().pow(2))
    return [T.coef(k) for k in range(-2, T.prec)]


# ------------------------------------------------------------------ interpolation in eps over F_p
def vander_inv(nodes, p):
    """inverse of V[t][k] = nodes[t]^k (Gauss-Jordan over F_p)"""
    n = len(nodes)
    A = [[pow(nodes[t], k, p) for k in range(n)] + [1 if j == t else 0 for j in range(n)] for t in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c])
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p); A[c] = [x * inv % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def interp_fp(vals, Vinv, p):
    """coefficients c_k of the polynomial with values vals[t] at the nodes of Vinv: c = Vinv . vals"""
    n = len(vals)
    return [sum(Vinv[k][t] * vals[t] for t in range(n)) % p for k in range(n)]


def interp_ser(sers, Vinv):
    """the same for Ser-valued data (coefficient series P_k)"""
    n = len(sers); out = []
    for k in range(n):
        acc = None
        for t in range(n):
            c = Vinv[k][t]
            if not c: continue
            term = sers[t].scale(c); acc = term if acc is None else acc.add(term)
        out.append(acc if acc is not None else Ser(sers[0].R, sers[0].prec, []))
    return out


def point_frame_check(S, H, graphs, kmax=10):
    """(a): s^{-10} x hatted sums at M = 1 equal the SRows sums (fq_msum.GraphSummer(S, identity)) through s^kmax"""
    R = S.R; I = SympM.identity(R)
    rs = GraphSummer(S, I).all_graphs(graphs); rh = hat_sums(H, I, graphs)
    bad = []
    for G in rs:
        d = rs[G].sub(rh[G].shift(-10)).truncate(kmax + 1)
        if d.c: bad.append((G, d.fmt(kmax)))
    return not bad, bad
