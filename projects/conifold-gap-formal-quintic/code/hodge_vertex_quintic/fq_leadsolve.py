# -*- coding: utf-8 -*-
"""fq_leadsolve.py -- the LEADING-ORDER system for a q-dependent factorisation R = R' M_0(zhat) of the formal-quintic R-matrix
at genus 2 (genus2_symbol programme, step 3), M_0(zhat) = exp(mu_1 zhat + mu_2 zhat^2 + mu_3 zhat^3) constant in zhat = z/s^4.

Substitution: the constant-M fit (fq_leadfit, M = exp(Sum m_k z^k)) at m_k = mu_k s^{-4k} turns a weight-w monomial, sitting at
s^{-10+4w} (SYMPLECTIC_SCAN s.2), into an s^{-10} term.  Hence the s^{-10} coefficient of the colour-summed graph Gamma is
   L_Gamma(mu) = c_Gamma + [weight 1 of s^-6](mu_1) + [weight 2 of s^-2](mu_1 mu_1, mu_2) + [weight 3 of s^+2](mu_1^3, mu_1 mu_2, mu_3),
a polynomial of weight <= 3 in the 40 entries of mu (lead_poly).  It must equal the genus-2 coefficient of the stable graph Gamma
in the FROZEN theory: rows R^ M_0^{-1} (R^ the conifold symbol, fq_symbol) and vertex CohFT M_0.omega^, omega^ the constant TFT
with Delta^_i^{1/2} = zeta^{4i}: that graph sum is fq_msum.GraphSummer run on constant series (FrozenRows: Q[a][i][n] = R^_{a,n}^{(i)},
x = 1, precision O(s^1)); frozen_sums(sym, M, graphs) returns its coefficients.  Checking L_Gamma(mu) = frozen sum at random mu
is the consistency check of the substitution (every valuation bookkeeping step of SYMPLECTIC_SCAN s.2 enters).
Structure: L_Gamma is cubic in mu_1 and LINEAR in (mu_2, mu_3) jointly (the only monomials containing mu_2 or mu_3 are m2, m1 m2, m3),
so for fixed mu_1 the system {L_Gamma = 0, Gamma = 2..7} is a 6 x 25 linear system (affine_in_23): its generic rank r gives the
dimension 15 + 25 - r of the solution variety on the dense open set where the rank is r, and random fibre points for the Jacobian.
Tautological point: mu_k = [zhat^k] log R~^ (colour basis), i.e. M_0 = R~^ + O(zhat^4): then R' = 1 + O(zhat^4) kills every edge
graph at genus 2 (z-orders <= 3) and L_1 = F_2 of M_0.omega^ = 1/1200.
Sub-families (restrict_poly with a linear parametrisation mu = A t): DIAGONAL (t = x_0..x_4, y_0..y_4: mu_1 = diag x, mu_3 = diag y;
mu_2 diagonal antisymmetric = 0, so 10 of the nominal 15 parameters), COVARIANT DIAGONAL (a_1, a_3: x_i = zeta^{-i} a_1,
y_i = zeta^{-3i} a_3), CIRCULANT = Z/5-covariant full M_0 ((mu_k)_{i+1,j+1} = zeta^{-k} (mu_k)_{ij}: (mu_k)_{ij} = zeta^{-ki} c_k(j-i),
c_1(0), c_1(1), c_1(2), c_2(1), c_2(2), c_3(0), c_3(1), c_3(2): 8 parameters; symmetry forces c_1(-d) = zeta^d c_1(d),
c_2(-d) = -zeta^{2d} c_2(d), c_3(-d) = zeta^{3d} c_3(d)), and the Fourier modes of the diagonal class
(x_i = (1/5) Sum_a zeta^{-ai} X_a, y_i likewise with Y_a).  Polynomials in t are dicts {sorted exponent tuple: coeff}.
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_sser import SympM
from fq_msum import GraphSummer
from fq_leadfit import MONOS3, MONO3_W, MONO3_INDEX, M1V, M2V, M3V, weight_part, poly_eval3
from fq_mfit import PAIRS1, PAIRS2, fp_rref, fp_nullspace, fp_solve

LEAD_ORDERS = {0: -10, 1: -6, 2: -2, 3: 2}


def lead_poly(polys, G):
    """L_Gamma as {mono idx: coeff} from the fitted s-coefficients (weight w taken from s^{-10+4w})"""
    out = {}
    for w, k in LEAD_ORDERS.items():
        for m, c in weight_part(polys[(G, k)], w).items(): out[m] = c
    return out


# ------------------------------------------------------------------ the frozen theory
class FrozenRows(object):
    """duck-typed fq_sser.SRows for GraphSummer: constant rows R^_{a,n}^{(i)}, x = 1, precision O(s^1)"""
    def __init__(self, sym):
        R = sym.R; self.R = R; self.Nz = sym.Nz; self.eta_inv = sym.eta_inv
        self.Q = [[[Ser.const(R, sym.Rh[a][i][n], 1) for n in range(sym.Nz + 1)] for i in range(5)] for a in range(5)]
        self.x = Ser.const(R, 1, 1); self.xinv = self.x


def frozen_sums(frozen, M, graphs, colmask=None):
    """{idx: F_p coefficient} of the frozen graph sum with rows R^ M^-1 and vertex M.omega^ (colmask: {idx: tuple} or None)"""
    GS = GraphSummer(frozen, M); out = {}
    for G in graphs:
        s = GS.graph(G, None if colmask is None else colmask[G['idx']])
        out[G['idx']] = s.coef(0) if s.prec > 0 else 0
    return out, GS


def mu_params(mu1, mu2, mu3, p):
    return ([mu1[i][j] % p for (i, j) in PAIRS1] + [mu2[i][j] % p for (i, j) in PAIRS2] + [mu3[i][j] % p for (i, j) in PAIRS1])


def params_to_mats(par, p):
    mu1 = [[0] * 5 for _ in range(5)]; mu2 = [[0] * 5 for _ in range(5)]; mu3 = [[0] * 5 for _ in range(5)]
    for v, (i, j) in enumerate(PAIRS1): mu1[i][j] = mu1[j][i] = par[v] % p; mu3[i][j] = mu3[j][i] = par[25 + v] % p
    for v, (i, j) in enumerate(PAIRS2): mu2[i][j] = par[15 + v] % p; mu2[j][i] = (-par[15 + v]) % p
    return mu1, mu2, mu3


def sympm_of(R, par):
    mu1, mu2, mu3 = params_to_mats(par, R.p)
    return SympM(R, mu1, mu2, mu3)


# ------------------------------------------------------------------ affine structure in (mu_2, mu_3), Jacobian
def affine_in_23(poly, par, p):
    """L(mu) = A + Sum_{v in M2V+M3V} B_v mu_v at fixed mu_1 = par[0:15]; returns (A, B as list over M2V + M3V)"""
    A = 0; B = [0] * 25
    for k, c in poly.items():
        mu = MONOS3[k]
        lin = [v for v in mu if v >= 15]
        if not lin:
            r = c
            for v in mu: r = r * par[v] % p
            A = (A + r) % p
        else:
            assert len(lin) == 1, "not linear in (mu_2, mu_3): %s" % (mu,)
            r = c
            for v in mu:
                if v < 15: r = r * par[v] % p
            B[lin[0] - 15] = (B[lin[0] - 15] + r) % p
    return A, B


def jacobian(polys_list, par, p):
    """rows: d L / d mu_v (v = 0..39) at par"""
    J = []
    for poly in polys_list:
        row = [0] * 40
        for k, c in poly.items():
            mu = MONOS3[k]
            for v in set(mu):
                e = mu.count(v); r = c * e % p
                rest = list(mu); rest.remove(v)
                for w in rest: r = r * par[w] % p
                row[v] = (row[v] + r) % p
        J.append(row)
    return J


def rank_fp(rows, ncols, p):
    return len(fp_rref(rows, ncols, p)[0]) if rows else 0


def fibre_point(polys_list, mu1par, p, rng):
    """random (mu_2, mu_3) solving the 6 x 25 affine system at mu_1 = mu1par: (par40 or None, rank, consistent)"""
    A = []; b = []
    for poly in polys_list:
        a, B = affine_in_23(poly, mu1par, p); A.append(B); b.append((-a) % p)
    x, N, rank = fp_solve(A, b, p)
    if x is None: return None, rank, False
    for n in N:
        t = rng.randrange(p); x = [(xi + t * ni) % p for xi, ni in zip(x, n)]
    return list(mu1par[:15]) + x, rank, True


# ------------------------------------------------------------------ restriction to parametrised sub-families
def restrict_poly(poly, A, p):
    """A: list over the 40 variables of linear forms {t_index: coeff} (mu_v = Sum_j A[v][j] t_j); returns {t-monomial: coeff}"""
    out = {}
    for k, c in poly.items():
        terms = {(): c}
        for v in MONOS3[k]:
            lf = A[v]
            if not lf: terms = {}; break
            new = {}
            for mono, cc in terms.items():
                for j, a in lf.items():
                    m2 = tuple(sorted(mono + (j,))); new[m2] = (new.get(m2, 0) + cc * a) % p
            terms = new
        for mono, cc in terms.items():
            if cc: out[mono] = (out.get(mono, 0) + cc) % p
    return {m: c for m, c in out.items() if c}


def class_diagonal(R):
    """t = (x_0..x_4, y_0..y_4); returns (A, names)"""
    A = [dict() for _ in range(40)]
    for i in range(5):
        A[PAIRS1.index((i, i))] = {i: 1}; A[25 + PAIRS1.index((i, i))] = {5 + i: 1}
    return A, ["x%d" % i for i in range(5)] + ["y%d" % i for i in range(5)]


def class_covariant_diagonal(R):
    """t = (a_1, a_3): x_i = zeta^{-i} a_1, y_i = zeta^{-3i} a_3"""
    A = [dict() for _ in range(40)]
    for i in range(5):
        A[PAIRS1.index((i, i))] = {0: R.z(-i)}; A[25 + PAIRS1.index((i, i))] = {1: R.z(-3 * i)}
    return A, ["a1", "a3"]


def class_circulant(R):
    """t = (c1_0, c1_1, c1_2, c2_1, c2_2, c3_0, c3_1, c3_2): (mu_k)_{ij} = zeta^{-ki} c_k(j-i)"""
    p = R.p
    def ck(k, d, base):           # returns {t: coeff} for c_k(d), d = 0..4
        if k == 1: red = {0: (0, 1), 1: (1, 1), 2: (2, 1), 3: (2, R.z(2)), 4: (1, R.z(1))}
        elif k == 2: red = {1: (0, 1), 2: (1, 1), 3: (1, (-R.z(4)) % p), 4: (0, (-R.z(2)) % p)}
        else: red = {0: (0, 1), 1: (1, 1), 2: (2, 1), 3: (2, R.z(1)), 4: (1, R.z(3))}
        if d not in red: return {}
        j, c = red[d]; return {base + j: c % p}
    A = [dict() for _ in range(40)]
    for v, (i, j) in enumerate(PAIRS1):
        A[v] = {t: c * R.z(-i) % p for t, c in ck(1, (j - i) % 5, 0).items()}
        A[25 + v] = {t: c * R.z(-3 * i) % p for t, c in ck(3, (j - i) % 5, 5).items()}
    for v, (i, j) in enumerate(PAIRS2):
        A[15 + v] = {t: c * R.z(-2 * i) % p for t, c in ck(2, (j - i) % 5, 3).items()}
    return A, ["c1_0", "c1_1", "c1_2", "c2_1", "c2_2", "c3_0", "c3_1", "c3_2"]


def class_params_to_mu(A, t, p):
    return [sum(c * t[j] for j, c in A[v].items()) % p for v in range(40)]


def modes_diagonal(R):
    """t = (X_0..X_4, Y_0..Y_4): x_i = (1/5) Sum_a zeta^{-ai} X_a, y_i = (1/5) Sum_a zeta^{-ai} Y_a"""
    p = R.p; inv5 = R.inv(5); A = [dict() for _ in range(40)]
    for i in range(5):
        A[PAIRS1.index((i, i))] = {a: inv5 * R.z(-a * i) % p for a in range(5)}
        A[25 + PAIRS1.index((i, i))] = {5 + a: inv5 * R.z(-a * i) % p for a in range(5)}
    return A, ["X%d" % a for a in range(5)] + ["Y%d" % a for a in range(5)]


def tpoly_eval(tp, t, p):
    r = 0
    for mono, c in tp.items():
        x = c
        for j in mono: x = x * t[j] % p
        r = (r + x) % p
    return r


def tpoly_fmt(tp, names, R):
    if not tp: return "0"
    return " + ".join("(%s)%s" % (R.fmt(c), ("" if not m else " " + "*".join(names[j] for j in m))) for m, c in sorted(tp.items(), key=lambda kv: (len(kv[0]), kv[0])))
