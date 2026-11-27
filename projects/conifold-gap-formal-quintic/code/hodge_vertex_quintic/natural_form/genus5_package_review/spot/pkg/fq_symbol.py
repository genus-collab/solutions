# -*- coding: utf-8 -*-
"""fq_symbol.py -- the CONIFOLD SYMBOL R^(zhat) of the formal-quintic R-matrix (genus2_symbol programme, step 1).

Normalised rows (fq_sser.SRows, conifold frame, s = 1/L, u = s^5):  Q_a^{(i)}(z) = Sum_n Q^{(i)}_{a,n} z^n with
val_s Q^{(i)}_{a,n} = -4n + c_a - 1, c = (0,-1,3,2,1), colour-independent.  In zhat := z/s^4 the rows have a limit:
   Q_a^{(i)}(z) = s^{c_a-1} ( R^_a^{(i)}(zhat) + O(s^5) ),    R^^{(i)}_{a,n} := [s^{-4n+c_a-1}] Q^{(i)}_{a,n}   (one F_p number per (a,n,i)).
The prefactors s^{c_a-1} cancel in every eta-pairing (eta^{ab} = (1/5) delta_{a+b=3 mod 5}, c_a + c_b = 2 on the pairs), so
R^ is a constant symplectic series:  Sum_ab R^_a^{(i)}(zhat) eta^{ab} R^_b^{(j)}(-zhat) = delta_ij.
Z/5 covariance: P_a^{(i)}(z) = zeta^{ia} P_a^{(0)}(zeta^{-i} z) and Q = P/Delta^{1/2}, Delta_i^{1/2} = zeta^{4i} x give
   Q^{(i)}_{a,n} = zeta^{i(a+1-n)} Q^{(0)}_{a,n},  hence  R^_a^{(i)}(zhat) = zeta^{i(a+1)} R^_a^{(0)}(zeta^{-i} zhat)
(the extra zeta^i relative to the P-rows comes from Delta_i^{-1/2}; the s-scaling carries no phase).  Both phases are tested.
Teleman's colour-colour form:  R~_{ji}(zhat) = eta(e~_j, R^ e~_i) = Sum_ab R^^{(j)}_{b,0} eta^{ba} R^_a^{(i)}(zhat)
(the frozen e~_j in flat components are R^^{(j)}_{.,0} = eta(H^., e~_j) at z = 0); R~(0) = 1 by symplecticity.
Colour shift: R~_{j+1,i+1}(zhat) = R~_{ji}(zeta^{-1} zhat), i.e. the coefficient matrices satisfy C_n[j+1][i+1] = zeta^{-n} C_n[j][i].
Matrix log:  log R~(zhat) = Sum_k mu_k zhat^k;  symplecticity <=> mu_k^T = (-1)^{k+1} mu_k.
Galois: sigma_m (zeta -> zeta^m) acts by colour relabelling i -> m i (flat indices fixed): sigma_m(C_n[j][i]) = C_n[mj][mi],
   used with fq_scalar_lib.qzeta_coords to reconstruct entries in Q(zeta).
Hodge/Gamma test: a scalar series g(zhat) with g(0) = 1 is exp(odd) iff log g has no even coefficients; then with
   b_k := [zhat^{2k-1}] log g,  p_{2k-1} := -b_k 2k(2k-1)/B_{2k}  is the power sum Sum_a v_a^{1-2k} if g = h_v (QRR constant of
   the weight multiset v, fq_hodge.hv convention h_v = exp(-Sum_k B_{2k}/(2k(2k-1)) zhat^{2k-1} p_{2k-1}(1/v))).
API: Symbol(S) with .Rh[a][i][n], .val_ok, .val_bad, .symplectic_check(nmax) -> failures, .covariance() -> dict,
     .colour_matrix() -> [C_0..C_Nz], .is_diagonal(C), .circulant_check(C), .diag_series(C, i), .identity_row();
     mat_series_log(C, R, kmax) -> [mu_1..]; ser_log(coeffs, R) -> log coefficients; exp_odd_test(coeffs, R) -> (even, odd);
     power_sums_from_b(R, b) -> {2k-1: p}; galois_entry(R, F, j, i) -> Q(zeta) coords; h_series(R, i, nmax, sign).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import bernoulli
from fq_hodge import hv, wt
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt

CEXP = (0, -1, 3, 2, 1)


class Symbol(object):
    def __init__(self, S):
        R = S.R; p = R.p; self.R = R; self.p = p; self.Nz = S.Nz; self.eta_inv = S.eta_inv
        N = self.Nz
        self.Rh = [[[0] * (N + 1) for i in range(5)] for a in range(5)]
        self.val_bad = []
        for a in range(5):
            for i in range(5):
                for n in range(N + 1):
                    Q = S.Q[a][i][n]; exp_v = -4 * n + CEXP[a] - 1
                    if not Q.c or Q.v != exp_v or Q.prec <= exp_v:
                        self.val_bad.append((a, i, n, Q.v if Q.c else None, Q.prec))
                    else:
                        self.Rh[a][i][n] = Q.c[0] % p
        self.val_ok = not self.val_bad

    # ---- checks
    def symplectic_check(self, nmax=None):
        """failures of Sum_ab Sum_m R^_{a,m}^{(i)} eta^{ab} (-1)^{n-m} R^_{b,n-m}^{(j)} = delta_ij delta_n0, n <= nmax"""
        p = self.p; N = self.Nz if nmax is None else nmax; bad = []
        for i in range(5):
            for j in range(5):
                for n in range(N + 1):
                    acc = 0
                    for (a, b), e in self.eta_inv.items():
                        for m in range(n + 1):
                            t = self.Rh[a][i][m] * self.Rh[b][j][n - m] * e
                            acc += t if (n - m) % 2 == 0 else -t
                    acc %= p
                    if acc != (1 if (i == j and n == 0) else 0): bad.append((i, j, n, acc))
        return bad

    def covariance(self):
        """for every (a, n): the exponents k (mod 5) with R^_{a,n}^{(i)} = zeta^{k i} R^_{a,n}^{(0)} for all i, or None;
        returns dict (a, n) -> k, plus the two hypotheses tested: 'a+1-n' and 'a-n'"""
        R = self.R; p = self.p; out = {}
        for a in range(5):
            for n in range(self.Nz + 1):
                base = self.Rh[a][0][n]; found = None
                for k in range(5):
                    if all(self.Rh[a][i][n] == base * R.z(k * i) % p for i in range(5)): found = k; break
                out[(a, n)] = found
        hyp1 = all(out[(a, n)] == (a + 1 - n) % 5 for (a, n) in out)
        hyp2 = all(out[(a, n)] == (a - n) % 5 for (a, n) in out)
        return dict(k=out, phase_a_plus_1_minus_n=hyp1, phase_a_minus_n=hyp2)

    # ---- colour-colour form
    def colour_matrix(self):
        """C[n][j][i] = Sum_ab R^_{b,0}^{(j)} eta^{ba} R^_{a,n}^{(i)}  (R~_{ji}(zhat) = Sum_n C[n][j][i] zhat^n)"""
        p = self.p; out = []
        for n in range(self.Nz + 1):
            C = [[0] * 5 for _ in range(5)]
            for j in range(5):
                for i in range(5):
                    C[j][i] = sum(self.Rh[b][j][0] * e * self.Rh[a][i][n] for (a, b), e in self.eta_inv.items()) % p
            out.append(C)
        return out

    def is_diagonal(self, C):
        return all(C[n][j][i] == 0 for n in range(len(C)) for j in range(5) for i in range(5) if i != j)

    def circulant_check(self, C):
        """C_n[j+1][i+1] = zeta^{-n} C_n[j][i] for all n, i, j"""
        R = self.R; p = self.p
        return all(C[n][(j + 1) % 5][(i + 1) % 5] == C[n][j][i] * R.z(-n) % p for n in range(len(C)) for j in range(5) for i in range(5))

    def diag_series(self, C, i):
        return [C[n][i][i] for n in range(len(C))]

    def identity_row(self):
        """f(zhat) = R^_0^{(0)}(zhat): coefficient list"""
        return [self.Rh[0][0][n] for n in range(self.Nz + 1)]

    def galois_entry(self, C, n, j, i):
        """Q(zeta) coordinates of C_n[j][i] from its conjugates C_n[mj][mi]"""
        return qzeta_coords(self.R, [C[n][(m * j) % 5][(m * i) % 5] for m in range(1, 5)])


# ------------------------------------------------------------------ series helpers over F_p
def _mat_mul(A, B, p):
    return [[sum(A[i][k] * B[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
def _mat_add(A, B, p, c=1):
    return [[(A[i][j] + c * B[i][j]) % p for j in range(5)] for i in range(5)]
def _mat_zero(): return [[0] * 5 for _ in range(5)]


def mat_series_mul(X, Y, p, kmax):
    out = [_mat_zero() for _ in range(kmax + 1)]
    for a in range(kmax + 1):
        for b in range(kmax + 1 - a):
            if any(x for r in X[a] for x in r) and any(y for r in Y[b] for y in r):
                out[a + b] = _mat_add(out[a + b], _mat_mul(X[a], Y[b], p), p)
    return out


def mat_series_log(C, R, kmax):
    """log of the matrix series C[0] = 1, C[1], ...: returns [mu_0 = 0, mu_1, ..., mu_kmax] (log(1+X) = Sum (-1)^{m+1} X^m/m)"""
    p = R.p
    assert all(C[0][i][j] == (1 if i == j else 0) for i in range(5) for j in range(5)), "constant term must be the identity"
    X = [_mat_zero()] + [C[n] for n in range(1, kmax + 1)]
    out = [_mat_zero() for _ in range(kmax + 1)]; power = X
    for m in range(1, kmax + 1):
        c = R.inv(m) * (1 if m % 2 == 1 else p - 1) % p
        for k in range(kmax + 1): out[k] = _mat_add(out[k], power[k], p, c)
        power = mat_series_mul(power, X, p, kmax)
    return out


def mat_series_exp(mu, R, kmax):
    """exp of the matrix series mu (mu[0] = 0): [z^0..z^kmax]"""
    p = R.p
    out = [_mat_zero() for _ in range(kmax + 1)]; out[0] = [[1 if i == j else 0 for j in range(5)] for i in range(5)]
    power = [list(map(list, out[0]))] + [_mat_zero() for _ in range(kmax)]; fact = 1
    for m in range(1, kmax + 1):
        power = mat_series_mul(power, mu, p, kmax); fact *= m; c = R.inv(fact)
        for k in range(kmax + 1): out[k] = _mat_add(out[k], power[k], p, c)
    return out


def ser_log(c, R):
    """log of the scalar series with coefficients c (c[0] = 1): coefficient list of the same length"""
    p = R.p; N = len(c)
    assert c[0] % p == 1, "log needs constant term 1"
    # theta f / f integrated: l_n = c_n - (1/n) Sum_{k=1}^{n-1} k l_k c_{n-k}
    l = [0] * N
    for n in range(1, N):
        s = sum(k * l[k] * c[n - k] for k in range(1, n)) % p
        l[n] = (c[n] - s * R.inv(n)) % p
    return l


def ser_mul(a, b, p):
    N = min(len(a), len(b))
    return [sum(a[k] * b[n - k] for k in range(n + 1)) % p for n in range(N)]


def exp_odd_test(c, R):
    """(even coefficients of log c (n >= 2, should vanish), odd coefficients {2k-1: b_k})"""
    l = ser_log(c, R)
    even = {n: l[n] for n in range(2, len(l), 2)}
    odd = {n: l[n] for n in range(1, len(l), 2)}
    return even, odd


def power_sums_from_b(R, odd):
    """p_{2k-1} = -b_k 2k(2k-1)/B_{2k} from odd = {2k-1: b_k}"""
    out = {}
    for n, b in odd.items():
        k = (n + 1) // 2
        out[n] = (-b * R.of(Fr(2 * k * (2 * k - 1)) / bernoulli(2 * k))) % R.p
    return out


def power_sums_of(R, v, ns):
    """{n: Sum_a v_a^{-n}} for the weight multiset v (field elements)"""
    return {n: sum(R.pow(R.inv(x), n) for x in v) % R.p for n in ns}


def h_series(R, i, nmax, sign=1):
    """h_{wt_i}(zhat)^sign coefficients (fq_hodge.hv, the QRR constants of the formal quintic at the zeta point)"""
    lams = [R.z(k) for k in range(5)]
    return hv(R, wt(R, lams, i), nmax, sign)


def fmt_list(R, c):
    return "[" + ", ".join(R.fmt(x) for x in c) + "]"
