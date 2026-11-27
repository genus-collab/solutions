# -*- coding: utf-8 -*-
"""fq_sser.py -- the formal-quintic rows in the conifold frame as Laurent series in s (s = 1/L, u = s^5), in the
NORMALISED canonical basis e~_i = Delta_i^{1/2} e_i (eta(e~_i, e~_j) = delta_ij), and the constant symplectic
factorisation R = R' M (genus2_symplectic programme, step 1 and the row part of step 3).

Series (fq_series.Ser over fq_fp.FpRing), Cf applied once to the rows of fq_rows.Rows (gsum_q.to_sseries, alpha = beta = 0):
  x := Cf(C_0) / L = s f_0(s^5)            (f_0 = varpi^cf_0 = 1 + (2/625) u^3 + ...),
  Delta_i^{1/2} = C_0 zeta^{3i/2}/L with zeta^{1/2} := zeta^3, i.e.  Delta_i^{1/2} = zeta^{4i} x,  Delta_i = zeta^{3i} x^2
  (checked against Cf(fq_rows.Rows.Delta(i))),
  Q_a^{(i)}(z) = eta(H^a, R(z) e~_i) = P_a^{(i)}(z)/Delta_i^{1/2},   Q[a][i][n] = [z^n],  val Q_{a,n} = -4n + (0,-1,3,2,1)_a - 1,
  symplecticity  Sum_ab Q_a^{(i)}(z) eta^{ab} Q_b^{(j)}(-z) = delta_ij  (eta^{ab} = (1/5) delta_{a+b=3 mod 5}).
Precision: the rows are truncated at the absolute precision prec_abs (default None = keep the Cf precision, O(s^(90-4n))).
  A product of factors whose nominal valuations add to -10 (every genus-2 decoration, weight 0) is correct to
  O(s^(prec_abs - 12)) at least (the worst factor has nominal valuation +2: Q_{2,0}; P_{2,0} has +3); the drivers choose prec_abs
  accordingly (prec_abs = 24 gives the terms through s^10).
M (class SympM): M(z) = exp(m_1 z + m_2 z^2 + m_3 z^3), constant matrices in the normalised basis, m_1^T = m_1,
  m_2^T = -m_2, m_3^T = m_3 (so M(z)^T M(-z) = 1, checked to z^3); .M[k], .Minv[k] = [z^k] M(z), M(z)^{-1} (k <= 3).
R' = R M^{-1}:  Q'^{(i)}_{a,n} = Sum_{k<=3} Sum_j Q^{(j)}_{a,n-k} (M^{-1}_k)_{ji}          (rows_prime)
Outer edge  E'_{ij}(psi, psi') = (delta_ij - Sum_ab Q'^{(i)}_a(-psi) eta^{ab} Q'^{(j)}_b(-psi'))/(psi + psi'):
  [psi^m psi'^n] = Sum_{k=0}^n (-1)^k G'_{m+1+k, n-k},  G'_{MN} = delta_ij delta_{M0} delta_{N0} - (-1)^{M+N} Sum_ab eta^{ab} Q'_{a,M} Q'_{b,N}.
Outer dilaton leaf  T'(z) = z (1 - R'^{-1}(z) 1),  (R'^{-1}(z) 1)_i = Q'^{(i)}_0(-z)  [R'^{-1}(z) = R'(-z)^*, eta(1, .) = eta(H^0, .)],
  so  T'^{(i)}_d = [z^d] T'_i = (-1)^d Q'^{(i)}_{0,d-1}  (d >= 2);  at M = 1 this is Dil_i(d)/Delta_i^{1/2} (fq_rows.Dil).
  Alternative form  T'^{(i)}_d = -Sum_j Sum_k (M_k)_{ij} (-1)^{d-1-k} Q^{(j)}_{0,d-1-k}  (R'^{-1} = M R^{-1}; checked equal).
API: SRows(rows, C, prec_abs); .Q, .x, .zh (zeta^{4i}), .val_table(), .symplectic_check(nmax); SympM(R, m1, m2, m3) and the
     constructors SympM.identity / random / diagonal;  rows_prime(S, M) -> Q';  EdgeLeaf(S, Qp, M, nmax) with .E(m, n)
     (5x5 list of Ser), .T(d) (list of 5 Ser), .T_alt(d).
"""
import os, sys, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from gsum_q import to_sseries, cf_images_u, u_to_s


class SRows(object):
    def __init__(self, rows, C, prec_abs=None, Lam=0):
        R = rows.R; ring = rows.ring; self.R = R; self.rows = rows; self.Nz = rows.Nz
        imu = cf_images_u(C, R, Lam); self.imu = imu
        f0 = imu['c0']
        self.x = u_to_s(f0).shift(1)                       # Cf(C0)/L = s f0(s^5)
        self.xinv = self.x.inv()
        self.zh = [R.z(4 * i) for i in range(5)]           # Delta_i^{1/2} = zh[i] x
        self.prec_abs = prec_abs
        self.P = [[[None] * (self.Nz + 1) for i in range(5)] for a in range(5)]
        self.Q = [[[None] * (self.Nz + 1) for i in range(5)] for a in range(5)]
        for a in range(5):
            for i in range(5):
                for n in range(self.Nz + 1):
                    Pser = to_sseries(rows.P[a][i][n], ring, imu)
                    if prec_abs is not None: Pser = Pser.truncate(prec_abs + 4)
                    self.P[a][i][n] = Pser
                    Q = Pser.mul(self.xinv).scale(R.inv(self.zh[i]))
                    if prec_abs is not None: Q = Q.truncate(prec_abs)
                    self.Q[a][i][n] = Q
        self.eta_inv = rows.eta_inv
        # Delta check
        self.delta_ok = True
        for i in range(5):
            D = to_sseries(rows.Delta(i), ring, imu)
            D2 = self.x.mul(self.x).scale(R.z(3 * i))
            if not D.eq(D2, min(D.prec, D2.prec) - 1): self.delta_ok = False

    def val_table(self, nmax=4):
        """{(a, n): set of valuations over the colours}"""
        return {(a, n): sorted(set(self.Q[a][i][n].v if self.Q[a][i][n].c else None for i in range(5)))
                for a in range(5) for n in range(nmax + 1)}

    def symplectic_check(self, nmax=4):
        """max |deviation| of Sum_ab Q_a^{(i)}(z) eta^{ab} Q_b^{(j)}(-z) from delta_ij, orders z^0..z^nmax; returns list of failures"""
        R = self.R; bad = []
        for i in range(5):
            for j in range(5):
                for n in range(nmax + 1):
                    acc = None
                    for (a, b), e in self.eta_inv.items():
                        for m in range(n + 1):
                            t = self.Q[a][i][m].mul(self.Q[b][j][n - m]).scale(e * (1 if (n - m) % 2 == 0 else -1))
                            acc = t if acc is None else acc.add(t)
                    target = Ser.const(R, 1 if (i == j and n == 0) else 0, acc.prec)
                    if not acc.eq(target, acc.prec - 1): bad.append((i, j, n, acc.fmt(acc.v + 3)))
        return bad


# ------------------------------------------------------------------ matrices over F_p
def mat_mul(A, B, p):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) % p for j in range(n)] for i in range(n)]
def mat_add(A, B, p, c=1):
    return [[(A[i][j] + c * B[i][j]) % p for j in range(len(A))] for i in range(len(A))]
def mat_scale(A, c, p): return [[A[i][j] * c % p for j in range(len(A))] for i in range(len(A))]
def mat_id(n, p): return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
def mat_zero(n): return [[0] * n for _ in range(n)]
def mat_T(A): return [list(r) for r in zip(*A)]
def mat_eq(A, B, p): return all((A[i][j] - B[i][j]) % p == 0 for i in range(len(A)) for j in range(len(A)))


def exp_zmatrix(ms, p, kmax=3):
    """[z^0..z^kmax] of exp(Sum_k ms[k] z^k) (ms[0] unused/zero), truncated: 1 + m + m^2/2 + m^3/6 + ..."""
    n = len(ms[1])
    m = [mat_zero(n)] + [ms[k] for k in range(1, kmax + 1)]
    out = [mat_id(n, p)] + [mat_zero(n) for _ in range(kmax)]
    power = [mat_id(n, p)] + [mat_zero(n) for _ in range(kmax)]     # m^j as z-series
    fact = 1
    for j in range(1, kmax + 1):
        # power <- power * m (z-series product)
        new = [mat_zero(n) for _ in range(kmax + 1)]
        for a in range(kmax + 1):
            for b in range(1, kmax + 1 - a):
                new[a + b] = mat_add(new[a + b], mat_mul(power[a], m[b], p), p)
        power = new; fact *= j
        inv = pow(fact, p - 2, p)
        for k in range(kmax + 1): out[k] = mat_add(out[k], power[k], p, inv)
    return out


class SympM(object):
    """M(z) = exp(m1 z + m2 z^2 + m3 z^3), constant, symplectic in the normalised basis (m1, m3 symmetric, m2 antisymmetric)"""
    def __init__(self, R, m1, m2, m3, check=True):
        p = R.p; self.R = R; self.p = p
        self.m1 = [[x % p for x in r] for r in m1]; self.m2 = [[x % p for x in r] for r in m2]; self.m3 = [[x % p for x in r] for r in m3]
        assert mat_eq(self.m1, mat_T(self.m1), p) and mat_eq(self.m3, mat_T(self.m3), p), "m1, m3 must be symmetric"
        assert mat_eq(self.m2, mat_scale(mat_T(self.m2), p - 1, p), p), "m2 must be antisymmetric"
        self.M = exp_zmatrix([None, self.m1, self.m2, self.m3], p)
        self.Minv = exp_zmatrix([None, mat_scale(self.m1, p - 1, p), mat_scale(self.m2, p - 1, p), mat_scale(self.m3, p - 1, p)], p)
        if check:
            # M(z)^T M(-z) = 1 and M M^-1 = 1 to z^3
            for k in range(4):
                S1 = mat_zero(5); S2 = mat_zero(5)
                for a in range(k + 1):
                    S1 = mat_add(S1, mat_mul(mat_T(self.M[a]), self.M[k - a], p), p, (1 if (k - a) % 2 == 0 else p - 1))
                    S2 = mat_add(S2, mat_mul(self.M[a], self.Minv[k - a], p), p)
                target = mat_id(5, p) if k == 0 else mat_zero(5)
                assert mat_eq(S1, target, p), "M not symplectic at z^%d" % k
                assert mat_eq(S2, target, p), "M M^-1 != 1 at z^%d" % k
        self.is_identity = all(not any(x for r in m for x in r) for m in (self.m1, self.m2, self.m3))

    @staticmethod
    def identity(R): return SympM(R, mat_zero(5), mat_zero(5), mat_zero(5))

    @staticmethod
    def diagonal(R, a1, a2=None):
        """M = diag(exp(a1_i z + a2_i z^3)); a1, a2 lists of 5 field elements (a2 None: zero)"""
        p = R.p
        m1 = [[R.of(a1[i]) if i == j else 0 for j in range(5)] for i in range(5)]
        m3 = [[(R.of(a2[i]) if a2 is not None else 0) if i == j else 0 for j in range(5)] for i in range(5)]
        return SympM(R, m1, mat_zero(5), m3)

    @staticmethod
    def random(R, rng, with_m3=True, with_m2=True):
        p = R.p
        m1 = mat_zero(5); m2 = mat_zero(5); m3 = mat_zero(5)
        for i in range(5):
            for j in range(i, 5):
                m1[i][j] = m1[j][i] = rng.randrange(p)
                if with_m3: m3[i][j] = m3[j][i] = rng.randrange(p)
                if with_m2 and i != j:
                    x = rng.randrange(p); m2[i][j] = x; m2[j][i] = (-x) % p
        return SympM(R, m1, m2, m3)

    def params(self):
        """the 25 (+15) coordinates: m1 upper triangle (15, row-major i <= j), m2 strict upper (10), m3 upper (15)"""
        v1 = [self.m1[i][j] for i in range(5) for j in range(i, 5)]
        v2 = [self.m2[i][j] for i in range(5) for j in range(i + 1, 5)]
        v3 = [self.m3[i][j] for i in range(5) for j in range(i, 5)]
        return v1, v2, v3


def rows_prime(S, M):
    """Q'^{(i)}_{a,n} = Sum_k Sum_j Q^{(j)}_{a,n-k} (M^-1_k)_{ji}; returns Qp[a][i][n] (Ser), n <= S.Nz"""
    R = S.R; p = R.p
    if M.is_identity: return S.Q
    Qp = [[[None] * (S.Nz + 1) for i in range(5)] for a in range(5)]
    for a in range(5):
        for i in range(5):
            for n in range(S.Nz + 1):
                acc = None
                for k in range(0, min(n, 3) + 1):
                    for j in range(5):
                        c = M.Minv[k][j][i]
                        if not c: continue
                        t = S.Q[a][j][n - k].scale(c)
                        acc = t if acc is None else acc.add(t)
                Qp[a][i][n] = acc
    return Qp


class EdgeLeaf(object):
    """outer edge tables E'(m, n)[i][j] and leaf vectors T'(d)[i] from the rows Q' (see module docstring); cached"""
    def __init__(self, S, Qp, M):
        self.S = S; self.Qp = Qp; self.M = M; self.R = S.R; self._E = {}; self._T = {}

    def G(self, M_, N_, i, j):
        R = self.R
        acc = None
        for (a, b), e in self.S.eta_inv.items():
            t = self.Qp[a][i][M_].mul(self.Qp[b][j][N_]).scale(e)
            acc = t if acc is None else acc.add(t)
        acc = acc.scale(-((-1) ** (M_ + N_)))
        if i == j and M_ == 0 and N_ == 0: acc = acc.add(Ser.const(R, 1, acc.prec))
        return acc

    def E(self, m, n):
        key = (m, n); tab = self._E.get(key)
        if tab is None:
            if m + n + 1 > self.S.Nz: raise ValueError("z-order %d exceeds the rows" % (m + n + 1))
            tab = [[None] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    acc = None
                    for k in range(n + 1):
                        t = self.G(m + 1 + k, n - k, i, j)
                        if k % 2: t = t.neg()
                        acc = t if acc is None else acc.add(t)
                    tab[i][j] = acc
            self._E[key] = tab
        return tab

    def T(self, d):
        """T'^{(i)}_d = (-1)^d Q'^{(i)}_{0,d-1}, d >= 2"""
        r = self._T.get(d)
        if r is None:
            r = [self.Qp[0][i][d - 1].scale((-1) ** d) for i in range(5)]; self._T[d] = r
        return r

    def T_alt(self, d):
        """-Sum_j Sum_k (M_k)_{ij} (-1)^{d-1-k} Q^{(j)}_{0,d-1-k}  (from R'^-1 = M R^-1)"""
        R = self.R; out = []
        for i in range(5):
            acc = None
            for k in range(0, min(d - 1, 3) + 1):
                for j in range(5):
                    c = self.M.M[k][i][j]
                    if not c: continue
                    t = self.S.Q[0][j][d - 1 - k].scale(c * (-1) ** (d - 1 - k) * (-1))
                    acc = t if acc is None else acc.add(t)
            out.append(acc)
        return out
