# -*- coding: utf-8 -*-
"""rf_core.py -- referee's independent model data for the frozen theory V = M_0.omega^ (LG model x^4 + v/x, form dx, Z/5 line).

Nothing here is imported from the worker (allgenus/model/ms_*).  Derivation (written out in REFEREE.md s.2):
  flat basis H^a (a in Z/5), c = (0,-1,3,2,1), N = (1,-1,-1,1,1), metric eta(H^a,H^b) = 5 N_a N_b [a+b = 3 mod 5],
  product H * H^b = beta_b H^(b+1), beta = (1, -1/v, 1, -1, -1).  g^a := N_a v^(-c_a/5) H^a satisfy g^a g^b = g^(a+b): the group algebra
  of Z/5; idempotents e_i = (1/5) Sum_a zeta^(-ia) g^a, eta(e_i,e_i) = v^(-2/5) zeta^(-3i); normalised e~_i = v^(1/5) zeta^(4i) e_i.
  At v = 1:  (Psi^-1 H^a)_i = N_a zeta^(i(a+1)),  unit 1 = Sum_i zeta^(-4i) e~_i,  sqrt(Delta_i) = zeta^(4i),  U = diag(-(5/4) zeta^i).
  R-matrix: S = Psi R(z) e^(U/z) solves z d/dv S = (H *) S; homogeneity R_k(v) = r_k v^(-4k/5) gives, with Z = diag(zeta^i),
  V = Psi^-1 mu Psi, V_ij = (1/5) Sum_a mu_a zeta^((i-j)(a+1)), mu_a = (c_a - 1)/4:
      (r_(k+1))_ij = -(4/5) ((V + k) r_k)_ij / (zeta^j - zeta^i)   (i != j),     (k+1) (r_(k+1))_ii = - Sum_j V_ij (r_(k+1))_ji.
  The sign convention of z is NOT assumed: rf_regress.py tests both R(z) and R(-z) against the numbers on record.
psi-class intersection numbers (WK): own implementation of the string, dilaton and Virasoro recursions (exact rationals).
Arithmetic mod p (p = 2^61-1 or 2^89-1; both are = 1 mod 5).  Run as a script for self-checks.
"""
import sys
from fractions import Fraction as Fr
from functools import lru_cache

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1
C_EXP = (0, -1, 3, 2, 1)
N_SIGN = (1, -1, -1, 1, 1)


def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r


@lru_cache(maxsize=None)
def wk(g, ds):
    """<tau_d1 ... tau_dn>_g (psi-class intersection numbers), ds a sorted tuple"""
    n = len(ds)
    if g < 0 or 2 * g - 2 + n <= 0: return Fr(0)
    if sum(ds) != 3 * g - 3 + n or min(ds) < 0: return Fr(0)
    if g == 0 and n == 3: return Fr(1)
    if g == 1 and n == 1: return Fr(1, 24)
    ds = list(ds)
    if 0 in ds:                                   # string equation
        ds.remove(0); tot = Fr(0)
        for j in range(len(ds)):
            if ds[j] > 0:
                e = ds[:]; e[j] -= 1; tot += wk(g, tuple(sorted(e)))
        return tot
    if 1 in ds:                                   # dilaton equation
        ds.remove(1)
        return (2 * g - 2 + len(ds)) * wk(g, tuple(sorted(ds)))
    k = ds[-1] - 1; rest = ds[:-1]; m = len(rest); tot = Fr(0)   # Virasoro recursion on the largest insertion tau_(k+1), k >= 1
    for j in range(m):
        e = rest[:]; e[j] += k
        tot += Fr(dfact(2 * k + 2 * rest[j] + 1), dfact(2 * rest[j] - 1)) * wk(g, tuple(sorted(e)))
    for r in range(k):
        s = k - 1 - r; c = Fr(dfact(2 * r + 1) * dfact(2 * s + 1), 2)
        tot += c * wk(g - 1, tuple(sorted(rest + [r, s])))
        for mask in range(1 << m):
            I = [rest[i] for i in range(m) if mask >> i & 1]; J = [rest[i] for i in range(m) if not mask >> i & 1]
            for g1 in range(g + 1):
                a = wk(g1, tuple(sorted(I + [r])))
                if a: tot += c * a * wk(g - g1, tuple(sorted(J + [s])))
    return tot / dfact(2 * k + 3)


class Model(object):
    """frozen model at v = 1 mod p; sign = +1: R(z) from the recursion, sign = -1: R(-z)"""

    def __init__(self, p, K, sign=1):
        self.p = p; self.K = K; self.sign = sign
        z = None
        for a in range(2, 200):
            x = pow(a, (p - 1) // 5, p)
            if x != 1: z = x; break
        self.zeta = z
        self.zp = [pow(z, k, p) for k in range(5)]
        inv = lambda x: pow(x % p, p - 2, p)
        self.inv = inv
        mu = [Fr(C_EXP[a] - 1, 4) for a in range(5)]
        V = [[sum(self.fr(mu[a]) * self.zp[((i - j) * (a + 1)) % 5] for a in range(5)) * inv(5) % p for j in range(5)] for i in range(5)]
        self.V = V
        r = [[[1 if i == j else 0 for j in range(5)] for i in range(5)]]
        for k in range(K):
            A = [[(sum(V[i][l] * r[k][l][j] for l in range(5)) + k * r[k][i][j]) % p for j in range(5)] for i in range(5)]
            nr = [[0] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    if i != j: nr[i][j] = (-4) * A[i][j] * inv(5 * (self.zp[j] - self.zp[i])) % p
            for i in range(5):
                nr[i][i] = (-sum(V[i][j] * nr[j][i] for j in range(5) if j != i)) * inv(k + 1) % p
            r.append(nr)
        if sign == -1: r = [[[(-1) ** k * x % p for x in row] for row in r[k]] for k in range(K + 1)]
        self.R = r
        self.Rinv = self.series_inverse(r)
        self.unit = [self.zp[(-4 * i) % 5] for i in range(5)]
        self.D = [self.zp[(4 * i) % 5] for i in range(5)]                 # sqrt(Delta_i)
        self.flat = [[N_SIGN[a] * self.zp[(i * (a + 1)) % 5] % p for i in range(5)] for a in range(5)]   # (Psi^-1 H^a)_i
        # leaves T_e = -R^-1_(e-1) 1 (e >= 2), edges E[(k,l)] = Sum_t (-1)^t G[(k+1+t, l-t)], G[(M,N)] = -A_M A_N^T
        self.T = {e: [(-sum(self.Rinv[e - 1][i][j] * self.unit[j] for j in range(5))) % p for i in range(5)] for e in range(2, K + 2)}
        G = {}
        for M in range(K + 1):
            for N in range(K + 1 - M):
                if M == 0 and N == 0: continue
                G[(M, N)] = [[(-sum(self.Rinv[M][i][l] * self.Rinv[N][j][l] for l in range(5))) % p for j in range(5)] for i in range(5)]
        self.E = {}
        for k in range(K):
            for l in range(K - k):
                tab = [[0] * 5 for _ in range(5)]
                for t in range(l + 1):
                    Gt = G[(k + 1 + t, l - t)]; s = 1 if t % 2 == 0 else p - 1
                    for i in range(5):
                        for j in range(5): tab[i][j] = (tab[i][j] + s * Gt[i][j]) % p
                self.E[(k, l)] = tab

    def fr(self, x):
        x = Fr(x); return x.numerator % self.p * pow(x.denominator % self.p, self.p - 2, self.p) % self.p

    def series_inverse(self, r):
        p = self.p; K = len(r) - 1; out = [[[1 if i == j else 0 for j in range(5)] for i in range(5)]]
        for k in range(1, K + 1):
            acc = [[0] * 5 for _ in range(5)]
            for j in range(1, k + 1):
                for a in range(5):
                    for b in range(5):
                        acc[a][b] = (acc[a][b] - sum(r[j][a][c] * out[k - j][c][b] for c in range(5))) % p
            out.append(acc)
        return out

    def symplectic_defect(self):
        """max over k of the z^k coefficient of R(-z)^T R(z) - 1 (normalised canonical frame, eta = 1); 0 means symplectic"""
        p = self.p; bad = 0
        for k in range(1, self.K + 1):
            for i in range(5):
                for j in range(5):
                    s = sum((-1) ** a * self.R[a][l][i] * self.R[k - a][l][j] for a in range(k + 1) for l in range(5)) % p
                    if s: bad += 1
        return bad


if __name__ == "__main__":
    assert wk(1, (1,)) == Fr(1, 24) and wk(2, (4,)) == Fr(1, 1152) and wk(3, (7,)) == Fr(1, 82944)
    assert wk(2, (2, 3)) == Fr(29, 5760) and wk(1, (1, 1)) == Fr(1, 24) and wk(0, (0, 0, 0, 1)) == 1
    assert wk(4, (10,)) == Fr(1, 7962624) and wk(2, (2, 2, 2)) == Fr(7, 240)
    print("WK self-checks passed")
    for p in (P61, P89):
        m = Model(p, 12)
        print("p = 2^%d-1: zeta^5 = %d, symplectic defect of R through z^12: %d" % (p.bit_length(), pow(m.zeta, 5, p), m.symplectic_defect()))
