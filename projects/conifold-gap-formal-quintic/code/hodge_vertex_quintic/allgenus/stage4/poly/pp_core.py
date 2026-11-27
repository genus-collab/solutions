# -*- coding: utf-8 -*-
"""pp_core.py -- Stage 4 / poly: exact (over Q) hatted WKB polynomials ph_n(u) of the formal quintic, colour 0.

Independent re-derivation (no project module imported).  D = q d/dq = -(1-u) d/du, 3125 q = 1 - u, s = u^(1/5), z = s^4 w
(w = zhat).  On F = s^e Fh(w, u):  (z D + L) F = s^(e-1) A_e Fh,  A_e = 1 + w (1-u) ((4/5) theta_w - e/5 - theta_u)   (colour 0, L = 1/s).
psi = s^-1 Phi, Phi_j = A_(-j) Phi_(j-1); PF  u M^5 - 1 - (1-u)(3 z M^4 + 17/5 z^2 M^3 + 9/5 z^3 M^2 + 274/625 z^4 M + 24/625 z^5)
gives (HPF):  Phi_5 - Phi_0 - (1-u)[3 w Phi_4 + 17/5 u w^2 Phi_3 + 9/5 u^2 w^3 Phi_2 + 274/625 u^3 w^4 Phi_1 + 24/625 u^4 w^5 Phi_0] = 0.
Colour i: Phi^(i)(w, u) = Phi^(0)(w/zeta^i, u) (HPF depends on w/zeta^i only; h_(w_i)(z) = h_(w_0)(z/zeta^i)).
[w^(n+1)] HPF = (1-u) [ (4n - 5 theta_u) ph_n + rest_n ]:  ph_n polynomial, free coefficient [u^(4m)] ph_(5m).
Polynomials: lists of Fractions, index = power of u.
"""
from fractions import Fraction as Fr
from math import comb

E5 = (Fr(3), Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625))   # e_k(1/5,..,5/5), k = 1..5


# ------------------------------------------------------------------ polynomial helpers
def pclean(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a
def padd(a, b, c=1):
    n = max(len(a), len(b)); out = [Fr(0)] * n
    for k, x in enumerate(a): out[k] += x
    for k, x in enumerate(b): out[k] += c * x
    return pclean(out)
def pscale(a, c): return pclean([x * c for x in a]) if c else []
def pshift(a, m): return ([Fr(0)] * m + list(a)) if a else []
def ptheta(a): return pclean([k * x for k, x in enumerate(a)])
def pone_minus(a): return padd(a, pshift(a, 1), -1)                   # (1-u) a
def pdiv_one_minus(a):
    """a/(1-u) for a divisible by (1-u)"""
    out = []; s = Fr(0)
    for x in a[:-1] if a else []:
        s += x; out.append(s)
    assert (s + (a[-1] if a else 0)) == 0, "not divisible by 1-u"
    return pclean(out)
def pval(a, x):
    r = Fr(0)
    for c in reversed(a): r = r * x + c
    return r
def pmul(a, b):
    if not a or not b: return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] += x * y
    return pclean(out)


# ------------------------------------------------------------------ Bernoulli numbers and the QRR series h_(w_0)
def bernoulli_list(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


def cyclo_power_sums(M):
    """S[m] = Sum_(j=1..4) (1 - zeta^j)^(-m), m = 0..M (Newton; y_j = 1/(1-zeta^j) roots of y^4 - 2y^3 + 2y^2 - y + 1/5)"""
    e = [Fr(1), Fr(2), Fr(2), Fr(1), Fr(1, 5)]
    S = [Fr(4)] + [Fr(0)] * M
    for m in range(1, M + 1):
        if m <= 4: S[m] = sum((-1) ** (k - 1) * e[k] * S[m - k] for k in range(1, m)) + (-1) ** (m - 1) * m * e[m]
        else: S[m] = sum((-1) ** (k - 1) * e[k] * S[m - k] for k in range(1, 5))
    return S


def hcoef(nmax, B=None):
    """[z^n] h_(w_0)(z), w_0 = {1 - zeta^j, j=1..4} U {-5}:  h = exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_(a in w0) a^(1-2k))"""
    if B is None: B = bernoulli_list(nmax + 2)
    S = cyclo_power_sums(nmax + 1)
    lg = [Fr(0)] * (nmax + 1)
    for k in range(1, nmax // 2 + 2):
        m = 2 * k - 1
        if m > nmax: break
        ps = S[m] + Fr(-5) ** (-m)
        lg[m] = -B[2 * k] / (2 * k * (2 * k - 1)) * ps
    out = [Fr(1)] + [Fr(0)] * nmax
    for n in range(1, nmax + 1):
        out[n] = sum(k * lg[k] * out[n - k] for k in range(1, n + 1)) / n
    return out, lg


# ------------------------------------------------------------------ the WKB recursion
class WKB(object):
    """ph[n] (polynomials in u), colour 0.  free: 'qrr' (ph_n(1) = h_n), 'zero' (free coefficient 0), or a dict m -> value of
    [u^(4m)] ph_(5m)."""
    def __init__(self, NZ, free='qrr'):
        self.NZ = NZ; self.h, self.lg = hcoef(NZ + 1)
        self.ph = [[Fr(1)]]; self.Ph = [[[Fr(1)]] for _ in range(6)]   # Ph[j][m] = [w^m] Phi_j
        self.free = {}; self.h_auto = []; self.rest = [None]
        for n in range(1, NZ + 1): self._step(n, free)

    def _phij(self, j, m, ph_m):
        """[w^m] Phi_j given [w^m] Phi_(j-1) (= ph_m if j = 1 chain) and stored [w^(m-1)] Phi_(j-1)"""
        prev = self.Ph[j - 1][m - 1]
        c = Fr(4 * (m - 1) + j, 5)                                     # e = -j: (4(m-1) - e)/5
        return padd(ph_m, pone_minus(padd(pscale(prev, c), ptheta(prev), -1)))

    def _row(self, n, phn):
        """[w^n] Phi_j for j = 0..5 given ph_n (rows up to w^(n-1) stored at Ph[j][0..n-1])"""
        out = [phn]
        for j in range(1, 6): out.append(self._phij(j, n, out[-1]))
        return out

    def _row_prov(self, n):
        return self._row(n, [])

    def _lhs_coef(self, n):
        """[w^(n+1)] of HPF divided by (1-u), with ph_(n+1) := 0 (it cancels) and the stored rows up to w^n"""
        P5 = []                                                         # [w^(n+1)] Phi_5 with ph_(n+1) = 0
        cur = []
        for j in range(1, 6):
            prev = self.Ph[j - 1][n]
            c = Fr(4 * n + j, 5)
            cur = padd(cur, pone_minus(padd(pscale(prev, c), ptheta(prev), -1)))
        P5 = cur
        br = []
        for mm, j in ((1, 4), (2, 3), (3, 2), (4, 1), (5, 0)):
            if n + 1 - mm >= 0:
                br = padd(br, pshift(pscale(self.Ph[j][n + 1 - mm], E5[mm - 1]), mm - 1))
        tot = padd(P5, pone_minus(br), -1)
        return pdiv_one_minus(tot)

    def _step(self, n, free):
        row0 = self._row_prov(n)                                       # provisional rows with ph_n = 0
        for j in range(6): self.Ph[j].append(row0[j])
        r0 = self._lhs_coef(n)                                         # = (4n - 5 theta) ph_n + rest; here ph_n = 0
        rest = r0; self.rest.append(rest)
        res = (4 * n // 5) if n % 5 == 0 else None
        c = [Fr(0)] * max(len(rest), (res + 1) if res is not None else 0)
        for k, x in enumerate(rest):
            if k == res: assert x == 0, "resonance obstruction n=%d" % n
            else: c[k] = -x / (4 * n - 5 * k)
        if res is not None:
            if free == 'qrr': c[res] = self.h[n] - sum(c)
            elif free == 'zero': c[res] = Fr(0)
            else: c[res] = Fr(free.get(n // 5, 0))
            self.free[n // 5] = c[res]
        else:
            self.h_auto.append((n, sum(c) == self.h[n]))
        c = pclean(c); self.ph.append(c)
        row = self._row(n, c)
        for j in range(6): self.Ph[j][n] = row[j]
        assert not self._lhs_coef_check(n), "HPF fails at w^%d" % (n + 1)

    def _lhs_coef_check(self, n):
        return self._lhs_coef(n)


if __name__ == "__main__":
    W = WKB(12)
    for n, p in enumerate(W.ph): print(n, len(p) - 1, p[:3], "...", p[-1] if p else 0)
    print("free", W.free); print("h_auto", W.h_auto)
