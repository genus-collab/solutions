# -*- coding: utf-8 -*-
"""gm_qrr.py -- colour-0 WKB polynomials ph_n(u) of the formal quintic, QRR-normalised at q = 0 (over Q or F_p).

Independent re-implementation (no project import) of the recursion in natural_form/referee/ref_rows.py, colour i = 0:
psi_0 = L e^(u_0/z) Phi(zhat, u), Phi = Sum_n zhat^n ph_n(u), ph_0 = 1, A_e = 1 + zhat (1-u)((4/5) theta_zhat - e/5 - theta_u),
Phi_j = A_(-j) Phi_(j-1), (HPF) Phi_5 - Phi_0 - (1-u)[3 zhat Phi_4 + (17/5) u zhat^2 Phi_3 + (9/5) u^2 zhat^3 Phi_2
+ (274/625) u^3 zhat^4 Phi_1 + (24/625) u^4 zhat^5 Phi_0] = 0.  Its zhat^(n+1) coefficient, divided by (1-u) (the factor is
manifest below), reads (4n - 5 theta_u) ph_n + rho_n(u) = 0.  For 5 | n the coefficient of u^(4n/5) is free; it is fixed by
ph_n(1) = h_n (QRR) unless `free` prescribes it.  h(z) = exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_{a in w_0} a^(1-2k)),
w_0 = {1 - zeta^j : j = 1..4} U {-5}; the cyclotomic power sums are computed exactly (cyclo_powersums).
API: hseries(R, N) -> [h_0..h_N];  wkb0(R, N, free=None) -> (ph, info);  R_of_L(R, ph, L) -> [R_n(L) = L^(4n) ph_n(L^-5)].
"""
from fractions import Fraction as Fr
from math import comb
from gm_series import bernoulli_list


def cyclo_powersums(M, N=5):
    """S[m] = Sum_{j=1}^{N-1} (1 - zeta_N^j)^(-m), m = 0..M (exact): from Sum_j 1/(1 + t - zeta^j) = g(t),
    g(t) = N (1+t)^(N-1)/((1+t)^N - 1) - 1/t, so S[m+1] = (-1)^m [t^m] g."""
    K = M + 3
    num = [Fr(N * comb(N - 1, k)) if k <= N - 1 else Fr(0) for k in range(K + 2)]
    den = [Fr(comb(N, k + 1)) if k + 1 <= N else Fr(0) for k in range(K + 2)]
    q = [Fr(0)] * (K + 2)
    for k in range(K + 2):
        q[k] = (num[k] - sum(q[j] * den[k - j] for j in range(k))) / den[0]
    q[0] -= 1
    g = q[1:]
    return [Fr(N - 1)] + [g[m - 1] * (-1) ** (m - 1) for m in range(1, M + 1)]


def hseries(R, N):
    B = bernoulli_list(N + 2); S = cyclo_powersums(N + 2)
    lg = [Fr(0)] * (N + 1)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        lg[m] = -B[2 * k] / (2 * k * (2 * k - 1)) * (S[m] + Fr(-5) ** (-m))
    lg = [R.of(x) for x in lg]; out = [R.zero()] * (N + 1); out[0] = R.one()
    for n in range(1, N + 1):
        out[n] = R.norm(sum(k * lg[k] * out[n - k] for k in range(1, n + 1)) * R.inv(n))
    return out


def _padd(R, a, b, c=1):
    n = max(len(a), len(b)); out = [R.zero()] * n
    for k, x in enumerate(a): out[k] += x
    for k, x in enumerate(b): out[k] += c * x
    out = [R.norm(x) for x in out]
    while out and out[-1] == 0: out.pop()
    return out
def _shift(a, m): return ([0] * m + list(a)) if a else []
def _theta(R, a): return [R.norm(k * x) for k, x in enumerate(a)]
def _scale(R, a, c): return [R.norm(x * c) for x in a] if c else []
def _om(R, a): return _padd(R, a, _shift(a, 1), -1)            # (1 - u) a


def _step(R, prev, nm1, l):
    """(1-u)((4(n-1)+l)/5 - theta) applied to Phi_(l-1)[n-1]"""
    return _om(R, _padd(R, _scale(R, prev, R.of(Fr(4 * nm1 + l, 5))), _theta(R, prev), -1))


def wkb0(R, N, free=None):
    h = hseries(R, N); Ph = [[None] * (N + 1) for _ in range(6)]; ph = []
    info = {'deg': [], 'free': {}, 'h_auto': [], 'res_ok': []}
    for n in range(N + 1):
        if n == 0:
            c = [R.one()]
        else:
            T = [[]]
            for l in range(1, 5): T.append(_padd(R, T[l - 1], _step(R, Ph[l - 1][n - 1], n - 1, l)))
            rho = []
            for l in range(1, 6):
                rho = _padd(R, rho, _padd(R, _scale(R, T[l - 1], R.of(Fr(4 * n + l, 5))), _theta(R, T[l - 1]), -1))
            rho = _padd(R, rho, _scale(R, T[4], R.of(-3)))
            for (j, cf, up) in ((3, Fr(17, 5), 1), (2, Fr(9, 5), 2), (1, Fr(274, 625), 3), (0, Fr(24, 625), 4)):
                if n - up >= 0: rho = _padd(R, rho, _shift(_scale(R, Ph[j][n - up], R.of(-cf)), up))
            fr = 4 * n // 5 if n % 5 == 0 else None
            c = [R.zero()] * max(len(rho), (fr + 1) if fr is not None else 0)
            for k, x in enumerate(rho):
                if k == fr: info['res_ok'].append(x == 0); continue
                c[k] = R.norm(x * R.inv(5 * k - 4 * n))
            if fr is not None:
                c[fr] = R.of(free[n // 5]) if (free is not None and n // 5 in free) else R.norm(h[n] - sum(c))
                info['free'][n // 5] = c[fr]
            else:
                info['h_auto'].append(R.norm(sum(c) - h[n]) == 0)
            while c and c[-1] == 0: c.pop()
        ph.append(c); info['deg'].append(len(c) - 1); Ph[0][n] = c
        for l in range(1, 6):
            Ph[l][n] = Ph[l - 1][n] if n == 0 else _padd(R, Ph[l - 1][n], _step(R, Ph[l - 1][n - 1], n - 1, l))
    return ph, info


def R_of_L(R, ph, L):
    """R_n(L) = L^(4n) ph_n(L^-5) = Sum_k ph_{n,k} L^(4n - 5k)"""
    Lr = R.of(L); Li = R.inv(Lr); out = []
    def pw(e):
        if hasattr(R, 'p'): return pow(Lr, e, R.p) if e >= 0 else pow(Li, -e, R.p)
        return Fr(L) ** e
    for n, P in enumerate(ph):
        out.append(R.norm(sum(x * pw(4 * n - 5 * k) for k, x in enumerate(P))))
    return out
