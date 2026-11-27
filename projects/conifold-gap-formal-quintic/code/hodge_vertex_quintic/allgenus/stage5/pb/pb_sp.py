# -*- coding: utf-8 -*-
"""pb_sp.py -- Stage-5 'pb' audit: the full-Stirling stationary phase, computed in the FLAT variable sigma = z s at a numerical
saddle L in F_p (p = 1 mod 5, explicit zeta in F_p), then interpolated in u = L^-5.  Imports only pb_core.

Summand (POLY s.5 form): q^s Gamma(1+5s) / prod_(j=0..4) Gamma(1+s-zeta^j/z), every Gamma by Stirling:
  log Gamma(1+y) = (y+1/2) log y - y + (1/2) log 2pi + sum_k b_k y^(1-2k),  b_k = B_2k/(2k(2k-1)).
Phase S(sigma)/z, S = sigma log q + 5 sigma log(5 sigma) - sum_j (sigma - zeta^j) log(sigma - zeta^j); only S^(k)(L), k >= 2 enter:
  S^(k)(L) = (-1)^k (k-2)! [5/L^(k-1) - sum_j (L - zeta^j)^(1-k)].
Amplitude a = sqrt(5 sigma/(sigma^5-1)) exp(sum_k b_k z^(2k-1) [(5 sigma)^(1-2k) - sum_j (sigma - zeta^j)^(1-2k)]), normalised by its
z^0 value at the saddle.  Formal SP by Lagrange-Buermann: with F = S(L+t) - S(L) = (S2/2) t^2 Q(t),
  c(z;L) = sum_m (2m-1)!! (z/(-S2))^m [t^2m] ( g(t,z) Q(t)^(-m-1/2) ).
Hatted dictionary (to be TESTED, not assumed, against the gauge): ph_n(u) = L^(-4n) c_n(L), u = L^(-5).
Variants for controls: drop5 (no Stirling term for Gamma(1+5s)), flip (b_k -> -b_k).
(GAMMA's summand Gamma(1+5s) prod_j Gamma(zeta^j/z - s) (-q)^s has the same S^(k), k >= 2, and the same normalised amplitude:
d^2/dsigma^2 of (c - sigma) log(c - sigma) and of -(sigma - c) log(sigma - c) are both 1/(c - sigma), and the Bernoulli terms
b_k ((c - sigma)/z)^(1-2k) = -b_k ((sigma - c)/z)^(1-2k) agree; so it is the same computation and is not coded twice.)
"""
from fractions import Fraction as Fr
import pb_core as C

P61 = 2 ** 61 - 1
P89 = 2 ** 89 - 1


def zeta5(p):
    for a in range(2, 100):
        z = pow(a, (p - 1) // 5, p)
        if z != 1:
            assert pow(z, 5, p) == 1
            return z
    raise RuntimeError


def ser_mul(a, b, n, p):
    res = [0] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]):
                res[i + j] += x * y
    return [x % p for x in res]


def ser_pow_unit(a, r, n, p):
    """a^r for a[0] = 1 (r a Fraction), first n coefficients."""
    rr = r.numerator * pow(r.denominator, p - 2, p) % p
    b = [1]
    for m in range(1, n):
        s = 0
        for k in range(1, min(m, len(a) - 1) + 1):
            s += (rr * k - (m - k)) * a[k] * b[m - k]
        b.append(s % p * pow(m, p - 2, p) % p)
    return b


def binom_series(x, e, n, p):
    """(1 + x t)^e, e an integer, first n coefficients."""
    c = [1]
    for i in range(n - 1):
        c.append(c[-1] * ((e - i) % p) % p * pow(i + 1, p - 2, p) % p * x % p)
    return c


def sp_flat(L, nmax, p, bern, variant=''):
    """c_0..c_nmax of the normalised full-Stirling saddle series at sigma = L (L in F_p)."""
    zt = zeta5(p); Z = [pow(zt, j, p) for j in range(5)]
    T = 2 * nmax + 1
    inv = lambda x: pow(x % p, p - 2, p)
    dj = [(L - Z[j]) % p for j in range(5)]
    # Taylor coefficients S^(k)(L)/k!, k = 2..T+1
    Sk = {}
    fact = 1
    for k in range(2, T + 3):
        fact = fact * k % p
        sgn = 1 if k % 2 == 0 else -1
        fk2 = 1
        for i in range(2, k - 1):
            fk2 = fk2 * i % p
        val = 5 * inv(pow(L, k - 1, p)) - sum(inv(pow(d, k - 1, p)) for d in dj)
        Sk[k] = sgn * fk2 * val % p * inv(fact) % p
    S2 = 2 * Sk[2] % p                      # S''(L)
    a = (-S2) % p
    Qser = [Sk[k + 2] * inv(Sk[2]) % p for k in range(T + 1)]  # Q(t) = F/(S''t^2/2)
    Qm = ser_pow_unit(Qser, Fr(-1, 2), T + 1, p)
    Qinv = ser_pow_unit(Qser, Fr(-1), T + 1, p)
    # amplitude z^0 part, normalised: sqrt((L+t)/L) * sqrt((L^5-1)/((L+t)^5-1))
    g0 = binom_half([1, inv(L)], T + 1, p)
    den = [1] + [0] * T
    poly = [(pow(L, 5, p) - 1) % p, 5 * pow(L, 4, p), 10 * pow(L, 3, p), 10 * L * L, 5 * L, 1]
    i0 = inv(poly[0])
    poly = [x * i0 % p for x in poly]
    g0 = ser_mul(g0, ser_pow_unit(poly, Fr(-1, 2), T + 1, p), T + 1, p)
    # Stirling exponent G_r(t), r = 2k-1
    G = {}
    for k in range(1, nmax // 2 + 2):
        r = 2 * k - 1
        if r > nmax:
            break
        bk = bern[2 * k] / (2 * k * (2 * k - 1))
        bkp = bk.numerator % p * inv(bk.denominator) % p
        if 'flip' in variant:
            bkp = (-bkp) % p
        term = [0] * (T + 1)
        if 'drop5' not in variant:
            c5 = inv(pow(5 * L % p, 2 * k - 1, p))                 # (5L)^(1-2k)
            for i, x in enumerate(binom_series(inv(L), 1 - 2 * k, T + 1, p)):
                term[i] += c5 * x
        for j in range(5):
            cj = inv(pow(dj[j], 2 * k - 1, p))                     # (L - zeta^j)^(1-2k)
            for i, x in enumerate(binom_series(inv(dj[j]), 1 - 2 * k, T + 1, p)):
                term[i] -= cj * x
        G[r] = [bkp * x % p for x in term]
    # E = exp(sum_r G_r z^r): r E_r = sum_(odd j <= r) j G_j E_(r-j)
    E = [g0]
    for r in range(1, nmax + 1):
        acc = [0] * (T + 1)
        for j in range(1, r + 1, 2):
            if j in G:
                prod = ser_mul(G[j], E[r - j], T + 1, p)
                for i in range(T + 1):
                    acc[i] += j * prod[i]
        ir = inv(r)
        E.append([x % p * ir % p for x in acc])
    # Q^(-m-1/2), m = 0..nmax
    Qp = [Qm]
    for m in range(1, nmax + 1):
        Qp.append(ser_mul(Qp[-1], Qinv, T + 1, p))
    out = []
    ia = inv(a)
    for n in range(nmax + 1):
        s = 0; dfact = 1; apow = 1
        for m in range(0, n + 1):
            if m >= 1:
                dfact = dfact * (2 * m - 1) % p; apow = apow * ia % p
            e = E[n - m]; q = Qp[m]
            inner = sum(e[i] * q[2 * m - i] for i in range(2 * m + 1)) % p
            s += dfact * apow % p * inner
        out.append(s % p)
    return out


def binom_half(lin, n, p):
    """(1 + c t)^(1/2) for lin = [1, c]."""
    return ser_pow_unit(lin + [0] * (n - 2), Fr(1, 2), n, p)


def hatted_values(L, nmax, p, bern, variant=''):
    """ph_n(u) at u = L^-5 from the flat series: ph_n = L^(-4n) c_n(L)."""
    c = sp_flat(L, nmax, p, bern, variant)
    iL4 = pow(pow(L, 4, p), p - 2, p)
    out = []; f = 1
    for n in range(nmax + 1):
        out.append(c[n] * f % p); f = f * iL4 % p
    return (pow(pow(L, 5, p), p - 2, p), out)


def interpolate(xs, ys, p):
    """coefficients of the unique polynomial of degree < len(xs) through the points (Newton form, then expanded)."""
    n = len(xs); dd = list(ys)
    coef = [dd[0]]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) * pow((xs[i] - xs[i - j]) % p, p - 2, p) % p
        coef.append(dd[j])
    poly = [0] * n
    for j in range(n - 1, -1, -1):
        # poly = poly * (x - xs[j]) + coef[j]
        new = [0] * n
        for i in range(n - 1):
            new[i + 1] = (new[i + 1] + poly[i]) % p
            new[i] = (new[i] - xs[j] * poly[i]) % p
        new[0] = (new[0] + coef[j]) % p
        poly = new
    return poly


def h_qrr(nmax, p, bern):
    """h_(w_0)(z) mod p, w_0 = {1 - zeta^j, j = 1..4} U {-5}: [z^0..z^nmax]."""
    zt = zeta5(p); W = [(1 - pow(zt, j, p)) % p for j in range(1, 5)] + [(-5) % p]
    ex = [0] * (nmax + 1)
    for k in range(1, nmax // 2 + 2):
        r = 2 * k - 1
        if r > nmax:
            break
        bk = bern[2 * k] / (2 * k * (2 * k - 1))
        s = sum(pow(pow(a, r, p), p - 2, p) for a in W) % p
        ex[r] = (-(bk.numerator % p) * pow(bk.denominator, p - 2, p) % p * s) % p
    out = [1] + [0] * nmax
    for n in range(1, nmax + 1):
        out[n] = sum(j * ex[j] * out[n - j] for j in range(1, n + 1)) % p * pow(n, p - 2, p) % p
    return out
