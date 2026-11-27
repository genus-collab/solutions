# -*- coding: utf-8 -*-
"""ref_rows.py -- REFEREE re-derivation of the hatted rows of the formal quintic DIRECTLY in the variables (zhat, u), with no use of
the ring N, of the s-series or of any fq_* module.

Derivation (referee).  D = q d/dq = -(1-u) d/du, L = 1/s, u = s^5, z = s^4 zhat.  On F = s^e Fh(zhat, u):
      (z D + zeta^i L) F = s^(e-1) A_e Fh,      A_e := zeta^i + zhat (1-u) ( (4/5) theta_zhat - e/5 - theta_u ).
WKB function psi_i = L Phi_i(zhat, u), Phi_i = Sum_n zhat^n ph_n(u), ph_0 = 1.  With Phi_j := A_{-j} Phi_{j-1} (Phi_0 = Phi) the
Picard-Fuchs equation  [(M^5 - 1) - q (3125 M^5 + 9375 z M^4 + 10625 z^2 M^3 + 5625 z^3 M^2 + 1370 z^4 M + 120 z^5)] psi = 0,
M = zD + zeta^i L, 3125 q = 1 - u, becomes after division by u/s
   (HPF)   Phi_5 - Phi_0 - (1-u) [ 3 zhat Phi_4 + (17/5) u zhat^2 Phi_3 + (9/5) u^2 zhat^3 Phi_2 + (274/625) u^3 zhat^4 Phi_1
                                   + (24/625) u^4 zhat^5 Phi_0 ] = 0 .
Its zhat^(n+1) coefficient reads  zeta^(-i) (1-u) (4n - 5 theta_u) ph_n + (terms in ph_0..ph_(n-1)) = 0: ph_n is a POLYNOMIAL in u,
unique unless 5 | n, where the coefficient of u^(4n/5) is free (needs a vanishing residue, checked) and is fixed by
ph_n(u = 1) = [z^n] h_{w_i}(z), h_w(z) = exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_{a in w} a^(1-2k)), w_i = {zeta^i - zeta^j (j != i)} U {-5 zeta^i}
(the q = 0 normalisation; for 5 not | n the identity ph_n(1) = h_n is a CHECK of the whole WKB and of the h-convention).
Rows (normalised canonical basis, Q = P/Delta^(1/2) = zeta^i X_a):  Qh_0 = (zeta^i / C0) Phi_i,
      Qh_(a+1) = kh_a A_(e_a) Qh_a,   e = c - 1 = (-1, -2, 2, 1, 0),   kh = (1/C1, C0^2 C1^2, 1/C1, 1/C0),
C0 = varpi_0 = f_0 + alpha g + beta f_2 (conifold frame), C1 = D t_c, t_c = g/varpi_0, g = the vanishing period (coefficient of
log u in the logarithmic solution, g = u + gamma u^2 + ...), f_k = u^k + O(u^3) the analytic solutions of the PF operator
u D^4 - (1-u)(2 D^3 + (7/5) D^2 + (2/5) D + 24/625).
Teleman form:  Rt_{ji}(zhat; u) = (1/5) Sum_{a+b = 3 mod 5} Qh_{b,0}^{(j)} Qh_{a}^{(i)}(zhat).
API: periods(p, N) -> (f0, f1, f2, g, h) lists;  hcoef(p, zeta, i, nmax);  wkb(p, zeta, i, NZ) -> (list of polys, info);
     HatRows(p, NZ, K, alpha, beta, zeta).Q[a][i][n] (u-series, K terms), .teleman() -> {(n, k): 5x5}, .checks
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import inv, of, find_zeta, smul, sinv, stheta, sscale, sadd, sone_minus, sderiv, bernoulli, mzero, bclean


# ------------------------------------------------------------------ conifold periods by dual numbers (u^eps = 1 + eps log u)
def _dmul(a, b, p): return (a[0] * b[0] % p, (a[0] * b[1] + a[1] * b[0]) % p)
def _dadd(a, b, p): return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def _pf_mono(n, p):
    """PF applied to u^(n+eps): dict {k: dual} of the coefficients of u^(k+eps)"""
    def D(c):
        out = {}
        for k, x in c.items():
            for kk, m, sgn in ((k, k, 1), (k - 1, k, -1)):         # D u^(k+eps) = (k+eps)(u^(k+eps) - u^(k-1+eps))
                t = _dmul(((m * sgn) % p, sgn % p), x, p); out[kk] = _dadd(out.get(kk, (0, 0)), t, p)
        return out
    jets = [{n: (1, 0)}]
    for _ in range(4): jets.append(D(jets[-1]))
    out = {}
    def acc(c, coef, shift):
        for k, x in c.items():
            out[k + shift] = _dadd(out.get(k + shift, (0, 0)), (x[0] * coef % p, x[1] * coef % p), p)
    acc(jets[4], 1, 1)                                               # u D^4
    for j, cf in ((3, Fr(2)), (2, Fr(7, 5)), (1, Fr(2, 5)), (0, Fr(24, 625))):
        c = of(cf, p); acc(jets[j], (-c) % p, 0); acc(jets[j], c, 1)  # -(1-u) cf D^j
    return out


def periods(p, N):
    """f0, f1, f2 (analytic, f_k = u^k + O(u^3)), g (vanishing period u + gamma u^2 + ..), h (g log u + h solves PF): lists mod u^N"""
    mono = {n: _pf_mono(n, p) for n in range(0, N + 1)}
    def solve(b):                                                    # b: dict n -> dual for n = 0, 1, 2
        for k in range(-3, 0):
            r = (0, 0)
            for n in range(0, 3): r = _dadd(r, _dmul(b[n], mono[n].get(k, (0, 0)), p), p)
            assert r == (0, 0), ("obstruction at u^%d" % k, r)
        for n in range(3, N):
            k = n - 3; r = (0, 0)
            for m in range(max(0, k - 1), n): r = _dadd(r, _dmul(b[m], mono[m].get(k, (0, 0)), p), p)
            ind = mono[n][k]
            assert ind[0] == n * (n - 1) ** 2 * (n - 2) % p
            i0 = inv(ind[0], p); iinv = (i0, (-ind[1] * i0 * i0) % p)
            b[n] = _dmul(((-r[0]) % p, (-r[1]) % p), iinv, p)
        return b
    fs = []
    for k in range(3):
        b = solve({n: (0, 1 if n == k else 0) for n in range(3)})
        assert all(b[n][0] == 0 for n in range(N)); fs.append([b[n][1] for n in range(N)])
    x = mono[1][-1]; y = mono[2][-1]
    assert x[0] == 0 and y == (0, 2), (x, y)
    gamma = (-x[1] * inv(2, p)) % p
    b = solve({0: (0, 0), 1: (1, 0), 2: (gamma, 0)})
    return fs[0], fs[1], fs[2], [b[n][0] for n in range(N)], [b[n][1] for n in range(N)]


# ------------------------------------------------------------------ QRR constants
def hcoef(p, zeta, i, nmax):
    """[z^n] h_{w_i}(z), n = 0..nmax"""
    zi = pow(zeta, i, p)
    w = [(zi - pow(zeta, j, p)) % p for j in range(5) if j != i] + [(-5 * zi) % p]
    lg = [0] * (nmax + 1)
    for k in range(1, nmax // 2 + 2):
        if 2 * k - 1 > nmax: break
        ps = sum(pow(inv(a, p), 2 * k - 1, p) for a in w) % p
        lg[2 * k - 1] = (-of(bernoulli(2 * k) / (2 * k * (2 * k - 1)), p) * ps) % p
    out = [0] * (nmax + 1); out[0] = 1                               # exp: n out_n = Sum k lg_k out_(n-k)
    for n in range(1, nmax + 1):
        out[n] = sum(k * lg[k] * out[n - k] for k in range(1, n + 1)) % p * inv(n, p) % p
    return out


# ------------------------------------------------------------------ polynomials in u and the hatted WKB
def _padd(a, b, p, c=1):
    n = max(len(a), len(b)); out = [0] * n
    for k, x in enumerate(a): out[k] = x
    for k, x in enumerate(b): out[k] = (out[k] + c * x) % p
    while out and out[-1] == 0: out.pop()
    return out
def _pshift(a, m): return ([0] * m + list(a)) if a else []
def _pone_minus(a, p): return _padd(a, _pshift(a, 1), p, p - 1)
def _ptheta(a, p): return [k * x % p for k, x in enumerate(a)]
def _pscale(a, c, p): return [x * c % p for x in a] if c % p else []


def _A(e, F, zi, p, i5, theta, one_minus, add, scale):
    """A_e on a zhat-series F (list over n): out[n] = zeta^i F[n] + (1-u)((4(n-1) - e)/5 - theta_u) F[n-1]"""
    out = []
    for n in range(len(F)):
        t = scale(F[n], zi, p)
        if n >= 1:
            c = (4 * (n - 1) - e) * i5 % p
            t = add(t, one_minus(add(scale(F[n - 1], c, p), theta(F[n - 1], p), p, p - 1), p), p)
        out.append(t)
    return out


def hpf(Phi, zi, p):
    """the left-hand side of (HPF) for the zhat-series of polynomials Phi (list over n); returns the list over n"""
    i5 = inv(5, p); N = len(Phi); Ph = [Phi]
    for j in range(1, 6): Ph.append(_A(-j, Ph[-1], zi, p, i5, _ptheta, _pone_minus, _padd, _pscale))
    out = []
    for n in range(N):
        t = _padd(Ph[5][n], Ph[0][n], p, p - 1); br = []
        for m, cf, upow, j in ((1, Fr(3), 0, 4), (2, Fr(17, 5), 1, 3), (3, Fr(9, 5), 2, 2), (4, Fr(274, 625), 3, 1), (5, Fr(24, 625), 4, 0)):
            if n - m >= 0: br = _padd(br, _pshift(_pscale(Ph[j][n - m], of(cf, p), p), upow), p)
        out.append(_padd(t, _pone_minus(br, p), p, p - 1))
    return out


def wkb(p, zeta, i, NZ):
    """ph_n(u), n = 0..NZ, for the colour i; info: degrees, free coefficients, check ph_n(1) = h_n"""
    zi = pow(zeta, i, p); h = hcoef(p, zeta, i, NZ); Phi = [[1]]; info = dict(deg=[0], free=[], h_auto=[])
    for n in range(1, NZ + 1):
        r = hpf(Phi + [[], []], zi, p)[n + 1]                        # zhat^(n+1) coefficient with ph_n = 0 (ph_(n+1) does not enter)
        assert sum(r) % p == 0, "not divisible by 1 - u at order %d" % n
        rho = []; s = 0
        for x in r[:-1]: s = (s + x) % p; rho.append(s * zi % p)     # rho = zeta^i r/(1-u)
        free = 4 * n // 5 if n % 5 == 0 else None
        c = [0] * max(len(rho), 1, (free + 1) if free is not None else 0)
        for k, x in enumerate(rho):
            if k == free: assert x == 0, "resonance obstruction at n = %d" % n
            else: c[k] = x * inv(5 * k - 4 * n, p) % p
        if free is not None:
            c[free] = (h[n] - sum(c)) % p; info['free'].append((n, free))
        else: info['h_auto'].append(sum(c) % p == h[n])
        while c and c[-1] == 0: c.pop()
        Phi.append(c); info['deg'].append(len(c) - 1)
        chk = hpf(Phi + [[]], zi, p)[n + 1]
        assert not chk, "HPF fails at zhat^%d" % (n + 1)
    return Phi, info
