"""lf_lmi2 N: Theorem LMI in Method A (sigma = s* + eps + eps tau, z = eps w, over Q[eps]/(eps^(N+1))), exact over Q.

M1  the quintic as LMI data at s* = 1 (alpha0 = -1, M0 = 0; g, F = the five regular factors): c_n has no pole at
    L = 1 and c(z; 1) = h (h from its definition in Q(zeta_5)), n <= N.
M2  the quintic as LMI data at s* = 0 (alpha0 = 5, M0 = 1; g, F = the five factors Gamma((zeta^j - sigma)/z)):
    no pole at L = 0 and c(z; 0) = Ginf, n <= N.
M3  non-polynomial data: g_k = r^k/k! (g entire), log F with geometric coefficients and z-dependence,
    M0 in {-1, 0, 1, 2}: no pole, value = the LMI formula.
M4  beyond the theorem (RA_2 s.5 'alignment'): singular factor Gamma(h + a0/z) with NON-INTEGER h, Stirling series
    by DLMF 5.11.8 (Bernoulli polynomials): no pole and value = LMI formula with delta = -h z/alpha0
    (h = 1/5, 1/3, -1/2, 7/4). Also: for integer h = M the Bernoulli-polynomial amplitude IS the Poch form.
Failing controls: singular b_1 + 1/7 (M1); M0 = 0 in place of M0 = 1 (M2, value); B_k(0) in place of B_k(h) (M4).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
from lf_core import bern_poly, bk, gbinom, ser_exp, TS
from lf_wick import wick_1
from lf_models import q_chart1, q_chart0, lmi_scaled, lmi_predict, _sing_amp
from lf_log import Log, first_diff
from lf_quintic import h_definition, ginf


def method_a_data(G, LF, alpha0, M0, n, **ctrl):
    ph3, logA, C = lmi_scaled(G, LF, alpha0, M0, n, **ctrl)
    cp = wick_1(ph3, logA, C, n, TS([]), TS.const(1))
    poles = [(m, k) for m in range(n + 1) for k in range(m) if cp[m].c[k] != 0]
    return poles, [cp[m].c[m] for m in range(n + 1)]


def bernoulli_amp(alpha0, h, n, bern_at=None):
    """log of a0^(h - 1/2) exp(Sum_(k>=2) (-1)^k B_k(h) (z/a0)^(k-1)/(k(k-1))) in tau, a0 = alpha0 (1 + tau)."""
    hb = h if bern_at is None else bern_at
    out = _sing_amp(alpha0, 0, Fr(1), n, stirling=False, poch=False, expo=h - Fr(1, 2))
    for k in range(2, n + 2):
        p = k - 1
        c = Fr((-1) ** k) * bern_poly(k, hb) / (k * (k - 1))
        for d in range(0, 2 * n - 2 * p + 1):
            out[(d, p)] = out.get((d, p), 0) + c * Fr(alpha0) ** (-p) * gbinom(-p, d)
    return out


def method_a_h(G, LF, alpha0, h, n, bern_at=None):
    ph3, logA, C = lmi_scaled(G, LF, alpha0, 0, n, stirling=False, poch=False, expo=Fr(-1, 2))
    base = _sing_amp(alpha0, 0, Fr(1), n, stirling=False, poch=False, expo=Fr(-1, 2))
    for key, v in base.items():                       # remove the default singular part, insert the h-form
        logA[key] = logA[key] - TS.const(v)
    for key, v in bernoulli_amp(alpha0, h, n, bern_at).items():
        logA[key] = logA.get(key, TS([])) + TS.const(v)
    cp = wick_1(ph3, logA, C, n, TS([]), TS.const(1))
    poles = [(m, k) for m in range(n + 1) for k in range(m) if cp[m].c[k] != 0]
    return poles, [cp[m].c[m] for m in range(n + 1)]


def quintic_data(which, n):
    phs, las, phm, lam, H = (q_chart1(Fr(1), n) if which == 1 else q_chart0(Fr(0), n))
    G = dict(phs)
    G[2] = H[0] / 2
    return G, dict(las)


def main():
    N = int(sys.argv[1])
    log = Log("lf_lmi2_%d" % N)
    # M1
    G, LF = quintic_data(1, N)
    poles, val = method_a_data(G, LF, -1, 0, N)
    log.check("M1 quintic at s* = 1: c_n(1 + eps) has no pole, n <= %d" % N, not poles)
    log.check("M1 quintic at s* = 1: c(z; 1) = h_(w_0) (definition), n <= %d" % N, val == h_definition(N)[:N + 1])
    px, _ = method_a_data(G, LF, -1, 0, N, b1_shift=Fr(1, 7))
    log.control("M1 singular b_1 + 1/7: regular", not px, "first pole (n, eps^k) = %s" % (px[0],) if px else "")
    # M2
    G0, LF0 = quintic_data(0, N)
    poles, val = method_a_data(G0, LF0, 5, 1, N)
    log.check("M2 quintic at s* = 0: c_n(eps) has no pole, n <= %d" % N, not poles)
    log.check("M2 quintic at s* = 0: c(z; 0) = Ginf, n <= %d" % N, val == ginf(N))
    p0, v0 = method_a_data(G0, LF0, 5, 0, N)
    log.control("M2 with M0 = 0 (Gamma(5 s) for Gamma(1 + 5 s)): value Ginf", v0 == ginf(N),
                "regular: %s; first diff z^%s" % (not p0, first_diff(v0, ginf(N))))
    # M3
    for r, alpha0, M0 in [(Fr(1, 2), -1, 0), (Fr(-2, 3), 2, 1), (Fr(3), Fr(-1, 3), -1), (Fr(1, 5), Fr(5, 2), 2)]:
        G3 = {k: r ** k / factorial(k) for k in range(2, N + 3)}
        LF3 = {}
        for d in range(1, N + 2):
            LF3[(d, 0)] = Fr((-1) ** (d + 1), d * 3 ** d)                 # log(1 + x/3)
            for m in (1, 2):
                LF3[(d - 1, m)] = Fr(m, 2) * Fr(-1, 2) ** (d - 1)          # (m/2) z^m/(1 + x/2)
        poles, val = method_a_data(G3, LF3, alpha0, M0, N)
        pred = lmi_predict(G3, LF3, alpha0, M0, N)
        log.check("M3 entire g (r = %s), alpha0 = %s, M0 = %s: no pole and value = LMI formula, n <= %d"
                  % (r, alpha0, M0, N), not poles and val == pred)
    # M4
    ok = True
    for M in range(-2, 4):
        bpa = bernoulli_amp(Fr(-3, 2), Fr(M), N)
        pa = _sing_amp(Fr(-3, 2), M, Fr(1), N)
        keys = set(bpa) | set(pa)
        ok &= all(bpa.get(k, 0) == pa.get(k, 0) for k in keys)
    log.check("M4 integer h = M: Bernoulli-polynomial amplitude = a0^(M - 1/2) Poch_M exp(Stirling), M = -2..3", ok)
    G4 = {2: Fr(3, 7), 3: Fr(-1, 2), 4: Fr(2, 9)}
    LF4 = {(1, 0): Fr(1, 3), (2, 0): Fr(-2, 5), (0, 1): Fr(1, 4), (1, 1): Fr(3, 2), (0, 2): Fr(-1, 6)}
    for h, alpha0 in [(Fr(1, 5), 5), (Fr(1, 3), -1), (Fr(-1, 2), Fr(2, 3)), (Fr(7, 4), Fr(-5, 2))]:
        poles, val = method_a_h(G4, LF4, alpha0, h, N)
        pred = lmi_predict(G4, LF4, alpha0, h, N)
        log.check("M4 h = %s, alpha0 = %s: no pole and value = formula with delta = -h z/alpha0, n <= %d (beyond LMI)"
                  % (h, alpha0, N), not poles and val == pred)
    px, vx = method_a_h(G4, LF4, 5, Fr(1, 5), N, bern_at=Fr(0))
    log.control("M4 B_k(0) for B_k(1/5): regular", not px, "first pole (n, eps^k) = %s" % (px[0],) if px else "")
    log.close()


if __name__ == "__main__":
    main()
