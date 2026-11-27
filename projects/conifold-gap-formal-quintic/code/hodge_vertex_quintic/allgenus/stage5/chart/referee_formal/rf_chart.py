# -*- coding: utf-8 -*-
"""rf_chart.py -- referee: the chart integrand Omega_ch (CHART_PROOF s.3.4) and the 7-variable Omega_7 (s.3.1) against c(z; L).

Taylor data built here from the definitions (no closed forms for the critical point are assumed: linear terms are computed
and must vanish; K = q y^5/prod x_j is computed from q = (1 - L^-5)/3125). Covariance C = -H^-1 by Gauss-Jordan, then the
referee's closed form C = diag(C0) + c u u^T (c = L(L^5-1)/5; chart u = (-5,1,1,1,1), C0 = (5L, zeta^j - L); Omega_7
u = (-1,-5,1,1,1,1,1), C0 = (0, 5L, zeta^j - L)) is asserted entrywise and used by rf_wick.make_expect1.
Usage: py -I -B rf_chart.py NQZ NQZ7 NP NP7      (log logs/rf_chart_<args>.log)
  NQZ: chart order, exact Q(zeta) (L = 1/2 and L = 1); NQZ7: Omega_7 order exact (L = 1/2);
  NP: chart order mod two new primes (L = -3, 5/3, 1); NP7: Omega_7 order mod the primes (L = -3).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rf_core import QZ, FP, prime_1mod10, h_series
from rf_wick import wickN, make_expect1, make_expect_do, mat_inv, hessian, uni, tmul, padd
from rf_sigma import c_series


def logser(coef, a, D, one):
    """[x^k] of coef * log(a + x) - (a + x) for k = 0..D (constant term dropped)"""
    s = [one - one] * (D + 1)
    for k in range(1, D + 1): s[k] = coef * Fr((-1) ** (k + 1), k) / _p(a, k, one)
    s[1] = s[1] - 1
    return s


def _p(x, k, one):
    r = one
    for _ in range(k): r = r * x
    return r


def geom(nv, i, a, D, one, e=-1):
    """(1 + x_i/a)^e as a dict polynomial to degree D"""
    return uni(nv, i, [one * Fr(comb(-e + k - 1, k) * (-1) ** k if e < 0 else comb(e, k)) / _p(a, k, one) for k in range(D + 1)])


def chart_data(L, N, of, zeta, variant=None, qsign=1):
    one = of(1); Lf = of(L); D = 2 * N + 2; nv = 5
    q = (1 - 1 / _p(Lf, 5, one)) * Fr(qsign, 3125); a = [zeta(j) - Lf for j in range(5)]
    Phi = {}
    padd(Phi, uni(nv, 0, logser(of(5), Lf * 5, D, one)))
    for j in range(1, 5): padd(Phi, uni(nv, j, logser(zeta(j) - 1, a[j], D, one)))
    K = q * _p(Lf * 5, 5, one) / (a[1] * a[2] * a[3] * a[4])
    cp = geom(nv, 0, Lf * 5, D, one, 5)
    for j in range(1, 5): cp = tmul(cp, geom(nv, j, a[j], D, one), D)
    cp[(0,) * nv] = cp[(0,) * nv] - 1; padd(Phi, cp, K)
    amp = {(0,) * nv: one}
    for j in range(1, 5): amp = tmul(amp, geom(nv, j, a[j], 2 * N, one), 2 * N)
    if variant == 'dyy': amp = tmul(amp, geom(nv, 0, Lf * 5, 2 * N, one), 2 * N)
    if variant == 'nojac': amp = geom(nv, 0, Lf * 5, 2 * N, one, -5)
    return Phi, [amp] + [{} for _ in range(N)], K


def omega7_data(L, N, of, zeta, variant=None):
    one = of(1); Lf = of(L); D = 2 * N + 2; nv = 7
    a = [zeta(j) - Lf for j in range(5)]; Phi = {}
    padd(Phi, uni(nv, 1, logser(Lf * 5, Lf * 5, D, one)))                # 5 sigma log y - y, sigma = L part
    s1 = logser(of(5), Lf * 5, D, one); s1[1] = s1[1] + 1                   # 5 s log y, s = sigma - L
    padd(Phi, {tuple([1, k] + [0] * 5): c for k, c in enumerate(s1) if k >= 1 and k + 1 <= D})
    for j in range(5):
        padd(Phi, uni(nv, 2 + j, logser(a[j], a[j], D, one)))              # (zeta^j - L) log x_j - x_j
        sj = logser(-one, a[j], D, one); sj[1] = sj[1] + 1                  # -s log x_j
        padd(Phi, {tuple([1] + [0] * (1 + j) + [k] + [0] * (4 - j)): c for k, c in enumerate(sj) if k >= 1 and k + 1 <= D})
    amp = {(0,) * nv: one}
    for j in range(5): amp = tmul(amp, geom(nv, 2 + j, a[j], 2 * N, one), 2 * N)
    if variant == 'dyy': amp = tmul(amp, geom(nv, 1, Lf * 5, 2 * N, one), 2 * N)
    return Phi, [amp] + [{} for _ in range(N)]


def cov_closed(L, of, zeta, which):
    one = of(1); Lf = of(L); a = [zeta(j) - Lf for j in range(5)]; c = Lf * (_p(Lf, 5, one) - 1) * Fr(1, 5)
    if which == 'chart': return [Lf * 5] + a[1:], c, [of(-5)] + [one] * 4
    return [one - one, Lf * 5] + a, c, [of(-1), of(-5)] + [one] * 5


def sp_of(Phi, Amp, N, nv, L, of, zeta, which, use_do=False):
    one = of(1); H = hessian(Phi, nv, one); C = [[-x for x in r] for r in mat_inv(H, one)]
    C0, c, u = cov_closed(L, of, zeta, which)
    assert all(C[i][j] == (C0[i] if i == j else 0) + c * u[i] * u[j] for i in range(nv) for j in range(nv)), 'covariance'
    ex = make_expect_do(C, one) if use_do else make_expect1(C0, c, u, one)
    Nn = wickN(Phi, Amp, N, nv, ex, one); return [x / Nn[0] for x in Nn]


def c_field(L, N, of, zeta, variant=None):
    one = of(1); Lf = of(L); ys = [1 / (zeta(j) - Lf) for j in range(5)]; pw = [one] * 5; P = [of(5)]
    for k in range(1, 2 * N + 5):
        pw = [x * y for x, y in zip(pw, ys)]; P.append(pw[0] + pw[1] + pw[2] + pw[3] + pw[4])
    return c_series(Lf, P, None, N, one, variant)


if __name__ == '__main__':
    A = [int(x) for x in sys.argv[1:5]] if len(sys.argv) > 4 else [2, 1, 3, 2]
    NQZ, NQZ7, NP, NP7 = A
    LOGREL = 'logs/rf_chart_%d_%d_%d_%d.log' % tuple(A); LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
    FAILS = []; T0 = time.time()

    def out(s): print(s); LOG.write(s + '\n'); LOG.flush()

    def check(ok, msg, ctrl=False):
        good = (not ok) if ctrl else ok
        out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
        if not good: FAILS.append(msg)

    out('rf_chart: exact Q(zeta) chart z^%d, Omega_7 z^%d; mod p chart z^%d, Omega_7 z^%d; log %s' % (NQZ, NQZ7, NP, NP7, LOGREL))
    of, zt = QZ.of, QZ.zeta
    Phi, Amp, K = chart_data(Fr(1, 2), 2, of, zt)
    n1 = sp_of(Phi, Amp, 2, 5, Fr(1, 2), of, zt, 'chart'); n2 = sp_of(Phi, Amp, 2, 5, Fr(1, 2), of, zt, 'chart', True)
    check(n1 == n2, 'expectation algorithms agree (closed-form rank-one vs differential operator), chart at L = 1/2 through z^2')
    for L in [Fr(1, 2), Fr(1)]:
        t1 = time.time(); Phi, Amp, K = chart_data(L, NQZ, of, zt)
        check(K == of(L - 1), 'L = %s: K = q y^5/prod_(j>=1) x_j = L - 1 at p_ch(L) (computed from q)' % L)
        ch = sp_of(Phi, Amp, NQZ, 5, L, of, zt, 'chart')
        ref = [of(x) for x in h_series(NQZ)] if L == 1 else c_field(L, NQZ, of, zt)
        check(ch == ref, '[Q(zeta)] L = %s: chart SP = %s exactly through z^%d   %.1f s'
              % (L, 'h_(w_0)' if L == 1 else 'c(z; L)', NQZ, time.time() - t1))
    t1 = time.time(); Phi, Amp = omega7_data(Fr(1, 2), NQZ7, of, zt)
    o7 = sp_of(Phi, Amp, NQZ7, 7, Fr(1, 2), of, zt, 'omega7')
    check(o7 == c_field(Fr(1, 2), NQZ7, of, zt), '[Q(zeta)] L = 1/2: 7-variable SP of Omega_7 = c(z; L) exactly through z^%d   %.1f s'
          % (NQZ7, time.time() - t1))
    for p in [prime_1mod10(7 * 10 ** 18), prime_1mod10(29 * 10 ** 23)]:
        FP.setup(p); of, zt = FP.of, FP.zeta; t1 = time.time()
        for L in [Fr(-3), Fr(5, 3), Fr(1)]:
            Phi, Amp, K = chart_data(L, NP, of, zt); ch = sp_of(Phi, Amp, NP, 5, L, of, zt, 'chart')
            ref = [of(x) for x in h_series(NP)] if L == 1 else c_field(L, NP, of, zt)
            check(ch == ref, '[F_p, p = %d] L = %s: chart SP = %s through z^%d' % (p, L, 'h' if L == 1 else 'c(z; L)', NP))
        Phi, Amp = omega7_data(Fr(-3), NP7, of, zt)
        check(sp_of(Phi, Amp, NP7, 7, Fr(-3), of, zt, 'omega7') == c_field(Fr(-3), NP7, of, zt),
              '[F_p, p = %d] L = -3: 7-variable SP of Omega_7 = c(z; L) through z^%d   (prime block %.0f s)' % (p, NP7, time.time() - t1))
    # controls (mod the last prime)
    NCt = min(NP, 3); L = Fr(-3); ref = c_field(L, NCt, of, zt)
    for var, txt in [('nojac', 'chart without the Jacobian X_0 of (L2) (amplitude y^-5)'), ('dyy', 'chart with measure dy/y')]:
        Phi, Amp, K = chart_data(L, NCt, of, zt, var); ch = sp_of(Phi, Amp, NCt, 5, L, of, zt, 'chart')
        bad = [n for n in range(NCt + 1) if ch[n] != ref[n]]
        check(not bad, '%s differs from c(z; L), first at z^%s' % (txt, bad[0] if bad else '-'), ctrl=True)
    try:
        Phi, Amp, K = chart_data(L, 1, of, zt, None, -1); sp_of(Phi, Amp, 1, 5, L, of, zt, 'chart'); moved = False
    except AssertionError:
        moved = True
    check(not moved, 'q -> -q in the chart phase: p_ch(L) is no longer critical (linear terms non-zero)', ctrl=True)
    N7c = 2; Phi, Amp = omega7_data(L, N7c, of, zt, 'dyy'); o7 = sp_of(Phi, Amp, N7c, 7, L, of, zt, 'omega7')
    check(o7 == c_field(L, N7c, of, zt), 'Omega_7 with dy/y equals c(z; L) through z^2 (it first differs at z^2)', ctrl=True)
    check(o7 == c_field(L, N7c, of, zt, 'g5s'), 'Omega_7 with dy/y equals the Gamma(5s) sigma-form through z^2 (consistent pair)')
    out('RF_CHART %s: %s   (%.0f s)' % (' '.join(map(str, A)), 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
    LOG.close()
