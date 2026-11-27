# -*- coding: utf-8 -*-
"""rf_hess.py -- referee check of CHART_PROOF s.3.2 (critical point, every Hessian and determinant, both orders), exact in Q(zeta).
Hessians are read off the referee's own Taylor data (rf_chart.omega7_data / chart_data, built from the definitions), then
compared with the closed forms of the note. Polynomial identities in L are decided by evaluation at more points than the degree.
Usage: py -I -B rf_hess.py      (log logs/rf_hess.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from rf_core import QZ
from rf_wick import hessian, det, mat_inv
from rf_chart import omega7_data, chart_data, _p

LOGREL = 'logs/rf_hess.log'; LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s): print(s); LOG.write(s + '\n'); LOG.flush()


def check(ok, msg, ctrl=False):
    good = (not ok) if ctrl else ok
    out(('PASS ' if good else 'FAIL ') + ('control: ' if ctrl else '') + msg)
    if not good: FAILS.append(msg)


of, zt = QZ.of, QZ.zeta; one = of(1); zero = of(0)
Ls = [Fr(2), Fr(-2, 3), Fr(1, 2), Fr(5, 3), Fr(-3), Fr(7, 5), Fr(11, 10), Fr(9, 10)]
out('rf_hess: exact over Q(zeta_5) at L = %s; log %s' % ([str(x) for x in Ls], LOGREL))
ok = {k: True for k in ['crit', 'H7', 'det7', 'a', 'b', 'bschur', 'dch', 'log', 'ctlq', 'ctdet']}
for L in Ls:
    Lf = of(L); a = [zt(j) - Lf for j in range(5)]; q = (1 - 1 / _p(Lf, 5, one)) * Fr(1, 3125)
    P = a[0] * a[1] * a[2] * a[3] * a[4]; Lam1 = a[1] * a[2] * a[3] * a[4]; uL5 = 1 - _p(Lf, 5, one)
    ok['crit'] &= (-q * _p(Lf * 5, 5, one) == P) and P == uL5
    ok['ctlq'] &= (q * _p(Lf * 5, 5, one) == P)                       # control: +q
    Phi, _ = omega7_data(L, 0, of, zt); H = hessian(Phi, 7, one)
    exp = [[zero] * 7 for _ in range(7)]
    exp[0][1] = exp[1][0] = 1 / Lf; exp[1][1] = -1 / (Lf * 5)
    for j in range(5): exp[0][2 + j] = exp[2 + j][0] = -1 / a[j]; exp[2 + j][2 + j] = -1 / a[j]
    ok['H7'] &= H == exp
    d7 = det(H, one); ok['det7'] &= d7 == 1 / (Lf * Lf * uL5 * uL5); ok['ctdet'] &= d7 == 1 / (Lf * Lf * uL5)
    # order (a): v = (y, x_0..x_4), w = sigma
    Hvv = [r[1:] for r in H[1:]]; dvv = det(Hvv, one)
    diag = all(Hvv[i][j] == 0 for i in range(6) for j in range(6) if i != j)
    iv = mat_inv(Hvv, one); hs = H[0][1:]
    S = H[0][0] - sum((hs[i] * iv[i][j] * hs[j] for i in range(6) for j in range(6)), zero)
    ok['a'] &= diag and dvv == 1 / (Lf * 5 * uL5) and S == 5 / (Lf * uL5) and dvv * S == d7
    # order (b): v = (sigma, x_0) = indices (0, 2), w = (y, x_1..x_4) = (1, 3, 4, 5, 6)
    vi, wi = [0, 2], [1, 3, 4, 5, 6]
    Hb = [[H[i][j] for j in vi] for i in vi]; db = det(Hb, one); ib = mat_inv(Hb, one)
    Sb = [[H[i][j] - sum((H[i][vi[k]] * ib[k][l] * H[vi[l]][j] for k in range(2) for l in range(2)), zero) for j in wi] for i in wi]
    ok['b'] &= db == -1 / ((1 - Lf) * (1 - Lf)) and Hb[0][0] == 0
    Pc, _, K = chart_data(L, 0, of, zt); Hc = hessian(Pc, 5, one)
    ok['bschur'] &= Sb == Hc
    dc = det(Hc, one); ok['dch'] &= dc == -1 / (Lf * Lf * Lam1 * Lam1) and db * dc == d7
    Dg = [Lf * 5] + a[1:]; v = [of(5)] + [of(-1)] * 4; g = [-(Lf * 5)] + [-x for x in a[1:]]
    ok['log'] &= all(Dg[i] * Hc[i][j] * Dg[j] == (g[i] if i == j else 0) + (Lf - 1) * v[i] * v[j] for i in range(5) for j in range(5))
check(ok['crit'], 'sigma-equation exponentiated: -q (5L)^5 = prod_j (zeta^j - L) = 1 - L^5 with q = (1 - L^-5)/5^5 (u L^5 = 1)')
check(ok['H7'], 'H_7 from the Taylor data = the closed forms of s.3.2 entrywise (H_ss = 0, H_sy = 1/L, H_sx_j = H_x_jx_j = -1/(zeta^j - L), H_yy = -1/(5L))')
check(ok['det7'], 'det H_7 = 1/(L^2 (1 - L^5)^2)')
check(ok['a'], 'order (a): 6x6 block diagonal, det = 1/(5L(1 - L^5)), Schur complement in sigma = 5/(L(1 - L^5)) = S_0\'\'(L), product = det H_7')
check(ok['b'], 'order (b): H_(sigma,x0) = [[0, -1/(1-L)], [-1/(1-L), -1/(1-L)]], det = -1/(1 - L)^2')
check(ok['bschur'], 'order (b): the 5x5 Schur complement EQUALS the chart Hessian (from the chart Taylor data) entrywise')
check(ok['dch'], 'det H_ch = -1/(L^2 Lam1(L)^2), Lam1 = 1+L+..+L^4, and det H_(sigma,x0) det H_ch = det H_7')
check(ok['log'], 'log coordinates: D H_ch D = diag(-5L, -(zeta^j - L)) + (L - 1) v v^T, D = diag(5L, zeta^j - L), v = (5,-1,-1,-1,-1)')
# polynomial identities by evaluation (degree bounds: scaled H_7 has one row of degree <= 5, six of degree <= 1: det degree <= 11;
# G(L) = diag(-5L, -(zeta^j - L)) + (L-1) v v^T has entries of degree <= 1: det degree <= 5)
okp = True; okg = True
for Li in range(2, 42):
    Lf = of(Li); a = [zt(j) - Lf for j in range(5)]; uL5 = 1 - _p(Lf, 5, one)
    Phi, _ = omega7_data(Fr(Li), 0, of, zt); H = hessian(Phi, 7, one)
    sc = [Lf * uL5, Lf * 5] + a; M = [[H[i][j] * sc[i] for j in range(7)] for i in range(7)]
    okp &= det(M, one) == of(5)
    v = [of(5)] + [of(-1)] * 4; g = [-(Lf * 5)] + [-x for x in a[1:]]
    G = [[(g[i] if i == j else zero) + (Lf - 1) * v[i] * v[j] for j in range(5)] for i in range(5)]
    okg &= det(G, one) == of(-25)
check(okp, 'polynomial identity: det of H_7 with rows scaled by L(1-L^5), 5L, zeta^j - L is the constant 5 (40 points > degree 11)')
check(okg, 'polynomial identity: det(diag(-5L, -(zeta^j - L)) + (L-1) v v^T) = -25 identically (40 points > degree 5)')
Pc, _, K = chart_data(Fr(1), 0, of, zt); Hc = hessian(Pc, 5, one)
expc = [[zero] * 5 for _ in range(5)]; expc[0][0] = of(Fr(-1, 5))
for j in range(1, 5): expc[j][j] = -1 / (zt(j) - 1)
check(K == 0 and Hc == expc and det(Hc, one) == of(Fr(-1, 25)), 'q = 0 (L = 1): K = 0, chart Hessian = diag(-1/5, -1/(zeta^j - 1)), det = -1/25')
Pm, _, _ = chart_data(Fr(2), 0, of, zt, None, -1)
lin = [c for e, c in Pm.items() if sum(e) == 1 and not (c == 0)]
check(not lin, 'q -> -q in the chart phase leaves p_ch(L) critical (L = 2)', ctrl=True)
check(ok['ctlq'], 'sigma-equation with +q in place of -q', ctrl=True)
check(ok['ctdet'], 'wrong closed form det H_7 = 1/(L^2 (1 - L^5))', ctrl=True)
out('RF_HESS: %s   (%.0f s)' % ('ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
LOG.close()
