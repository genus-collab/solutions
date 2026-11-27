# -*- coding: utf-8 -*-
"""ch_hess.py -- Stage 5 / chart, critical points and Hessians of CHART_PROOF s.3.2, exact.

Part 1 (exact Q(zeta_5), at rational L != 0, 1): the Taylor expansions of the 7-variable phase at p(L) = (L, 5L, zeta^j - L) and of
the chart phase at p_ch(L) = (5L, zeta^j - L) are built by series arithmetic (ch_wick); checked: the exponentiated sigma-equation
-q (5L)^5 = prod_j (zeta^j - L) with q = (1 - L^-5)/5^5, vanishing gradients, K = q y^5/prod x_j = L - 1 at p_ch, the hand formulas
for all Hessian entries, det H_7 directly, by order (a) (6x6 block x Schur complement in sigma) and by order (b) ((sigma, x_0) block
x 5x5 Schur complement), Schur complement of order (b) = chart Hessian ENTRYWISE, and the closed forms
   det H_7 = 1/(L^2 (1-L^5)^2),  S_sigma = 5/(L(1-L^5)),  det H_(sigma,x0) = -1/(1-L)^2,  det H_ch = -(1-L)^2/(L^2 (1-L^5)^2).
Part 2 (symbolic in L, Q(zeta)[L], Bareiss): with rows of H_7 scaled by L(1-L^5), 5L, (zeta^j - L) the determinant is the
constant 5; the chart Hessian in log coordinates, diag(-5L, -(zeta^j - L)) + (L-1) v v^T, v = (5,-1,-1,-1,-1), has determinant
the constant -25; prod_j (zeta^j - L) = 1 - L^5 and Sum_j prod_(i!=j) (zeta^i - L) = 5 L^4.
Part 3 (q = 0, L = 1): chart Hessian = diag(-1/5, -1/(zeta^j - 1)), determinant -1/25.
Controls (must FAIL): q -> -q (critical equation); a wrong closed form for det H_7.
Usage: py -I -B ch_hess.py   (log logs/ch_hess.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import QZ, log1p_lin, geom
from ch_wick import guni, gadd, gmul, gconst, gtrunc, hessian, linear_is_zero, matinv, det

LOGREL = 'logs/ch_hess.log'
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()
F = QZ(); one = F.one; Z = [F.zeta(j) for j in range(5)]


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def phase7(L, D=3):
    """Taylor polynomial (degrees 1..D) of Phi at p(L); variables (s, eta_y, xi_0..xi_4); the s-linear term with transcendental
    coefficient lambda + 5 log 5L - Sum log(zeta^j - L) is omitted (it vanishes: critical equation, checked separately)"""
    n = 7; Lq = F.of(L); c = [F.sub(Z[j], Lq) for j in range(5)]
    lg = log1p_lin(F, F.smul(5, Lq), D)
    P = gadd(F, guni(F, n, 1, [F.zero] + [F.smul(5 * L, x) for x in lg[1:]]), guni(F, n, 1, [F.zero, F.neg(one)]))
    P = gadd(F, P, gmul(F, guni(F, n, 0, [F.zero, F.of(5)]), guni(F, n, 1, lg), D))
    for j in range(5):
        lj = log1p_lin(F, c[j], D)
        P = gadd(F, P, guni(F, n, 2 + j, [F.zero, F.neg(one)] + [F.zero] * (D - 1)))
        P = gadd(F, P, guni(F, n, 2 + j, [F.mul(c[j], x) for x in lj]))
        P = gadd(F, P, gmul(F, guni(F, n, 0, [F.zero, F.neg(one)]), guni(F, n, 2 + j, lj), D))
    return P


def phasech(L, D=3, qsign=1):
    """Taylor polynomial (degrees 1..D) of psi = 5 log y - y + Sum_(j>=1) (-x_j + (zeta^j - 1) log x_j) + q y^5/prod x_j at p_ch(L)"""
    n = 5; Lq = F.of(L); c = [F.sub(Z[j], Lq) for j in range(5)]
    q = Fr(qsign) * (1 - Fr(1) / Fr(L) ** 5) / 3125
    K = F.mul(F.of(q * (5 * Fr(L)) ** 5), F.inv(F.mul(F.mul(c[1], c[2]), F.mul(c[3], c[4]))))
    P = guni(F, n, 0, [F.zero] + [F.smul(5, x) for x in log1p_lin(F, F.smul(5, Lq), D)[1:]])      # 5 log(1 + eta/5L)
    P = gadd(F, P, guni(F, n, 0, [F.zero, F.neg(one)]))                                         # - eta
    prod = guni(F, n, 0, geom(F, F.smul(5, Lq), D, 5))                                          # (1 + eta/5L)^5
    for j in range(1, 5):
        P = gadd(F, P, guni(F, n, j, [F.zero, F.neg(one)] + [F.zero] * (D - 1)))
        P = gadd(F, P, guni(F, n, j, [F.mul(F.sub(Z[j], one), x) for x in log1p_lin(F, c[j], D)]))
        prod = gmul(F, prod, guni(F, n, j, geom(F, c[j], D, -1)), D)
    prod = {d: m for d, m in prod.items() if d >= 1}
    P = gadd(F, P, prod, K)
    return P, K


def schur(Hm, I, J):
    """Schur complement of the block I in Hm, on the indices J"""
    A = [[Hm[i][k] for k in I] for i in I]; Ai = matinv(F, A)
    S = []
    for a in J:
        row = []
        for b in J:
            v = Hm[a][b]
            for x in range(len(I)):
                for y in range(len(I)):
                    v = F.sub(v, F.mul(F.mul(Hm[a][I[x]], Ai[x][y]), Hm[I[y]][b]))
            row.append(v)
        S.append(row)
    return S


def q_of(L): return F.of((1 - Fr(1) / Fr(L) ** 5) / 3125)


out('ch_hess: exact over Q(zeta_5); log %s' % LOGREL)
LS = [Fr(2), Fr(3, 2), Fr(-2, 3), Fr(7, 5), Fr(1, 3), Fr(-3), Fr(11, 10), Fr(9, 10)]
ok = {k: True for k in ('crit', 'grad', 'entries', 'det7', 'orda', 'ordb', 'schur=chart', 'closed', 'K', 'log25')}
for L in LS:
    Lq = F.of(L); c = [F.sub(Z[j], Lq) for j in range(5)]; prodc = F.mul(F.mul(F.mul(c[0], c[1]), F.mul(c[2], c[3])), c[4])
    ok['crit'] &= F.eq(F.mul(F.neg(q_of(L)), F.of((5 * L) ** 5)), prodc)
    P7 = phase7(L); ok['grad'] &= linear_is_zero(F, P7); H = hessian(F, P7, 7)
    hand = [[F.zero] * 7 for _ in range(7)]
    hand[0][1] = hand[1][0] = F.of(1 / L); hand[1][1] = F.of(-1 / (5 * L))
    for j in range(5):
        hand[0][2 + j] = hand[2 + j][0] = F.neg(F.inv(c[j])); hand[2 + j][2 + j] = F.neg(F.inv(c[j]))
    ok['entries'] &= all(F.eq(H[i][k], hand[i][k]) for i in range(7) for k in range(7))
    d7 = det(F, H); e7 = F.of(1 / (L ** 2 * (1 - L ** 5) ** 2)); ok['det7'] &= F.eq(d7, e7)
    W = list(range(1, 7)); Ss = schur(H, W, [0])[0][0]; dW = det(F, [[H[i][k] for k in W] for i in W])
    ok['orda'] &= F.eq(F.mul(dW, Ss), d7) and F.eq(Ss, F.of(5 / (L * (1 - L ** 5))))
    V = [0, 2]; Wb = [1, 3, 4, 5, 6]; dV = det(F, [[H[i][k] for k in V] for i in V]); Sb = schur(H, V, Wb)
    ok['ordb'] &= F.eq(F.mul(dV, det(F, Sb)), d7) and F.eq(dV, F.of(-1 / (1 - L) ** 2))
    Pc, K = phasech(L); ok['K'] &= F.eq(K, F.of(L - 1)); ok['grad'] &= linear_is_zero(F, Pc); Hc = hessian(F, Pc, 5)
    ok['schur=chart'] &= all(F.eq(Hc[i][k], Sb[i][k]) for i in range(5) for k in range(5))
    ok['closed'] &= F.eq(det(F, Hc), F.of(-(1 - L) ** 2 / (L ** 2 * (1 - L ** 5) ** 2)))
    Dg = [F.smul(5, Lq)] + c[1:]
    ok['log25'] &= F.eq(det(F, [[F.mul(F.mul(Dg[i], Dg[k]), Hc[i][k]) for k in range(5)] for i in range(5)]), F.of(-25))
LSs = ', '.join(str(x) for x in LS)
check('exponentiated sigma-equation -q (5L)^5 = prod_j (zeta^j - L), q = (1 - L^-5)/3125, at L = %s' % LSs, ok['crit'])
check('gradients of Phi at p(L) and of psi at p_ch(L) vanish (Taylor code)', ok['grad'])
check('K = q y^5/prod_(j>=1) x_j = L - 1 at p_ch(L)', ok['K'])
check('Hessian H_7 from the Taylor code equals the hand formulas (s.3.2) entrywise', ok['entries'])
check('det H_7 (direct elimination) = 1/(L^2 (1-L^5)^2)', ok['det7'])
check('order (a): det(6x6 block, diagonal) x Schur complement in sigma = det H_7, Schur = S_0\'\'(L) = 5/(L(1-L^5))', ok['orda'])
check('order (b): det H_(sigma,x0) = -1/(1-L)^2 and det H_(sigma,x0) x det(5x5 Schur complement) = det H_7', ok['ordb'])
check('order (b): the 5x5 Schur complement EQUALS the chart Hessian entrywise', ok['schur=chart'])
check('det H_ch = -(1-L)^2/(L^2 (1-L^5)^2); in log coordinates det = -25', ok['closed'] and ok['log25'])
from ch_hess_sym import symbolic_part, controls_part
symbolic_part(F, check, out)
controls_part(F, check, out, phasech, q_of)
out('CH_HESS: %s   (%.0f s)' % ('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
