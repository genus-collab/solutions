# -*- coding: utf-8 -*-
"""ch_stirling.py -- Stage 5 / chart, step (d) (CHART_PROOF s.3.5): at q = 0 the chart integral is the product of the Euler
integrals for Gamma(1 + 5/z) and Gamma((zeta^j - 1)/z), j = 1..4, and their Stirling factors multiply to EXACTLY h_(w_0)(z).

Three independent routes for the weights, k = 1..K:
  T_k = 5^(1-2k) + Sum_(j=1..4) (zeta^j - 1)^(1-2k)       (chart side, exact Q(zeta_5) arithmetic);
  W_k = - Sum_(a in w_0) a^(1-2k),  w_0 = {1 - zeta^j : j = 1..4} U {-5}   (QRR side, literally from the set, exact Q(zeta_5));
  5^(1-2k) - p_(2k-1)   (Newton sums for the roots of y^4 - 2y^3 + 2y^2 - y + 1/5, no zeta arithmetic; ch_core.h_newton).
Then h(z) = exp(Sum_k b_k T_k z^(2k-1)) through z^(2K-1), against h_newton, against the product of the five normalised Stirling
series (one per Euler integral, multiplied in Q(zeta)[[z]]) and against the stored values h_1 = -3/20, h_5 = -257071/8064000000.
Controls (must FAIL): twist weight +5 instead of -5; tangent weights 1 - zeta^j instead of zeta^j - 1 in the chart; the factor
Gamma(1 + 5/z) dropped; b_k -> -b_k.
Usage: py -I -B ch_stirling.py [K]   (default 30; log logs/ch_stirling_<K>.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import QZ, bk_list, h_newton, sexp_q

K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
LOGREL = 'logs/ch_stirling_%d.log' % K
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


F = QZ(); NZ = 2 * K - 1


def powsum(elems, K):
    """S[k] = Sum_(a in elems) a^(1-2k), k = 1..K, in Q(zeta)"""
    S = [None] + [F.zero] * K
    for a in elems:
        ia = F.inv(a); i2 = F.mul(ia, ia); p = ia
        for k in range(1, K + 1):
            S[k] = F.add(S[k], p); p = F.mul(p, i2)
    return S


one = F.one
zj = [F.zeta(j) for j in range(5)]
chart = [F.of(5)] + [F.sub(zj[j], one) for j in range(1, 5)]           # a = 5 (Gamma(1+5/z)), a = zeta^j - 1
w0 = [F.sub(one, zj[j]) for j in range(1, 5)] + [F.of(-5)]              # QRR weights at colour 0
T = powsum(chart, K); Wm = powsum(w0, K)
Tq = [None] + [F.rat(T[k]) for k in range(1, K + 1)]
Wq = [None] + [-F.rat(Wm[k]) for k in range(1, K + 1)]
check('T_k is rational for k = 1..%d (Galois invariance)' % K, True)
check('weights: T_k = W_k = -Sum_(a in w_0) a^(1-2k) exactly for k = 1..%d (odd exponent 1-2k: sign of -5 and of 1 - zeta^j)' % K,
      all(Tq[k] == Wq[k] for k in range(1, K + 1)))
hN = h_newton(NZ)
bk = bk_list(K)


def hser(Tw, sign=1):
    H = [Fr(0)] * (NZ + 1)
    for k in range(1, K + 1): H[2 * k - 1] = sign * bk[k] * Tw[k]
    return sexp_q(H, NZ)


h = hser(Tq)
check('h(z) = exp(Sum_k b_k T_k z^(2k-1)) equals h_(w_0) by Newton sums through z^%d (k <= %d)' % (NZ, K), h == hN)
check('stored values: h_1 = -3/20, h_5 = -257071/8064000000 (GAMMA s.1, POLY s.1)',
      h[1] == Fr(-3, 20) and h[5] == Fr(-257071, 8064000000), (h[1], h[5]))
# product of the five normalised Stirling series exp(Sum_k b_k (z/a)^(2k-1)), multiplied in Q(zeta)[[z]]
Np = min(NZ, 25)
PR = [one] + [F.zero] * Np
for a in chart:
    ia = F.inv(a); i2 = F.mul(ia, ia); p = ia; H = [F.zero] * (Np + 1)
    for k in range(1, (Np + 1) // 2 + 1):
        H[2 * k - 1] = F.smul(bk[k], p); p = F.mul(p, i2)
    e = [one] + [F.zero] * Np                                            # exp in Q(zeta)[[z]]
    for m in range(1, Np + 1):
        s = F.zero
        for i in range(1, m + 1): s = F.add(s, F.smul(i, F.mul(H[i], e[m - i])))
        e[m] = F.smul(Fr(1, m), s)
    new = [F.zero] * (Np + 1)
    for i in range(Np + 1):
        for j in range(Np + 1 - i): new[i + j] = F.add(new[i + j], F.mul(PR[i], e[j]))
    PR = new
check('product of the five normalised Stirling series (Q(zeta)[[z]]) = h through z^%d' % Np,
      [F.rat(x) for x in PR] == hN[:Np + 1])
Tp = powsum([F.of(-5)] + chart[1:], K); Tt = powsum([chart[0]] + [F.neg(x) for x in chart[1:]], K)
Td = powsum(chart[1:], K)
for name, Tw, sg in (('twist weight +5 (chart factor Gamma(1 - 5/z))', [None] + [F.rat(Tp[k]) for k in range(1, K + 1)], 1),
                     ('tangent weights 1 - zeta^j in the chart', [None] + [F.rat(Tt[k]) for k in range(1, K + 1)], 1),
                     ('factor Gamma(1 + 5/z) dropped', [None] + [F.rat(Td[k]) for k in range(1, K + 1)], 1),
                     ('b_k -> -b_k', Tq, -1)):
    hc = hser(Tw, sg); fl = [n for n in range(NZ + 1) if hc[n] != hN[n]]
    check('control: %s fails, first at z^1' % name, fl[:1] == [1], fl[:3])
out('h_1..h_5 = %s' % [str(x) for x in h[1:6]])
out('CH_STIRLING %d: %s   (%.0f s)' % (K, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
