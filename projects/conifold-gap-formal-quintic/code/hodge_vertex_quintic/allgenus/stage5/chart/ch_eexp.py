# -*- coding: utf-8 -*-
"""ch_eexp.py -- Stage 5 / chart, step (e) (CHART_PROOF s.3.7): the full-Stirling saddle series, re-expanded at q = 0, is h(z) times
the colour-0 I-function, coefficientwise in q:

    L e^(Dpsi(q)/z) c(z; L(q)) = h(z) Sum_d q^d I_d(z)     in Q((z))[[q]],
    I_d(z) = z^-d prod_(k=1..5d) (5 + k z) / ( d! prod_(i=1..d) Lam1(1 + i z) ),   Lam1(X) = X^4 + X^3 + X^2 + X + 1,
    L(q) = (1 - 3125 q)^(-1/5),  Dpsi(q) = 5 log L - 4 (L - 1) - Sum_(m>=2) ((L-1)^m/m) Sum_(j=1..4) (zeta^j - 1)^(1-m).

Inputs, all independent of each other: c_n(L) from the sigma-form (ch_phases.sigma_data + ch_wick.sp1, exact over Q at rational L),
interpolated as a Laurent polynomial with exponents 4n - 5k (k in a padded range, fit confirmed on EXTRA further points, padding
coefficients must vanish); h by Newton sums; I_d from the product formula (FORMAL_QUINTIC_genus2 s.1 with lambda_i = 1).
Also D_q Dpsi = L - 1 is checked. Controls (must FAIL): the factor L dropped; Dpsi dropped; h with twist weight +5 (h_1 -> -11/60).
Usage: py -I -B ch_eexp.py [DQ] [NZ] [EXTRA]   (default 4 6 3; log logs/ch_eexp_<DQ>_<NZ>.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import QF, QZ, h_newton, sexp_q
from ch_wick import sp1, normalise
from ch_phases import sigma_data

DQ = int(sys.argv[1]) if len(sys.argv) > 1 else 4
NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 6
EXTRA = int(sys.argv[3]) if len(sys.argv) > 3 else 3
LOGREL = 'logs/ch_eexp_%d_%d.log' % (DQ, NZ)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()
Q = QF(); NC = NZ + DQ; PAD = 2


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def solve(M, y):
    n = len(y); A = [r[:] + [y[i]] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [x - f * w for x, w in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]


def smul(a, b, K):
    return [sum(a[i] * b[k - i] for i in range(k + 1) if i < len(a) and k - i < len(b)) for k in range(K + 1)]


# 1. c_n(L) as Laurent polynomials: exponents 4n - 5k, k = -PAD .. floor(4n/5) + PAD
pts = sorted({Fr(a, b) for b in range(1, 6) for a in range(-9, 10)} - {Fr(0), Fr(1), Fr(-1)}, key=lambda x: (abs(x.numerator) + x.denominator, x))
need = (4 * NC) // 5 + 1 + 2 * PAD + EXTRA; pts = pts[:need]
QZf = QZ(); vals = {L: [QZf.rat(x) for x in normalise(QZf, sp1(QZf, *sigma_data(QZf, L, NC), NC))] for L in pts}
coef = {}; fitok = padok = True
for n in range(NC + 1):
    ks = list(range(-PAD, (4 * n) // 5 + PAD + 1)); m = len(ks)
    a = solve([[L ** (4 * n - 5 * k) for k in ks] for L in pts[:m]], [vals[L][n] for L in pts[:m]])
    fitok &= all(sum(x * L ** (4 * n - 5 * k) for x, k in zip(a, ks)) == vals[L][n] for L in pts[m:])
    padok &= all(x == 0 for x, k in zip(a, ks) if k < 0 or k > (4 * n) // 5)
    coef[n] = {k: x for x, k in zip(a, ks) if x}
check('c_n(L), n <= %d: Laurent fit with exponents 4n-5k confirmed on %d extra points each' % (NC, EXTRA), fitok)
check('padding coefficients vanish: c_n in Q[L], exponents 4n - 5k, 0 <= k <= floor(4n/5) (n <= %d)' % NC, padok)
# 2. q-series: L(q), c_n(L(q)), Dpsi(q)
D = DQ + 1
Lq = [Fr(1)]
for k in range(1, D): Lq.append(Lq[-1] * (Fr(-1, 5) - (k - 1)) / k * (-3125))           # (1 - 3125 q)^(-1/5)


def spow_int(a, e, K):                    # a^e for an integer e (a[0] != 0)
    if e >= 0:
        r = [Fr(1)] + [Fr(0)] * K
        for _ in range(e): r = smul(r, a, K)
        return r
    ia = [1 / a[0]] + [Fr(0)] * K
    for m in range(1, K + 1): ia[m] = -sum(a[i] * ia[m - i] for i in range(1, m + 1) if i < len(a)) / a[0]
    return spow_int(ia, -e, K)


cq = {n: [sum(x * spow_int(Lq, 4 * n - 5 * k, DQ)[d] for k, x in coef[n].items()) for d in range(D)] for n in range(NC + 1)}
F = QZ(); tw = []                                                                         # T'_m = Sum_(j>=1) (zeta^j - 1)^(1-m)
for m in range(0, D + 1):
    s = F.zero
    for j in range(1, 5):
        b = F.sub(F.zeta(j), F.one); ib = F.inv(b); p = F.one
        for _ in range(m - 1): p = F.mul(p, ib)
        s = F.add(s, p if m >= 1 else F.mul(b, F.one))
    tw.append(F.rat(s))
lm1 = [Fr(0)] + Lq[1:]                                                                    # L - 1
lg = [Fr(0)] + [Fr(3125) ** k / k for k in range(1, D)]                                   # 5 log L = -log(1 - 3125 q)
Dp = [lg[d] - 4 * lm1[d] for d in range(D)]; pw = lm1[:]
for m in range(2, D):
    pw = smul(pw, lm1, DQ)
    Dp = [Dp[d] - pw[d] * tw[m] / m for d in range(D)]
dDp = [d * Dp[d] for d in range(D)]
check('D_q Dpsi = L - 1 (critical value of the chart phase, D_q = q d/dq), through q^%d' % DQ, dDp == lm1)
# 3. both sides in Q((z))[[q]]: dictionary (d, n) -> coefficient, n >= -d
def lhs(withL=True, withD=True):
    R = {}; E = [Fr(1)] + [Fr(0)] * DQ; Lf = Lq if withL else [Fr(1)] + [Fr(0)] * DQ
    for j in range(DQ + 1 if withD else 1):
        for n in range(-j, NZ + 1):
            t = smul(smul(Lf, E, DQ), cq[n + j], DQ)
            for d in range(D): R[(d, n)] = R.get((d, n), 0) + t[d]
        E = [x / (j + 1) for x in smul(E, Dp, DQ)]
    return R


def rhs(h):
    R = {}
    for d in range(D):
        num = [Fr(1)]
        for k in range(1, 5 * d + 1): num = smul(num, [Fr(5), Fr(k)], NZ + d)
        den = [Fr(1)]
        for i in range(1, d + 1): den = smul(den, [Fr(5), Fr(10 * i), Fr(10 * i * i), Fr(5 * i ** 3), Fr(i ** 4)], NZ + d)
        fact = Fr(1)
        for i in range(1, d + 1): fact *= i
        I = smul(num, spow_int(den, -1, NZ + d), NZ + d)                                  # z^d I_d(z) as a power series
        t = smul(h, I, NZ + d)
        for n in range(-d, NZ + 1): R[(d, n)] = t[n + d] / fact
    return R


h = h_newton(NZ + DQ); A = lhs(); B = rhs(h)
check('step (e): L e^(Dpsi/z) c(z; L(q)) = h(z) Sum_d q^d I_d(z), all coefficients q^d z^n, d <= %d, -d <= n <= %d' % (DQ, NZ),
      all(A[k] == B[k] for k in B), sum(1 for k in B if B[k] != 0))
A1 = lhs(withL=False); fl = [k for k in sorted(B) if A1.get(k, 0) != B[k]]
check('control: factor L dropped fails (first at q^1)', fl[:1] and fl[0][0] == 1, fl[:2])
A2 = lhs(withD=False); fd = [k for k in sorted(B) if A2.get(k, 0) != B[k]]
check('control: e^(Dpsi/z) dropped fails (first at q^1)', fd[:1] and fd[0][0] == 1, fd[:2])
hw = smul(h, sexp_q([Fr(0), Fr(-1, 30)] + [Fr(0)] * (NZ + DQ - 1), NZ + DQ), NZ + DQ)   # log-coefficient -3/20 -> -11/60
Bw = rhs(hw); fh = [k for k in sorted(B) if A[k] != Bw[k]]
check('control: h with the z^1 log-coefficient of twist +5 (-11/60 for -3/20) fails at q^0 z^1', fh[:1] == [(0, 1)], fh[:2])
out('CH_EEXP %d %d: %s   (%.0f s)' % (DQ, NZ, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
