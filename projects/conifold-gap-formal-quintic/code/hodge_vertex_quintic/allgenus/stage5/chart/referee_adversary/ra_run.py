# -*- coding: utf-8 -*-
"""ra_run.py -- stage5/chart/referee_adversary: adversarial tests of the chart argument (CHART_PROOF s.3) by Laurent
expansion at the colour point (ra_sp.py).  Usage: py -I -B ra_run.py NEX NP   (exact order NEX over Q, order NP mod PR1, PR2)
Sections: E engine cross-checks; C the identity (Gamma) itself; P predictions of the argument for modified integrands to which
it applies (must hold); I modifications of a REGULAR factor that (L3) does not cover (the residue picture predicts them);
X falsifiers: modifications of the SINGULAR factor or of its Euler structure (the argument must fail, and the identity too).
Log: logs/ra_run_NEX_NP.log (relative paths only).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from ra_core import FQ, FP, RSym, RNum, PR1, PR2, stirling_b, bernoulli
from ra_sp import laurent, predict, sexp
from ra_models import quintic, colour_i, rational_weights, core_lemma, orbifold, XPsumM, XPsumZeta, single
from ra_sp import Group, Sing

A = [int(x) for x in sys.argv[1:]] + [None, None]
NEX = A[0] if A[0] is not None else 6
NP = A[1] if A[1] is not None else 8
LOGREL = 'logs/ra_run_%d_%d.log' % (NEX, NP)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()
B = stirling_b(40)


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def h_closed(F, N, route, d=5, colour_twist=None):
    """exp(Sum_k b_k z^(2k-1) [d^(1-2k) + Sum_j e_j^(1-2k)]), e_j = zeta^j - 1, by the given power-sum route"""
    xp = XPsumM() if route == 'm' else XPsumZeta(F)
    L = [F.zero] * (N + 1)
    for k in range(1, N // 2 + 2):
        if 2 * k - 1 > N:
            break
        T = xp.ps(F, 2 * k - 1)
        if d:
            T = F.red(T + F.c(Fr(d) ** (1 - 2 * k)))
        L[2 * k - 1] = F.red(F.c(B[k]) * T)
    return sexp(F, L)


def mulser(F, a, b):
    N = len(a) - 1
    return [F.red(sum(a[i] * b[n - i] for i in range(n + 1))) for n in range(N + 1)]


def analyse(spec, F, N):
    """Laurent data at sigma*: (poles = {n: pole order}, values c_n(sigma*), prediction)"""
    Nn = laurent(spec, RSym(F, N + 1), N)
    poles = {}
    for n in range(N + 1):
        nz = [k for k in range(n) if not F.iszero(Nn[n][k])]
        if nz:
            poles[n] = n - min(nz)
    vals = [Nn[n][n] for n in range(N + 1)]
    return poles, vals, predict(spec, F, N), Nn


def same(F, a, b):
    return all(F.iszero(F.red(x - y)) for x, y in zip(a, b))


def case(name, spec, F, N, closed=None, expect=True, pspec=None):
    t = time.time(); poles, vals, pr, Nn = analyse(spec, F, N)
    if pspec is not None:
        pr = predict(pspec, F, N)
    okval = same(F, vals, pr) and not poles
    okc = True if closed is None else same(F, pr, closed)
    if expect:
        check('%s [%s, z^%d]: regular at the colour point, value = residue prediction%s' % (
            name, F.name, N, '' if closed is None else ' = closed form'), okval and okc,
            '%.1f s' % (time.time() - t) if okval and okc else 'poles %s, first mismatch %s, closed ok %s' % (
                poles, next((n for n in range(N + 1) if not F.iszero(F.red(vals[n] - pr[n]))), None), okc))
    else:
        n1 = min(poles) if poles else None
        mm = next((n for n in range(N + 1) if not F.iszero(F.red(vals[n] - pr[n]))), None)
        check('%s [%s, z^%d]: identity FAILS as expected' % (name, F.name, N), bool(poles) or mm is not None,
              'first pole at z^%s (order %s; poles %s)' % (n1, poles.get(n1), poles) if poles else 'no pole; first value mismatch z^%s' % mm)
    return vals, pr


FQ_ = FQ(); P1, P2 = FP(PR1), FP(PR2)
out('ra_run: exact Q through z^%d, mod PR1 = %d and PR2 = %d through z^%d; log %s' % (NEX, PR1, PR2, NP, LOGREL))

# ---- E: engine cross-checks
c2 = laurent(quintic(), RNum(FQ_, 1), 5)
check('E1 c_n(2), n <= 5, by evaluation at L = 2 = worker/Stage-4 stored values',
      [str(x) for x in c2[1:]] == ['-12/5', '72/25', '354278/625', '-73040834/3125', '11495493197/31500'], [str(x) for x in c2[1:]])
cm = laurent(quintic(), RNum(FQ_, Fr(-5, 3)), 3); cm = [x / Fr(-5, 3) ** n for n, x in enumerate(cm)]   # c_n = eps^-n N_n
check('E1 c_n(-2/3), n <= 3, by evaluation at L = -2/3 = worker values', [str(x) for x in cm[1:]] ==
      ['-4/135', '8/18225', '4067314/110716875'], [str(x) for x in cm[1:]])
K = 7; Nn = laurent(quintic(), RSym(FQ_, K), 2)
ok = all(Nn[1][1 + e] == Fr(-3, 20) * comb(4, e) for e in range(K - 1)) and \
    all(Nn[2][2 + e] == Fr(9, 800) * comb(8, e) for e in range(K - 2))
check('E2 Laurent data at L = 1: c_1 = -(3/20)(1+eps)^4, c_2 = (9/800)(1+eps)^8 (all available eps-orders)', ok)
for Fp in (P1, P2):
    a = laurent(quintic('m'), RSym(Fp, NP + 1), NP); b = laurent(quintic('zeta', Fp), RSym(Fp, NP + 1), NP)
    check('E3 power-sum routes m and zeta give identical Laurent data [%s, z^%d]' % (Fp.name, NP), a == b)

# ---- C: the identity (Gamma) itself
hQ = h_closed(FQ_, NEX, 'm')
check('C0 h_1..h_5 = stored values (GAMMA s.1, POLY s.1)', [str(x) for x in hQ[1:6]] ==
      ['-3/20', '9/800', '-269/80000', '2823/6400000', '-257071/8064000000'] if NEX >= 5 else True)
case('C1 (Gamma): quintic, colour 0', quintic(), FQ_, NEX, closed=hQ)
for Fp in (P1, P2):
    case('C2 (Gamma): quintic, colour 0, zeta route', quintic('zeta', Fp), Fp, NP, closed=h_closed(Fp, NP, 'zeta'))
    hz = h_closed(Fp, NP, 'zeta'); zi = Fp.zeta5()
    for i in (1, 2):
        rot = [Fp.red(hz[n] * pow(Fp.inv(pow(zi, i, Fp.p)), n, Fp.p)) for n in range(NP + 1)]
        case('C3 colour %d at sigma* = zeta^%d: value h(zeta^-%d z)' % (i, i, i), colour_i(Fp, i), Fp, NP, closed=rot)

# ---- P: predictions where the chart argument applies
N = NEX; z = [FQ_.zero] * (N + 1)


def poly(coeffs):
    s = [FQ_.zero] * (N + 1)
    for i, c in enumerate(coeffs[:N + 1]):
        s[i] = Fr(c)
    return s


def poch(e, M):
    """Poch_M(e z) as a series"""
    s = poly([1])
    if M >= 0:
        for i in range(M):
            s = mulser(FQ_, s, poly([1, Fr(i) * e]))
    else:
        for i in range(1, -M + 1):
            s = mulser(FQ_, s, [Fr(e * i) ** n for n in range(N + 1)])
    return s


for My in (0, 2, 3):
    case('P1 Gamma(%d + 5s) in place of Gamma(1 + 5s)' % My, quintic(My=My), FQ_, N, closed=mulser(FQ_, poch(Fr(1, 5), My), hQ))
mz = poly([5, 10, 10, 5, 1]); mzm = poly([5, -10, 10, -5, 1])        # m(z), m(-z)


def _inv(s):
    o = [Fr(0)] * (N + 1); o[0] = 1 / s[0]
    for n in range(1, N + 1):
        o[n] = -sum(s[i] * o[n - i] for i in range(1, n + 1)) / s[0]
    return o


P6 = poch(Fr(1, 5), 6)
case('P2 Gamma(1 + (1 - sigma)/z) for the j = 0 factor (M0 = 1)', quintic(M0=1), FQ_, N,
     closed=mulser(FQ_, mulser(FQ_, P6, _inv([x / 5 for x in mz])), hQ))
case('P2 Gamma(-1 + (1 - sigma)/z) for the j = 0 factor (M0 = -1)', quintic(M0=-1), FQ_, N,
     closed=mulser(FQ_, poch(Fr(1, 5), -4), hQ))
case('P2 Gamma(2 + (1 - sigma)/z) for the j = 0 factor (M0 = 2)', quintic(M0=2), FQ_, N)
case('P3 Gamma(1 + (zeta^j - sigma)/z), j = 1..4 (Mx = 1)', quintic(Mx=1), FQ_, N, closed=hQ)
case('P3 Gamma(2 + (zeta^j - sigma)/z), j = 1..4 (Mx = 2)', quintic(Mx=2), FQ_, N,
     closed=mulser(FQ_, [x / 5 for x in mzm], hQ))
case('P3 Gamma(-1 + (zeta^j - sigma)/z), j = 1..4 (Mx = -1)', quintic(Mx=-1), FQ_, N,
     closed=mulser(FQ_, _inv([x / 5 for x in mz]), hQ))
for d in (0, 3, 7):
    case('P4 degree %d twist Gamma(1 + %ds) (fiat Stirling; d = 0 is equivariant P^4)' % (d, d), quintic(d=d), FQ_, N,
         closed=h_closed(FQ_, N, 'm', d=d))
case('P5 amplitude x e^(2 sigma)', quintic(fextra={(0, 1): 2}), FQ_, N, closed=hQ)
pe, ps = single(1, 4)
case('P5 amplitude x 1/(3 + sigma)', quintic(extra=[Group('r', 1, pe, ps, phase=False, bk=None, expo=-1)]), FQ_, N, closed=hQ)
pe, ps = single(1, 3)
case('P5 amplitude x (2 + sigma)', quintic(extra=[Group('r', 1, pe, ps, phase=False, bk=None, expo=1)]), FQ_, N, closed=hQ)
case('P5 amplitude x exp(z sigma^2)', quintic(fextra={(1, 0): 1, (1, 1): 2, (1, 2): 1}), FQ_, N,
     closed=mulser(FQ_, sexp(FQ_, poly([0, 1])), hQ))
for g in ({2: Fr(3, 7)}, {3: Fr(-2, 5)}, {2: 1, 3: -1, 4: 2}):
    case('P6 extra phase g = %s in powers of (sigma - 1)' % {k: str(v) for k, v in g.items()}, quintic(gstar=g), FQ_, N, closed=hQ)
case('P7 rational weights (1; 3, -1/2, 5/2, -2), d = 5', rational_weights([3, Fr(-1, 2), Fr(5, 2), -2]), FQ_, N)
case('P7 rational weights (1; 2, -3), d = 2', rational_weights([2, -3], d=2), FQ_, N)
gR = {2: Fr(1, 3), 3: -2, 5: Fr(1, 2)}; fR = {(0, 1): 1, (0, 2): Fr(-1, 2), (1, 0): 2, (1, 1): Fr(1, 5), (2, 3): -3}
for a0, M0 in ((-1, 0), (-1, 1), (-1, -1), (2, 0), (2, 1), (Fr(-1, 3), 0), (Fr(-1, 3), 2)):
    case('P8 core lemma: one singular factor (alpha0 %s, M0 %d), random g, F' % (a0, M0), core_lemma(a0, M0, gR, fR), FQ_, N)
case('P9 quintic with Gamma(2(1 - sigma)/z) for the j = 0 factor', quintic(alpha0=-2), FQ_, N, closed=hQ)

# ---- I: regular factors outside (L3) (residue picture)
case('I1 y-factor b_1 + 1/7', quintic(y_bk={1: Fr(1, 7)}), FQ_, N)
case('I2 x-factors b_2 + 3', quintic(x_bk={2: 3}), FQ_, N)
case('I3 x-factors: Stirling series dropped', quintic(x_bk=None), FQ_, N)
case('I4 y-factor exponent +1/2 for -1/2', quintic(y_expo=Fr(1, 2)), FQ_, N)

# ---- X: falsifiers (the identity must fail)
for nm, kw in (('X1 singular factor: Stirling series dropped', dict(sing_bk=None)),
               ('X2 singular factor: b_1 + 1/7', dict(sing_bk={1: Fr(1, 7)})),
               ('X3 all six factors: b_1 + 1/7 (prediction with modified regular b_1)',
                dict(sing_bk={1: Fr(1, 7)}, y_bk={1: Fr(1, 7)}, x_bk={1: Fr(1, 7)})),
               ('X4 singular factor squared: Gamma((1 - sigma)/z)^2', dict(count=2)),
               ('X5 singular factor: exponent 0 for -1/2 (no Gamma function)', dict(sing_expo=Fr(0))),
               ('X7 singular factor: b_3 + 1/7 (visible from z^5)', dict(sing_bk={3: Fr(1, 7)})),
               ('X8 singular factor: b_2 + 1/7', dict(sing_bk={2: Fr(1, 7)}))):
    if 'b_3' in nm and N < 5:
        continue
    case(nm, quintic(**kw), FQ_, N, expect=False)
# two 'falsifiers' of the first draft that are NOT falsifiers: identities of the Gamma function turn them into covered cases
case('P2b singular exponent +1/2 is Gamma(1 + w) = w Gamma(w): same as M0 = 1', quintic(sing_expo=Fr(1, 2)), FQ_, N,
     pspec=quintic(M0=1))
q10 = quintic(); q10['sing'] = Sing(alpha0=1, M0=1)
case('P10 1/Gamma((1 - sigma)/z) = Gamma(1 + (sigma - 1)/z) x e^(linear)/const (reflection): residue at sigma = 1 - z',
     quintic(count=-1), FQ_, N, pspec=q10)
# ---- O: the orbifold point L = 0 (q = oo): singular factor Gamma(1 + 5 sigma/z), residue at sigma = -z/5
Bn = bernoulli(5 * (N // 5) + 2)
def bpoly(n, x):
    return sum(comb(n, k) * Bn[k] * Fr(x) ** (n - k) for k in range(n + 1))
GL = [Fr(0)] * (N + 1)
for l in range(1, N // 5 + 1):
    GL[5 * l] = Fr((-1) ** (l + 1)) * bpoly(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1))
case('O1 orbifold point L = 0: c(z; 0) = Ginf(z) = exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1))) (Stage-4 closed form)',
     orbifold(), FQ_, N, closed=sexp(FQ_, GL))
out('RA_RUN %d %d: %s   (%.0f s)' % (NEX, NP, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
