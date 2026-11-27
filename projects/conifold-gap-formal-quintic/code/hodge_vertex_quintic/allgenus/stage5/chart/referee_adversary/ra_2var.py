# -*- coding: utf-8 -*-
"""ra_2var.py -- stage5/chart/referee_adversary: the TWO-variable chart (referee's simplification, REFEREE.md s.3), tested.

Only the singular factor Gamma_St(M0 + a0/z), a0 = alpha0 (sigma - sigma*), is written as an Euler integral, x = e^nu; the
regular factors stay in Stirling form, as a phase g(sigma) and an amplitude F(sigma, z).  After nu = mu - lambda/alpha0:
    Phi_2(sigma, mu) = g(sigma) + alpha0 (sigma - sigma*) mu - Q e^mu,      amplitude e^(M0 mu) F(sigma, z),
with critical point (L, mu_c), Q e^(mu_c) = alpha0 (L - sigma*).  In the scaled variables sigma = L + eps tau, z = eps z'
(eps = L - sigma*), mu = mu_c + m, the integrand is
    exp([Sum_(k>=2) eps^(k-1) g_k tau^k + alpha0 tau m - alpha0 Sum_(k>=2) m^k/k!]/z') e^(M0 m) F(L + eps tau, eps z').
Claim tested (Fubini (L1b) + (L3) for the ONE Euler integral): its normalised 2-variable Wick series equals the normalised
one-variable series of the full Stirling summand (ra_sp.laurent), order by order in z' and in eps.
The 2-variable Wick expectation uses <eta1^a eta2^b> for the 2x2 covariance C = -H^-1 directly (no Fubini inside the code).
Usage: py -I -B ra_2var.py N   (log logs/ra_2var_N.log)
"""
import os, sys, time, copy
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb, factorial
from ra_core import FQ, FP, RSym, RNum, PR1
from ra_sp import build, bexp, laurent, Sing
from ra_models import quintic, orbifold, core_lemma

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
LOGREL = 'logs/ra_2var_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def pmul(R, a, b, D):
    """product of r-series with eta-polynomial coefficients: a[k] = list over eta-degree; truncate r^D"""
    out = [[] for _ in range(D + 1)]
    for k1, p in enumerate(a):
        for k2, q in enumerate(b):
            if k1 + k2 > D or not p or not q:
                continue
            o = out[k1 + k2]
            need = len(p) + len(q) - 1
            if len(o) < need:
                o.extend(R.zero() for _ in range(need - len(o)))
            for i, x in enumerate(p):
                if R.iszero(x):
                    continue
                for j, y in enumerate(q):
                    o[i + j] = R.add(o[i + j], R.mul(x, y))
    return out


def pexp(R, V, D):
    """exp of an r-series V (V[0] = 0, each V[k] an eta-polynomial), truncated at r^D: E = Sum V^n/n!"""
    E = [[R.const(1)]] + [[] for _ in range(D)]
    term = [[R.const(1)]] + [[] for _ in range(D)]
    for n in range(1, D + 1):
        term = pmul(R, term, V, D)
        term = [[R.scal(Fr(1, n), x) for x in p] for p in term]
        if all(not p for p in term):
            break
        for k in range(D + 1):
            p = term[k]
            if p:
                if len(E[k]) < len(p):
                    E[k].extend(R.zero() for _ in range(len(p) - len(E[k])))
                for i, x in enumerate(p):
                    E[k][i] = R.add(E[k][i], x)
    return E


def moments(R, C11, C12, C22, amax, bmax):
    """<eta1^a eta2^b> for the 2x2 covariance (Isserlis)"""
    M = {}
    for a in range(amax + 1):
        for b in range(bmax + 1):
            s = R.zero()
            for c in range(0, min(a, b) + 1):
                if (a - c) % 2 or (b - c) % 2:
                    continue
                i, j = (a - c) // 2, (b - c) // 2
                cnt = Fr(factorial(a) * factorial(b), factorial(c) * factorial(i) * factorial(j) * 2 ** (i + j))
                t = R.scal(cnt, R.const(1))
                for _ in range(c):
                    t = R.mul(t, C12)
                for _ in range(i):
                    t = R.mul(t, C11)
                for _ in range(j):
                    t = R.mul(t, C22)
                s = R.add(s, t)
            M[(a, b)] = s
    return M


def twovar(spec, R, N):
    S = spec['sing']; a0, M0 = S.alpha0, S.M0
    reg = copy.copy(spec); reg['sing'] = Sing(alpha0=a0, count=0)
    Fk, LA = build(reg, R, N); A = bexp(R, LA, N); D = 2 * N
    # tau-part: exp(Sum_(k>=3) Fk r^(k-2) eta^k) * Sum A[m][j] r^(2m+j) eta^j
    V1 = [[] for _ in range(D + 1)]
    for k in range(3, len(Fk)):
        if k - 2 <= D:
            V1[k - 2] = [R.zero()] * k + [Fk[k]]
    Am = [[] for _ in range(D + 1)]
    for m in range(N + 1):
        for j, x in enumerate(A[m]):
            if 2 * m + j <= D:
                p = Am[2 * m + j]
                if len(p) < j + 1:
                    p.extend(R.zero() for _ in range(j + 1 - len(p)))
                p[j] = R.add(p[j], x)
    T = pmul(R, pexp(R, V1, D), Am, D)
    # mu-part: exp(-alpha0 Sum_(k>=3) r^(k-2) eta^k / k!) * exp(M0 r eta)
    V2 = [[] for _ in range(D + 1)]
    for k in range(3, D + 3):
        V2[k - 2] = [R.zero()] * k + [R.const(-a0 / factorial(k))]
    if M0:                                          # the source term M0 r eta (from e^(M0 m)); V2[1] holds the cubic term
        V2[1] = V2[1] + [R.zero()] * max(0, 2 - len(V2[1]))
        V2[1][1] = R.add(V2[1][1], R.const(M0))
    Mp = pexp(R, V2, D)
    # covariance C = -H^-1, H = [[2 F2, alpha0], [alpha0, -alpha0]]
    h11, h12, h22 = R.scal(2, Fk[2]), R.const(a0), R.const(-a0)
    det = R.sub(R.mul(h11, h22), R.mul(h12, h12)); idet = R.inv(det)
    C11 = R.scal(-1, R.mul(h22, idet)); C22 = R.scal(-1, R.mul(h11, idet)); C12 = R.mul(h12, idet)
    amax = max((len(p) for p in T), default=1); bmax = max((len(p) for p in Mp), default=1)
    Mo = moments(R, C11, C12, C22, amax, bmax)
    res = []
    for n in range(N + 1):
        s = R.zero()
        for k1 in range(0, 2 * n + 1):
            k2 = 2 * n - k1
            for a, x in enumerate(T[k1]):
                if R.iszero(x):
                    continue
                for b, y in enumerate(Mp[k2]):
                    if (a + b) % 2 == 0 and not R.iszero(y):
                        s = R.add(s, R.mul(R.mul(x, y), Mo[(a, b)]))
        res.append(s)
    i0 = R.inv(res[0])
    return [R.mul(x, i0) for x in res]


FQ_ = FQ()
out('ra_2var: two-variable chart (one Euler integral) against the one-variable full-Stirling series, through z^%d; log %s'
    % (N, LOGREL))
gR = {2: Fr(1, 3), 3: -2, 5: Fr(1, 2)}; fR = {(0, 1): 1, (0, 2): Fr(-1, 2), (1, 0): 2, (1, 1): Fr(1, 5), (2, 3): -3}
cases = [('quintic, colour 0', quintic()), ('quintic with Gamma(1 + (1 - sigma)/z) (M0 = 1)', quintic(M0=1)),
         ('orbifold point L = 0 (singular factor Gamma(1 + 5 sigma/z))', orbifold()),
         ('core lemma, random g, F (alpha0 = -1, M0 = 0)', core_lemma(-1, 0, gR, fR)),
         ('core lemma, random g, F (alpha0 = 2, M0 = 1)', core_lemma(2, 1, gR, fR))]
for nm, spec in cases:
    t = time.time()
    for R, tag in ((RNum(FQ_, Fr(1, 2)), 'eps = 1/2'), (RNum(FQ_, Fr(-3, 7)), 'eps = -3/7'), (RSym(FQ_, N + 2), 'eps-series')):
        a = twovar(spec, R, N); b = laurent(spec, R, N)
        b = [R.mul(x, R.inv(b[0])) for x in b]
        check('%s [%s, %s]: 2-variable chart = 1-variable full-Stirling series through z^%d' % (nm, FQ_.name, tag, N),
              all(R.iszero(R.sub(x, y)) for x, y in zip(a, b)))
    out('     (%.1f s)' % (time.time() - t))
# control: the singular factor's Stirling series dropped in the 1-variable form must NOT match the 2-variable chart
spec = quintic(); R = RNum(FQ_, Fr(1, 2))
a = twovar(spec, R, N); b = laurent(quintic(sing_bk=None), R, N)
check('control: 1-variable form without the singular Stirling series differs from the chart (eps = 1/2)',
      not all(R.iszero(R.sub(x, y)) for x, y in zip(a, b)))
out('RA_2VAR %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not FAILS else 'FAILURES: %s' % FAILS, time.time() - T0))
