# -*- coding: utf-8 -*-
"""lm_run_x.py N -- referee lmi_adversary: integrands that violate ONE hypothesis of Theorem LMI (exact over Q).
Each check PASSes when the analogue of LMI fails in the PREDICTED way (pole order, leading coefficient, log-dependence),
or, for the positive controls marked (+), when the identity holds although the integrand looks dangerous."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import random
from fractions import Fraction as Fr
from math import factorial
from lm_core import LS, lsc, F, fr, fpow, setmod
from lm_wick import wick1, beta
from lm_data import Sing, cseries, lmi_value
from lm_models import rand_data
from lm_fals import nonlinear_data, lam_extra, qdep_extra
from lm_harness import Run, regular, poles, values, first_pole

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
setmod(None)
R = Run('lm_run_x_%d' % N)
rng = random.Random(1411)


def lead(cs, n):
    x = cs[n]
    return (-x.v, x.c[0]) if (x.c and x.v < 0) else (0, None)


def dfact(k):
    r = 1
    for i in range(1, 2 * k, 2):
        r *= i
    return r


# X1-X3: two singular factors at the same point (hypothesis: g holomorphic at s*)
G, f = rand_data(rng, 6, 5, 3)
for al in (-1, 2):
    cs = cseries(N, [Sing(F(al), 0), Sing(F(al), 0)], G, f)
    kap = [Fr(dfact(n) ** 2, factorial(n) * 16 ** n) * Fr(-1) ** n for n in range(N + 1)]   # alpha0 = -1 gives GEN's
    ok = all(lead(cs, n) == (n, F(kap[n] * Fr(1, al) ** n)) for n in range(1, N + 1))
    R.check('X1 Gamma(a/z)^2, alpha0=%d: pole of order n at every z^n, leading ((2n-1)!!)^2/(n! 16^n) (-1/alpha0)^n' % al,
            ok, 'poles %s' % poles(cs))
cs = cseries(N, [Sing(F(-1), 0), Sing(F(1), 0)], G, f)
R.check('X2 Gamma(a/z) Gamma(-a/z) (phases cancel): poles (observed order 2n)', not regular(cs), 'poles %s' % poles(cs))
cs = cseries(N, [Sing(F(-1), 0), Sing(F(-2), 0)], G, f)
R.check('X3 Gamma(a/z) Gamma(2a/z): poles', not regular(cs), 'poles %s' % poles(cs))
cs = cseries(N, [Sing(F(-1), 0), Sing(F(-1), 1)], G, f)
R.check('X3b Gamma(a/z) Gamma(1 + a/z) = (a/z) Gamma(a/z)^2: poles', not regular(cs), 'poles %s' % poles(cs))

# X4: nonlinear singular argument a = alpha0 x + alpha1 x^2 (hypothesis: a0 linear); ell = log a(L) is a free parameter
for M0 in (0, 1):
    res = {}
    for ell in (F(0), F(1), F(Fr(5, 2))):
        p, a, C = nonlinear_data(N, F(-1), F(Fr(1, 3)), M0, G, f, ell)
        res[ell] = wick1(p, a, N, C)
    r0, r1, r2 = res[F(0)], res[F(1)], res[F(Fr(5, 2))]
    reg = all(regular(r) for r in res.values())
    same0 = values(r0) == values(r1) == values(r2)
    naive = lmi_value(N, F(-1), F(M0), G, f)
    if M0 == 0:
        dep = [n for n in range(1, N + 1) if r0[n].coef(1) != r1[n].coef(1) or r0[n].coef(1) != r2[n].coef(1)]
        R.check('X4 nonlinear a, M0=0: regular at fixed ell, [e^1] c_n depends on ell = log a(L): c_n has e log e terms, '
                'NOT holomorphic at s*', reg and len(dep) > 0, 'ell-dependent at e^1 for n in %s' % dep)
        R.check('X4 nonlinear a, M0=0 (observation of draft 1, now checked): the limit exists, is ell-free and equals '
                'F(s*,z)/F(s*,0)', same0 and values(r0) == naive)
    else:
        dep = [n for n in range(1, N + 1) if r0[n].coef(0) != r1[n].coef(0)]
        R.check('X4 nonlinear a, M0=1: [e^0] c_n depends on ell = log a(L): c_n diverges like log(L - s*), no limit',
                reg and len(dep) > 0 and not same0, 'ell-dependent at e^0 for n in %s' % dep)

# X5: amplitude depending on lambda = -g'(L) - alpha0 log a0(L) (data not holomorphic in Q)
v = {}
for ell in (F(0), F(1)):
    cs = cseries(N, [Sing(F(-1), 0)], G, f, extraH=lam_extra(G, F(-1), ell, F(Fr(2, 3)), [F(1), F(1)], N))
    v[ell] = cs
diff = fr(v[F(1)][1].coef(0) - v[F(0)][1].coef(0))
R.check('X5 F = exp((2/3) z lambda (1 + x)): c_1(s*) depends on ell: d c_1/d ell = (2/3)(-alpha0) = 2/3, log-divergent',
        diff == F(Fr(2, 3)), 'difference %s' % diff)
# X5+: amplitude depending holomorphically on Q (positive control)
for M0 in (0, 1, -1):
    f1 = [F(rng.randint(-5, 5)) for _ in range(4)]
    cs = cseries(N, [Sing(F(-1), M0)], G, f, extraH=qdep_extra(G, F(-1), f1, N, N + 8))
    R.check('X5+ (+) F_Q = exp(z Q f1(sigma)), Q-holomorphic, M0=%d: no pole, value = formula with F_Q at Q = 0' % M0,
            regular(cs) and values(cs) == lmi_value(N, F(-1), F(M0), G, f))

# X6-X7: the singular factor is not an exact Stirling form (hypothesis: A = Stirling form of Gamma(M0 + a0/z))
for al in (-1, 2):
    for n0, eta in ((1, Fr(1, 7)), (2, Fr(1, 7)), (3, Fr(-2, 3))):
        if n0 > N:
            continue
        cs = cseries(N, [Sing(F(al), 0, stir=lambda n, n0=n0, eta=eta: beta(n, 0) + (eta if n == n0 else 0))], G, f)
        fp = first_pole(cs)
        R.check('X6 alpha0=%d: (z/a)^%d coefficient + %s: first pole at z^%d, order %d, coefficient eta/alpha0^%d'
                % (al, n0, eta, n0, n0, n0), fp == (n0, n0, F(eta * Fr(1, al) ** n0)), 'first pole %s' % (fp,))
    cs = cseries(N, [Sing(F(al), Fr(1, 3), stir=lambda n: beta(n, 0))], G, f)
    fp = first_pole(cs)
    R.check('X6 alpha0=%d: elementary part of Gamma(1/3 + w), Bernoulli part of Gamma(w): pole at z^1, coef 1/(9 alpha0)'
            % al, fp == (1, 1, F(Fr(1, 9 * al))), 'first pole %s' % (fp,))
    for M, E in ((0, 0), (1, 1), (Fr(1, 3), Fr(1, 3)), (0, Fr(1, 2))):
        cs = cseries(N, [Sing(F(al), M, expo=E)], G, f)
        pred = (Fr(1, 24) + beta(1, Fr(M)) - Fr(E) ** 2 / 2) / al
        fp = first_pole(cs)
        if pred != 0:
            R.check('X7 alpha0=%d: Gamma(%s + w) with exponent a^%s: pole at z^1, coef (1/24 + beta_1 - E^2/2)/alpha0 = %s'
                    % (al, M, E, pred), fp is not None and fp[0] == 1 and lead(cs, 1) == (1, F(pred)), 'first pole %s' % (fp,))
        else:
            R.check('X7 (+) alpha0=%d: exponent a^(+1/2) with M = 0 is Gamma(1 + w) = w Gamma(w): regular, LMI value of M0=1'
                    % al, regular(cs) and values(cs) == lmi_value(N, F(al), F(1), G, f))

# X8: amplitude with a pole at s* (hypothesis: F holomorphic); observed n! (-1/alpha0)^n (E_1-type divergent series)
for al in (-1, 2):
    cs = cseries(N, [Sing(F(al), 0)], G, f, pows=[-1])
    ok = all(lead(cs, n) == (n, F(factorial(n) * Fr(-1, al) ** n)) for n in range(1, N + 1))
    R.check('X8 alpha0=%d, F = 1/(sigma - s*): pole of order n at z^n, leading n! (-1/alpha0)^n' % al, ok,
            'leading %s' % [lead(cs, n)[1] for n in range(1, min(N, 6) + 1)])

# X9: F(s*, 0) = 0 (hypothesis F(s*, 0) != 0): F = (sigma - s*) F~
for al, M0 in ((-1, 0), (3, 0)):
    cs = cseries(N, [Sing(F(al), M0)], G, f, pows=[1])
    R.check('X9 (+) alpha0=%d M0=0, F = (sigma - s*) F~: regular, value = LMI with M0 = 1 (shift rule)' % al,
            regular(cs) and values(cs) == lmi_value(N, F(al), F(1), G, f))
for al, M0 in ((-1, 1), (3, 2), (2, -1)):
    cs = cseries(N, [Sing(F(al), M0)], G, f, pows=[1])
    V = lmi_value(N, F(al), F(M0), G, f)
    k = fr(-F(M0) * F(Fr(1, al)))
    simple = all((not x.c) or x.v >= -1 for x in cs)
    resid = all(cs[n].coef(-1) == fr(k * V[n - 1]) for n in range(1, N + 1))
    R.check('X9 alpha0=%d M0=%d, F = (sigma - s*) F~: only simple poles, [e^-1] c_n = (-M0/alpha0) V_(n-1) (chart + (L2\'))'
            % (al, M0), simple and resid and not regular(cs), 'poles %s' % poles(cs))

# X10: value controls (regular, wrong value)
cs = cseries(N, [Sing(F(-1), 1)], G, f)
R.check('X10 M0 = 1 (measure dx) read with the M0 = 0 value (delta = 0): value differs',
        regular(cs) and values(cs) != lmi_value(N, F(-1), F(0), G, f))
R.check('X10 M0 = 1: delta with the wrong sign (+M0 z/alpha0): value differs',
        values(cs) != lmi_value(N, F(-1), F(-1), G, f))
R.done()
