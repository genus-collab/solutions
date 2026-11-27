# -*- coding: utf-8 -*-
"""lm_run_eps.py N FIELD -- referee lmi_adversary: engine checks (E), LMI predictions (P), the proof steps (S).
FIELD = Q (exact rationals) or p1 / p2 (two primes = 1 mod 20, found by lm_harness.prime_below).
E: stored quintic numbers; (L3) for integer AND non-integer M; the one-Gamma identity.
P: Theorem LMI on random data (no pole, value = RA's formula), integer and non-integer M0, g''(s*) = 0, complex alpha0.
S: step 1+2 (the 1-variable series equals the two-variable chart at every L = s* + e, coefficientwise), step 3
   (chart at e = 0 equals (L2')), and two controls that must break exactly one step each."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import random
from fractions import Fraction as Fr
from lm_core import LS, lsc, lsm, F, fr, fpow, finv, setmod
from lm_wick import wick1, gexp, beta, bk
from lm_data import Sing, cseries, lmi_value, zexp
from lm_models import quintic_fixed, quintic_s1, rand_data, h_quintic
from lm_chart import chart_series
from lm_harness import Run, ls_equal, regular, poles, values, prime_below, sqrt_m1

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
FIELD = sys.argv[2] if len(sys.argv) > 2 else 'Q'
PR = {'Q': None, 'p1': prime_below(2 ** 62), 'p2': prime_below(2 ** 70)}[FIELD]
setmod(PR)
R = Run('lm_run_eps_%d_%s' % (N, FIELD))
R.say('field %s%s; z-order N = %d' % (FIELD, '' if PR is None else ' (p = %d)' % PR, N))

# ---- E: engine ------------------------------------------------------------------------------------------------------
st = [Fr(-12, 5), Fr(72, 25), Fr(354278, 625), Fr(-73040834, 3125), Fr(11495493197, 31500)]
p, a, C = quintic_fixed(2, 5)
cs = wick1(p, a, 5, C)
R.check('E1 quintic c_1..c_5(2) = stored values (SYN4 s.3 E3, RA E1)', values(cs)[1:] == [F(x) for x in st])
G, f = quintic_s1(3 * 6 + 12)
cs = cseries(6, [Sing(F(-1), 0)], G, f, exact=False)
e1 = LS(0, [F(1)] * 1)
ref1, ref2 = lsc(Fr(-3, 20)), lsc(Fr(9, 800))
for _ in range(4):
    ref1 = ref1 * LS(0, [F(1), F(1)])
for _ in range(8):
    ref2 = ref2 * LS(0, [F(1), F(1)])
ok1, n1 = ls_equal(cs[1], ref1)
ok2, n2 = ls_equal(cs[2], ref2)
R.check('E2 quintic c_1(1+e) = -(3/20)(1+e)^4, c_2 = (9/800)(1+e)^8', ok1 and ok2, '(%d, %d coefficients)' % (n1, n2))
R.check('E2b h_1..h_6 = -3/20, 9/800, -269/80000, 2823/6400000, -257071/8064000000, h_6',
        h_quintic(6)[1:6] == [F(x) for x in [Fr(-3, 20), Fr(9, 800), Fr(-269, 80000), Fr(2823, 6400000),
                                             Fr(-257071, 8064000000)]])


def euler_N(M, NN, fieldM=False):
    """N[log t - t, t^(M-1); 1] (L3 at a = 1, measure dt): phase p_k = (-1)^(k+1)/k, amplitude (1+x)^(M-1)"""
    pp = {k: lsc(Fr((-1) ** (k + 1), k)) for k in range(3, 2 * NN + 3)}
    m1 = fr(M - F(1)) if fieldM else F(Fr(M) - 1)
    H = {(0, j): lsc(Fr((-1) ** (j + 1), j)).scal(m1) for j in range(1, 2 * NN + 1)}
    return values(wick1(pp, gexp(H, 2 * NN), NN, lsc(1)))


def stir_beta(M, NN, fieldM=False):
    from lm_data import beta_field
    h = [F(0)] + [beta_field(n, M) if fieldM else F(beta(n, Fr(M))) for n in range(1, NN + 1)]
    return zexp(h, NN)


def poch_stir(M, NN, sign=1):
    s = [F(0)] * (NN + 1)
    for k in range(1, NN // 2 + 2):
        if 2 * k - 1 <= NN:
            s[2 * k - 1] = F(bk(k))
    E = zexp(s, NN)
    fac = [F(1)] + [F(0)] * NN
    rng_ = range(0, M) if M >= 0 else range(1, -M + 1)
    for i in rng_:
        if M >= 0:
            lin = [F(1), F(sign * i)]
        else:
            lin = [F(sign * i) ** 0] + [fpow(F(sign * i), k) for k in range(1, NN + 1)]   # 1/(1 - i z)
        fac = [fr(sum(fac[k - j] * lin[j] for j in range(0, min(k, len(lin) - 1) + 1))) for k in range(NN + 1)]
    return [fr(sum(fac[k - j] * E[j] for j in range(k + 1))) for k in range(NN + 1)]


NE = max(N, 10)
for M in range(-3, 5):
    v = euler_N(M, NE)
    R.check('E3 (L3) M = %d: SP of Euler = Poch_M exp(Sum b_k z^(2k-1)) = exp(Sum beta_n(M) z^n), z^%d' % (M, NE),
            v == poch_stir(M, NE) and v == stir_beta(M, NE))
for M in [Fr(1, 3), Fr(-7, 5), Fr(5, 2), Fr(2, 7)]:
    R.check('E3 (L3) non-integer M = %s: SP of Euler = exp(Sum beta_n(M) z^n), z^%d' % (M, NE),
            euler_N(M, NE) == stir_beta(M, NE))
R.check('E3 control: M = 1/3 against the M = 0 Bernoulli series differs', euler_N(Fr(1, 3), NE) != stir_beta(0, NE))
R.check('E3 control: M = 2 against Poch with (1 - z) differs', euler_N(2, NE) != poch_stir(2, NE, sign=-1))
if PR is not None:
    rg = random.Random(7)
    for _ in range(2):
        Mf = rg.randrange(PR)
        R.check('E3 (L3) generic M in F_p: SP of Euler = exp(Sum beta_n(M) z^n), z^%d' % NE,
                euler_N(Mf, NE, True) == stir_beta(Mf, NE, True))
for al, M in [(1, 0), (-1, 0), (3, 2), (Fr(-1, 2), -1), (2, Fr(1, 3))]:
    cs = cseries(N, [Sing(F(al), M)], [F(0)] * 3, {})
    R.check('E4 one-Gamma identity alpha=%s M=%s: c_n = 0 identically (n = 1..%d)' % (al, M, N),
            all(not x.c and x.prec >= 1 for x in cs[1:]), 'prec %s' % [x.prec for x in cs[1:]])

# ---- P: LMI on random data -----------------------------------------------------------------------------------------
rng = random.Random(20261114)
cases = [(-1, 0), (-1, 1), (-1, -2), (2, 0), (3, 1), (Fr(-1, 3), 2), (5, 1), (Fr(1, 2), -1), (Fr(-7, 4), 3),
         (-1, Fr(1, 3)), (Fr(2, 3), Fr(-2, 5)), (-3, Fr(7, 2))]
for al, M in cases:
    G, f = rand_data(rng, 6, 5, 3)
    cs = cseries(N, [Sing(F(al), M)], G, f)
    pred = lmi_value(N, F(al), F(M), G, f)
    R.check('P LMI alpha0=%s M0=%s random g (deg 6), F (z^3, x^5): no pole, value = formula' % (al, M),
            regular(cs) and values(cs) == pred, 'min prec %d' % min(x.prec for x in cs))
G, f = rand_data(rng, 6, 5, 3)
G[2] = F(0)
cs = cseries(N, [Sing(F(-2), 1)], G, f)
R.check('P LMI with g\'\'(s*) = 0 (GEN SKETCH needs g\'\' != 0)', regular(cs) and values(cs) == lmi_value(N, F(-2), F(1), G, f))
if PR is not None:
    i_ = sqrt_m1(PR)
    for al, Mf in [(i_, F(0)), (i_, F(2)), (rng.randrange(1, PR), rng.randrange(PR)), (rng.randrange(1, PR), rng.randrange(PR))]:
        G, f = rand_data(rng, 6, 5, 3)
        cs = cseries(N, [Sing(al, Mf, field_M=True)], G, f)
        R.check('P LMI over F_p: alpha0 %s, M0 %s generic' % ('= sqrt(-1)' if al == i_ else 'random', 'random' if Mf > 2 else Mf),
                regular(cs) and values(cs) == lmi_value(N, al, Mf, G, f))

# ---- S: the steps of the proof ----------------------------------------------------------------------------------------
for al, M in [(-1, 0), (2, 1), (Fr(-1, 3), -2), (-1, Fr(1, 3))]:
    G, f = rand_data(rng, 6, 5, 3)
    one = cseries(N, [Sing(F(al), M)], G, f)
    ch = chart_series(N, F(al), F(M), G, f)
    oks = [ls_equal(x, y) for x, y in zip(one, ch)]
    R.check('S1 steps 1+2: 1-variable c_n(s*+e) = 2-variable chart N_n(e), alpha0=%s M0=%s' % (al, M),
            all(o for o, _ in oks), 'coefficients compared per n: %s' % [c for _, c in oks])
    R.check('S2 chart data regular at e = 0 (no negative powers) and chart(0) = (L2\') value',
            regular(ch) and values(ch) == lmi_value(N, F(al), F(M), G, f))
G, f = rand_data(rng, 6, 5, 3)
ch = chart_series(N, F(-1), F(1), G, f, mcube=F(Fr(2, 7)))
R.check('S3 control: a pure m^3 term (2/7) m^3 breaks (L2\'): chart(0) != formula',
        values(ch) != lmi_value(N, F(-1), F(1), G, f))
one = cseries(N, [Sing(F(-1), 1, stir=lambda n: beta(n, 1) + (Fr(1, 7) if n == 1 else 0))], G, f)
ch = chart_series(N, F(-1), F(1), G, f)
R.check('S4 control: singular b_1 + 1/7 breaks step 1 (1-variable != chart) and gives a pole',
        not ls_equal(one[1], ch[1])[0] and not regular(one), 'poles %s' % poles(one))
R.done()
