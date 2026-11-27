# -*- coding: utf-8 -*-
"""lm_run_q.py N FIELD -- referee lmi_adversary: the formal quintic through Theorem LMI (s.3.4 of RA), with controls.
FIELD = Q or p1 / p2.  Every series c_n(s* + e) is computed by the plain Wick sum at L = s* + e with Laurent-series
coefficients (lm_core); the regular data are Taylor series at s*, truncated with tracked precision.
Q1 colour 0 at L = 1: no pole, value h.   Q2 Gamma(5s) for Gamma(1 + 5s): value still h (blind), c_n changes off L = 1.
Q3 singular factor Gamma(1 + (1 - sigma)/z) (measure dx_0): no pole, value = LMI with delta = z, NOT h.
Q4 q = oo (L = 0): Gamma(1 + 5 sigma/z) singular (alpha0 = 5, M0 = 1): value Ginf (LMI formula AND Bernoulli closed form).
Q5 q = oo with Gamma(5s) (M0 = 0): value exp(5 Sum_(5 | 2k-1) b_k z^(2k-1)), NOT Ginf.
Q6 (F_p only) colours i = 1, 2 with an explicit zeta: value h(zeta^-i z)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from lm_core import F, fr, fpow, setmod
from lm_wick import bk
from lm_data import Sing, cseries, lmi_value, zexp
from lm_models import quintic_s1, quintic_s0, quintic_col, h_quintic, ginf
from lm_harness import Run, ls_equal, regular, poles, values, prime_below, root_of_unity

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
FIELD = sys.argv[2] if len(sys.argv) > 2 else 'Q'
PR = {'Q': None, 'p1': prime_below(2 ** 62), 'p2': prime_below(2 ** 70)}[FIELD]
setmod(PR)
R = Run('lm_run_q_%d_%s' % (N, FIELD))
R.say('field %s%s; z-order N = %d' % (FIELD, '' if PR is None else ' (p = %d)' % PR, N))
K = 3 * N + 12
h = h_quintic(N)

G, f = quintic_s1(K)
c1 = cseries(N, [Sing(F(-1), 0)], G, f, exact=False)
R.check('Q1 quintic colour 0: no pole at L = 1 and c_n(1) = h_n, n <= %d' % N, regular(c1) and values(c1) == h,
        'min prec %d' % min(x.prec for x in c1))
R.check('Q1 LMI formula F(1,z)/F(1,0) for the five regular factors = h', lmi_value(N, F(-1), F(0), G, f) == h)
f5 = {m: list(r) for m, r in f.items()}
for i in range(1, len(f5[0])):
    f5[0][i] = fr(f5[0][i] - F(Fr((-1) ** (i - 1), i)))         # (5 sigma)^(+1/2) -> (5 sigma)^(-1/2)
c5 = cseries(N, [Sing(F(-1), 0)], G, f5, exact=False)
R.check('Q2 Gamma(5s) for Gamma(1 + 5s): still no pole and value h at L = 1 (the twist measure is invisible there)',
        regular(c5) and values(c5) == h)
dn = [n for n in range(N + 1) if not ls_equal(c1[n], c5[n])[0]]
R.check('Q2 ... but c_n(L) itself changes for L != 1 (orders n where the two series differ)', len(dn) > 0, 'n in %s' % dn)
cm = cseries(N, [Sing(F(-1), 1)], G, f, exact=False)
vm = lmi_value(N, F(-1), F(1), G, f)
R.check('Q3 singular factor Gamma(1 + (1 - sigma)/z): no pole, value = LMI with delta = z, and that value != h',
        regular(cm) and values(cm) == vm and vm != h, 'first differing order %s' %
        next((n for n in range(N + 1) if vm[n] != h[n]), None))

G0, f0 = quintic_s0(K)
c0 = cseries(N, [Sing(F(5), 1)], G0, f0, exact=False)
gi = ginf(N)
R.check('Q4 q = oo: no pole at L = 0 and c_n(0) = [z^n] Ginf (Bernoulli closed form), n <= %d' % N,
        regular(c0) and values(c0) == gi)
R.check('Q4 LMI formula at s* = 0 (delta = -z/5) = Ginf', lmi_value(N, F(5), F(1), G0, f0) == gi)
cz = cseries(N, [Sing(F(5), 0)], G0, f0, exact=False)
hz = [F(0)] * (N + 1)
for k in range(1, N // 2 + 2):
    if 2 * k - 1 <= N and (2 * k - 1) % 5 == 0:
        hz[2 * k - 1] = F(5 * bk(k))
vz = zexp(hz, N)
R.check('Q5 q = oo with Gamma(5s) (measure dy/y): no pole, value exp(5 Sum_(5|2k-1) b_k z^(2k-1)) = LMI(M0 = 0)',
        regular(cz) and values(cz) == vz == lmi_value(N, F(5), F(0), G0, f0))
R.check('Q5 ... and this value differs from Ginf (the shift 1/5 comes from the measure dy)', N < 5 or vz != gi)

if PR is not None:
    zeta = root_of_unity(PR, 5)
    for i0 in (1, 2):
        Gc, fc = quintic_col(K, i0, zeta)
        cc = cseries(N, [Sing(F(-1), 0)], Gc, fc, exact=False)
        hi = [fr(h[n] * fpow(zeta, -i0 * n)) for n in range(N + 1)]
        R.check('Q6 colour %d at L = zeta^%d (explicit zeta mod p): no pole, value h(zeta^-%d z)' % (i0, i0, i0),
                regular(cc) and values(cc) == hi)
        R.check('Q6 control: colour %d value != h(z) (the zeta-twist is seen)' % i0, N < 1 or values(cc) != h)
R.done()
