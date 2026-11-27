# -*- coding: utf-8 -*-
"""ch_zpow.py -- Stage 5 / chart, step (a) (CHART_PROOF s.3.1): the powers of z cancel IDENTICALLY in sigma.

Each Gamma factor Gamma(w), w = c + (alpha + beta sigma)/z, becomes an Euler integral after t = x/z, which produces z^(-w).
The exponents are elements of Q(zeta_5)[sigma, 1/z], stored as {(sigma-degree, 1/z-degree): coefficient}; the identities below are
identities in that ring (symbolic in sigma), not evaluations. Checked:
  (i)   Sum over the six factors of -w = -1 (so Omega_7 carries z^-1, independent of sigma);
  (ii)  order (a): the six inner Stirling factors z^(a/z + m) (y: a = 5 sigma, m = 1; x_j: a = zeta^j - sigma, m = 0) give z^(+1);
  (iii) the log z terms of the six Stirling forms: coefficient of log z equals 2 (the elementary factor z^2 of GAMMA s.2);
  (iv)  q = 0 chart (a = 5, m = 1; a = zeta^j - 1, m = 0): z^(+1); step (e), the q^d term: (q/z)^d times y^(5d) and x_j^(-d): z^(+1).
Controls (must FAIL (i)): the j = 0 factor dropped; Gamma(1 + 4s) in place of Gamma(1 + 5s). Gamma(5s) in place of Gamma(1 + 5s)
(measure dy/y) keeps (i) sigma-free but changes the constant power (z^0): printed, not a failure.
Usage: py -I -B ch_zpow.py   (log logs/ch_zpow.log)
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from ch_core import QZ

LOGREL = 'logs/ch_zpow.log'
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []
F = QZ()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def padd(A, B, s=1):
    R = dict(A)
    for k, v in B.items():
        R[k] = F.add(R.get(k, F.zero), F.smul(s, v))
    return {k: v for k, v in R.items() if not F.iszero(v)}


def w_of(c, alpha, beta):
    """w = c + (alpha + beta sigma)/z as {(i, j): coeff} (sigma^i z^-j)"""
    return padd(padd({(0, 0): F.of(c)}, {(0, 1): alpha}), {(1, 1): F.of(beta)})


def show(P):
    return ' + '.join('(%s) sigma^%d z^-%d' % ([str(x) for x in v], i, j) for (i, j), v in sorted(P.items())) or '0'


zj = [F.zeta(j) for j in range(5)]
FACT = [w_of(1, F.zero, 5)] + [w_of(0, zj[j], -1) for j in range(5)]      # Gamma(1 + 5 sigma/z), Gamma((zeta^j - sigma)/z)
E = {}
for w in FACT: E = padd(E, w, -1)
check('(i) Sum of the Euler-substitution exponents -w over the six factors = -1, identically in sigma', E == {(0, 0): F.of(-1)}, show(E))
# (ii) inner Stirling factors z^(a/z + m)
A2 = padd(padd({(0, 0): F.of(1)}, {(1, 1): F.of(5)}), {})
for j in range(5): A2 = padd(A2, padd({(0, 1): zj[j]}, {(1, 1): F.of(-1)}))
check('(ii) order (a): Sum (a/z + m) over the six inner Euler integrals = +1, identically in sigma (cancels z^-1)',
      A2 == {(0, 0): F.of(1)}, show(A2))
# (iii) log z coefficient of log Gamma_St(w), w = num/z: -(w - 1/2); Gamma_St(1 + 5 sigma/z) = (5 sigma/z) Gamma_St(5 sigma/z)
LZ = padd({(0, 0): F.of(-1)}, padd({(0, 0): F.of(Fr(1, 2))}, {(1, 1): F.of(-5)}))
for j in range(5): LZ = padd(LZ, padd({(0, 0): F.of(Fr(1, 2))}, padd({(0, 1): F.neg(zj[j])}, {(1, 1): F.of(1)})))
check('(iii) coefficient of log z in the six Stirling forms = 2, identically in sigma (elementary factor z^2)',
      LZ == {(0, 0): F.of(2)}, show(LZ))
# (iv) q = 0 chart and the q^d term of step (e)
c0 = padd({(0, 0): F.of(1)}, {(0, 1): F.of(5)})
for j in range(1, 5): c0 = padd(c0, {(0, 1): F.sub(zj[j], F.one)})
check('(iv) q = 0 chart: Sum (a/z + m) = +1 (a = 5, m = 1; a = zeta^j - 1, m = 0)', c0 == {(0, 0): F.of(1)}, show(c0))
ok = True
for d in range(0, 9):
    cd = padd({(0, 0): F.of(1 + 5 * d - d)}, {(0, 1): F.of(5)})          # (q/z)^d: -d; y: m = 1 + 5d
    for j in range(1, 5): cd = padd(cd, padd({(0, 1): F.sub(zj[j], F.one)}, {(0, 0): F.of(-d)}))
    ok = ok and cd == {(0, 0): F.of(1)}
check('(iv) step (e): the q^d term has the same power z^(+1) as the q^0 term, d = 0..8', ok)
for name, fac in (('j = 0 factor dropped', [FACT[0]] + FACT[2:]), ('Gamma(1 + 4s) in place of Gamma(1 + 5s)', [w_of(1, F.zero, 4)] + FACT[1:])):
    Ec = {}
    for w in fac: Ec = padd(Ec, w, -1)
    check('control: %s gives a sigma-dependent power of z' % name, any(i > 0 for (i, j) in Ec), show(Ec))
Eg = {}
for w in [w_of(0, F.zero, 5)] + FACT[1:]: Eg = padd(Eg, w, -1)
out('note: Gamma(5s) (measure dy/y) gives the exponent %s: still sigma-free, only the constant power changes' % show(Eg))
out('CH_ZPOW: %s' % ('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS))
LOG.close()
