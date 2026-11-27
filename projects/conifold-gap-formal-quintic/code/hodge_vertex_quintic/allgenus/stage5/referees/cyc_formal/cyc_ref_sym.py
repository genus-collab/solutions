# -*- coding: utf-8 -*-
"""cyc_ref_sym.py K -- referee cyc_formal: (CYC) for the GENERIC conifold frame, alpha and beta indeterminates.

y_0 = f_0 + alpha g + beta f_2 (fq_conifold's parametrisation), t = g/y_0, coefficients in Q[alpha, beta] (cyc_ref_core.P2).
  Y1  B(y_0, g) = 0 identically in (alpha, beta) (bilinearity and C4 on the basis).
  Y2  M := (1/y_0) D (1/Dt) D (u y_0^2 Dt^2) D (1/Dt) D (1/y_0), expanded as an operator, equals PF4 coefficientwise below
      u^(K-4): every u-coefficient of M - PF4 is the zero POLYNOMIAL in alpha, beta.
  Y3  the operator form of the Corollary, O o (1/y_0) = (1/y_0) o D o PF4, identically in (alpha, beta), below u^(K-5).
  Control: y_1 = f_1 (gamma = 0) instead of g: some u-coefficient of M - PF4 is a non-zero polynomial (the obstruction is
  visible already at alpha = beta = 0).
Usage: py -I -B cyc_ref_sym.py K    (stored as logs/cyc_ref_sym_K.log)
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U, P2
from cyc_ref_ops import Op, REl
from cyc_ref_frob import PF4, basis, bmat, pairing

K = int(sys.argv[1]) if len(sys.argv) > 1 else 10
NF = [0]
DOP = Op.Dk(1)


def check(name, ok, control=False):
    if control:
        print(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        NF[0] += bool(ok)
    else:
        print(('PASS ' if ok else 'FAIL ') + name)
        NF[0] += (not ok)


def Mop(y0, t):
    e = t.D(); ie = e.inv(); iy = y0.inv()
    mid = U * y0 * y0 * e * e
    return Op.fn(iy) * (DOP * (Op.fn(ie) * (DOP * (Op.fn(mid) * (DOP * (Op.fn(ie) * (DOP * Op.fn(iy))))))))


def prec_of(op):
    return min(s.N for s in op.A.values())


t0 = time.time()
bas = basis(K + 8)
B = bmat(bas['Y'])[0]
f = [x.trunc(K) for x in bas['f']]; g = bas['g'].trunc(K)
al, be = P2({(1, 0): 1}), P2({(0, 1): 1})
v0 = [Fr(1), al, al * bas['gamma'] + be, Fr(0)]                # y_0 = f_0 + alpha (f_1 + gamma f_2) + beta f_2
Bsym = sum((v0[i] * B[i][j] * [0, 1, bas['gamma'], 0][j] for i in range(4) for j in range(4)), P2())
check('Y1 B(y_0, g) = 0 identically in (alpha, beta)', Bsym == 0)
y0 = f[0] + g.sc(al) + f[2].sc(be)
t = g * y0.inv()
dt = max((max(i + j for (i, j) in x.d) for x in t.c if isinstance(x, P2) and x.d), default=0)
check('Y0 the symbols are live: t = u + (7/10 - alpha) u^2 + ..., coefficients of t of total degree up to %d' % dt,
      t.coef(2) == P2({(0, 0): Fr(7, 10), (1, 0): -1}) and dt >= K - 3)
M = Mop(y0, t)
p = prec_of(M)
nz = sum(1 for k, s in M.A.items() for x in s.c if not x == 0)
deg = max((max(i + j for (i, j) in x.d) for s in M.A.values() for x in s.c if isinstance(x, P2) and x.d), default=0)
print('   M has %d non-zero coefficient entries below u^%d (top total degree in alpha, beta seen: %d)' % (nz, p, deg))
check('Y2 M = PF4 coefficientwise below u^%d, identically in (alpha, beta)' % p, M.diff(PF4, p))
ie = t.D().inv(); iy0 = y0.inv()
ic = [ie, U * y0 * y0 * t.D() * t.D(), ie, iy0, iy0]
O = Op.fn(ic[0]) * DOP
for a in (1, 2, 3, 4):
    O = Op.fn(ic[a]) * (DOP * O)
O = O * Op.fn(iy0); R = Op.fn(iy0) * (DOP * PF4)
po = min(prec_of(O), prec_of(R))
check('Y3 O o (1/y_0) = (1/y_0) o D o PF4 below u^%d, identically in (alpha, beta)' % po, O.diff(R, po))
Mc = Mop(y0, f[1] * y0.inv())
check('Y2 control: y_1 = f_1 (non-isotropic): M = PF4 identically', Mc.diff(PF4, prec_of(Mc)), control=True)
print('CYC_REF_SYM %d: %s   (%.1f s)' % (K, 'ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
