# -*- coding: utf-8 -*-
"""cyc_ref_xproj.py -- referee cyc_formal: is CYC's PF4 the project's?  Cross-check against the project's own code.

Imports SCRATCH COPIES xproj/fq_conifold.py, xproj/fq_fp.py, xproj/fq_series.py (byte-identical copies of the project files
fq_conifold.py, fq_fp.py, fq_series.py in the package folder; nothing in the project is run or written) and compares
Conifold(QRing, Nu=20) (exact over Q) with this referee's own Frobenius basis (cyc_ref_frob):
  X1  the analytic solutions f_0, f_1, f_2 (all coefficients below u^28), gamma, the log companion h (u^3..u^27);
  X2  fq_conifold's pairing matrix in its basis (f_0, g, f_2) against the referee's B, and its flags Bconst_ok, pairing_vacuous;
  X3  varpi0 and tc in the frames (alpha, beta) = (0, 0), (1/3, -2/7) against cyc_ref_frob.frame(.., 'project').
Usage: py -I -B cyc_ref_xproj.py     (stored as logs/cyc_ref_xproj.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'xproj'))
from fractions import Fraction as Fr
from fq_fp import QRing
from fq_conifold import Conifold
from cyc_ref_core import Ser
from cyc_ref_ops import REl
from cyc_ref_frob import basis, bmat, frame

NF = [0]


def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    NF[0] += (not ok)


t0 = time.time()
C = Conifold(QRing(), Nu=20)
bas = basis(40)
M = 28
check('X1 project f_0, f_1, f_2 = referee f_0, f_1, f_2 below u^%d' % M,
      all(C.f[i].coef(k) == bas['f'][i].coef(k) for i in range(3) for k in range(M)))
check('X1 project gamma = referee gamma = 7/10', C.gamma == bas['gamma'] == Fr(7, 10))
check('X1 project log companion h = referee h (u^0..u^%d), log obstruction zero' % (M - 1),
      all(C.h.coef(k) == bas['h'].coef(k) for k in range(M)) and all(x == 0 for x in C.log_obstruction))
gl = [REl.of(bas['f'][0]), REl.of(bas['g']), REl.of(bas['f'][2])]
B3 = bmat(gl)[0]
print('   fq_conifold Bmat (basis f_0, g, f_2) = %s' % [[str(x) for x in r] for r in C.Bmat])
check('X2 fq_conifold Bmat (basis f_0, g, f_2) = referee B in that basis; Bconst_ok and pairing_vacuous are True',
      [[Fr(x) for x in r] for r in C.Bmat] == B3 and C.Bconst_ok and C.pairing_vacuous)
for (a, b) in ((0, 0), (Fr(1, 3), Fr(-2, 7))):
    Cab = Conifold(QRing(), Nu=20, alpha=a, beta=b)
    y0, y1, t = frame(bas, Fr(a), Fr(b))
    check('X3 frame (%s, %s): project varpi0 and tc = referee y_0 and t = g/y_0 below u^%d' % (a, b, M),
          all(Cab.varpi0.coef(k) == y0.coef(k) and Cab.tc.coef(k) == t.coef(k) for k in range(M)))
print('CYC_REF_XPROJ: %s   (%.1f s)' % ('ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
