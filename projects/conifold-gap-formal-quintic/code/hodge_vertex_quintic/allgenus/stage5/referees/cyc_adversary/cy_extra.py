# -*- coding: utf-8 -*-
"""cy_extra.py N -- referee cyc_adversary, by-product: is Psi := varpi0 (Y^C_Q)_0 (the quintic vanishing solution of GI Lemma Q,
built from the cyclic system alone, times varpi0) independent of the conifold frame?  PB s.2.7 states V = C0 (Y^C_Q)_0 "in every
frame" with V the frame-free residue of the MB summand, which predicts this.  Exact over Q, z^-n coefficients n <= 16, u^(N-1).
Usage: py -I -B cy_extra.py 40"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cy_core import Ser, ONE_M
from cy_frob import pf4_poly, frobenius
from cy_log import Log

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LG = Log('cy_extra_%d' % N)
LG('cy_extra: exact over Q, u^%d, z^-n for n <= 16' % (N - 1))
FR = frobenius(pf4_poly(), N)
f0, f1, f2 = [Ser(0, x, N) for x in FR['f']]; g = Ser(0, FR['g'], N); U = Ser.poly({1: 1})


def lin(*ps):
    acc = Ser(0, [])
    for x, s in ps: acc = acc.add(s.scale(x))
    return acc


def psi_vanishing(y0, y1, nmax=16):
    Dt = y1.mul(y0.inv()).D(); c1 = U.mul(y0).mul(y0).mul(Dt).mul(Dt).inv()
    b = [x.mul(ONE_M.inv(cap=N)).neg() for x in (Dt, c1, Dt, y0, y0)]
    y = [[Ser(0, []), Ser(0, [-1]), Ser(0, []), Ser(0, []), Ser(0, [])]]
    for n in range(nmax): y.append([b[a].mul(y[-1][(a + 1) % 5]).integ0() for a in range(5)])
    return [y0.mul(yn[0]) for yn in y]


ref = psi_vanishing(f0, g)
same = lambda P, Q: all(P[n].sub(Q[n]).is_zero() for n in range(len(P)))
for lab, y0, y1 in (('programme (1/3,-2/7)', lin((1, f0), (Fr(1, 3), g), (Fr(-2, 7), f2)), g),
                    ('programme (7/11,-3/5)', lin((1, f0), (Fr(7, 11), g), (Fr(-3, 5), f2)), g),
                    ('programme (-2,0)', lin((1, f0), (-2, g)), g),
                    ('critic-convention (1/3,-2/7) = f0 + f1/3 - 2f2/7', lin((1, f0), (Fr(1, 3), f1), (Fr(-2, 7), f2)), g)):
    LG.check('(X1) varpi0 (Y^C_Q)_0 in frame %s equals the frame-(0,0) one, n <= 16' % lab, same(psi_vanishing(y0, y1), ref))
P2 = psi_vanishing(lin((2, f0), (Fr(2, 3), g), (Fr(-4, 7), f2)), g.scale(3))
LG.check('(X2) rescaled frame (2 varpi0, 3 g): Psi = 3 Psi(0,0) (the vanishing solution carries the scale of g, not of varpi0)',
         all(P2[n].sub(ref[n].scale(3)).is_zero() for n in range(17)))
LG.check('(X3) control: wrong gamma (y1 = f1 + f2/2) gives the same Psi', same(psi_vanishing(f0, lin((1, f1), (Fr(1, 2), f2))), ref),
         control=True)
LG('    Psi_1 = %s ; Psi_6 = %s' % (ref[1].head(3), ref[6].head(3)))
LG.close()
