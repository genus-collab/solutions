# -*- coding: utf-8 -*-
"""xy_explicit.py -- the explicit low (g, n) formulas (code convention, = (-1)^n ABDKS = EO07 kernel with y -> -y):
   omega_{0,3} = Sum_i D_{x_i}[ dz_1 dz_2 dz_3 / prod_{j != i} (z_i - z_j)^2 ]
   omega_{1,1} = -D_x^2 (h dz),                         h = z^2/2 + v/(12 z^3)
   omega_{1,2} = -D_1^3[h(z_1) B] - D_2^3[h(z_2) B] + (1/4)(D_1^2 + 2 D_1 D_2 + D_2^2)[B^(2)],  B = dz1 dz2/(z1-z2)^2, B^(2) = dz1 dz2/(z1-z2)^4
   omega_{2,1} = -D_x^4 (k dz) + (1/2) D_x^5 (h^2 dz),   k = 1/80 + v/(80 z^5)
with D_x omega = d(omega/dx), dx = (4 z^3 - v/z^2) dz.  Each is checked against own direct TR (exact, Q(zeta_5)) at rational
points; omega_{1,1}, omega_{2,1} are printed as N(z, v)/(4 z^5 - v)^K.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from xy_field import QZ
from xy_tr import TR
from xy_closed import Dx, common, pmul, ev_rat
from xy_jets import jmul, uni, edge_jet, apply_D, taylor_invxp, taylor_laurent

V = 4


def two_point_12(b, v):
    """the omega_{1,2} formula above at z = b (jets)"""
    n = 2; tot = Fr(0)
    h = {2: Fr(1, 2), -3: Fr(v, 12)}
    for (i, r, m, cf) in [(0, 3, 2, -1), (1, 3, 2, -1)]:
        R = [r if q == i else 0 for q in range(n)]
        jet = jmul(uni(i, n, taylor_laurent(h, b[i], R[i])), edge_jet(0, 1, n, m, b, R), R)
        jet, _ = apply_D(jet, i, r, taylor_invxp(b[i], v, r), R)
        tot += cf * jet.get((0, 0), 0)
    for (r0, r1, cf) in [(2, 0, Fr(1, 4)), (1, 1, Fr(1, 2)), (0, 2, Fr(1, 4))]:
        R = [r0, r1]; jet = edge_jet(0, 1, n, 4, b, R); Rc = list(R)
        jet, Rc = apply_D(jet, 0, r0, taylor_invxp(b[0], v, r0), Rc)
        jet, Rc = apply_D(jet, 1, r1, taylor_invxp(b[1], v, r1), Rc)
        tot += cf * jet.get((0, 0), 0)
    return tot


def three_point_03(b, v):
    n = 3; tot = Fr(0)
    for i in range(3):
        R = [1 if q == i else 0 for q in range(n)]; jet = {(0, 0, 0): Fr(1)}
        for j in range(3):
            if j != i: jet = jmul(jet, edge_jet(i, j, n, 2, b, R), R)
        jet, _ = apply_D(jet, i, 1, taylor_invxp(b[i], v, 1), R)
        tot += jet.get((0, 0, 0), 0)
    return tot


def main():
    tr = TR(QZ, 2, 4)
    for (g, n) in [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1)]: tr.step(g, n)
    h = {(2, 0): Fr(1, 2), (-3, 1): Fr(1, 12)}
    N, K = h, 0
    for _ in range(2): N, K = Dx(N, K)
    w11 = ({kk: -c for kk, c in N.items()}, K)
    k5 = {(0, 0): Fr(1, 80), (-5, 1): Fr(1, 80)}
    A, KA = k5, 0
    for _ in range(4): A, KA = Dx(A, KA)
    Bq, KB = pmul(h, h), 0
    for _ in range(5): Bq, KB = Dx(Bq, KB)
    w21 = common([({kk: -c for kk, c in A.items()}, KA), ({kk: c / 2 for kk, c in Bq.items()}, KB)])
    for name, (N, K), g in [("omega_{1,1}", w11, 1), ("omega_{2,1}", w21, 2)]:
        ok = all((tr.value(g, 1, (z,)) - QZ.frac(ev_rat(N, K, z, V))).is_zero() for z in (Fr(2), Fr(-3, 7), Fr(5, 3)))
        print("CHECK %s = explicit D_x formula (TR, exact, 3 points): %s" % (name, ok))
        terms = sorted(N.items(), key=lambda t: (-t[0][1], t[0][0]))
        print("   %s = [%s] / (4 z^5 - v)^%d dz" % (name, " + ".join("(%s) z^%d v^%d" % (c, ez, ev) for (ez, ev), c in terms), K))
    for b in [(Fr(2), Fr(-1, 3)), (Fr(-5, 2), Fr(7, 4))]:
        print("CHECK omega_{1,2} formula = TR at %s: %s" % (b, (tr.value(1, 2, b) - QZ.frac(two_point_12(b, V))).is_zero()))
    for b in [(Fr(2), Fr(-1, 3), Fr(3, 2)), (Fr(-5, 2), Fr(7, 4), Fr(1, 3))]:
        print("CHECK omega_{0,3} formula = TR at %s: %s" % (b, (tr.value(0, 3, b) - QZ.frac(three_point_03(b, V))).is_zero()))


if __name__ == "__main__":
    main()
