# -*- coding: utf-8 -*-
"""cy_proj.py NU -- referee cyc_adversary, task item (a): are CYC's objects the programme's objects?
Imports the project modules READ-ONLY (fq_conifold, fq_ring, fq_rows, fq_fp; nothing is written by them) and compares:
 (P1) fq_conifold.Conifold (exact over Q internally): f0, f1, f2, gamma, g, h, varpi0, t_c, the pairing table, in the frames
      (alpha, beta) = (0, 0) and (1/3, -2/7), with this referee's independent Frobenius data (cy_core);
 (P2) the ACTUAL row coefficients of fq_rows.Rows (kinv[a] = 1/(C1 kappa_a), a = 0..3, elements of the ring N), pushed through the
      project's conifold map Cf (Conifold.cf_images, two values of Lam), against 1/c_a, c = (Dt, L^5/(varpi0^2 Dt^2), Dt, varpi0) of
      CYC s.0, as s-series mod p = 2^89-1 (u = s^5); also C1 = W/c0^2 -> Dt and c0 -> varpi0;
 (P3) UR s.1's b-vector against b_a = -c_a/(1-u) from CYC's c (exact over Q).   Usage: py -I -B cy_proj.py 20"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
HVQ = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
sys.path.insert(0, HERE); sys.path.insert(1, HVQ)
from fractions import Fraction as Fr
from cy_core import Ser, ONE_M
from cy_frob import pf4_poly, frobenius
from cy_log import Log
import fq_conifold, fq_ring, fq_rows
from fq_fp import FpRing

NU = int(sys.argv[1]) if len(sys.argv) > 1 else 20
LG = Log('cy_proj_%d' % NU)
LG('cy_proj: project modules imported read-only from the package folder (relative: ../../../../); NU = %d' % NU)
N = NU + 8
FR = frobenius(pf4_poly(), N)
mine_f = [Ser(0, x, N) for x in FR['f']]; mine_g = Ser(0, FR['g'], N); mine_h = Ser(0, FR['h'], N)
R = FpRing(); p = R.p


def same_q(proj, mine, upto):
    """project Ser over Q vs own Ser, coefficients u^0..u^(upto-1)"""
    return all(Fr(proj.coef(k)) == mine.coef(k) for k in range(upto))


def modp(x): x = Fr(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def to_s_modp(f, ns):
    """own u-series -> {s-exponent: value mod p}, exponents < ns"""
    return {5 * (f.v + i): modp(x) for i, x in enumerate(f.c) if x and 5 * (f.v + i) < ns}


def same_s(projser, own, lo, hi):
    lo = min([lo, projser.v] + [k for k in own if k < hi])
    return all(projser.coef(k) % p == own.get(k, 0) for k in range(lo, hi))


for al, be in ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7))):
    C = fq_conifold.Conifold(R, Nu=NU, alpha=al, beta=be)
    up = NU
    okf = all(same_q(C.f[k], mine_f[k], up) for k in range(3))
    okg = C.gamma == FR['gamma'] and same_q(C.g, mine_g, up) and same_q(C.h, mine_h, up)
    y0 = mine_f[0].add(mine_g.scale(al)).add(mine_f[2].scale(be)); tc = mine_g.mul(y0.inv())
    okv = same_q(C.varpi0, y0, up) and same_q(C.tc, tc, up)
    LG.check('(P1) frame (%s, %s): fq_conifold f0, f1, f2 = own closed-form Frobenius basis through u^%d' % (al, be, up - 1), okf)
    LG.check('(P1) frame (%s, %s): gamma = %s, g, h (h_0 = h_1 = h_2 = 0) agree' % (al, be, C.gamma), okg)
    LG.check('(P1) frame (%s, %s): varpi0 = f0 + alpha g + beta f2 and t_c = g/varpi0 agree (programme convention)' % (al, be), okv)
    LG.check('(P1) frame (%s, %s): pairing table in the basis (f0, g, f2) = [[0,0,2],[0,0,0],[-2,0,0]], constant, pairing_vacuous'
             % (al, be), C.Bconst_ok and C.pairing_vacuous and [[Fr(x) for x in r] for r in C.Bmat] == [[0, 0, 2], [0, 0, 0], [-2, 0, 0]],
             info='project Bmat = %s' % [[str(x) for x in r] for r in C.Bmat])
    rows = fq_rows.Rows(R, Nz=1, reduce=False, colours=[0])
    ring = rows.ring
    Dt = tc.D(); c = [Dt, Ser.poly({1: 1}).mul(y0).mul(y0).mul(Dt).mul(Dt).inv(), Dt, y0]
    for Lam in (0, 3):
        img = C.cf_images(Lam)
        ns = min(v.prec for v in img.values()) - 12
        res = []
        for a in range(4):
            ps = rows.kinv[a].subst(img, ring); own = to_s_modp(c[a].inv(), ns)
            res.append(same_s(ps, own, ps.v, min(ps.prec, ns)))
        c1p = ring.C1().subst(img, ring); c0p = ring.gen('c0').subst(img, ring)
        okc1 = same_s(c1p, to_s_modp(Dt, ns), c1p.v, min(c1p.prec, ns)); okc0 = same_s(c0p, to_s_modp(y0, ns), c0p.v, min(c0p.prec, ns))
        LG.check('(P2) frame (%s, %s), Lam = %d: Cf(kinv[a]) = 1/c_a for a = 0, 1, 2, 3 (fq_rows order = CYC order)' % (al, be, Lam),
                 all(res), info='per a: %s; s-orders < %d, mod 2^89-1' % (res, ns))
        LG.check('(P2) frame (%s, %s), Lam = %d: Cf(C1) = Cf(W)/Cf(c0)^2 = D t_c and Cf(c0) = varpi0' % (al, be, Lam), okc1 and okc0)
    wrong = to_s_modp(c[1].scale(5).inv(), ns)
    LG.check('(P2) control, frame (%s, %s): Cf(kinv[1]) = 1/(5 c_1) (a Yukawa factor 5 in the row coefficient)' % (al, be),
             same_s(rows.kinv[1].subst(img, ring), wrong, -20, ns), control=True)   # img: the Lam = 3 images
    # (P3) UR s.1: b = (t_c', -1/(u varpi0^2 (1-u)^3 t_c'^2), t_c', -varpi0/(1-u), -varpi0/(1-u)), ' = d/du
    tp = tc.deriv(); im = ONE_M.inv(cap=N)
    bUR = [tp, Ser.poly({1: 1}).mul(y0).mul(y0).mul(ONE_M).mul(ONE_M).mul(ONE_M).mul(tp).mul(tp).inv().neg(), tp,
           y0.mul(im).neg(), y0.mul(im).neg()]
    cfull = c + [y0]
    LG.check('(P3) frame (%s, %s): UR s.1 b_a = -c_a/(1-u) for all five a (incl. the closing coefficient c_4 = varpi0)' % (al, be),
             all(bUR[a].sub(cfull[a].mul(im).neg()).is_zero() for a in range(5)))
LG.close()
