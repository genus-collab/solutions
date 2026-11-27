# -*- coding: utf-8 -*-
"""pb_cyc.py N -- Stage-5 'pb' audit: the hypothesis (CYC) used by Theorem P(b), exact over Q.

(CYC): the cyclic system z D Y_a = c_a Y_(a+1), c = (Dt, L^5/(varpi0^2 (Dt)^2), Dt, varpi0, varpi0), D = q d/dq, L^5 = 1/u, is the
formal quintic's Picard-Fuchs equation PF5 = z^5 D PF4 - 1 for Psi = varpi0 Y_0.  Since c_4 = varpi0 this is equivalent to the
operator identity (PB_AUDIT s.2.6)
    PF4 = (1/varpi0) D (1/Dt) D (u varpi0^2 (Dt)^2) D (1/Dt) D (1/varpi0),
PF4 = u D^4 - 2(1-u) D^3 - (7/5)(1-u) D^2 - (2/5)(1-u) D - (24/625)(1-u).  Two fourth-order operators agree iff they agree on
u^0..u^4 (non-zero Wronskian); tested on u^0..u^6 to O(u^N), frames (0,0), (1/3,-2/7), (7/11,-3/5).
Controls: t = f2/f0 (a non-isotropic pair: f2 is not the vanishing period) and the vanishing period with gamma = 1/2.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C
from pb_core import LS
from fractions import Fraction as Fr

LOG = []
def log(s=''):
    print(s); LOG.append(s)


def Dop(f):
    one_m = LS.from_list([1, -1], 0, f.prec + 2)
    return one_m.mul(f.deriv()).neg()


def rhs(f, w0, t):
    tp = Dop(t); itp = tp.inv(); iw0 = w0.inv()
    mid = w0.mul(w0).mul(tp).mul(tp).shift(1)
    g = f.mul(iw0); g = Dop(g).mul(itp); g = Dop(g).mul(mid); g = Dop(g).mul(itp); g = Dop(g).mul(iw0)
    return g


def lhs(f):
    zero = LS(f.v - 8, [Fr(0)] * (f.prec - f.v + 8), f.prec)
    return C.PF4_pair(f, zero)[0]


PRECS = []


def check(w0, t, N):
    worst = None
    for n in range(7):
        f = LS.from_list([1], n, N)
        d = lhs(f).sub(rhs(f, w0, t))
        PRECS.append(d.prec)
        d = d.norm()
        if d.c:
            k = d.v
            worst = k if worst is None else min(worst, k)
    return worst


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    C.set_field(0); t0 = time.time()
    log('pb_cyc.py %d   (exact over Q; operator identity tested on u^0..u^6 to O(u^%d))' % (N, N - 8))
    K = C.Conifold(N + 10)
    ok = True
    for fr in ((0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5))):
        K.set_frame(*fr)
        bad = check(K.varpi0.trunc(N + 4), K.tc.trunc(N + 4), N)
        log('  frame (%s,%s): PF4 = (1/varpi0) D (1/Dt) D (u varpi0^2 Dt^2) D (1/Dt) D (1/varpi0) on u^0..u^6, '
            'every coefficient below u^%d: %s' % (fr[0], fr[1], min(PRECS[-7:]), 'PASS' if bad is None else 'FAIL first at u^%d' % bad))
        ok = ok and bad is None
    f0, f1, f2 = K.f
    bad = check(f0.trunc(N + 4), f2.div(f0).trunc(N + 4), N)
    log('  control [t = f2/f0, non-isotropic pair]: %s'
        % ('first mismatch at u^%d: PASS (fails as it must)' % bad if bad is not None else 'CONTROL DID NOT BEHAVE'))
    ok = ok and bad is not None
    gw = f1.add(f2.scale(Fr(1, 2)))
    bad = check(f0.trunc(N + 4), gw.div(f0).trunc(N + 4), N)
    log('  control [vanishing period with gamma = 1/2]: %s'
        % ('first mismatch at u^%d: PASS (fails as it must)' % bad if bad is not None else 'CONTROL DID NOT BEHAVE'))
    ok = ok and bad is not None
    log('total %.1f s' % (time.time() - t0))
    log('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    with open(os.path.join(HERE, 'logs', 'pb_cyc_%d.log' % N), 'w') as f:
        f.write('\n'.join(LOG) + '\n')
    print('log written to logs/pb_cyc_%d.log' % N)


if __name__ == '__main__':
    main()
