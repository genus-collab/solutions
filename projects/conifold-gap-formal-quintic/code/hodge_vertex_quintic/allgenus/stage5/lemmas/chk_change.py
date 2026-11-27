# -*- coding: utf-8 -*-
"""chk_change.py -- Lemma L1a (formal change of variables) and Lemma FL1 (exact forms), exact over Q.

L1a: SP^norm_x[f; a] = SP^norm_y[f o phi; (a o phi) det Dphi / det J] for a random polynomial map
     phi(y) = J y + (quadratic) + (cubic), J = Dphi(0) invertible, random phase f and amplitude a(x; z).
FL1: SP^norm[f; sum_i (d_i h_i + z^-1 (d_i f) h_i)] = 0 for random h_i(x; z).
Controls (must fail): L1a without det Dphi; L1a with f o phi replaced by f o J (phase not pulled back);
FL1 without the z^-1 (d_i f) h_i term.
Usage: py -I -B chk_change.py [K] [NINST]     (through z^K, default 8; NINST instances per dimension, default 2)
"""
import os, sys, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from sp_core import Ctx, rand_poly, zseries
from sp_iter import total_sp
from sp_rand import rand_phase, rand_amp, rand_sym, det_q
from chk_util import Log, first_diff

KTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NI = int(sys.argv[2]) if len(sys.argv) > 2 else 2
log = Log(HERE, 'chk_change_%d_%d' % (KTOP, NI))
log.out('chk_change: Lemmas L1a and FL1 through z^%d, %d instances per dimension; log %s' % (KTOP, NI, log.rel))


def rand_map(cb, rng, xs, ys):
    """phi_i(y) = (J y)_i + random quadratic + cubic terms; returns (subs dict, J)"""
    n = len(xs)
    while True:
        J = [[rng.randint(-2, 2) for _ in range(n)] for _ in range(n)]
        if det_q(J) != 0:
            break
    subs = {}
    for a, i in enumerate(xs):
        s = {}
        for b, v in enumerate(ys):
            s = cb.add(s, cb.var(v, J[a][b]))
        subs[i] = cb.add(s, rand_poly(cb, rng, ys, 2, 3, 2, 0.6))
    return subs, J


for n, K in [(1, KTOP), (2, KTOP), (3, min(KTOP, 6))]:
    xs, ys = list(range(n)), list(range(n, 2 * n))
    c = Ctx(2 * n, [1] * (2 * n), 2 * K); cb = Ctx(2 * n, [1] * (2 * n), 2 * K + 2)
    for inst in range(NI):
        rng = random.Random(104729 * n + inst); t = time.time()
        H, f = rand_phase(cb, rng, xs, 5 if n < 3 else 4); a = rand_amp(c, rng, xs)
        lhs = zseries(c, total_sp(c, cb, f, a, xs))
        subs, J = rand_map(cb, rng, xs, ys)
        fphi = cb.compose(f, subs); aphi = cb.compose(a, subs)
        Dphi = [[cb.deriv(subs[i], v) for v in ys] for i in xs]
        jac = c.scal(c.trunc(cb.det(Dphi)), 1 / det_q(J))
        rhs = zseries(c, total_sp(c, cb, fphi, c.mul(c.trunc(aphi), jac), ys))
        fd = first_diff(lhs, rhs, K)
        log.check('L1a n = %d inst %d: SP_x[f; a] = SP_y[f o phi; (a o phi) det Dphi / det J] through z^%d (J = %s)'
                  % (n, inst, K, J), fd is None, 'z^1: %s  (%.1f s)' % (lhs.get(1, 0), time.time() - t))
        if inst == 0:
            r1 = zseries(c, total_sp(c, cb, fphi, c.trunc(aphi), ys))
            log.control('L1a n = %d: without det Dphi' % n, first_diff(lhs, r1, K) is not None,
                        'first differs at z^%s' % first_diff(lhs, r1, K))
            lin = {i: {k: v for k, v in subs[i].items() if sum(k[0]) == 1} for i in xs}
            r2 = zseries(c, total_sp(c, cb, cb.compose(f, lin), c.mul(c.trunc(aphi), jac), ys))
            log.control('L1a n = %d: phase pulled back by J only' % n, first_diff(lhs, r2, K) is not None,
                        'first differs at z^%s' % first_diff(lhs, r2, K))
        # FL1: an exact form d( exp(f/z) sum_i h_i dx^(i) ) has SP zero
        hs = [c.add(rand_poly(c, rng, xs, 0, 3, 3, 0.7), c.zp(rand_poly(c, rng, xs, 0, 2, 3, 0.7), 1)) for _ in xs]
        ex, ex_bad = {}, {}
        for i, h in zip(xs, hs):
            term = c.add(c.deriv(h, i), c.zp(c.mul(c.deriv(c.trunc(f), i), h), -1))
            ex = c.add(ex, term); ex_bad = c.add(ex_bad, c.deriv(h, i))
        v = zseries(c, total_sp(c, cb, f, ex, xs))
        log.check('FL1 n = %d inst %d: SP of an exact form vanishes through z^%d' % (n, inst, K),
                  all(v.get(k, 0) == 0 for k in range(-2, K + 1)), 'nonzero orders: %s' % sorted(v)[:3])
        if inst == 0:
            vb = zseries(c, total_sp(c, cb, f, ex_bad, xs))
            log.control('FL1 n = %d: without the z^-1 (d_i f) h_i term' % n, any(vb.get(k, 0) for k in range(-2, K + 1)),
                        'first nonzero order z^%s' % (min(vb) if vb else None))
log.finish('CHK_CHANGE %d %d' % (KTOP, NI))
