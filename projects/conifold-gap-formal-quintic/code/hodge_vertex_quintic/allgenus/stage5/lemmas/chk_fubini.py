# -*- coding: utf-8 -*-
"""chk_fubini.py -- Lemma L1b (formal Fubini) on random polynomial phases in 2 and 3 variables, exact over Q.

Total side: SP^norm over all variables at the critical point 0 (Wick, constant covariance).
Iterated side (L1b): inner critical point xi(y) (formal IFT), inner Taylor data at xi(y), inner covariance
-K(y)^-1 with K(y) the x-Hessian at xi(y) (y-dependent Wick moments), inner series b(y; z); then the outer SP of
G(y) = F(xi(y), y) with amplitude (det K(y)/det K(0))^(-1/2) b(y; z).
Controls (must fail): no determinant factor; inner covariance frozen at y = 0; inner expansion at x = 0 with the
linear term dropped (i.e. ignoring the motion of the inner critical point).
Usage: py -I -B chk_fubini.py [K] [NINST]      (top order z^K, default 8; NINST instances per shape, default 2)
"""
import os, sys, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from sp_core import Ctx, rand_poly, zseries, inv_q
from sp_iter import total_sp, iterated_sp
from chk_util import Log, first_diff

KTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NI = int(sys.argv[2]) if len(sys.argv) > 2 else 2
log = Log(HERE, 'chk_fubini_%d_%d' % (KTOP, NI))
log.out('chk_fubini: Lemma L1b through z^%d, %d random instances per shape (nx, ny); log %s' % (KTOP, NI, log.rel))


from sp_rand import rand_phase, rand_amp

# shape (nx, ny) and top order; the 4-variable shape (2, 2) runs to z^min(K, 6)
for (nx, ny, K) in [(1, 1, KTOP), (1, 2, KTOP), (2, 1, KTOP), (2, 2, min(KTOP, 6))]:
    n = nx + ny; xs, ys = list(range(nx)), list(range(nx, n))
    c = Ctx(n, [1] * n, 2 * K); cb = Ctx(n, [1] * n, 2 * K + 2)
    dmax = 5 if n == 2 else 4
    for inst in range(NI):
        rng = random.Random(7919 * (10 * nx + ny) + inst)
        H, F = rand_phase(cb, rng, list(range(n)), dmax, nx); amp = rand_amp(c, rng, list(range(n)))
        t = time.time()
        tot = zseries(c, total_sp(c, cb, F, amp, list(range(n))))
        G0, it = iterated_sp(c, cb, F, lambda s: cb.compose(amp, s), xs, ys)
        it = zseries(c, it)
        fd = first_diff(tot, it, K)
        log.check('(nx, ny) = (%d, %d) inst %d: total SP = iterated SP through z^%d (H = %s, %d phase terms)'
                  % (nx, ny, inst, K, H, len(F)), fd is None and G0 == 0,
                  'z^1: %s, z^2: %s  (%.1f s)' % (tot.get(1, 0), tot.get(2, 0), time.time() - t))
        if inst == 0:
            for name, kw in [('no factor (det K(y)/det K(0))^(-1/2)', {'nodet': True}),
                             ('inner covariance frozen at y = 0', {'freeze_cov': True}),
                             ('inner expansion at x = 0, linear term dropped', {'at_zero': True})]:
                _, ct = iterated_sp(c, cb, F, lambda s: cb.compose(amp, s), xs, ys, **kw)
                fdc = first_diff(tot, zseries(c, ct), K)
                log.control('(nx, ny) = (%d, %d): %s' % (nx, ny, name), fdc is not None, 'first differs at z^%s' % fdc)
log.finish('CHK_FUBINI %d %d' % (KTOP, NI))
