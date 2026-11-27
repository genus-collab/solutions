# -*- coding: utf-8 -*-
"""chk_delta.py -- Lemma L2 (formal delta lemma), with and without outer variables, exact over Q.

Phase Phi = sigma G(t, w) + f(t, w), amplitude a(t, w; z) free of sigma, critical point (sigma*, 0, 0) with
G(0,0) = 0, G_t(0,0) = g_t != 0 (so the (sigma, t)-block of the Hessian has determinant -g_t^2).
(i)  no outer variables: SP^norm over (sigma, t) equals a(0; z) EXACTLY (all higher Wick terms cancel);
(ii) outer variables w (1 or 2): SP^norm over (sigma, t, w) equals SP^norm_w[ f(tau(w), w) ;
     a(tau(w), w; z) g_t / G_t(tau(w), w) ], tau(w) the root of G(tau(w), w) = 0 (L2 combined with L1b).
Controls (must fail): amplitude depending on sigma; (ii) without the Jacobian factor g_t / G_t(tau(w), w).
Usage: py -I -B chk_delta.py [K] [NINST]     (through z^K, default 8; NINST instances per case, default 2)
"""
import os, sys, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from sp_core import Ctx, rand_poly, zseries
from sp_iter import total_sp
from sp_rand import rand_amp
from chk_util import Log, first_diff

KTOP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NI = int(sys.argv[2]) if len(sys.argv) > 2 else 2
log = Log(HERE, 'chk_delta_%d_%d' % (KTOP, NI))
log.out('chk_delta: Lemma L2 through z^%d, %d instances per case; log %s' % (KTOP, NI, log.rel))


def setup(cb, rng, nw):
    """variables: 0 = s = sigma - sigma*, 1 = t, 2.. = w.  Returns sigma*, G, f, Phi (expanded at 0)."""
    tw = [1] + list(range(2, 2 + nw))
    gt = rng.choice([-3, -2, -1, 1, 2, 3]); sig = Fr(rng.randint(-5, 5), rng.randint(1, 4))
    lin = cb.var(1, gt)
    for v in tw[1:]:
        lin = cb.add(lin, cb.var(v, rng.randint(-2, 2)))
    G = cb.add(lin, rand_poly(cb, rng, tw, 2, 4, 3, 0.6))
    f = cb.add(cb.scal(lin, -sig), rand_poly(cb, rng, tw, 2, 4, 3, 0.7))
    Phi = cb.add(cb.add(cb.scal(G, sig), cb.mul(cb.var(0), G)), f)
    return gt, sig, G, f, Phi


for nw, K in [(0, KTOP), (1, KTOP), (2, min(KTOP, 6))]:
    n = 2 + nw; allv = list(range(n)); wv = list(range(2, n))
    c = Ctx(n, [1] * n, 2 * K); cb = Ctx(n, [1] * n, 2 * K + 2)
    for inst in range(NI):
        rng = random.Random(15485863 * (nw + 1) + inst); t0 = time.time()
        gt, sig, G, f, Phi = setup(cb, rng, nw)
        a = rand_amp(c, rng, [1] + wv)
        tot = zseries(c, total_sp(c, cb, Phi, a, allv))
        if nw == 0:
            a0 = {k[1]: v for k, v in a.items() if not any(k[0])}
            fd = first_diff(tot, a0, K)
            log.check('L2 (i) inst %d: SP over (sigma, t) = a(0; z) exactly through z^%d (sigma* = %s, g_t = %d)'
                      % (inst, K, sig, gt), fd is None, 'a(0; z) = %s  (%.1f s)' % (a0, time.time() - t0))
        else:
            tau = {}
            for _ in range(cb.W2 + 3):                      # G(tau(w), w) = 0, fixed-point iteration
                new = cb.trunc(cb.add(tau, cb.compose(G, {1: tau}), Fr(-1, gt)))
                if new == tau:
                    break
                tau = new
            assert not cb.compose(G, {1: tau}), "root of G not converged"
            fres = cb.compose(f, {1: tau})
            Gt = cb.compose(cb.deriv(G, 1), {1: tau})       # = g_t (1 + eps)
            jac = c.binom_pow(c.trunc(cb.add(cb.scal(Gt, Fr(1, gt)), cb.const(1), -1)), -1)
            ares = c.trunc(cb.compose(a, {1: tau}))
            rhs = zseries(c, total_sp(c, cb, fres, c.mul(ares, jac), wv))
            fd = first_diff(tot, rhs, K)
            log.check('L2 (ii) nw = %d inst %d: SP over (sigma, t, w) = SP_w[f(tau(w), w); a(tau(w), w) g_t/G_t(tau(w), w)]'
                      ' through z^%d' % (nw, inst, K), fd is None, 'z^1: %s  (%.1f s)' % (tot.get(1, 0), time.time() - t0))
            if inst == 0:
                r0 = zseries(c, total_sp(c, cb, fres, ares, wv))
                log.control('L2 (ii) nw = %d: without the Jacobian g_t / G_t(tau(w), w)' % nw,
                            first_diff(tot, r0, K) is not None, 'first differs at z^%s' % first_diff(tot, r0, K))
        if inst == 0:
            asig = c.add(a, c.mul(c.var(0), c.add(c.const(1), c.var(1))))   # amplitude now involves sigma
            ts = zseries(c, total_sp(c, cb, Phi, asig, allv))
            log.control('L2 nw = %d: amplitude depending on sigma (a + s(1 + t))' % nw, first_diff(tot, ts, K) is not None,
                        'first differs at z^%s' % first_diff(tot, ts, K))
log.finish('CHK_DELTA %d %d' % (KTOP, NI))
