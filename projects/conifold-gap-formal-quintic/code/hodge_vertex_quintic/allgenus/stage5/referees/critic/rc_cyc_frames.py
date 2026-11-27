# -*- coding: utf-8 -*-
"""rc_cyc_frames.py -- part E of rc_cyc.py (executed inside it; shares its globals f, g, A, XS, N, check, say, ...).
Frames, the factorised operator M, the Corollary in operator form, unit shapes, controls; then the log is written."""


def pf4_tw(x):                    # u^-x PF4 u^x
    out, gg = None, ONE
    for k in range(5):
        t = A[k] * gg
        out = t if out is None else out + t
        gg = gg.Dx(x)
    return out


PF4TW = {x: pf4_tw(x) for x in XS}


def frame_data(y0, y1, mid=None):
    t = y1 / y0
    Dt = t.D()
    c = [Dt, U.inv() / (y0 * y0 * Dt * Dt), Dt, y0, y0]          # L^5 = 1/u
    middle = (U if mid is None else mid) * y0 * y0 * Dt * Dt
    chain = [('m', y0.inv()), ('D',), ('m', Dt.inv()), ('D',), ('m', middle), ('D',), ('m', Dt.inv()), ('D',),
             ('m', y0.inv())]
    return t, Dt, c, chain


def agree(chain, target):         # twisted images at all x agree; returns (ok, precision reached)
    ok, worst = True, BIG
    for x in XS:
        d = twisted(chain, x) - target(x)
        worst = min(worst, d.prec)
        ok &= (d.c == {})
    return ok, worst


def cor_chain(c, y0):
    ch = []
    for a in (4, 3, 2, 1, 0):
        ch += [('m', c[a].inv()), ('D',)]
    return ch + [('m', y0.inv())]


def cor_target(y0):
    return lambda x: y0.inv() * PF4TW[x].Dx(x)


f0, f1, f2 = f
fr = [('programme (alpha,beta) = (0,0)', f0, g, (0, -1, 0, 0, 0)),
      ('programme (1/3,-2/7)', f0 + g * Fr(1, 3) - f2 * Fr(2, 7), g, (0, -1, 0, 0, 0)),
      ('programme (7/11,-3/5)', f0 + g * Fr(7, 11) - f2 * Fr(3, 5), g, (0, -1, 0, 0, 0)),
      ('critic (a,b) = (1/3,-2/7): y0 = f0 + f1/3 - 2f2/7', f0 + f1 * Fr(1, 3) - f2 * Fr(2, 7), g, (0, -1, 0, 0, 0)),
      ('critic (7/11,-3/5)', f0 + f1 * Fr(7, 11) - f2 * Fr(3, 5), g, (0, -1, 0, 0, 0)),
      ('rescaled (2 f0, 3 g)', f0 * 2, g * 3, (0, -1, 0, 0, 0)),
      ('K2: y0 = f1 (y0(0) = 0)', f1, g, (0, -3, 0, 1, 1)),
      ('K2: y0 = f2', f2, g, (-2, -1, -2, 2, 2)),
      ('K2: y0 = f1 - 2 f2', f1 - f2 * 2, g, (0, -3, 0, 1, 1)),
      ('K2: Lagrangian pair (f1, f2)', f1, f2, (0, -3, 0, 1, 1))]
d1 = fr[3][1] - (f0 + g * Fr(1, 3) - f2 * Fr(109, 210))
d2 = fr[4][1] - (f0 + g * Fr(7, 11) - f2 * Fr(23, 22))
d3 = fr[1][1] - (f0 + f1 * Fr(1, 3) - f2 * Fr(11, 210))
check('E0 relabelling (E3 = K1): critic (1/3,-2/7) = programme (1/3,-109/210); critic (7/11,-3/5) = programme (7/11,-23/22); '
      'programme (1/3,-2/7) = f0 + f1/3 - (11/210) f2', d1.c == {} and d2.c == {} and d3.c == {})
check('E0 u^2 coefficients: programme (1/3,-2/7) -> -11/210, critic (1/3,-2/7) -> -2/7',
      fr[1][1].coef(2) == Fr(-11, 210) and fr[3][1].coef(2) == Fr(-2, 7))
t00, Dt00, c00, _ = frame_data(f0, g)
check('E0 frame (0,0): t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4 (FORMAL_QUINTIC_genus2 s.1), varpi0 = 1 + (2/625)u^3 + ...',
      [t00.coef(k) for k in range(5)] == [0, 1, Fr(7, 10), Fr(41, 75), Fr(9, 20)] and f0.coef(3) == Fr(2, 625))
tp = t00.der()
urb = [tp, -(U * f0 * f0 * (ONE - U) * (ONE - U) * (ONE - U) * tp * tp).inv(), tp, -f0 / (ONE - U), -f0 / (ONE - U)]
check('E0 UR s.0 item 1 b-vector = -c/(1-u) in frame (0,0), all five entries (below u^%d)' % (N - 2),
      all(((c00[a] * (-1) / (ONE - U)) - urb[a]).zero_below(N - 2) for a in range(5)))
for name, y0, y1, shape in fr:
    t, Dt, c, chain = frame_data(y0, y1)
    ok, pr = agree(chain, lambda x: PF4TW[x])
    check('E1 M(y0, y1) = PF4 on u^x (six x) below u^%d, %s' % (pr, name), ok and pr >= N - 8)
    ok2, pr2 = agree(cor_chain(c, y0), cor_target(y0))
    check('E2 Corollary (1/c4)D(1/c3)D(1/c2)D(1/c1)D(1/c0)D o (1/y0) = (1/y0) D PF4 below u^%d, %s' % (pr2, name),
          ok2 and pr2 >= N - 12)
    check('E3 unit shape val(c) = %s, %s' % (shape, name), tuple(x.val() for x in c) == shape)
say('   controls (each must fail):')
for name, y0, y1 in [('non-isotropic pair (f0, f1)', f0, f1), ('wrong gamma: (f0, f1 + f2/2)', f0, f1 + f2 * Fr(1, 2)),
                     ('wrong gamma: (f0, f1 + (701/1000) f2)', f0, f1 + f2 * Fr(701, 1000))]:
    t, Dt, c, chain = frame_data(y0, y1)
    check('E4 M = PF4 for %s' % name, agree(chain, lambda x: PF4TW[x])[0], control=True)
    check('E4 the Theorem-U unit shape (0,-1,0,0,0) still holds for %s (K3: the shape cannot see the frame)' % name,
          tuple(x.val() for x in c) == (0, -1, 0, 0, 0))
t, Dt, c, chain = frame_data(f0, g, mid=ONE)
check('E4 M = PF4 with middle factor y0^2 Dt^2 (no u), frame (0,0)', agree(chain, lambda x: PF4TW[x])[0], control=True)
t, Dt, c, chain = frame_data(f0, g)
csw = [c[0], c[2], c[1], c[3], c[4]]
check('E4 Corollary with c_1, c_2 swapped, frame (0,0)', agree(cor_chain(csw, f0), cor_target(f0))[0], control=True)
c5 = [c[0], c[1] * 5, c[2], c[3], c[4]]
check('E4 Corollary with 5 c_1 (Yukawa factor 5 in the row coefficient, K6), frame (0,0)',
      agree(cor_chain(c5, f0), cor_target(f0))[0], control=True)
say('rc_cyc %d: %s   (%.1f s)' % (N, 'ALL CHECKS PASSED' if NFAIL[0] == 0 else '%d FAILURES' % NFAIL[0], time.time() - T0))
with open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8') as fh:
    fh.write('# %s (py -I -B rc_cyc.py %d)\n' % (LOGREL, N))
    fh.write('\n'.join(OUT[:-1]) + '\n' + OUT[-1].split('   (')[0] + '\n')
print('log written: ' + LOGREL)
