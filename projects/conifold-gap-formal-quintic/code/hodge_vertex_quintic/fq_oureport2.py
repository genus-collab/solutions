# -*- coding: utf-8 -*-
"""fq_oureport2.py -- report blocks of genus2_ou.py, steps 4 (solving J n = -(g1_2, g1_4) on the classes, right and left placement),
5 (verification by eps -> u of one solution per solvable class) and the frame-comparison summary of step 6."""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_rref, VNAME
from fq_scalar_lib import qzeta_coords, qzeta_fmt
from gsum_q import is_u_series
from fq_hat import ucoef, uprec
from fq_ou import (GRAPHS, class_columns, restricted, solve_restricted, expand_solution, N_of_params, verify_solution)
from fq_ourec import circ_qzeta, circ_conjugates, fmt_rat, fmt_qz_fp, CIRC_NAMES
from fq_oureport import F2_U, GNAME

CLASSES = [('cov', "(i) covariant scalar: nu_1 = b_1 diag(zeta^-i), nu_3 = b_3 diag(zeta^-3i)", 2),
           ('diag', "(ii) diagonal, Fourier modes X_a, Y_a (x_i = (1/5) Sum_a zeta^-ai X_a)", 10),
           ('circ', "(iii) circulant (Z/5-covariant): c_1(0..2), c_2(1..2), c_3(0..2)", 8),
           ('full', "(iv) full sp: 40 coordinates", 40)]


def step4(J, g1, M0s, log, check, R, rec, tag):
    p = R.p; t0 = time.time(); b = [(-g1[2]) % p, (-g1[4]) % p]; sols = {}; cov_t = {}
    log("   system: J[2] . n = -g1_2 = %s,  J[4] . n = -g1_4 = %s" % (fmt_rat(R, b[0]), fmt_rat(R, b[1])))
    for cls, desc, k in CLASSES:
        colsA, names = class_columns(R, M0s.M, M0s.Minv, cls, 'right')
        spans = {}
        for pl in ('right', 'left'):
            cols = colsA if pl == 'right' else class_columns(R, M0s.M, M0s.Minv, cls, 'left')[0]
            spans[pl] = cols
            Jr = restricted(J, cols, p); x, null, rank = solve_restricted(Jr, b, p); cons = x is not None
            dim = len(null) if cons else None
            key = "%s:%s:%s" % (tag, cls, pl)
            log("   %s, %s placement [M = %s]: restricted 2 x %d matrix rank %d, consistent %s, solution dimension %s" % (desc, pl.upper(), "M_0 exp(uN)" if pl == 'right' else "exp(uN') M_0", k, rank, cons, dim))
            if cls == 'cov':
                for g, G in enumerate((2, 4)):
                    log("      equation graph %d: (%s) b_1 + (%s) b_3 = %s" % (G, fmt_rat(R, Jr[g][0]), fmt_rat(R, Jr[g][1]), fmt_rat(R, b[g])))
                    for j in range(2): rec.rat("%s:J%d:%s" % (key, G, names[j]), Jr[g][j])
                check("%s: 2 x 2 matrix invertible (rank 2, the two equations independent), unique solution; RATIONAL entries and solution" % key,
                      rank == 2 and cons and all(R.ratrec(v) is not None for row in Jr for v in row) and (p < 10 ** 20 or all(R.ratrec(v) is not None for v in x)))
                if cons:
                    log("      solution (b_1, b_3) = (%s, %s)" % (fmt_rat(R, x[0]), fmt_rat(R, x[1])))
                    rec.rat("%s:b1" % key, x[0]); rec.rat("%s:b3" % key, x[1]); cov_t[pl] = x
            elif cls == 'diag':
                for g, G in enumerate((2, 4)):
                    log("      graph %d row (X_0..X_4 | Y_0..Y_4): [%s]" % (G, ", ".join(fmt_rat(R, v) for v in Jr[g])))
                    for j in range(k): rec.rat("%s:J%d:%s" % (key, G, names[j]), Jr[g][j])
                check("%s: restricted matrix RATIONAL in the Fourier modes; rank %d; consistent %s" % (key, rank, cons), cons and all(R.ratrec(v) is not None for row in Jr for v in row))
                if cons:
                    log("      particular solution (RREF, free modes 0): %s; nullspace dimension %d" % (", ".join("%s = %s" % (names[j], fmt_rat(R, x[j])) for j in range(k) if x[j]), dim))
                    for j in range(k): rec.rat("%s:sol:%s" % (key, names[j]), x[j])
                    xi = [sum(R.z(-a * i) * x[a] for a in range(5)) * R.inv(5) % p for i in range(5)]; yi = [sum(R.z(-a * i) * x[5 + a] for a in range(5)) * R.inv(5) % p for i in range(5)]
                    log("      i.e. x_i = [%s], y_i = [%s] (Q(zeta) via sigma_m x_i = x_{mi})" % (", ".join(qzeta_fmt(R, qzeta_coords(R, [xi[(m * i) % 5] for m in range(1, 5)])) for i in range(5)),
                                                                                                     ", ".join(qzeta_fmt(R, qzeta_coords(R, [yi[(m * i) % 5] for m in range(1, 5)])) for i in range(5))))
            elif cls == 'circ':
                for g, G in enumerate((2, 4)):
                    qz = circ_qzeta(R, Jr[g]); conj = circ_conjugates(R, Jr[g])
                    log("      graph %d row (%s): [%s]" % (G, ", ".join(CIRC_NAMES), ", ".join(fmt_qz_fp(R, c) for c in qz)))
                    for j in range(k): rec.qz("%s:J%d:%s" % (key, G, names[j]), [conj[m][j] for m in range(4)])
                emb = [0] * k
                if pl in cov_t: emb[0] = cov_t[pl][0]; emb[5] = cov_t[pl][1]
                ok_emb = pl in cov_t and all(sum(Jr[g][j] * emb[j] for j in range(k)) % p == b[g] for g in range(2))
                check("%s: consistent, rank %d, solution set of dimension %s; the covariant-scalar solution (c_1(0) = b_1, c_3(0) = b_3) is a particular solution" % (key, rank, dim), cons and ok_emb)
                if cons: x = emb if ok_emb else x
            else:
                emb = expand_solution(colsA, [0] * 40, p)
                if pl in cov_t:
                    cc, _ = class_columns(R, M0s.M, M0s.Minv, 'cov', 'right'); emb = expand_solution(cc, cov_t[pl], p)
                ok_emb = pl in cov_t and all(sum(Jr[g][j] * emb[j] for j in range(40)) % p == b[g] for g in range(2))
                check("%s: rank 2, consistent, affine solution space of dimension 38; the covariant-scalar solution is a particular solution" % key, rank == 2 and cons and dim == 38 and ok_emb)
                if cons: x = emb if ok_emb else x
            if cons:
                Ncls = expand_solution(colsA, x, p)                     # the class element (N for right, N' for left)
                sols[(cls, pl)] = dict(t=x, N=Ncls, null=null, rank=rank, dim=dim, cols=colsA)
            else:
                sols[(cls, pl)] = dict(t=None, N=None, null=null, rank=rank, dim=None, cols=colsA)
        both = spans['right'] + spans['left']; rk_both = len(fp_rref([list(c) for c in both], 40, p)[0])
        same_space = rk_both == k
        if cls == 'circ':
            check("%s:circ: the circulant class is Ad(M_0)-stable (M_0 is circulant): right and left placement give the SAME set of N (combined rank %d = 8)" % (tag, rk_both), same_space)
        elif cls == 'full':
            check("%s:full: left = right as sets (Ad(M_0) invertible, combined rank %d = 40)" % (tag, rk_both), same_space)
        else:
            log("      right and left placement span the same N-space: %s (combined rank %d)" % (same_space, rk_both))
    log("   step 4 in %.1f s" % (time.time() - t0))
    return sols


def step5(H, M0s, sols, graphs, Vinv, log, check, R, rec, tag, tot1, c1, keys=None, random_full=True):
    """tot1: the point-frame total [u^0..u^4] of this frame (M-independent), c1 = [u^2] t_c"""
    p = R.p; t0 = time.time(); out = {}; g1tot = (-2 * c1 * R.of(Fr(1, 1200))) % p
    keys = keys or [k for k in sols if sols[k]['N'] is not None]
    for key in keys:
        s = sols[key]
        if s['N'] is None: log("   %s: no solution to verify" % (key,)); continue
        cls, pl = key; N = N_of_params(s['N'], p)
        vs = verify_solution(H, M0s.M, M0s.Minv, N, graphs, R, Vinv, pl)
        name = "%s:%s:%s" % (tag, cls, pl); uc = vs['uc']
        log("   %s: eps-polynomial P_Gamma(eps; u) from 10 hatted sums with %s, degree in eps %d; eps = u substituted; u-coefficients through u^4:" % (name, "M_0 exp(eps N)" if pl == 'right' else "exp(eps N') M_0", vs['maxdeg']))
        for G in GRAPHS: log("      graph %d: [%s]" % (G, ", ".join(fmt_rat(R, c) for c in uc[G])))
        check("%s: eps-degree <= 3; eps = u: ALL SIX edge graphs are O(u^2) (REGULAR at the conifold point)" % name, vs['maxdeg'] <= 3 and vs['edge_regular'])
        check("%s: vertex graph [u^0] = 1/1200, [u^1] = -2 c_1/1200 = %s (1/t_c coefficient 0): the whole polar part (1/1200) t_c^-2 sits in the vertex graph" % (name, fmt_rat(R, g1tot)),
              uc[1][0] == R.of(Fr(1, 1200)) and uc[1][1] == g1tot)
        check("%s: total through u^4 = the point-frame total F_2^cf s^10 of this frame (M-independence holds for u-dependent M); all series pure u-series" % name,
              all(vs['total'][k] == tot1[k] for k in range(min(5, len(tot1)))) and all(is_u_series(vs['Gu'][G]) for G in GRAPHS))
        log("      vertex graph u^2, u^3, u^4 coefficients: %s, %s, %s;  edge graphs' u^2 coefficients: %s" % (fmt_rat(R, uc[1][2]), fmt_rat(R, uc[1][3]), fmt_rat(R, uc[1][4]), ", ".join(fmt_rat(R, uc[G][2]) for G in GRAPHS if G != 1)))
        for G in GRAPHS:
            for k in range(2, 5): rec.rat("%s:ver:g%d:u%d" % (name, G, k), uc[G][k])
        out[key] = uc
    if random_full and sols.get(('full', 'right'), {}).get('N') is not None:
        rng = random.Random(p % 1009 + 5); s = sols[('full', 'right')]
        n = list(s['N'])
        for nv in s['null']:
            t = rng.randrange(p); n = [(a + t * c) % p for a, c in zip(n, nv)]
        vs = verify_solution(H, M0s.M, M0s.Minv, N_of_params(n, p), graphs, R, Vinv, 'right')
        check("%s: a RANDOM point of the 38-dimensional solution space (particular + random nullspace combination) also regularises all six edge graphs (eps-degree %d)" % (tag, vs['maxdeg']), vs['edge_regular'] and vs['maxdeg'] <= 3)
        log("      (vertex graph there: [%s])" % ", ".join(R.fmt(c) for c in vs['uc'][1]))
    log("   step 5 in %.1f s" % (time.time() - t0))
    return out


def frame_summary(dA, dB, log, check, R, rec):
    """what changed between (alpha, beta) = (0, 0) [dA] and (1/3, -2/7) [dB]; both dicts hold g1, A (t_c data), sols"""
    p = R.p
    log("   (alpha, beta) = (0, 0) -> (1/3, -2/7): c_1 of t_c %s -> %s" % (R.fmt(dA['c1']), R.fmt(dB['c1'])))
    log("   g1_Gamma: %s -> %s" % ({G: fmt_rat(R, dA['g1'][G]) for G in GRAPHS}, {G: fmt_rat(R, dB['g1'][G]) for G in GRAPHS}))
    log("   g1_2, g1_4 changed with (alpha, beta): %s; g1_1 shifted by 2 alpha/1200: %s; at the second frame g1_3 = g1_5 = g1_6 = g1_7 = 0 and Sum g1 = -2 c_1/1200: %s"
        % (any(dA['g1'][G] != dB['g1'][G] for G in (2, 4)), (dB['g1'][1] - dA['g1'][1]) % p == (2 * dB['alpha'] * R.of(Fr(1, 1200))) % p,
           all(dB['g1'][G] == 0 for G in (3, 5, 6, 7)) and sum(dB['g1'].values()) % p == (-2 * dB['c1'] * R.of(Fr(1, 1200))) % p))
    log("   per-graph 1/t_c coefficients with M_0: %s -> %s" % ({G: fmt_rat(R, dA['A'][G][1]) for G in GRAPHS}, {G: fmt_rat(R, dB['A'][G][1]) for G in GRAPHS}))
    log("   per-graph 1/t_c coefficients at M = 1 (point frame): %s -> %s" % ({G: fmt_rat(R, dA['A1'][G][1]) for G in GRAPHS}, {G: fmt_rat(R, dB['A1'][G][1]) for G in GRAPHS}))
    for pl in ('right', 'left'):
        a = dA['sols'][('cov', pl)]['t']; bb = dB['sols'][('cov', pl)]['t']
        log("   covariant-scalar solution (b_1, b_3), %s placement: (%s, %s) -> (%s, %s)" % (pl, fmt_rat(R, a[0]), fmt_rat(R, a[1]), fmt_rat(R, bb[0]), fmt_rat(R, bb[1])))
    log("   covariant-scalar solutions changed with (alpha, beta): %s" % {pl: dA['sols'][('cov', pl)]['t'] != dB['sols'][('cov', pl)]['t'] for pl in ('right', 'left')})
    check("frame dependence: the totals t_c^-2 = 1/1200 and t_c^-1 = 0 persist at (1/3, -2/7) with M_0 (per-graph values logged above)",
          sum(dB['A'][G][0] for G in GRAPHS) % p == R.of(Fr(1, 1200)) and sum(dB['A'][G][1] for G in GRAPHS) % p == 0)
