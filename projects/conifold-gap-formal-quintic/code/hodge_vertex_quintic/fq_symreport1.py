# -*- coding: utf-8 -*-
"""fq_symreport1.py -- report blocks of genus2_symbol.py, steps 1 (the conifold symbol) and 2 (the weight-<= 3 fit).
Each block takes (log, check) callables and returns data for the cross-prime comparison (block 4 of the driver)."""
import random, time, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_symbol import (Symbol, mat_series_log, mat_series_exp, ser_log, ser_mul, exp_odd_test, power_sums_from_b, h_series, fmt_list)
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt
from fq_wsearch import named_candidates, power_sums, search_blocks, wt_weights
from fq_leadfit import (MONOS3, MONO3_W, mono_name3, expected_weights, weights_of, poly_eval3, fit3, collect_samples, poly_qzeta3,
                        shift_check, M3V)

F_EXPECTED = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000)]


def ratbig(R, a, B=10 ** 13):
    """rational reconstruction with a larger bound (2^89-1 only; spurious rate ~ B^2/p ~ 1e-1 -- cross-checked at the second prime)"""
    r = R.ratrec(a, B)
    return str(r) if r is not None else "[%d mod p]" % (a % R.p)


def qz(R, vals4):
    return qzeta_fmt(R, qzeta_coords(R, vals4))


def step1(S, log, check, R):
    p = R.p; t0 = time.time(); sym = Symbol(S); N = sym.Nz; data = {}
    check("symbol: val_s Q^{(i)}_{a,n} = -4n + c_a - 1 exactly (non-zero leading coefficient) for all 5 x 5 x %d (a, i, n), n <= %d" % (N + 1, N), sym.val_ok, sym.val_bad[:5])
    bad = sym.symplectic_check(); check("symbol: R^ symplectic as a constant series: Sum_ab R^_a^{(i)}(zhat) eta^{ab} R^_b^{(j)}(-zhat) = delta_ij through zhat^%d" % N, not bad, bad[:3])
    cov = sym.covariance()
    check("symbol: Z/5 covariance R^_a^{(i)}(zhat) = zeta^{i(a+1)} R^_a^{(0)}(zeta^{-i} zhat) (phase a+1-n per coefficient)", cov['phase_a_plus_1_minus_n'], cov['k'])
    log("   (the phase zeta^{ia} of the P-rows, R^_a^{(i)}(zhat) = zeta^{ia} R^_a^{(0)}(zeta^{-i} zhat), holds: %s -- the extra zeta^i is Delta_i^{-1/2} = zeta^{-4i}/x)" % cov['phase_a_minus_n'])
    # flat-basis grading
    Rf = [[[sum(sym.Rh[a][i][n] * sym.Rh[b][i][0] for i in range(5)) % p for b in range(5)] for a in range(5)] for n in range(N + 1)]
    grade = all(Rf[n][a][b] == 0 for n in range(N + 1) for a in range(5) for b in range(5) if (a + b - n + 2) % 5 != 0)
    check("symbol in the flat basis: R^_{ab,n} := Sum_i R^_{a,n}^{(i)} R^_{b,0}^{(i)} vanishes unless a + b = n - 2 mod 5 (zhat^n shifts the flat degree by n mod 5)", grade)
    # (a) colour-colour form
    C = sym.colour_matrix(); diag = sym.is_diagonal(C); circ = sym.circulant_check(C)
    log("   (a) Teleman colour-colour form R~_{ji}(zhat) = Sum_ab R^_{b,0}^{(j)} eta^{ba} R^_a^{(i)}(zhat); coefficient matrices C_n (rows j, columns i), Q(zeta) via Galois (sigma_m: C_n[mj][mi]):")
    for n in range(N + 1):
        log("   C_%d:" % n)
        for j in range(5): log("      [%s]" % ", ".join(qz(R, [C[n][(m * j) % 5][(m * i) % 5] for m in range(1, 5)]) if n < 6 else ratbig(R, C[n][j][i]) for i in range(5)))
    check("symbol: C_0 = identity (symplecticity at zhat^0)", all(C[0][j][i] == (1 if i == j else 0) for i in range(5) for j in range(5)))
    log("   R~ diagonal in colours: %s;  Z/5-circulant with phase (C_n[j+1][i+1] = zeta^{-n} C_n[j][i]): %s" % (diag, circ))
    check("symbol: R~ is NOT diagonal in colours (every off-diagonal entry of C_1 is non-zero) and IS circulant with phase zeta^{-n}",
          (not diag) and circ and all(C[1][j][i] for i in range(5) for j in range(5) if i != j))
    mu = mat_series_log(C, R, N)
    sym_ok = all((mu[k][i][j] - (-1) ** (k + 1) * mu[k][j][i]) % p == 0 for k in range(1, N + 1) for i in range(5) for j in range(5))
    check("symbol: log R~ = Sum mu_k zhat^k has mu_k^T = (-1)^{k+1} mu_k (symplectic log), k <= %d" % N, sym_ok)
    log("   first rows (j = 0) of mu_1, mu_2, mu_3 in Q(zeta) [sigma_m(mu_k[0][i]) = mu_k[0][mi]]:")
    for k in (1, 2, 3): log("      mu_%d[0][.] = [%s]" % (k, ", ".join(qz(R, [mu[k][0][(m * i) % 5] for m in range(1, 5)]) for i in range(5))))
    check("symbol: mu_2 != 0 (an antisymmetric EVEN part is present: R~ is not exp(odd) even as a matrix)", any(x for r in mu[2] for x in r))
    # (c) diagonal part exponential test
    d0 = sym.diag_series(C, 0); ev, od = exp_odd_test(d0, R)
    log("   (c) diagonal entry R~_00(zhat) = %s" % fmt_list(R, d0))
    log("       log R~_00: even coefficients %s; odd %s" % ({k: R.fmt(v) for k, v in ev.items()}, {k: R.fmt(v) for k, v in od.items()}))
    check("symbol: diagonal entries R~_ii are exp(odd series)", all(v == 0 for v in ev.values()), "even coefficients non-zero: %s" % {k: R.fmt(v) for k, v in ev.items() if v})
    check("symbol: R~_ii(zhat) = R~_00(zeta^{-i} zhat) (covariance of the diagonal)", all(C[n][i][i] == C[n][0][0] * R.z(-n * i) % p for n in range(N + 1) for i in range(5)))
    # identity row f
    f = sym.identity_row(); log("   identity row f(zhat) = R^_0^{(0)}(zhat) = [%s]" % ", ".join(ratbig(R, x) for x in f))
    check("symbol: f coefficients = 1, -3/20, 9/800, 11571/80000, -2415609/6400000, 47735919/128000000 (zhat^0..5)", all(f[n] == R.of(F_EXPECTED[n]) for n in range(6)))
    evf, odf = exp_odd_test(f, R)
    log("   log f: even %s; odd (b_k at zhat^{2k-1}) %s" % ({k: ratbig(R, v) for k, v in evf.items()}, {k: ratbig(R, v) for k, v in odf.items()}))
    check("symbol: f is exp(odd series)", all(v == 0 for v in evf.values()), "even coefficients: %s" % {k: ratbig(R, v) for k, v in evf.items() if v})
    ps = power_sums_from_b(R, odf); log("   implied power sums p_{2k-1} = -b_k 2k(2k-1)/B_{2k}: %s" % {k: ratbig(R, v) for k, v in ps.items()})
    log("   (reference: wt_0 = {1 - zeta^j} U {-5} has p_1, p_3, p_5 = %s)" % [R.fmt(x) for x in power_sums(R, wt_weights(R, 0))])
    target = (ps[1], ps[3], ps[5])
    nc = named_candidates(R); hits = [name for name, v in nc if power_sums(R, v) == target]
    hits1 = [name for name, v in nc if power_sums(R, v)[0] == target[0]]
    log("   named candidates (c * wt_i, c in {+-1, +-5, +-1/5}, all sub-multisets: %d): full (p_1, p_3, p_5) match %s; p_1 alone matched by %s" % (len(nc), hits, hits1))
    t1 = time.time(); bl, nb = search_blocks(R, target, 6); ctrl, _ = search_blocks(R, power_sums(R, wt_weights(R, 0)), 6)
    log("   block search (%d Galois-orbit/rational blocks, multisets of size <= 6, %.1f s): matches for the symbol's (p_1, p_3, p_5): %s; control (wt_0's power sums): %s" % (nb, time.time() - t1, bl, ctrl))
    check("symbol: no weight multiset of size <= 6 from the candidate pool realises (p_1, p_3, p_5) of f (and f is not exp(odd) anyway)", not bl and not hits)
    # QRR relation
    h = h_series(R, 0, N); hinv = h_series(R, 0, N, -1); g = ser_mul(hinv, f, p)
    log("   h_{wt_0}(zhat) = [%s]" % ", ".join(ratbig(R, x) for x in h))
    log("   h_{wt_0}^{-1} f = [%s]" % ", ".join(ratbig(R, x) for x in g))
    lg = ser_log(g, R); log("   log(h^{-1} f) = [%s]" % ", ".join(ratbig(R, x) for x in lg))
    check("symbol: f = h_{wt_0}(zhat) (1 + (37/250) zhat^3 + O(zhat^4)): agreement with the QRR constant through zhat^2 only", g[1] == 0 and g[2] == 0 and g[3] == R.of(Fr(37, 250)))
    log("   step 1 in %.1f s" % (time.time() - t0))
    data.update(C=C, mu=mu, f=f, g=g, lg=lg, ps=ps, d0=d0, hits=hits, bl=bl)
    return sym, data


def step2(S, graphs, log, check, R, rng, nsamp=1010, nver=20, SympM=None, GraphSummer=None):
    p = R.p; t0 = time.time(); ORDERS = list(range(-10, 3)); keys = [(G['idx'], k) for G in graphs for k in ORDERS]
    log("   %d + %d random symplectic M (fq_sser.SympM.random), orders s^-10..s^2, %d monomials of weight <= 3" % (nsamp, nver, len(MONOS3)))
    samples = collect_samples(S, graphs, nsamp + nver, rng, ORDERS, SympM, GraphSummer, log=log)
    log("   evaluations in %.1f s" % (time.time() - t0)); t1 = time.time()
    polys, rank, nbad = fit3(samples[:nsamp], keys, p)
    check("fit: rank %d = 991 monomials, %d inconsistent surplus rows (%.1f s)" % (rank, nbad, time.time() - t1), rank == 991 and nbad == 0)
    check("fit verified on %d extra random M for all %d (graph, order) coefficients" % (nver, len(keys)), all(poly_eval3(polys[k], par, p) == vals[k] for (par, vals) in samples[nsamp:] for k in keys))
    allok = True
    for G in range(1, 8):
        parts = []
        for k in ORDERS:
            ws = weights_of(polys[(G, k)]); exp_w = expected_weights(k)
            if ws: parts.append("s^%d: w%s" % (k, ws))
            if not set(ws) <= set(exp_w): allok = False
        log("   graph %d: %s  (all other orders identically 0)" % (G, "  ".join(parts)))
    check("weight structure: s^-10, s^-5, s^0 weight 0; s^-6, s^-1 weight 1; s^-2 weight 2; s^+2 weight 3 ONLY (no weight-0 part); s^-9..-7, -4, -3, 1 zero",
          allok and all(weights_of(polys[(G, 2)]) == [3] and weights_of(polys[(G, -2)]) == [2] and weights_of(polys[(G, -6)]) == [1] for G in range(1, 8)))
    check("m3 never at a polar order; m3 (weight 3) present at s^+2 for graphs 1, 2, 4 only; m2 at s^-2 and s^+2 (via m1 m2) for graphs 1, 2, 4, 5",
          all(not any(v in M3V for k in range(-10, 0) for m in polys[(G, k)] for v in MONOS3[m]) for G in range(1, 8))
          and [G for G in range(1, 8) if any(v in M3V for m in polys[(G, 2)] for v in MONOS3[m])] == [1, 2, 4])
    nco = sum(len(v) for v in polys.values()); nsh = sum(shift_check(polys[k], R) for k in polys)
    check("colour-shift covariance c_{shift mu} = zeta^w c_mu on all %d fitted coefficients" % nco, nsh == 0, nsh)
    qzd = {k: poly_qzeta3(R, polys[k]) for k in polys}
    nrec = sum(1 for v in qzd.values() for cs in v.values() if qzeta_ratrec(R, cs) is None)
    log("   Galois reconstruction in Q(zeta) of the %d coefficients: %d not reconstructible by fq_fp.ratrec (bound ~1e10 at 2^89-1, ~1e7 at 2^61-1)" % (nco, nrec))
    for G in (1, 2, 4):
        ms = sorted((m for m in polys[(G, 2)] if any(v in M3V for v in MONOS3[m])), key=lambda m: MONOS3[m])
        log("   graph %d, s^+2, m3 coefficients: %s" % (G, "; ".join("%s: %s" % (mono_name3(MONOS3[m]), qzeta_fmt(R, qzd[(G, 2)][m])) for m in ms[:5]) + (" ..." if len(ms) > 5 else "")))
    log("   step 2 in %.1f s" % (time.time() - t0))
    return polys, qzd, samples[nsamp:]
