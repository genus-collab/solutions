# -*- coding: utf-8 -*-
"""fq_symreport3.py -- report block of genus2_symbol.py, step 3 (the leading-order system for M_0(zhat)).
step3(polys, sym, graphs, log, check, R, rng) returns data for the cross-prime comparison."""
import random, time, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_sser import SympM
from fq_msum import GraphSummer
from fq_mfit import params_of, fp_solve
from fq_symbol import mat_series_log
from fq_scalar_lib import qzeta_coords, qzeta_fmt, p_fmt
from fq_leadfit import MONOS3, MONO3_W, mono_name3, poly_eval3, fit3
from fq_leadsolve import (lead_poly, FrozenRows, frozen_sums, mu_params, affine_in_23, jacobian, rank_fp, fibre_point, restrict_poly,
                          class_diagonal, class_covariant_diagonal, class_circulant, modes_diagonal, tpoly_fmt, tpoly_eval)
from fq_leadclass import tp_vars, tp_to_uni, eliminate_linear, uni_gcd_all, uni_roots, tp_eval
from fq_leadmodes import diagonal_solve, random_W_points, tp_eval_partial
from fq_leadline import solve_on_plane, on_line_check, residual
from fq_leadcol import coloured_polys, monochrome_table, solve_percolour, theta_constants, mixed_structure, covariant_percoloured


def circ_qzeta(R, c):
    """Q(zeta) coordinates of a Galois-stable circulant parameter vector (c_0, c_1, c_2): sigma_m(c(e)) = c(m e), c(3) = zeta^2 c(2), c(4) = zeta c(1)"""
    p = R.p; full = [c[0], c[1], c[2], R.z(2) * c[2] % p, R.z(1) * c[1] % p]
    return [qzeta_fmt(R, qzeta_coords(R, [full[(m * e) % 5] for m in range(1, 5)])) for e in range(3)]


def step3(polys, sym, graphs, log, check, R, rng):
    p = R.p; t0 = time.time(); data = {}
    L = {G: lead_poly(polys, G) for G in range(1, 8)}; edge = [L[G] for G in range(2, 8)]
    log("   L_Gamma(mu) = c_Gamma + [w1 of s^-6] + [w2 of s^-2] + [w3 of s^+2]: monomial counts %s" % {G: len(L[G]) for G in L})
    fr = FrozenRows(sym)
    res0, _ = frozen_sums(fr, SympM.identity(R), graphs)
    check("frozen theory (rows R^, Delta^_i^{1/2} = zeta^{4i}) at M_0 = 1 reproduces every c_Gamma and the total 1/1200",
          all(res0[G] == L[G].get(0, 0) for G in L) and sum(res0.values()) % p == R.of(Fr(1, 1200)))
    bad = 0
    for _ in range(5):
        M = SympM.random(R, rng); res, _ = frozen_sums(fr, M, graphs); par = params_of(M)
        if any(poly_eval3(L[G], par, p) != res[G] for G in L): bad += 1
        if sum(res.values()) % p != R.of(Fr(1, 1200)): bad += 1
    check("KEY CHECK: L_Gamma(mu) (substitution m_k = mu_k s^-4k in the s-series fit) = frozen graph sum with rows R^ M_0^-1 and vertex M_0.omega^ at 5 random mu; frozen total = 1/1200", bad == 0, bad)
    t1 = time.time(); fs = []
    for t in range(1030):
        M = SympM.random(R, rng); res, _ = frozen_sums(fr, M, graphs); fs.append((params_of(M), res))
    fpolys, rank, nbad = fit3(fs[:1010], list(range(1, 8)), p)
    check("frozen-theory fit (1010 random mu, 991 monomials): rank %d, %d inconsistent; verified on 20 more; EVERY coefficient equals L_Gamma (%.1f s)" % (rank, nbad, time.time() - t1),
          rank == 991 and nbad == 0 and all(poly_eval3(fpolys[G], par, p) == v[G] for (par, v) in fs[1010:] for G in range(1, 8)) and all(fpolys[G] == L[G] for G in L))
    # tautological point
    C = sym.colour_matrix(); mu = mat_series_log(C, R, 3); tpar = mu_params(mu[1], mu[2], mu[3], p); Mt = SympM(R, mu[1], mu[2], mu[3])
    vals = [poly_eval3(L[G], tpar, p) for G in range(1, 8)]
    check("(i) tautological mu = log R~ (zhat^1..3): L_1 = 1/1200 and L_Gamma = 0 for the six edge graphs; same from the frozen graph sum",
          vals[0] == R.of(Fr(1, 1200)) and all(v == 0 for v in vals[1:]) and frozen_sums(fr, Mt, graphs)[0] == {G: (R.of(Fr(1, 1200)) if G == 1 else 0) for G in range(1, 8)})
    # (ii) structure of the full system
    J = jacobian(edge, tpar, p); rk = rank_fp(J, 40, p)
    log("   (ii) Jacobian of (L_2..L_7) at the tautological point: rank %d; rows non-zero for graphs %s (graphs with >= 2 edges vanish to second order: every edge factor is O(delta))" % (rk, [G for G, r in zip(range(2, 8), J) if any(r)]))
    check("tautological point is a SINGULAR point of the solution variety: Jacobian rank 2 (only the one-edge graphs 2, 4), Zariski tangent space dim 38", rk == 2)
    nu_free = [G for G in range(2, 8) if not any(v >= 15 for m in L[G] for v in MONOS3[m])]
    log("   graphs without any mu_2, mu_3 dependence: %s; L_Gamma is cubic in mu_1 and LINEAR in (mu_2, mu_3) jointly" % nu_free)
    ranks = []
    for _ in range(3):
        mu1 = [rng.randrange(p) for _ in range(15)]; pt, rB, cons = fibre_point(edge, mu1, p, rng); ranks.append((rB, cons))
    log("   at random mu_1: rank of the 6 x 25 (mu_2, mu_3)-coefficient matrix and consistency: %s" % ranks)
    check("(mu_2, mu_3)-coefficient matrix has rank 3 (rows 2, 4, 5) and the system is inconsistent at generic mu_1: mu_1 must lie on W = {L_3 = L_6 = L_7 = 0}", all(r == (3, False) for r in ranks))
    ptT, rBT, consT = fibre_point(edge, tpar[:15], p, rng)
    check("fibre over the tautological mu_1: coefficient rank %d, consistent: a 23-dimensional affine space of (mu_2, mu_3) completing mu_1^taut to solutions" % rBT, rBT == 2 and consT and all(poly_eval3(e, ptT, p) == 0 for e in edge))
    log("      Jacobian rank at a random point of that fibre: %d" % rank_fp(jacobian(edge, ptT, p), 40, p))
    t1 = time.time(); Wpts = random_W_points([L[3], L[6], L[7]], p, rng, ntries=8)
    log("   random points of W (random affine 3-planes in A^15, 3 cubics in 3 unknowns by resultants): %d found in %.1f s" % (len(Wpts), time.time() - t1))
    gen = []
    for mu1 in Wpts[:2]:
        pt, rB, cons = fibre_point(edge, mu1, p, rng)
        jr = rank_fp(jacobian(edge, pt, p), 40, p) if pt else None
        jW = rank_fp([r[:15] for r in jacobian([L[3], L[6], L[7]], mu1 + [0] * 25, p)], 15, p)
        gen.append((rB, cons, jr, jW))
    log("      (rank B, consistent, full Jacobian rank, Jacobian rank of (L_3, L_6, L_7) in mu_1): %s" % gen)
    check("generic point of W: W smooth of dim 12 there, fibre dim 22, full Jacobian rank 6: the solution variety has local dimension 34 at generic points", bool(gen) and all(g == (3, True, 6, 3) for g in gen))
    data['gen'] = gen; data['rk_taut'] = rk; data['rBT'] = rBT
    # classes: covariant diagonal
    A2, n2 = class_covariant_diagonal(R); tc = {G: restrict_poly(L[G], A2, p) for G in range(1, 8)}
    log("   COVARIANT DIAGONAL (a_1, a_3) [mu_2 = 0 forced: antisymmetric diagonal]:")
    for G in range(1, 8): log("      graph %d: %s" % (G, tpoly_fmt(tc[G], n2, R)))
    red, cpiv, _ = eliminate_linear([tc[G] for G in range(2, 8)], 1, 0, p)      # pivot graph 2 (a_3 coefficient -5/48)
    unis = [tp_to_uni(e, 0) for e in red if e]; g = uni_gcd_all(unis, p)
    log("      eliminating a_3 with graph 2: %d univariate polynomials in a_1 of degrees %s; gcd = %s" % (len(unis), [len(u) - 1 for u in unis], p_fmt(g, R, "a1")))
    check("COVARIANT DIAGONAL: no solution (gcd 1); colour-summed version of SCALAR_SCAN s.5", len(g) == 1)
    data['cov_gcd_deg'] = len(g) - 1
    # diagonal (non-covariant), Fourier modes
    Am, nm = modes_diagonal(R); tm = {G: restrict_poly(L[G], Am, p) for G in range(1, 8)}
    log("   DIAGONAL M_0 (x_i, y_i; 10 parameters) in Fourier modes X_a = Sum_i zeta^{ai} x_i, Y_a:")
    for G in range(1, 8): log("      graph %d: %s" % (G, tpoly_fmt(tm[G], nm, R)))
    ds = diagonal_solve(tm, R)
    check("DIAGONAL: every graph is c + alpha X_1 + beta Y_3 + delta X_0X_2 + eps X_1^2 + phi X_3X_4 + gamma C(X), C = X_0^2X_3 + 2X_0X_1X_2 + X_0X_4^2 + X_1^3/3 + 2X_1X_3X_4 + X_2^2X_4 + X_2X_3^2 (Z/5 grading; the cubic parts all proportional)", ds['structure_ok'])
    log("      affine in w = (Y_3, X_0X_2, X_3X_4, C) with constant 6 x 4 matrix of rank %d; %d left null vectors give conditions on X_1:" % (ds['rank'], len(ds['left_null'])))
    for c in ds['conditions']: log("         %s = 0" % p_fmt(c, R, "X1"))
    log("      gcd of the conditions: %s; roots in F_p: %s" % (p_fmt(ds['gcd'], R, "X1"), ds['roots']))
    check("DIAGONAL: no solution at all (the two conditions on X_1 are coprime linear polynomials); any (X_0X_2, X_3X_4, C) would be realisable, so this is the only obstruction", ds['rank'] == 4 and len(ds['gcd']) == 1)
    data['diag_conditions'] = ds['conditions']
    # circulant
    A3, n3 = class_circulant(R); tcc = {G: restrict_poly(L[G], A3, p) for G in range(1, 8)}
    log("   CIRCULANT (Z/5-covariant) M_0: 8 parameters c_1(0..2), c_2(1..2), c_3(0..2); variables per graph: %s" % {G: [n3[v] for v in tp_vars(tcc[G])] for G in range(1, 8)})
    tt = [mu[1][0][0], mu[1][0][1], mu[1][0][2], mu[2][0][1], mu[2][0][2], mu[3][0][0], mu[3][0][1], mu[3][0][2]]
    check("tautological point lies in the circulant class (log R~ is circulant) and solves its equations", all(tpoly_eval(tcc[G], tt, p) == (R.of(Fr(1, 1200)) if G == 1 else 0) for G in range(1, 8)))
    eqs = [tcc[3], tcc[6], tcc[7]]; P = tt[:3]
    planes = []
    for _ in range(3):
        a, b, c = ([rng.randrange(p) for _ in range(3)] for _ in range(3)); planes.append(solve_on_plane(eqs, a, b, c, p))
    log("      pure-c_1 equations (graphs 3, 6, 7) on 3 random planes: gcd degrees %s, points %s" % ([q['gcd_degree'] for q in planes], [len(q['points']) for q in planes]))
    Q = next((q['points'][0] for q in planes if q['points']), None)
    if Q is None: check("circulant: a point of the c_1 locus found on a random plane", False); return data
    d = [(Q[j] - P[j]) % p for j in range(3)]; dn = [x * R.inv(d[0]) % p for x in d] if d[0] else d
    check("circulant: the c_1 locus is a LINE ell through the tautological c_1 (gcd degree 1 on every plane; all three cubics vanish identically along P + t d; plane points collinear)",
          all(q['gcd_degree'] == 1 for q in planes) and all(on_line_check(eqs, P, d, p)) and all(all(((X[j] - P[j]) * d[k] - (X[k] - P[k]) * d[j]) % p == 0 for j in range(3) for k in range(3)) for q in planes for X in q['points']))
    log("      ell = c_1^taut + t d, direction (d_0 = 1) in Q(zeta) [sigma_m c(e) = c(me)]: %s;  c_1^taut = %s" % (circ_qzeta(R, dn), circ_qzeta(R, P)))
    cross = [(P[1] * dn[2] - P[2] * dn[1]) % p, (P[2] * dn[0] - P[0] * dn[2]) % p, (P[0] * dn[1] - P[1] * dn[0]) % p]
    log("      ell passes through the origin c_1 = 0: %s" % (cross == [0, 0, 0]))
    inv2 = R.inv(2)
    check("circulant: direction d(j) = (zeta^j + zeta^{3j})/2, i.e. the mu_1-matrix D_{ij} = zeta^{-i} d(j-i) = (zeta^{j-2i} + zeta^{i-2j})/2 (symmetrised rank-2 tensor of the flat vectors zeta^j and zeta^{-2i})",
          dn[1] == (R.z(1) + R.z(3)) * inv2 % p and dn[2] == (R.z(1) + R.z(2)) * inv2 % p)
    fib = []
    for t in [0, rng.randrange(p), rng.randrange(p)]:
        c1 = [(P[j] + t * d[j]) % p for j in range(3)]; A = []; bv = []
        for G in (2, 4, 5):
            tp = tp_eval_partial(tcc[G], {0: c1[0], 1: c1[1], 2: c1[2]}, p); A.append([tp.get((j,), 0) for j in range(3, 8)]); bv.append((-tp.get((), 0)) % p)
        x, N, rB = fp_solve(A, bv, p); fib.append((rB, x is not None))
    log("      fibres (c_2, c_3) over ell: at the tautological point and two random points: (rank of the 3 x 5 system, consistent) = %s" % fib)
    check("circulant: over ell the (c_2, c_3) system is consistent with rank 2 at the tautological point (fibre dim 3) and rank 3 generically (fibre dim 2): a 3-dimensional family of solutions through the tautological point", fib[0] == (2, True) and all(f == (3, True) for f in fib[1:]))
    offs = []
    for _ in range(2):
        e1, e2 = ([rng.randrange(p) for _ in range(3)] for _ in range(2)); rr = residual(eqs, P, d, e1, e2, p)
        offs.append((rr['orders'], rr[2]['result']['identically_zero'], len(rr[2]['result']['elim']) - 1, rr[2]['result']['dd'], [tuple(X) for X in rr[2]['off_line']]))
    log("      blow-up along ell (chart v = lam w, two random coordinate choices): (vanishing orders along ell, eliminant identically zero, eliminant degree, dd pattern, off-line points) = %s" % [(o[0], o[1], o[2], o[3], len(o[4])) for o in offs])
    X = offs[0][4][0] if offs[0][4] else None
    same = X is not None and all(o[4] == offs[0][4] for o in offs)
    check("circulant: exactly one further point X of the c_1 locus off ell (both coordinate choices agree; eliminant squarefree part has two F_p roots: w(X) and the exceptional w = 0)", same and all(o[3] == {1: 2} for o in offs))
    if X is not None:
        log("      X = %s (Q(zeta) via the circulant Galois rule); L_3, L_6, L_7 at X: %s" % (circ_qzeta(R, X), [tp_eval(E, {0: X[0], 1: X[1], 2: X[2]}, p) for E in eqs]))
        A = []; bv = []
        for G in (2, 4, 5):
            tp = tp_eval_partial(tcc[G], {0: X[0], 1: X[1], 2: X[2]}, p); A.append([tp.get((j,), 0) for j in range(3, 8)]); bv.append((-tp.get((), 0)) % p)
        x, N, rB = fp_solve(A, bv, p); log("      fibre over X: rank %d, consistent %s -> a %s-dimensional family of circulant solutions over X" % (rB, x is not None, 5 - rB if x is not None else None))
        e = [(X[j] - P[j]) % p for j in range(3)]
        check("circulant: X - c_1^taut = (2/5) e with e(j) = (1 + zeta^{-j})/2, i.e. the mu_1-matrix (zeta^{-i} + zeta^{-j})/2 x 2/5",
              e[0] == R.of(Fr(2, 5)) and e[1] == R.of(Fr(1, 5)) * (1 + R.z(-1)) % p and e[2] == R.of(Fr(1, 5)) * (1 + R.z(-2)) % p)
        data['X'] = circ_qzeta(R, X); data['X_fibre'] = (rB, x is not None)
    data['ell_dir'] = circ_qzeta(R, dn); data['fib'] = fib
    # per-coloured-graph version (diagonal classes)
    t1 = time.time(); cp, rank, nbad, drop_ok = coloured_polys(fr, graphs, R, rng, 75)
    check("per-coloured frozen fit (115 coloured graphs, 61 monomials of weight <= 3 in x_i, y_i, 75 diagonal M_0): rank %d, %d inconsistent (%.1f s); leafless foreign terms of the vertex graph exactly zeta^{3j}(x_j^3/10 + y_j/48)" % (rank, nbad, time.time() - t1), rank == 61 and nbad == 0 and drop_ok)
    tab, okv, okt = monochrome_table(cp, R)
    check("monochrome coloured graphs depend only on (x_c, y_c) and reproduce the SCALAR_SCAN leading columns (T0_-10, T1_-6, T2_-2, T3_2, T4_2) at colour 0 for all 7 graphs", okv and okt)
    th, okth = theta_constants(cp, R); mb = mixed_structure(cp, R)
    check("mixed theta colourings (20) are non-zero constants; mixed colourings of graphs 4, 5, 7 use only 1, x_c, x_c', x_c x_c'", okth and len(th) == 20 and not mb, mb)
    sp = solve_percolour(cp, R)
    log("      per colour c, monochrome equations in (x_c, y_c): (gcd degree all 7 / edge graphs only, F_p roots): %s" % {c: (len(sp[c]['gcd_all']) - 1, len(sp[c]['gcd_edge']) - 1, sp[c]['roots_all'], sp[c]['roots_edge']) for c in range(5)})
    check("per-coloured DIAGONAL system: inconsistent (mixed theta constants non-zero; and already the monochrome edge equations of each colour have gcd 1)", all(len(sp[c]['gcd_edge']) == 1 for c in range(5)))
    cv = covariant_percoloured(cp, R, True); cv2 = covariant_percoloured(cp, R, False)
    log("      per-coloured COVARIANT system (110 edge equations, %d distinct, a_3 eliminated with graph 2 colour 0): gcd degree %d; with the vertex graph (115, pivot graph 1): gcd degree %d" % (cv['n_distinct'], len(cv['gcd']) - 1, len(cv2['gcd']) - 1))
    check("per-coloured COVARIANT system: gcd 1 (reproduces SCALAR_SCAN s.5)", len(cv['gcd']) == 1 and len(cv2['gcd']) == 1)
    data['percolour'] = {c: (len(sp[c]['gcd_all']) - 1, len(sp[c]['gcd_edge']) - 1) for c in range(5)}
    log("   step 3 in %.1f s" % (time.time() - t0))
    data['L'] = L
    return data
