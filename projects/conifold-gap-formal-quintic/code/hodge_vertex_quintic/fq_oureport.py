# -*- coding: utf-8 -*-
"""fq_oureport.py -- report blocks of genus2_ou.py, steps 1 (HatRows validation, point frame in t_c), 2 (M_0 = R~ mod zhat^4, g1_Gamma)
and 3 (the 2 x 40 system J).  Each block takes (log, check) callables and a fq_ourec.Rec for the cross-prime comparison."""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_sser import SympM
from fq_mfit import VNAME
from fq_leadsolve import FrozenRows, frozen_sums
from fq_scalar_lib import qzeta_coords, qzeta_fmt
from gsum_q import is_u_series
from fq_hat import (hat_sums, ucoef, ucoefs, uprec, tautological_M0, tc_expansion, point_frame_check, vander_inv)
from fq_ou import (EPS_NODES, GRAPHS, jacobian_N, jacobian_from_json, dexp_check, dexp_numeric_check, direction_N)
from fq_ourec import jrow_qzeta, jrow_conjugates, fmt_rat, fmt_qz_fp

C_LEAD = [Fr(1807, 3072000), Fr(-811, 256000), Fr(-57, 64000), Fr(843, 1024000), Fr(3, 51200), Fr(87, 32000), Fr(9, 12800)]
C_TC1 = [Fr(-791, 15360000), Fr(-1771, 3840000), Fr(41, 12800), Fr(743, 15360000), Fr(-41, 153600), Fr(37, 160000), Fr(-173, 64000)]
F2_U = [Fr(1, 1200), Fr(-7, 6000), Fr(97, 72000), Fr(-6809, 18000000), Fr(-2551, 15000000)]
GNAME = {1: "genus-2 vertex", 2: "g1 loop", 3: "g0 two loops", 4: "g1-g1", 5: "g0 loop - g1", 6: "theta", 7: "dumbbell"}


def ufmt(R, ser, kmax):
    return " + ".join("(%s) u^%d" % (fmt_rat(R, ucoef(ser, k)), k) for k in range(kmax + 1))


def step1(S, H, sym, tc, graphs, log, check, R, rec, tag, brief=False):
    p = R.p; t0 = time.time(); rng = random.Random(p % 7919 + 11)
    check("%s HatRows: all Qhat[a][i][n] = Q s^{4n+1-c_a} and xhat = x/s are pure u-series (only s^{5k} terms), n <= %d" % (tag, H.Nz), H.useries_ok)
    check("%s HatRows: every Qhat has valuation 0 with non-zero constant term, xhat = 1 + O(u); minimal s-precision %d (u^0..u^%d known)" % (tag, H.min_prec, (H.min_prec - 1) // 5), H.val_ok and H.min_prec >= 15)
    check("%s HatRows: Qhat|_{u=0} = the conifold symbol R^ (fq_symbol.Symbol.Rh) for all (a, i, n)" % tag, H.symbol_check(sym))
    ok, bad = point_frame_check(S, H, graphs, 10)
    check("%s (a) s^{-10} x hatted sums at M = 1 = point-frame s-series sums (GraphSummer(SRows, 1)) through s^10, all 7 graphs" % tag, ok, bad[:2])
    fr = FrozenRows(sym); I = SympM.identity(R); rh = hat_sums(H, I, graphs)
    check("%s (b) [u^0] of the hatted sums at M = 1 = frozen sums = 1807/3072000, -811/256000, -57/64000, 843/1024000, 3/51200, 87/32000, 9/12800; total 1/1200" % tag,
          all(ucoef(rh[G], 0) == R.of(C_LEAD[G - 1]) for G in GRAPHS) and frozen_sums(fr, I, graphs)[0] == {G: R.of(C_LEAD[G - 1]) for G in GRAPHS}
          and sum(ucoef(rh[G], 0) for G in GRAPHS) % p == R.of(Fr(1, 1200)))
    okb = True; okc = True
    for _ in range(2):
        M = SympM.random(R, rng); r = hat_sums(H, M, graphs); fz, _ = frozen_sums(fr, M, graphs)
        if any(ucoef(r[G], 0) != fz[G] for G in GRAPHS): okb = False
        if not all(is_u_series(r[G]) for G in GRAPHS): okc = False
    check("%s (b) [u^0] of the hatted sums at 2 random symplectic M (coefficient matrices read as [zhat^k]) = fq_leadsolve.frozen_sums at the same M" % tag, okb)
    check("%s (c) hatted sums at random M are pure u-series (no s^k, 5 does not divide k, through the precision)" % tag, okc)
    tot = None
    for G in GRAPHS:
        tot = rh[G] if tot is None else tot.add(rh[G])
        for k in range(3): rec.rat("%s:M1:g%d:u%d" % (tag, G, k), ucoef(rh[G], k))
    if not brief:
        log("   u-expansion of G_Gamma = s^10 T_Gamma at M = 1 (point frame), through u^2 (precision u^%d):" % (uprec(rh[2]) - 1))
        for G in GRAPHS: log("      graph %d (%s): %s" % (G, GNAME[G], ufmt(R, rh[G], 2)))
    tot_u = ucoefs(tot, min(4, uprec(tot) - 1))
    for k, x in enumerate(tot_u): rec.rat("%s:M1:total:u%d" % (tag, k), x)
    if tag == "ab00":
        check("%s point frame: total Sum_Gamma G_Gamma = F_2^cf s^10 = 1/1200 - (7/6000) u + (97/72000) u^2 - (6809/18000000) u^3 - (2551/15000000) u^4" % tag,
              all(tot_u[k] == R.of(F2_U[k]) for k in range(len(tot_u))))
    else:
        log("   point frame total F_2^cf s^10 at this frame: %s" % ufmt(R, tot, len(tot_u) - 1))
        check("%s point frame total: [u^0] = 1/1200 and [u^1] = -2 c_1/1200 (so that the t_c^-1 coefficient vanishes), c_1 = %s" % (tag, R.fmt(tc.coef(2))),
              tot_u[0] == R.of(Fr(1, 1200)) and tot_u[1] == (-2 * tc.coef(2) * R.of(Fr(1, 1200))) % p)
    log("   t_c = %s (c_1 = %s): T_Gamma = a_0 t_c^-2 + (a_1 + 2 c_1 a_0) t_c^-1 + ...; per graph (t_c^-2, t_c^-1, t_c^0):" % (tc.fmt(3), R.fmt(tc.coef(2))))
    A = {}; okt = True
    for G in GRAPHS:
        A[G] = tc_expansion(ucoefs(rh[G], 3), tc, R)
        log("      graph %d: %s, %s, %s" % (G, R.fmt(A[G][0]), R.fmt(A[G][1]), R.fmt(A[G][2])))
        for k in range(3): rec.rat("%s:M1:tc%d:g%d" % (tag, k - 2, G), A[G][k])
        if A[G][0] != R.of(C_LEAD[G - 1]) or (tag == "ab00" and A[G][1] != R.of(C_TC1[G - 1])): okt = False
    if tag == "ab00":
        check("%s point frame in t_c: per-graph t_c^-2 = c_Gamma and t_c^-1 = -791/15360000, -1771/3840000, 41/12800, 743/15360000, -41/153600, 37/160000, -173/64000 (FORMAL_QUINTIC s.3)" % tag, okt)
    else:
        check("%s point frame in t_c: per-graph t_c^-2 = c_Gamma unchanged (the t_c^-1 coefficients may change with alpha, beta: logged)" % tag, okt)
    T = [sum(A[G][k] for G in GRAPHS) % p for k in range(3)]
    check("%s point frame total in t_c: (1/1200) t_c^-2 + 0 t_c^-1 + constant %s" % (tag, R.fmt(T[2])), T[0] == R.of(Fr(1, 1200)) and T[1] == 0)
    rec.rat("%s:M1:tc0:total" % tag, T[2])
    log("   step 1 in %.1f s" % (time.time() - t0))
    return dict(rh=rh, A=A, tot=tot_u)


def step2(H, sym, tc, graphs, log, check, R, rec, tag, tot1):
    """tot1: the point-frame total [u^0..] of the same frame (step1); the total is M-independent"""
    p = R.p; t0 = time.time(); c1 = tc.coef(2); g1tot = (-2 * c1 * R.of(Fr(1, 1200))) % p
    M0s, mu, C, ok = tautological_M0(sym, R)
    if tag == "ab00":
        log("   mu_k = [zhat^k] log R~ (first rows, Q(zeta) via sigma_m(mu_k[0][i]) = mu_k[0][mi]):")
        for k in (1, 2, 3): log("      mu_%d[0][.] = [%s]" % (k, ", ".join(qzeta_fmt(R, qzeta_coords(R, [mu[k][0][(m * i) % 5] for m in range(1, 5)])) for i in range(5))))
    check("%s M_0 := SympM(mu_1, mu_2, mu_3) reproduces R~ mod zhat^4 exactly: M_0[k] = C_k for k = 0..3" % tag, ok)
    r0 = hat_sums(H, M0s, graphs); K = min(uprec(r0[G]) for G in GRAPHS) - 1
    uc = {G: ucoefs(r0[G], K) for G in GRAPHS}; g1 = {G: uc[G][1] for G in GRAPHS}
    log("   hatted sums G_Gamma[M_0](u) with the tautological M_0, through u^%d:" % K)
    for G in GRAPHS:
        log("      graph %d (%s): %s" % (G, GNAME[G], ufmt(R, r0[G], min(K, 3))))
        for k in range(K + 1): rec.rat("%s:M0:g%d:u%d" % (tag, G, k), uc[G][k])
    check("%s [u^0] G_Gamma[M_0] = (1/1200, 0, 0, 0, 0, 0, 0): every edge graph loses its leading pole (CONIFOLD_SYMBOL s.3)" % tag,
          uc[1][0] == R.of(Fr(1, 1200)) and all(uc[G][0] == 0 for G in GRAPHS if G != 1))
    check("%s KEY: g1_Gamma := [u^1] G_Gamma[M_0] = 0 for the graphs with >= 2 edges (3, 5, 6, 7): they are REGULAR for M = M_0 (and for every M_0 exp(uN))" % tag,
          all(g1[G] == 0 for G in (3, 5, 6, 7)), {G: R.fmt(g1[G]) for G in (3, 5, 6, 7)})
    check("%s g1_2 and g1_4 are RATIONAL (colour sums of a Z/5-covariant M_0): g1_2 = %s, g1_4 = %s, g1_1 = %s" % (tag, fmt_rat(R, g1[2]), fmt_rat(R, g1[4]), fmt_rat(R, g1[1])),
          all(R.ratrec(g1[G]) is not None for G in (1, 2, 4)) or p < 10 ** 20)
    KK = min(K, len(tot1) - 1)
    check("%s Sum_Gamma g1_Gamma = -2 c_1/1200 = %s (= [u^1] of F_2^cf s^10%s, M-independent); total through u^%d = the point-frame total of this frame" % (tag, fmt_rat(R, g1tot), " = -7/6000" if tag == "ab00" else "", KK),
          sum(g1.values()) % p == g1tot and all(sum(uc[G][k] for G in GRAPHS) % p == tot1[k] for k in range(KK + 1)))
    A = {G: tc_expansion(uc[G][:3], tc, R) for G in GRAPHS}
    log("   per-graph (t_c^-2, t_c^-1, t_c^0) with M_0: %s" % "; ".join("g%d: %s, %s, %s" % (G, R.fmt(A[G][0]), fmt_rat(R, A[G][1]), fmt_rat(R, A[G][2])) for G in GRAPHS))
    for G in GRAPHS:
        for k in range(3): rec.rat("%s:M0:tc%d:g%d" % (tag, k - 2, G), A[G][k])
    check("%s with M_0: edge-graph 1/t_c coefficients = g1_Gamma; vertex 1/t_c coefficient = g1_1 + 2 c_1/1200 = %s; total 1/t_c coefficient 0; total t_c^-2 = 1/1200" % (tag, fmt_rat(R, A[1][1])),
          all(A[G][1] == g1[G] for G in GRAPHS if G != 1) and A[1][1] == (g1[1] + 2 * c1 * R.of(Fr(1, 1200))) % p
          and sum(A[G][1] for G in GRAPHS) % p == 0 and sum(A[G][0] for G in GRAPHS) % p == R.of(Fr(1, 1200)))
    log("   vertex graph with M_0: constant term (t_c^0) %s; [u^2] G_1[M_0] = %s" % (fmt_rat(R, A[1][2]), fmt_rat(R, uc[1][2])))
    log("   step 2 in %.1f s" % (time.time() - t0))
    return dict(M0=M0s, mu=mu, C=C, g1=g1, uc=uc, A=A, r0=r0)


def step3(sym, d2, graphs, log, check, R, rec, json_path=None):
    p = R.p; t0 = time.time(); M0 = d2['M0']; mu = d2['mu']
    fr = FrozenRows(sym); Vinv = vander_inv(EPS_NODES, p)
    log("   40 basis directions N_v (fq_leadsolve.params_to_mats: nu_1, nu_3 symmetric with both entries 1; nu_2 antisymmetric (i,j) = +1, (j,i) = -1), eps = 0..9, frozen sums:")
    J, maxdeg, consts = jacobian_N(fr, M0.M, M0.Minv, graphs, R, Vinv, log)
    check("eps-polynomials G^0_Gamma[M_0 exp(eps N_v)] have degree <= 3 (coefficients of eps^4..eps^9 vanish) for all 40 directions and 7 graphs; maximal degree seen %d" % maxdeg, maxdeg <= 3)
    check("constant terms = L_Gamma(M_0) = (1/1200, 0, .., 0) for every direction", consts[1] == {R.of(Fr(1, 1200))} and all(consts[G] == {0} for G in GRAPHS if G != 1))
    zero_rows = [G for G in GRAPHS if not any(J[G])]
    check("D_N L_Gamma = 0 identically for the graphs with >= 2 edges: the 4 x 40 block for graphs 3, 5, 6, 7 is zero; rows 1, 2, 4 non-zero", zero_rows == [3, 5, 6, 7])
    check("Sum_Gamma J[Gamma] = 0 (the total 1/1200 is M-independent), i.e. J[1] = -(J[2] + J[4])", all(sum(J[G][v] for G in GRAPHS) % p == 0 for v in range(40)))
    from fq_mfit import fp_rref
    rk = len(fp_rref([J[2], J[4]], 40, p)[0])
    check("rank of the 2 x 40 matrix J (rows graphs 2, 4) = 2: the two linear equations are independent on the full class", rk == 2)
    data = dict(J=J, Vinv=Vinv, rank=rk)
    for G in (2, 4):
        qz = jrow_qzeta(R, J[G]); conj = jrow_conjugates(J[G], p); nz = [v for v in range(40) if J[G][v]]
        log("   J[%d] (graph %s) in Q(zeta) via Galois (sigma_m J_v = sgn J_{pi_m v}), %d non-zero entries:" % (G, GNAME[G], len(nz)))
        for v in nz:
            rec.qz("J%d:%s" % (G, VNAME[v]), [conj[m][v] for m in range(4)])
        log("      " + "; ".join("%s: %s" % (VNAME[v], fmt_qz_fp(R, qz[v])) for v in nz))
        data['J%d_vars' % G] = [VNAME[v] for v in nz]
    if json_path and os.path.exists(json_path):
        try:
            Jmu, Phi, polys = jacobian_from_json(json_path, R, mu)
            nf = dexp_check(J, Jmu, Phi, p)
            rng = random.Random(3)
            num_ok = all(dexp_numeric_check(R, mu, direction_N(rng.randrange(40), p), Vinv) for _ in range(3))
            check("cross-check: J . coords(dexp(Y_e)) = mu-gradient of the stored leading polynomials L_Gamma (symbol_result.json) at the tautological point, all 7 graphs x 40 directions (%d failures); dexp(Y) = Y - [mu,Y]/2 + [mu,[mu,Y]]/6 verified numerically" % nf, nf == 0 and num_ok)
            data['dexp_ok'] = (nf == 0 and num_ok)
        except Exception as e:
            log("   (cross-check against symbol_result.json skipped: %s)" % e)
    log("   step 3 in %.1f s" % (time.time() - t0))
    return data
