# -*- coding: utf-8 -*-
"""genus2_diag.py -- formal quintic, genus 2: the DIAGONAL, COLOUR-DEPENDENT constant scalar factorisation
R = R' diag(phi_i(z)),  phi_i(z) = exp(a1_i z + a2_i z^3)  (ten independent constants, i = 0..4), vertex CohFT phi_i . 1.
Every coloured genus-2 graph term is a polynomial of a-weight <= 3 in the a's of its colours (fq_diag_lib: 5 or 12 monomials),
obtained by evaluating the unchanged graph sum (gsum_q.genus2_sums) at NRAND random points of F_p^10 and solving for the
monomial coefficients (first m points; the remaining points and the special points are consistency checks).
Run:  py -I genus2_diag.py   (writes logs/genus2_diag.log; both primes 2^89-1 and 2^61-1).  A failing check STOPS the run.
Steps per prime: (1) graph sums; (2) coefficient solve; (3) checks in N (zero point = F5, covariant Hodge/true patterns = F4/F0,
covariant random point = fq_phi.phi_frame, colour-shift phases, colour-summed total a-independent); (4) exact exponent classes,
Cf of every coefficient, valuation bounds; (5) tables; (6) the covariant scan recovered; (7) solves: FULL, WEAK, u-series.
Then (8) cross-prime comparison by residues.
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_graphs import Decorated
from fq_conifold import Conifold
from gsum_q import make_frames, genus2_sums, total, cf_images_u
from fq_phi import phi_frame
from fq_diag_lib import (diag_frame, covariant_point, zero_point, random_point, local_vars, local_monos, mono_weight, mono_global,
                         gmono_name, gmono_weight, solve_coeffs, eval_local, GNAMES)
from fq_diag_report import (cf_all, exact_classes, val_checks, zero_pattern, table_reps, table_all, leading_w3, parse_prev_log, compare_prev)
from fq_diag_solve import solve_full, solve_weak, solve_useries
from fq_scalar_lib import qzeta_ratrec

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_diag.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
class StopRun(Exception): pass
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: raise StopRun(name)

NRAND = 15; NU = 10
HODGE_PT = (Fr(-1, 12), Fr(1, 360)); TRUE_PT = (Fr(-3, 20), Fr(-7, 2500))


def sums(frame, deco, label):
    t0 = time.time(); S = genus2_sums(frame, deco)
    log("   %-40s %d coloured graphs summed in %.1f s" % (label, len(S), time.time() - t0)); return S


def run(p, results):
    R = FpRing(p); t0 = time.time()
    log("=" * 98); log("prime p = %d" % p)
    rows = Rows(R, Nz=6); ring = rows.ring; deco = Decorated(2, hodge=True)
    C = Conifold(R, Nu=NU); imu0 = cf_images_u(C, R, 0)
    log("rows, decorations, conifold periods (u-precision %d) in %.1f s" % (min(x.prec for x in imu0.values()), time.time() - t0))
    # ---- wrapper check: the diagonal frame at the covariant point coincides with fq_phi.phi_frame (rows and vertex tables)
    a1r, a2r = 7, 11
    Fd = diag_frame(rows, covariant_point(R, a1r, a2r), 'dcheck'); Fc = phi_frame(rows, a1r, a2r, 'ccheck')
    ok = all(Fd.P[a][i][n].eq(Fc.P[a][i][n], ring) for a in range(5) for i in range(5) if Fd.P[a][i] is not None for n in range(rows.Nz + 1))
    ds = [(g, d) for g in (1, 2) for d in [(), (0,), (1,), (2,), (0, 0), (1, 1), (0, 1), (4,), (2, 1), (0, 0, 0)]]
    ok = ok and all(Fd.vt[i].val(g, d) == Fc.vt[i].val(g, d) for i in range(5) for (g, d) in ds)
    check("diag_frame at a1_i = zeta^-i a1, a2_i = zeta^-3i a2 has the rows and vertex tables of fq_phi.phi_frame (independent per-colour values accepted)", ok)
    # ---- (1) graph sums
    log("-" * 98); log("(1) graph sums: frames F5, F4, F0; the zero point; %d random points of F_p^10; a covariant random point" % NRAND)
    fr = make_frames(rows, ['F0', 'F4', 'F5'])
    S = {n: sums(fr[n], deco, "frame %s:" % n) for n in ('F5', 'F4', 'F0')}
    rng = random.Random(p % 1000003 + 1)
    pts = [random_point(rng, p) for _ in range(NRAND)]
    log("   random points (seed p mod 1000003 + 1): first point %s" % (pts[0],))
    SP = [sums(diag_frame(rows, pt, 'd%d' % j), deco, "random point %d:" % j) for j, pt in enumerate(pts)]
    S0 = sums(diag_frame(rows, zero_point(), 'd0'), deco, "zero point:")
    acov = (rng.randrange(1, p), rng.randrange(1, p))
    Scov = sums(phi_frame(rows, acov[0], acov[1], 'cov'), deco, "fq_phi covariant frame at random (a1, a2):")
    keys = sorted(S['F5'])
    # ---- (2) the coefficient solve
    log("-" * 98); log("(2) solve for the monomial coefficients per coloured graph (first m points), consistency on the remaining points")
    t1 = time.time(); T = {}; ms = {}
    for key in keys: T[key], ms[key] = solve_coeffs([Sj[key] for Sj in SP], pts, key[1], ring)
    log("   %d coloured graphs: %d with 5 monomials, %d with 12; %.1f s" % (len(keys), sum(1 for k in keys if ms[k] == 5), sum(1 for k in keys if ms[k] == 12), time.time() - t1))
    bad = [(key, j) for key in keys for j in range(ms[key], NRAND) if not eval_local(T[key], key[1], pts[j], ring).eq(SP[j][key], ring)]
    check("the interpolant reproduces every unused random point exactly in N (>= 3 extra points per coloured graph: no missing monomial)", not bad, bad[:3])
    # ---- (3) checks in N
    log("-" * 98); log("(3) checks exactly in N")
    check("zero point: directly computed diagonal frame = point frame F5 for every coloured term", all(S0[k].eq(S['F5'][k], ring) for k in keys))
    c0 = {k: next(m for m in T[k] if sum(m) == 0) for k in keys}
    check("constant coefficient (node 0 of the interpolant) = F5 for every coloured term", all(T[k][c0[k]].eq(S['F5'][k], ring) for k in keys))
    cpH = covariant_point(R, *HODGE_PT); cpT = covariant_point(R, *TRUE_PT); cpR = covariant_point(R, *acov)
    check("interpolant at the covariant Hodge pattern a1_i = -zeta^-i/12, a2_i = zeta^-3i/360 = frame F4 for every coloured term",
          all(eval_local(T[k], k[1], cpH, ring).eq(S['F4'][k], ring) for k in keys))
    check("interpolant at the covariant true-class pattern a1_i = -3 zeta^-i/20, a2_i = -7 zeta^-3i/2500 = frame F0 for every coloured term",
          all(eval_local(T[k], k[1], cpT, ring).eq(S['F0'][k], ring) for k in keys))
    check("interpolant at the covariant random point = the fq_phi covariant frame (previous scan) for every coloured term",
          all(eval_local(T[k], k[1], cpR, ring).eq(Scov[k], ring) for k in keys))
    phase_ok = {}
    for sgn in (1, -1):
        phase_ok[sgn] = all(T[(idx, tuple((c + 1) % 5 for c in col))][m].eq(T[(idx, col)][m].scale(R.z(sgn * mono_weight(m, local_vars(col)))), ring)
                            for (idx, col) in keys for m in T[(idx, col)])
    log("   colour shift col -> col + 1: coefficient of a local monomial of a-weight W multiplied by zeta^W: %s; by zeta^-W: %s" % (phase_ok[1], phase_ok[-1]))
    check("colour-shift covariance coef_{col+1}(mono) = zeta^(W) coef_col(mono) for every coloured graph and monomial", phase_ok[1])
    G = {}
    for key in keys:
        vs = local_vars(key[1])
        for m, el in T[key].items():
            g = mono_global(m, vs); G[g] = el if g not in G else G[g].add(el, ring)
    F2 = total(S['F5'], ring)
    check("colour-summed total: constant part = F_2 = total(F5) = total(F0) exactly in N", G[(0,) * 10].eq(F2, ring) and total(S['F0'], ring).eq(F2, ring))
    nonconst = [g for g in G if sum(g)]
    check("colour-summed total independent of all ten parameters: the sum over coloured graphs of the coefficient of each of the %d non-constant global monomials is 0 exactly in N" % len(nonconst),
          all(G[g].iszero() for g in nonconst), [gmono_name(g) for g in nonconst if not G[g].iszero()][:5])
    log("   global monomials present: %d (weights %s)" % (len(G), sorted(set(gmono_weight(g) for g in G))))
    # ---- (4) exact classes, Cf, valuations, zero pattern
    log("-" * 98); log("(4) exact exponent classes from the L-powers; Cf (L -> 1/s, u = s^5) of every coefficient; valuation bounds")
    exact_classes(T, keys, log, check)
    t1 = time.time(); cf = cf_all(T, ring, imu0)
    log("   %d Cf substitutions in %.1f s; precision O(s^%d)" % (sum(len(v) for v in cf.values()), time.time() - t1, min(s.prec for v in cf.values() for s in v.values())))
    val_checks(cf, keys, R, log, check)
    pat = zero_pattern(T, keys, log)
    reps = [k for k in keys if k[1][0] == 0]
    check("23 colour-shift orbit representatives", len(reps) == 23, len(reps))
    # ---- (5) tables
    log("-" * 98); log("(5) tables of the polar coefficients (polar part = const_-10 s^-10 + const_-5 s^-5 + Sum_c a1_c (.. s^-6 + .. s^-1) + Sum_{c<=c'} a1_c a1_c' (.. s^-2))")
    data, nrec = table_reps(cf, reps, R, log, check)
    leading_w3(cf, reps, R, log)
    table_all(cf, keys, R, log)
    # ---- (6) the covariant scan recovered
    log("-" * 98); log("(6) comparison with the covariant scan genus2_scalar.log at this prime (leading coefficients of the 23 representatives)")
    compare_prev(cf, reps, R, parse_prev_log(os.path.join(HERE, "logs", "genus2_scalar.log"), p, R), log, check)
    # ---- (7) solves
    log("-" * 98); log("(7a) FULL polar system of the 110 coloured edge graphs")
    nfull = solve_full(cf, keys, R, log, check)
    log("-" * 98); log("(7b) WEAK system: the a-dependent polar coefficients only")
    weak = solve_weak(cf, keys, R, log, check)
    log("-" * 98); log("(7c) u-series statement (all orders to the series precision)")
    solve_useries(cf, keys, R, log, check)
    results[p] = dict(R=R, data=data, nrec=nrec, weak={k: v for k, v in weak.items()}, pat=pat, nfull=nfull,
                      vals=tuple((key, tuple((m, (S.v if S.c else None)) for m, S in sorted(cf[key].items()))) for key in keys))
    log("time for prime %d: %.1f s" % (p, time.time() - t0))


def cross_prime(results):
    log("=" * 98); log("(8) cross-prime comparison by residues (2^89-1 Q(zeta) coordinates reduced mod 2^61-1 vs the 2^61-1 residues)")
    A, B = results[P2], results[P1]; RA, RB = A['R'], B['R']
    bad = []; nrec = 0
    for k, (x, cs, rs) in A['data'].items():
        if rs is None: nrec += 1; continue
        val = sum(RB.mul(RB.of(r), RB.z(j)) for j, r in enumerate(rs)) % RB.p
        if val != B['data'][k][0]: bad.append(k)
    check("all %d table entries have Q(zeta) coordinates reconstructible at 2^89-1" % len(A['data']), nrec == 0, nrec)
    check("table entries: 2^89-1 coordinates re-assembled mod 2^61-1 equal the 2^61-1 residues", not bad, bad[:5])
    check("zero patterns, Cf valuations, full-system counts identical at both primes", A['pat'] == B['pat'] and A['vals'] == B['vals'] and A['nfull'] == B['nfull'])
    check("weak system: rank, pivots, kernel dimension and solution set identical at both primes", A['weak'] == B['weak'], (A['weak'], B['weak']))


if __name__ == "__main__":
    results = {}; T0 = time.time()
    try:
        run(P2, results); run(P1, results); cross_prime(results)
        log("total time %.1f s" % (time.time() - T0)); log("ALL DIAG CHECKS PASSED")
    except StopRun as e:
        log("total time %.1f s" % (time.time() - T0)); log("STOPPED AT FAILING CHECK: %s" % e)
    LOG.close()
