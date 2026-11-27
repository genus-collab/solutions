# -*- coding: utf-8 -*-
"""genus2_scalar.py -- PLAN_scalar_2026-11-04.md: genus-2 coloured graph terms of the formal quintic in the general constant
covariant scalar frame R = R' phi_i, phi_0 = exp(a1 z + a2 z^3), vertex class phi_i . 1 (fq_phi.py), as polynomials
T_{Gamma,c}(a) = T0 + T1 a1 + T2 a1^2 + T3 a1^3 + T4 a2 (T_k in N), their conifold-frame (Cf) expansions and the solves.
Run:  py -I genus2_scalar.py   (writes logs/genus2_scalar.log).  Two primes (2^89-1 first, then 2^61-1).
Steps: (1) graph sums (gsum_q.genus2_sums, unchanged) at the nodes (0,0), (1,0), (-1,0), (2,0), (0,1); Vandermonde interpolation
to T0..T4 per coloured graph (fq_scalar_lib);  (2) checks: node (0,0) = F5; interpolant at (-1/12, 1/360) = F4, at
(-3/20, -7/2500) = F0 and at a random point = the directly computed phi frame, per coloured term; colour-shift equality of each
T_k; Sum_k T_k = 0 (k >= 1), Sum T0 = F_2 exactly in N;  (3) Cf of T0..T4 for all 115 coloured graphs with the valuation bounds
(val T0 = -10, T1 >= -6, T2 >= -2, T3, T4 >= 2) and exponent classes mod 5 (0, 4, 3, 2, 2);  (4) tables per orbit
representative (23) and totals;  (5) the solves (full, weak);  (6) OPTIONAL q-dependent leading system;  (7) cross-prime equality.
"""
import os, sys, time, itertools, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_rows import Rows
from fq_graphs import Decorated
from fq_conifold import Conifold
from gsum_q import make_frames, genus2_sums, total, to_sseries, cf_images_u, ser_fmt
from fq_phi import phi_frame
from fq_scalar_lib import NODES, MONO_NAMES, vander_inverse, interpolate, eval_poly, ser_classes, ser_val
from fq_scalar_report import coef_table, solve_full, solve_weak, solve_qdep, COLNAMES

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_scalar.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

NU = 10
HODGE_PT = (Fr(-1, 12), Fr(1, 360)); TRUE_PT = (Fr(-3, 20), Fr(-7, 2500))
VAL_MIN = {0: -10, 1: -6, 2: -2, 3: 2, 4: 2}; CLASSES = {0: {0}, 1: {4}, 2: {3}, 3: {2}, 4: {2}}


def sums(frame, deco, label):
    t0 = time.time(); S = genus2_sums(frame, deco)
    log("   %-34s %d coloured graphs summed in %.1f s" % (label, len(S), time.time() - t0)); return S


def run(p, results):
    R = FpRing(p); t0 = time.time()
    log("=" * 98); log("prime p = %d" % p)
    rows = Rows(R, Nz=6); ring = rows.ring
    deco = Decorated(2, hodge=True)
    C = Conifold(R, Nu=NU); imu0 = cf_images_u(C, R, 0)
    log("rows, decorations, conifold periods (u-precision %d) in %.1f s" % (min(x.prec for x in imu0.values()), time.time() - t0))
    # ---- (1) graph sums at the nodes and the reference frames
    log("-" * 98); log("(1) graph sums: nodes (a1, a2) and reference frames")
    fr = make_frames(rows, ['F0', 'F4', 'F5'])
    S = {}
    for node in NODES: S[node] = sums(phi_frame(rows, node[0], node[1], 'phi%s' % (node,)), deco, "node (a1, a2) = (%s, %s):" % node)
    for n in ('F5', 'F4', 'F0'): S[n] = sums(fr[n], deco, "frame %s:" % n)
    rng = random.Random(p % 1000003)
    ar = (rng.randrange(1, p), rng.randrange(1, p))
    log("   random point a1 = %d, a2 = %d (mod p)" % ar)
    S['rand'] = sums(phi_frame(rows, ar[0], ar[1], 'phirand'), deco, "random point:")
    keys = sorted(S['F5'])
    Minv = vander_inverse()
    T = {key: interpolate({node: S[node][key] for node in NODES}, ring, Minv) for key in keys}
    log("   interpolated T0..T4 (monomials %s) for %d coloured graphs" % (MONO_NAMES, len(keys)))
    # ---- (2) checks in N
    log("-" * 98); log("(2) checks exactly in N")
    check("node (0,0) equals the point frame F5 for every coloured term", all(S[NODES[0]][k].eq(S['F5'][k], ring) for k in keys))
    check("interpolant at (-1/12, 1/360) equals frame F4 (Hodge_(zeta^i), R h^-1) for every coloured term",
          all(eval_poly(T[k], HODGE_PT[0], HODGE_PT[1], ring).eq(S['F4'][k], ring) for k in keys))
    check("interpolant at (-3/20, -7/2500) equals frame F0 (true class, R h_wt^-1) for every coloured term",
          all(eval_poly(T[k], TRUE_PT[0], TRUE_PT[1], ring).eq(S['F0'][k], ring) for k in keys))
    check("interpolant at the random point equals the directly computed phi frame for every coloured term (no hidden monomials)",
          all(eval_poly(T[k], ar[0], ar[1], ring).eq(S['rand'][k], ring) for k in keys))
    for k in range(5):
        check("T%d: colour-shift equality T(idx, c+1) = T(idx, c) for every coloured graph" % k,
              all(T[(idx, col)][k].eq(T[(idx, tuple((c + 1) % 5 for c in col))][k], ring) for (idx, col) in keys))
    F2 = total(S['F5'], ring)
    check("Sum over coloured graphs of T0 equals F_2 = total(F5) = total(F0) exactly in N",
          total({k: T[k][0] for k in keys}, ring).eq(F2, ring) and total(S['F0'], ring).eq(F2, ring))
    for k in range(1, 5):
        check("Sum over coloured graphs of T%d = 0 exactly in N" % k, total({kk: T[kk][k] for kk in keys}, ring).iszero())
    # ---- (3) Cf of T0..T4
    log("-" * 98); log("(3) conifold-frame expansions of T0..T4 (L -> 1/s, u = s^5), all 115 coloured graphs")
    t1 = time.time(); cf = {}
    for key in keys: cf[key] = [to_sseries(T[key][k], ring, imu0) for k in range(5)]
    log("   %d substitutions in %.1f s; series precision O(s^%d)" % (5 * len(keys), time.time() - t1, min(s.prec for v in cf.values() for s in v)))
    vals = {key: tuple(ser_val(s) for s in cf[key]) for key in keys}
    check("val T0 = -10 for all 115 coloured graphs", all(vals[k][0] == -10 for k in keys))
    for j in range(1, 5):
        bad = [k for k in keys if vals[k][j] is not None and vals[k][j] < VAL_MIN[j]]
        check("val T%d >= %d for all coloured graphs (zero series allowed)" % (j, VAL_MIN[j]), not bad, bad[:3])
        log("      T%d: zero series for %d of %d coloured graphs; min valuation %s" % (j, sum(1 for k in keys if vals[k][j] is None), len(keys),
            min([vals[k][j] for k in keys if vals[k][j] is not None], default=None)))
    for j in range(5):
        bad = [k for k in keys if not ser_classes(cf[k][j], R) <= CLASSES[j]]
        check("exponent classes mod 5 of T%d in %s (to the series precision)" % (j, sorted(CLASSES[j])), not bad, bad[:3])
    for i in range(5):
        log("   vertex graph, colour %d: T0: %s" % (i, ser_fmt(cf[(1, (i,))][0], R, kmax=0)))
        log("      T1: %s" % ser_fmt(cf[(1, (i,))][1], R, kmax=4)); log("      T2: %s" % ser_fmt(cf[(1, (i,))][2], R, kmax=8))
        log("      T3: %s" % ser_fmt(cf[(1, (i,))][3], R, kmax=7)); log("      T4: %s" % ser_fmt(cf[(1, (i,))][4], R, kmax=7))
    reps = [k for k in keys if k[1][0] == 0]
    check("23 colour-shift orbit representatives (col[0] = 0)", len(reps) == 23, len(reps))
    # ---- (4) tables
    log("-" * 98); log("(4) table per orbit representative: %s  (polar part = T0_-10 s^-10 + T0_-5 s^-5 + a1 (T1_-6 s^-6 + T1_-1 s^-1) + a1^2 T2_-2 s^-2)" % COLNAMES)
    data, tot, nrec, ngal = coef_table(cf, reps, R, log)
    check("Galois back-substitution of the Q(zeta) coordinates reproduces every table entry", ngal == 0, ngal)
    check("all Q(zeta) coordinates of the table entries reconstruct to rationals", nrec == 0, nrec)
    exp_tot = [R.of(Fr(1, 1200)), R.of(Fr(-7, 6000)), 0, 0, 0]
    check("totals row = (1/1200, -7/6000, 0, 0, 0)", tot == exp_tot, [R.fmt(x) for x in tot])
    # ---- (5) solves
    log("-" * 98); log("(5) the solve for a constant covariant scalar (a1, a2)")
    full_ok, nfull = solve_full(cf, R, log)
    check("full system inconsistent: T0_-10 != 0 for all 110 coloured edge graphs", full_ok, nfull)
    cnt, ser_nz = solve_weak(cf, reps, R, log)
    # ---- (6) optional q-dependent leading system
    log("-" * 98); log("(6) OPTIONAL (beyond the user's constant-scalar choice): q-dependent a_1 in L^4 Q[[u]], a_2 in L^12 Q[[u]];")
    log("    leading s^-10 system T0_-10 + al_1 T1_-6 + al_1^2 T2_-2 + al_1^3 T3_2 + al_2 T4_2 = 0 over the 23 representatives")
    qd = solve_qdep(cf, reps, R, log)
    results[p] = dict(data=data, vals=vals, cnt=cnt, ser_nz=ser_nz, nfull=nfull, qd_deg=qd.get('deg'), qd_roots=qd.get('roots_fmt', []),
                      tot=[R.fmt(x) for x in tot])
    log("time for prime %d: %.1f s" % (p, time.time() - t0))


if __name__ == "__main__":
    results = {}; T0 = time.time()
    run(P2, results); run(P1, results)
    log("=" * 98); log("(7) cross-prime comparison (2^89-1 vs 2^61-1)")
    a, b = results[P2], results[P1]
    check("table entries (Q(zeta) coordinates, reconstructed) identical at both primes",
          all(a['data'][k] == b['data'][k] for k in a['data'] if k not in ('gsums', 'tot')))
    check("per-graph colour sums and totals identical at both primes", a['data']['gsums'] == b['data']['gsums'] and a['data']['tot'] == b['data']['tot'])
    check("Cf valuations of T0..T4 identical at both primes", a['vals'] == b['vals'])
    check("weak-system counts and u-series counts identical at both primes", a['cnt'] == b['cnt'] and a['ser_nz'] == b['ser_nz'] and a['nfull'] == b['nfull'])
    check("optional q-dependent system: gcd degree and roots identical at both primes", a['qd_deg'] == b['qd_deg'] and a['qd_roots'] == b['qd_roots'],
          "(%s, %s) vs (%s, %s)" % (a['qd_deg'], a['qd_roots'], b['qd_deg'], b['qd_roots']))
    log("total time %.1f s" % (time.time() - T0))
    log("ALL SCALAR CHECKS PASSED" if not FAILED else "SCALAR FAILURES: %s" % FAILED)
    LOG.close()
