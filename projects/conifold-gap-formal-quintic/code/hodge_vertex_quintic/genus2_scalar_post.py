# -*- coding: utf-8 -*-
"""genus2_scalar_post.py -- post-run checks for the scalar scan (PLAN_scalar_2026-11-04.md; referee items of the write-up
SCALAR_SCAN_genus2.md).  Run:  py -I genus2_scalar_post.py   (writes logs/genus2_scalar_post.log).  Both primes.
No existing module is modified; genus2_scalar.py and its log stay the record of the scan itself.

(A) EXACT exponent classes and valuation bounds: the Cf image of a numerator key with L-power m is s^-m times a u-series
    (c_j, i_j, 1/c_0, 1/W -> u-series; L -> 1/s), so the exponents of Cf(T_k) lie in {-m mod 5 : m an L-power of T_k} and
    val Cf(T_k) >= -maxLdeg(T_k).  Checked on the interpolated NEls T0..T4 of all 115 coloured graphs (all orders, no series).
(B) Cross-prime comparison of the 23 x 5 table by RESIDUES: the Q(zeta) coordinates reconstructed at 2^89-1 (all 115 entries
    reconstruct there), reduced mod 2^61-1 and re-assembled with zeta mod 2^61-1, must equal the 2^61-1 entries.  This replaces
    the 'reconstructed' comparison of genus2_scalar.py, which cannot reconstruct denominators 15360000 at 2^61-1 (bound 10737418).
(C) The optional q-dependent leading system: gcd over F_p after eliminating al_2 with the vertex graph, for the sub-systems
    (all 22; without the four theta mixed colourings; monochrome only; non-theta mixed only) and the pairwise common factors.
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_graphs import Decorated
from fq_conifold import Conifold
from gsum_q import genus2_sums, to_sseries, cf_images_u
from fq_phi import phi_frame
from fq_scalar_lib import (NODES, vander_inverse, interpolate, ser_coef, qzeta_coords, qzeta_ratrec, galois_colouring,
                           p_gcd, p_sub, p_scale, p_trim, p_monic, p_fmt)
from fq_scalar_report import COLS, COLNAMES

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_scalar_post.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

CLASSES = {0: {0}, 1: {4}, 2: {3}, 3: {2}, 4: {2}}
MAXL = {0: 10, 1: 6, 2: 2, 3: -2, 4: -2}                       # expected max L-degree of T_k  (val >= -maxL)
LEAD = [(0, -10), (1, -6), (2, -2), (3, 2), (4, 2)]              # the leading system's coefficients


def run(p, store):
    R = FpRing(p); t0 = time.time()
    log("=" * 98); log("prime p = %d" % p)
    rows = Rows(R, Nz=6); ring = rows.ring; deco = Decorated(2, hodge=True)
    C = Conifold(R, Nu=10); imu0 = cf_images_u(C, R, 0)
    S = {}
    for node in NODES:
        t1 = time.time(); S[node] = genus2_sums(phi_frame(rows, node[0], node[1], 'phi%s' % (node,)), deco)
        log("   node (a1, a2) = (%s, %s): 115 coloured graphs summed in %.1f s" % (node[0], node[1], time.time() - t1))
    keys = sorted(S[NODES[0]]); Minv = vander_inverse()
    T = {key: interpolate({node: S[node][key] for node in NODES}, ring, Minv) for key in keys}
    # ---- (A) exact exponent classes and valuation bounds from the L-powers of the numerator keys
    log("-" * 98); log("(A) exact exponent classes and L-degree bounds of the NEls T0..T4 (no series truncation)")
    check("all T_k are Lam-free (no log q)", all(not T[key][k].has('Lam') for key in keys for k in range(5)))
    for k in range(5):
        cls = set(); maxL = None; minL = None
        for key in keys:
            el = T[key][k]
            if el.iszero(): continue
            cls |= set((-m) % 5 for (m, *_rest) in el.num)
            lo, hi = el.Ldeg(); maxL = hi if maxL is None else max(maxL, hi); minL = lo if minL is None else min(minL, lo)
        check("T%d: exponent classes {-m mod 5 : m an L-power} = %s (exact, all orders)" % (k, sorted(CLASSES[k])), cls == CLASSES[k], sorted(cls))
        check("T%d: max L-degree over the 115 coloured graphs = %d (so val Cf(T%d) >= %d exactly)" % (k, MAXL[k], k, -MAXL[k]), maxL == MAXL[k],
              "L-degree range (%s, %s)" % (minL, maxL))
    # ---- Cf images and the table entries
    log("-" * 98); log("(B) table entries (23 representatives x 5 columns) and their Q(zeta) coordinates")
    t1 = time.time(); cf = {key: [to_sseries(T[key][k], ring, imu0) for k in range(5)] for key in keys}
    log("   %d Cf substitutions in %.1f s" % (5 * len(keys), time.time() - t1))
    reps = [key for key in keys if key[1][0] == 0]
    ent = {}; coords = {}
    for key in reps:
        idx, col = key
        for (k, e) in COLS:
            ent[(key, k, e)] = ser_coef(cf[key][k], e)
            xs = [ser_coef(cf[(idx, galois_colouring(col, m))][k], e) for m in range(1, 5)]
            coords[(key, k, e)] = qzeta_coords(R, xs)
    lead = {key: [ser_coef(cf[key][k], e) for (k, e) in LEAD] for key in reps}
    store[p] = dict(R=R, ent=ent, coords=coords, lead=lead, reps=reps)
    # ---- (C) optional leading system: gcds of sub-systems and pairwise common factors
    log("-" * 98); log("(C) optional q-dependent leading system: sub-system gcds after eliminating al_2 with the vertex graph")
    cP = p_trim(list(lead[(1, (0,))][:4])); EP = lead[(1, (0,))][4]; Einv = pow(EP, p - 2, p)
    polys = {}
    for key in reps:
        if key == (1, (0,)): continue
        c = p_trim(list(lead[key][:4])); E = lead[key][4]
        polys[key] = p_sub(c, p_scale(cP, E * Einv % p, p), p)
    theta_mixed = [key for key in polys if key[0] == 6 and len(set(key[1])) > 1]
    mono = [key for key in polys if len(set(key[1])) == 1]
    mixed_nontheta = [key for key in polys if len(set(key[1])) > 1 and key[0] != 6]
    def gcd_of(sub):
        g = None
        for key in sub:
            f = polys[key]
            if not f: continue
            g = f if g is None else p_gcd(g, f, p)
        return p_monic(g, p) if g else g
    for label, sub in [("all 22 equations", list(polys)), ("without the 4 theta mixed colourings", [k for k in polys if k not in theta_mixed]),
                       ("monochrome equations only (graphs 2-7, colour 0)", mono), ("non-theta mixed colourings only", mixed_nontheta)]:
        g = gcd_of(sub)
        degs = sorted(set(len(polys[k]) - 1 for k in sub))
        log("   %-48s %2d equations, degrees %s, gcd degree %d: %s" % (label, len(sub), degs, len(g) - 1, p_fmt(g, R, "al_1")))
        check("gcd = 1 for the sub-system: %s" % label, len(g) == 1)
    check("the theta mixed colourings' equations are non-zero constants (T0_-10 alone)",
          all(len(polys[k]) == 1 for k in theta_mixed), [len(polys[k]) - 1 for k in theta_mixed])
    pairs = []
    pk = sorted(polys)
    for a, b in itertools.combinations(pk, 2):
        if not polys[a] or not polys[b]: continue
        g = p_gcd(polys[a], polys[b], p)
        if len(g) > 1: pairs.append((a, b, len(g) - 1))
    log("   pairs of equations with a common factor (degree): %d pairs; degree-3 pairs (proportional cubics): %s; degree-2 pairs: %s"
        % (len(pairs), [(a, b) for a, b, d in pairs if d == 3], [(a, b) for a, b, d in pairs if d == 2]))
    store[p]['pairs'] = [(a, b, d) for a, b, d in pairs]
    log("time for prime %d: %.1f s" % (p, time.time() - t0))


def cross_prime(store):
    log("=" * 98); log("(B) cross-prime comparison by residues (2^89-1 rationals reduced mod 2^61-1 vs the 2^61-1 residues)")
    A, B = store[P2], store[P1]; RA, RB = A['R'], B['R']
    nrec = 0; bad = []; bad_coords = []
    for key in A['reps']:
        for (k, e) in COLS:
            rs = qzeta_ratrec(RA, A['coords'][(key, k, e)])
            if rs is None: nrec += 1; continue
            # re-assemble at the small prime and compare with the entry there
            val = sum(RB.mul(RB.of(r), RB.z(j)) for j, r in enumerate(rs)) % RB.p
            if val != B['ent'][(key, k, e)]: bad.append((key, k, e))
            if tuple(RB.of(r) for r in rs) != tuple(B['coords'][(key, k, e)]): bad_coords.append((key, k, e))
    check("all 115 table entries have Q(zeta) coordinates reconstructible at 2^89-1", nrec == 0, nrec)
    check("115 table entries: 2^89-1 Q(zeta) coordinates re-assembled mod 2^61-1 equal the 2^61-1 entries", not bad, bad[:5])
    check("115 table entries: 2^89-1 coordinates reduced mod 2^61-1 equal the 2^61-1 coordinate residues", not bad_coords, bad_coords[:5])
    for key in [(1, (0,)), (4, (0, 0))]:
        rs = qzeta_ratrec(RA, A['coords'][(key, 0, -10)])
        log("   %s T0_-10 = %s at 2^89-1; mod 2^61-1: %d; 2^61-1 residue: %d" % (str(key), rs[0], RB.of(rs[0]), B['ent'][(key, 0, -10)]))
    # leading-system coefficients (T3_2, T4_2 included) and the pair structure
    badl = []
    for key in A['reps']:
        for j in range(5):
            x = RA.ratrec(A['lead'][key][j])
            if x is not None and RB.of(x) != B['lead'][key][j] and RA.fmt(A['lead'][key][j]) == str(x): badl.append((key, j))
    check("rational leading-system coefficients identical at both primes", not badl, badl[:5])
    check("pairwise common-factor structure of the optional system identical at both primes", A['pairs'] == B['pairs'])


if __name__ == "__main__":
    store = {}; T0 = time.time()
    run(P2, store); run(P1, store); cross_prime(store)
    log("total time %.1f s" % (time.time() - T0))
    log("ALL POST CHECKS PASSED" if not FAILED else "POST FAILURES: %s" % FAILED)
    LOG.close()
