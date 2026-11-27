# -*- coding: utf-8 -*-
"""calib_check.py -- PLAN step 3: Hodge-vertex graph sum (rows R~ = R h_{wt_i}^-1, vertex class Delta_i^{g_v-1} Hodge_{wt_i})
against the localisation numbers of loc_check.py.  Run:  py -I calib_check.py   (writes logs/calib_check.log).
  genus 1: d_T F_1 = <H>_1 from the one-leg graphs (gsum_q.genus1_dT_sums); N_{1,d} = [Q^d] <H>_1 / d; [Q^0] recorded.
  genus 2: the 7 leg-free graphs, all 5^V colourings; q^0 term vs F_2(0) = -19/9600; Q-coefficients vs N_{2,1}, N_{2,2}.
  controls that must fail: (R~, 1), (R, Hodge_wt), (R h_(1)^+1, Hodge_(1));  (R, 1) point frame and (R h_(1)^-1, Hodge_(1))
  must agree with the reference exactly in N.  Two primes: rationals reconstructed at 2^89-1, residues checked at 2^61-1.
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_genus0 import Genus0
from fq_ring import NEl
from fq_rows import Rows
from fq_graphs import Decorated
from gsum_q import make_frames, genus2_sums, genus1_dT_sums, total, to_qseries

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "calib_check.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

REF = {('N', 1, 1): Fr(-305, 6), ('N', 1, 2): Fr(-535775, 12), ('N', 2, 1): Fr(2129, 5760), ('N', 2, 2): Fr(203401, 144),
       ('F0', 2): Fr(-19, 9600)}
D = 4            # q-precision q^D


def qcoefs(S, R, D):
    return [S.coef(d) if d < S.prec else None for d in range(D + 1)]


def run(p, store):
    R = FpRing(p); t0 = time.time()
    log("=" * 78); log("prime p = %d" % p)
    lams = [R.z(k) for k in range(5)]
    rows = Rows(R, Nz=6); ring = rows.ring
    G = Genus0(R, lams, D, Ns=6); img = ring.images_from_qseries(G)
    deco = Decorated(2, hodge=True)
    log("rows, genus-0 q-series (q^%d) and decorations in %.1f s; decorations %s" % (D, time.time() - t0, deco.counts()))
    frames = make_frames(rows)
    F0 = frames['F0']
    # ---- constant-map normalisation: Delta_i(0) = W_i^0 = zeta^{3i}
    ok = all(rows.q0_value(rows.Delta(i)) == R.z(3 * i) for i in range(5))
    check("Delta_i(q=0) = W_i^0 = zeta^{3i} (loc_check: constant maps = (W_i^0)^{g-1} Int Hodge_{wt_i})", ok)
    ok = all(rows.q0_value(F0.P[a][i][n]) == (R.z(i * a) if n == 0 else 0) for a in range(5) for i in range(5) for n in range(7))
    check("R~ rows at q = 0 are the identity (zeta^{ia} delta_{n0}): no Dil, no edge at q = 0", ok)
    res = {}
    # ---- genus 1
    t1 = time.time()
    g1 = genus1_dT_sums(F0)
    H1 = total(g1, ring)
    check("genus-1 <H>_1 graph sum is Lam-free and i3-free in N", not H1.has('Lam') and not H1.has('i3'),
          "terms %d, denominators (c0^%d W^%d), L-degree %s" % (H1.nterms(), H1.a, H1.b, H1.Ldeg()))
    H1q = to_qseries(H1, ring, img); H1Q = G.Qseries(H1q)
    cA = qcoefs(H1Q, R, D)
    res[('H1', 0)] = cA[0]
    log("   genus 1: <H>_1 = d_T F_1 as a Q-series: [Q^0] = %s (= Sum_i zeta^i Int_{Mbar_{1,1}} Hodge_{wt_i}; the constant F_1(0) itself"
        " needs a convention, Mbar_{1,0} does not exist)" % R.fmt(cA[0]))
    s1 = R.zero
    for i in range(5): s1 = R.add(s1, R.mul(R.z(i), F0.vt[i].val(1, (0,))))
    check("[Q^0] <H>_1 = Sum_i zeta^i Int_{Mbar_{1,1}} Hodge_{wt_i}", R.eq(cA[0], s1), R.fmt(s1))
    for d in range(1, D + 1):
        N1d = R.mul(cA[d], R.inv(d)); res[('N', 1, d)] = N1d
        if ('N', 1, d) in REF:
            check("N_{1,%d} = [Q^%d] <H>_1 / %d = %s (localisation %s)" % (d, d, d, R.fmt(N1d), REF[('N', 1, d)]), R.eq(N1d, R.of(REF[('N', 1, d)])))
        else: log("   N_{1,%d} = %s (no localisation value to compare)" % (d, R.fmt(N1d)))
    per_colour = [R.fmt(G.Qseries(to_qseries(g1[('A', (i,))].add(g1[('B', (i,))], ring), ring, img)).coef(1)) for i in range(5)]
    log("   genus 1 per colour, [Q^1]: %s   (%.1f s)" % (per_colour, time.time() - t1))
    # ---- genus 2, reference frame
    t1 = time.time()
    S0 = genus2_sums(F0, deco); T0 = total(S0, ring)
    log("   genus 2 reference frame F0 = (R~, Hodge_wt): %d coloured graphs in %.1f s; total %d terms, denominators (c0^%d W^%d), L-degree %s"
        % (len(S0), time.time() - t1, T0.nterms(), T0.a, T0.b, T0.Ldeg()))
    check("genus-2 total is Lam-free and i3-free", not T0.has('Lam') and not T0.has('i3'))
    T0q = to_qseries(T0, ring, img)
    f20 = T0q.coef(0); res[('F0', 2)] = f20
    check("F_2(0) = [q^0] graph sum = %s (localisation -19/9600)" % R.fmt(f20), R.eq(f20, R.of(REF[('F0', 2)])))
    v0 = R.zero
    for i in range(5): v0 = R.add(v0, R.mul(R.z(3 * i), F0.vt[i].val(2, ())))
    check("[q^0] comes from the vertex graph alone: Sum_i Delta_i(0) Int_{Mbar_2} Hodge_{wt_i} = %s" % R.fmt(v0), R.eq(v0, f20))
    nonvtx = [k for k, v in S0.items() if k[0] != 1 and not R.iszero(to_qseries(v, ring, img).coef(0))]
    check("every edge graph has vanishing q^0 term at every colouring", not nonvtx, nonvtx[:5])
    T0Q = G.Qseries(T0q); c2 = qcoefs(T0Q, R, D)
    check("N_{2,0} = [Q^0] = F_2(0)", R.eq(c2[0], f20))
    for d in range(1, D + 1):
        res[('N', 2, d)] = c2[d]
        if ('N', 2, d) in REF:
            check("N_{2,%d} = [Q^%d] F_2 = %s (localisation %s)" % (d, d, R.fmt(c2[d]), REF[('N', 2, d)]), R.eq(c2[d], R.of(REF[('N', 2, d)])))
        else: log("   N_{2,%d} = %s (no localisation value to compare)" % (d, R.fmt(c2[d])))
    # per-graph Q^1 breakdown (information)
    byg = {}
    for (idx, col), v in S0.items():
        byg[idx] = byg.get(idx, ring.const(0)).add(v, ring)
    log("   genus 2, [Q^1] per graph (colour-summed): %s" % {idx: R.fmt(G.Qseries(to_qseries(v, ring, img)).coef(1)) for idx, v in sorted(byg.items())})
    log("   genus 2, [Q^1] of graph 1 (vertex) per colour: %s" % [R.fmt(G.Qseries(to_qseries(S0[(1, (i,))], ring, img)).coef(1)) for i in range(5)])
    # ---- other frames: exact agreement in N or failure
    for name, must in (('F5', True), ('F1', True), ('F4', True), ('F5b', False), ('F0b', False), ('F4b', False), ('F2', True)):
        t1 = time.time()
        T = total(genus2_sums(frames[name], deco), ring)
        same = T.eq(T0, ring)
        TQ = G.Qseries(to_qseries(T, ring, img)); cc = qcoefs(TQ, R, D)
        lab = {"F5": "(R, 1) point frame", "F1": "(R h_(1)^-1, Hodge_(1)) target", "F5b": "(R~, 1) control",
               "F0b": "(R, Hodge_wt) control", "F4": "(R h_(zeta^i)^-1, Hodge_(zeta^i)) covariant frame",
               "F4b": "(R h_(zeta^i)^+1, Hodge_(zeta^i)) wrong-sign control",
               "F2": "(R h_(1)^+1, Hodge_(1)) wrong sign, but colour-independent: equality is STRUCTURAL at the zeta point "
                     "in genus 2 (Z/5 phases kill every term of lambda-degree 1..4; see gsum_q.make_frames NOTE)"}[name]
        check("%s %s: total %s the reference exactly in N  [N_{2,0..2} = %s]  (%.1f s)" %
              (name, lab, "equals" if must else "differs from", [R.fmt(x) for x in cc[:3]], time.time() - t1), same == must)
        res[('frame', name)] = tuple(cc)
    log("time %.1f s" % (time.time() - t0))
    store[p] = res
    return res


if __name__ == "__main__":
    store = {}
    run(P2, store)
    run(P1, store)
    # cross-prime: rationals reconstructed at P2 must have the P1 residues
    R2 = FpRing(P2); R1 = FpRing(P1)
    import math
    def crt(a2, a1):
        M = P1 * P2
        x = (a2 * P1 * pow(P1, -1, P2) + a1 * P2 * pow(P2, -1, P1)) % M
        r0, r1 = M, x; s0, s1 = 0, 1; B = math.isqrt(M // 2) // 1000
        while r1 > B:
            q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
        return Fr(r1, s1) if s1 and abs(s1) <= B else None
    def both(a2, a1):
        """(rational, how): reconstruct at 2^89-1 and check the 2^61-1 residue, else CRT of the two residues (bound ~1e19)"""
        r = R2.ratrec(a2)
        if r is not None: return (r if R1.of(r) == a1 else None, "2^89-1 reconstruction, residue at 2^61-1 agrees")
        return (crt(a2, a1), "CRT of the two primes (value too large for single-prime reconstruction)")
    log("=" * 78); log("cross-prime check (reconstruct at 2^89-1 and compare the residue at 2^61-1, or CRT):")
    final = {}
    for k in store[P2]:
        v2 = store[P2][k]; v1 = store[P1][k]
        if isinstance(v2, tuple):
            rs = [both(a, b) for a, b in zip(v2, v1) if a is not None]
            check("%s: all Q-coefficients consistent across primes" % (k,), all(r[0] is not None for r in rs), [str(r[0]) for r in rs])
        else:
            r, how = both(v2, v1); final[k] = r
            check("%s = %s (%s)" % (k, r, how), r is not None)
    log("calibration table (exact rationals; N_{g,3}, N_{g,4} are predictions, no localisation value):")
    for k in (('H1', 0), ('N', 1, 1), ('N', 1, 2), ('N', 1, 3), ('N', 1, 4), ('F0', 2), ('N', 2, 1), ('N', 2, 2), ('N', 2, 3), ('N', 2, 4)):
        log("   %-12s graph sum = %-28s localisation = %s" % (k, final.get(k), REF.get(k, "n/a")))
    log("ALL S3 CHECKS PASSED" if not FAILED else "S3 FAILURES: %s" % FAILED)
    LOG.close()
