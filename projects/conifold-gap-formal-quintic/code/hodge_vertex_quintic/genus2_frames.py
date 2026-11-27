# -*- coding: utf-8 -*-
"""genus2_frames.py -- PLAN step 4: genus-2 coloured graph terms of the formal quintic in several Hodge frames, their
formal L-degree and their conifold-frame (Cf) Laurent expansion in s (L = 1/s, s^5 = u = 1 - 3125 q).
Run:  py -I genus2_frames.py   (writes logs/genus2_frames.log).
Frames (gsum_q.make_frames): F0 (R~, Hodge_wt) reference; F1 (R h_(1)^-1, Hodge_(1)) target; F2 (R h_(1)^+1, Hodge_(1))
control; F3 (R h_(5)^-1, Hodge_(5)); F4 (R h_(zeta^i)^-1, Hodge_(zeta^i)); F5 (R, 1) point frame; F0b (R, Hodge_wt) and
F5b (R~, 1) controls.  Per coloured graph (all 5^V colourings): sum over decorations -> NEl, reduced; (i) L-degree range;
(ii) Cf valuation in s (negative = singular), bound 5 Sum_v (2 g_v - 2)_+ in s-units; Lam-independence (Lam = 0, 1).
Vertex graph: s-expansion per colour and, when the series is in u = s^5, the polar part c/t_c^2 + c'/t_c.
Totals: exact equality in N across frames; Cf expansion; polar part in t_c.  (iii) sensitivity of F1 to (alpha, beta) of
varpi^cf_0.  Two primes.
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_rows import Rows
from fq_graphs import Decorated
from fq_conifold import Conifold
from gsum_q import (make_frames, genus2_sums, total, to_sseries, cf_images_u, is_u_series, s_to_u, polar_in_tc, ser_fmt)

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_frames.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

FRAMES = ['F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F0b', 'F5b', 'F4b']
COVARIANT = {'F0', 'F4', 'F5', 'F0b', 'F5b', 'F4b'}  # colour shift i -> i+1 is a symmetry of the coloured terms
# totals that must equal F0: F5 and F4 by the frame-change identity; F1, F2, F3 structurally (zeta point, genus 2: see
# gsum_q.make_frames NOTE); F0b, F5b, F4b must differ.
MUST_EQUAL = {'F1', 'F2', 'F3', 'F4', 'F5'}
NU = 10


def expand_info(S, R, tc_u, nterms=5):
    """text: s-expansion (first nterms non-zero terms) and, for a u-series, the t_c polar part"""
    if not S.c: return "0 + O(s^%d)" % S.prec, None
    terms = [(S.v + i, x) for i, x in enumerate(S.c) if not R.iszero(x)][:nterms]
    txt = " + ".join("(%s) s^%d" % (R.fmt(x), k) for k, x in terms) + " + ..."
    pol = None
    if is_u_series(S):
        U = s_to_u(S)
        T = polar_in_tc(U, tc_u)
        pol = (T.coef(-2) if T.v <= -2 else R.zero, T.coef(-1) if T.v <= -1 else R.zero, T.coef(0) if T.prec > 0 else None)
        txt += "   [u-series; in t_c: (%s)/t_c^2 + (%s)/t_c + (%s) + ...]" % (R.fmt(pol[0]), R.fmt(pol[1]), R.fmt(pol[2]) if pol[2] is not None else "?")
    return txt, pol


def frame_report(name, F, deco, ring, R, imu0, imu1, tc_u, ref_total):
    t0 = time.time()
    S = genus2_sums(F, deco); T = total(S, ring)
    tsum = time.time() - t0
    log("-" * 78); log("frame %s: %d coloured graphs summed in %.1f s" % (name, len(S), tsum))
    t0 = time.time()
    ser = {}; vals = {}
    lam_ok = True
    for key, el in S.items():
        s0 = to_sseries(el, ring, imu0); ser[key] = s0
        s1 = to_sseries(el, ring, imu1)
        if not s0.eq(s1, min(s0.prec, s1.prec) - 1): lam_ok = False
        vals[key] = (s0.v if s0.c else None, el.Ldeg())
    check("%s: Cf images independent of Lam (0 vs 1) for all coloured graphs" % name, lam_ok)
    log("   %-6s %-10s %-24s %6s %8s %8s %8s %8s" % ("graph", "g_v", "edges", "bound", "nsing", "minval", "maxLdeg", "minLdeg"))
    tab = {}
    for G in deco.graphs:
        idx = G['idx']; keys = [k for k in S if k[0] == idx]
        vs = [vals[k][0] for k in keys if vals[k][0] is not None]
        lds = [vals[k][1] for k in keys if vals[k][1] is not None]
        nsing = sum(1 for v in vs if v < 0)
        bound = 5 * G['bound']
        row = (bound, nsing, len(keys), min(vs) if vs else None, max(l[1] for l in lds) if lds else None, min(l[0] for l in lds) if lds else None)
        tab[idx] = row
        log("   %-6d %-10s %-24s %6d %8s %8s %8s %8s" % (idx, G['gv'], G['ed'], bound, "%d/%d" % (nsing, len(keys)), row[3], row[4], row[5]))
    # Z/5 symmetry
    zok = True; zval = True
    for (idx, col), el in S.items():
        col2 = tuple((c + 1) % 5 for c in col)
        if vals[(idx, col2)][0] != vals[(idx, col)][0]: zval = False
        if name in COVARIANT and not el.eq(S[(idx, col2)], ring): zok = False
    check("%s: Cf valuations constant along Z/5 colour orbits" % name, zval)
    if name in COVARIANT: check("%s: coloured terms exactly Z/5 invariant in N (colour shift)" % name, zok)
    else: log("   (%s is not Z/5 covariant: common weight; exact colour-shift invariance not expected)" % name)
    # vertex graph
    log("   vertex graph (idx 1) per colour:")
    for i in range(5):
        txt, pol = expand_info(ser[(1, (i,))], R, tc_u)
        log("      colour %d: Ldeg %s, val %s: %s" % (i, vals[(1, (i,))][1], vals[(1, (i,))][0], txt))
    # total
    Tser = None
    for s in ser.values(): Tser = s if Tser is None else Tser.add(s)
    txt, pol = expand_info(Tser, R, tc_u, 7)
    log("   total: %d terms, denominators (c0^%d W^%d), L-degree %s; Cf: val %s: %s" % (T.nterms(), T.a, T.b, T.Ldeg(), Tser.v if Tser.c else None, txt))
    Tdirect = to_sseries(T, ring, imu0)
    check("%s: Cf(total NEl) = sum of coloured Cf series" % name, Tdirect.eq(Tser, min(Tdirect.prec, Tser.prec) - 1))
    if ref_total is not None:
        same = T.eq(ref_total, ring)
        must = name in MUST_EQUAL
        note = " (structural at the zeta point, genus 2: Z/5 phases; not a test of the exponent sign)" if name in ('F1', 'F2', 'F3') else ""
        check("%s: total %s the F0 total exactly in N%s" % (name, "equals" if must else "differs from", note), same == must,
              "" if same == must else ("equal" if same else "different"))
    log("   (Cf substitutions %.1f s)" % (time.time() - t0))
    return S, T, tab, vals, pol, Tser


def run(p, results):
    R = FpRing(p); t0 = time.time()
    log("=" * 78); log("prime p = %d" % p)
    rows = Rows(R, Nz=6); ring = rows.ring
    deco = Decorated(2, hodge=True)
    C = Conifold(R, Nu=NU)
    imu0 = cf_images_u(C, R, 0); imu1 = cf_images_u(C, R, 1)
    tc_u = Ser(R, C.tc.v, [R.of(x) for x in C.tc.c])
    log("rows, decorations, conifold periods (u-precision %d) in %.1f s; t_c = %s" % (min(x.prec for x in imu0.values()), time.time() - t0, tc_u.fmt(3)))
    frames = make_frames(rows, FRAMES)
    ref = None; out = {}
    for name in FRAMES:
        S, T, tab, vals, pol, Tser = frame_report(name, frames[name], deco, ring, R, imu0, imu1, tc_u, ref)
        if name == 'F0': ref = T
        out[name] = dict(tab=tab, vals={k: v[0] for k, v in vals.items()}, pol=pol, S=S, Tser=Tser)
    results[p] = out
    # (iii) sensitivity of the F1 verdict to (alpha, beta)
    if p == P2:
        log("-" * 78); log("(iii) sensitivity of the F1 Cf verdict to varpi^cf_0 = f0 + alpha g + beta f2:")
        S1 = out['F1']['S']; base = out['F1']['vals']
        for (al, be) in ((Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(2, 3), Fr(-1, 7))):
            Cab = Conifold(R, Nu=NU, alpha=al, beta=be); im = cf_images_u(Cab, R, 0)
            changed = []; lead_changed = 0
            for key, el in S1.items():
                s = to_sseries(el, ring, im)
                v = s.v if s.c else None
                if v != base[key]: changed.append(key)
                s0 = to_sseries(el, ring, imu0)
                if s.c and s0.c and s.coef(s.v) != s0.coef(s0.v): lead_changed += 1
            log("   (alpha, beta) = (%s, %s): valuations changed for %d/%d coloured graphs %s; leading coefficient changed for %d"
                % (al, be, len(changed), len(S1), changed[:4], lead_changed))
    log("time for prime %d: %.1f s" % (p, time.time() - t0))


if __name__ == "__main__":
    results = {}
    T0 = time.time()
    run(P2, results)
    run(P1, results)
    log("=" * 78); log("cross-prime comparison of the valuation tables:")
    for name in FRAMES:
        same = results[P2][name]['vals'] == results[P1][name]['vals'] and results[P2][name]['tab'] == results[P1][name]['tab']
        check("%s: per-coloured-graph Cf valuations and tables identical at 2^89-1 and 2^61-1" % name, same)
    log("summary (p = 2^89-1): per frame, per graph: nsing/total, min valuation, max Ldeg  [bound in s-units]")
    for name in FRAMES:
        tab = results[P2][name]['tab']
        log("   %-4s " % name + "  ".join("g%d:%s/%s,v%s,L%s[%d]" % (idx, r[1], r[2], r[3], r[4], r[0]) for idx, r in sorted(tab.items())))
    log("total time %.1f s" % (time.time() - T0))
    log("ALL S4 CHECKS PASSED" if not FAILED else "S4 FAILURES: %s" % FAILED)
    LOG.close()
