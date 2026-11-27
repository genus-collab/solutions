# -*- coding: utf-8 -*-
"""genus3_step5.py -- u-dependent passes at genus 3 for the covariant-scalar class on the GENUS-3-ALONE [u^1] system (the combined genus-2+3
system has no covariant-scalar solution): N_1(t) = N_1^(0) + t N_1^(free) runs over the 1-dimensional solution family; for every edge graph
the [u^2] coefficient of G_Gamma[M_0 exp(u N_1(t))] is a polynomial of degree <= 2 in t (two edge factors O(u), each linear in N_1), the
[u^3] coefficient of degree <= 3.  Passes at t = 0, 1, 2, 3 (E <= 3 graphs; the full table at t = 0) give these polynomials exactly; the
2-edge graphs' [u^2] polynomials must have a COMMON root for the scalar class to survive at O(u^2) at genus 3 alone.
Also: the LEFT-placed family, and the first u-correction of the tautological factorisation (fq3_taut) as the control point where every
edge graph is O(u^2) by construction (the 1-edge graphs then need N_2).
Run:  py -I genus3_step5.py --prime 89 [--nproc 24] [--placement right|left]  -> appends to logs/genus3_ou_<p>.log, state g3_state_<p>.json
"""
import os, sys, time, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_hat import vander_inv, interp_fp, tc_series
from fq3_steps import frame_of, State, hat_pass, IDX_ALL, IDX_VERTEX, IDX_BY_E, GD, total_of, fmt_u, tc_fmt, log_table, solve_system
from fq3_m import sp_mats, sp_params, expand, class_columns, NSP, sp_add, sp_scale
from fq3_msum import tc_expansion3
from fq3_par import setup
from fq3_rec import Rec
from fq3_taut import taut_N1


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--nproc", type=int, default=24)
    ap.add_argument("--placement", default="right"); ap.add_argument("--NN", type=int, default=21); ap.add_argument("--taut", action="store_true"); ap.add_argument("--only_taut", action="store_true")
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1
    frame = frame_of(p, (0, 0), 34, 6); state = State(os.path.join(HERE, "g3_state_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_ou_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    st = setup(frame); R = st['R']; rec = Rec(R); M0 = st['M0']; tc = tc_series(st['C'], R, 10); T0 = time.time()
    tot_id = state.get('total_id'); J3 = {int(i): v for i, v in state.get('J').items()}; A0 = {int(i): v for i, v in state.get('A_M0').items()}
    log("-" * 110); log("(5) u-dependent passes, covariant-scalar class, genus-3-alone [u^1] system, placement %s  (%s)" % (a.placement, time.strftime("%Y-%m-%d %H:%M")))
    if a.only_taut:
        log("   (only the tautological control pass)")
        N1t, chk = taut_N1(st); log("   taut checks %s" % chk)
        A, info = hat_pass(frame, ('M0exp', [N1t], 'right'), IDX_ALL, a.nproc, log, a.NN, label="taut")
        log_table(log, R, rec, A, info, IDX_ALL); state.put('A_taut_N1', {str(i): A[i] for i in A})
        check("tautological N_1: every edge graph O(u^2); graphs with E >= 2 O(u^4) (regular)", all(A[i][k] == 0 for i in IDX_ALL if i != 1 for k in range(2)) and all(A[i][k] == 0 for i in IDX_ALL if GD[i]['E'] >= 2 for k in range(4)))
        check("tautological N_1: total = point-frame total", total_of(R, A) == tot_id)
        log("STEP 5 (taut only) ALL CHECKS PASSED" if not FAILED else "STEP 5 FAILURES: %s" % FAILED); LOG.close(); return
    cols, names = class_columns(R, 'cov', M0, a.placement)
    rows = [J3[2], J3[5]]; rhs = [(-A0[2][1]) % p, (-A0[5][1]) % p]
    Jr, x, ns, rank, cons = solve_system(R, rows, rhs, cols)
    assert cons and len(ns) == 1, (cons, len(ns))
    log("   family (b1, b3, b5)(t) = %s + t %s" % ([rec.fmt(c) for c in x], [rec.fmt(c) for c in ns[0]]))
    N0 = sp_mats(expand(cols, x, p), p); Nf = sp_mats(expand(cols, ns[0], p), p)
    tvals = [0, 1, 2, 3]; Vinv = vander_inv(tvals, p); tables = {}
    for t in tvals:
        Nt = sp_add(N0, sp_scale(Nf, t, p), p)
        idxs = IDX_ALL if t == 0 else (IDX_BY_E[0] + IDX_BY_E[1] + IDX_BY_E[2] + IDX_BY_E[3])
        t0 = time.time(); A, info = hat_pass(frame, ('M0exp', [Nt], a.placement), idxs, a.nproc, log, a.NN, label="t=%d" % t)
        log("   pass t = %d (%d graphs): %.0f s; max eps-degree %d" % (t, len(idxs), time.time() - t0, max(i['maxdeg'] or 0 for i in info.values())))
        tables[t] = A; state.put('A_cov5_%s_t%d' % (a.placement, t), {str(i): A[i] for i in A})
        if t == 0:
            log_table(log, R, rec, A, info, IDX_ALL)
            check("t = 0: total = point-frame total (M-independence)", total_of(R, A) == tot_id)
            check("t = 0: 1-edge graphs 2, 5 vanish through u^1 (the [u^1] equations are satisfied)", all(A[i][k] == 0 for i in (2, 5) for k in range(2)))
            check("t = 0: every graph is O(u^E)", all(A[i][k] == 0 for i in IDX_ALL for k in range(GD[i]['E'])))
            reg = [i for i in IDX_ALL if i != IDX_VERTEX and all(A[i][k] == 0 for k in range(4))]
            log("   t = 0: regular edge graphs %s; non-regular %s" % (reg, [i for i in IDX_ALL if i != IDX_VERTEX and i not in reg]))
    # polynomials in t of [u^2] (2-edge graphs) and [u^3] (3-edge graphs), degrees <= 2 and <= 3
    log("   obstruction polynomials along the family (coefficients of t^0, t^1, t^2, t^3):")
    common = None
    for i in IDX_BY_E[2]:
        c2 = interp_fp([tables[t][i][2] for t in tvals], Vinv, p); c3 = interp_fp([tables[t][i][3] for t in tvals], Vinv, p)
        check("2-edge graph %d: [u^2] is a polynomial of degree <= 2 in t" % i, c2[3] == 0)
        log("      graph %2d (E=2): [u^2](t) = %s ; [u^3](t) = %s" % (i, [rec.fmt(c) for c in c2[:3]], [rec.fmt(c) for c in c3]))
        state.put('obs_u2_%s_%d' % (a.placement, i), c2[:3])
        # roots of the quadratic over F_p (if any): t = (-b +- sqrt(b^2 - 4ac)) / 2a
        A2, B2, C2 = c2[2], c2[1], c2[0]
        if A2:
            disc = (B2 * B2 - 4 * A2 * C2) % p; r = pow(disc, (p + 1) // 4, p) if p % 4 == 3 else None
            roots = set()
            if r is not None and r * r % p == disc:
                roots = {(-B2 + r) * pow(2 * A2, p - 2, p) % p, (-B2 - r) * pow(2 * A2, p - 2, p) % p}
            log("         F_p roots of [u^2](t): %s" % ([rec.fmt(z) for z in roots] if roots else "none (discriminant not a square / no root)"))
        elif B2: roots = {(-C2) * pow(B2, p - 2, p) % p}; log("         linear: root t = %s" % rec.fmt(list(roots)[0]))
        else: roots = None if C2 else "all"; log("         constant %s" % rec.fmt(C2))
        if roots == "all": continue
        common = set(roots or set()) if common is None else (common & set(roots or set()))
    check("KEY: the five 2-edge graphs' [u^2] obstructions have a common zero on the covariant-scalar family (genus 3 alone)", bool(common), "common roots: %s" % ([rec.fmt(z) for z in common] if common else "none"))
    for i in IDX_BY_E[3]:
        c3 = interp_fp([tables[t][i][3] for t in tvals], Vinv, p)
        log("      graph %2d (E=3): [u^3](t) = %s" % (i, [rec.fmt(c) for c in c3]))
    state.put('obs_common_%s' % a.placement, sorted(common) if common else [])
    if a.taut:
        log("   control: tautological N_1 (every edge graph O(u^2) by construction), full table")
        N1t, chk = taut_N1(st); log("   taut checks %s" % chk)
        A, info = hat_pass(frame, ('M0exp', [N1t], 'right'), IDX_ALL, a.nproc, log, a.NN, label="taut")
        log_table(log, R, rec, A, info, IDX_ALL); state.put('A_taut_N1', {str(i): A[i] for i in A})
        check("tautological N_1: every edge graph O(u^2); graphs with E >= 2 O(u^4) (regular)", all(A[i][k] == 0 for i in IDX_ALL if i != 1 for k in range(2)) and all(A[i][k] == 0 for i in IDX_ALL if GD[i]['E'] >= 2 for k in range(4)))
        check("tautological N_1: total = point-frame total", total_of(R, A) == tot_id)
    log("step 5 time %.0f s" % (time.time() - T0)); log("STEP 5 ALL CHECKS PASSED" if not FAILED else "STEP 5 FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
