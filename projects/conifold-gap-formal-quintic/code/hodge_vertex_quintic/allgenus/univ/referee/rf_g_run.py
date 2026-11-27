# -*- coding: utf-8 -*-
"""rf_g_run.py -- REFEREE (univ): gauges of the worker's cases and of two NEW systems with the referee's own solver (rf_gauge), exact
over Q; regression of the stored-normalisation constants printed in the worker's e4.log; (G) and symplecticity in both normalisations.
Run: py -I -B rf_g_run.py  -> rf_g_run.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_pf import system
from rf_gauge import gauge, check_equation, regular, P11, times_c
from rf_gchk import check_G, symplectic_defects

F = Fr
CASES = {   # name: (F, Gm, (alpha, beta), worker's stored P_11(0,z) from e4.log or None)
    "quintic": ((1,), (1,), (0, 0), {5: F(-43397, 807975000), 10: F(-765368951400317, 1305647201250000000)}),
    "quintic_ab": ((1,), (1,), (F(1, 3), F(-2, 7)), {5: F(-43397, 807975000), 10: F(-765368951400317, 1305647201250000000)}),
    "c1": ((1, F(1, 3), F(-2, 7), F(5, 11)), (1, F(-3, 5), F(1, 2), F(2, 3)), (0, 0),
           {5: F(31764344287, 1163868750), 10: F(5761596837414242725693984281037, 41740938549382073203125000)}),
    "c2": ((F(3, 2), F(1, 3)), (1, F(-3, 5), F(1, 2)), (0, 0),
           {0: F(2, 3), 5: F(915238702031, 359952862500), 10: F(-382006806145740134909048581153, 221713447385390412375000000)}),
    "c3": ((1, F(-7, 4), 3, F(-1, 9)), (1, F(11, 6), F(-5, 3), F(4, 13)), (0, 0), {5: F(38700127538971973, 59156697600000)}),
    "c4": ((1, F(2, 5), F(-1, 3)), (1, F(1, 7), F(3, 4)), (F(1, 3), F(-2, 7)), {5: F(-171815241648683, 24387917400000)}),
    "r1": ((2, F(-1, 4), F(5, 3), F(1, 6)), (F(3, 2), F(2, 9), -1, F(1, 5)), (0, 0), None),
    "r2": ((F(1, 2), 3, -2), (1, F(-1, 3), 4), (F(1, 3), F(-2, 7)), None),
}


def build(name, KG=16, N=40, norm="canon"):
    Fp, Gp, ab, _ = CASES[name]
    S = system(N, Fp, Gp, *ab)
    p, info = gauge(S['b'], S['bt'], KG, norm)
    return S, p, info


def main():
    LOG = open(os.path.join(HERE, "rf_g_run.log"), "a", encoding="ascii")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    fails = []
    def check(name, ok, det=""):
        log(("PASS " if ok else "FAIL ") + name + (": %s" % (det,) if det else ""))
        if not ok: fails.append(name)
    log("=" * 100); log("rf_g_run: own gauge solver (direct Laurent undetermined coefficients), exact over Q   [%s]" % time.strftime("%Y-%m-%dT%H:%M:%S"))
    names = sys.argv[1].split(",") if len(sys.argv) > 1 else list(CASES)
    for name in names:
        t0 = time.time(); Fp, Gp, ab, ref = CASES[name]
        log("-" * 100); log("case %s: F = %s, Gm = %s, frame %s" % (name, [str(x) for x in Fp], [str(x) for x in Gp], tuple(str(x) for x in ab)))
        S, ps, info = build(name, norm="stored")
        Sc, pc, infoc = build(name, norm="canon")
        log("   resonances (z-order: index) %s; precision of p_k: %s" % ({k: r for k, r in info['res'].items() if r is not None}, info['prec']))
        check("   gauge equation satisfied to the tracked precision (stored and canonical)", not check_equation(ps, S['b'], S['bt']) and not check_equation(pc, Sc['b'], Sc['bt']))
        check("   every p_{k,b} (k <= 16) is a power series in u (Theorem U), both normalisations", not regular(ps) and not regular(pc))
        st = P11(ps); cn = P11(pc)
        log("   stored    P_11(0, z) = %s" % "  ".join("[z^%d] %s" % (k, x) for k, x in st.items()))
        log("   canonical P_11(0, z) = %s" % "  ".join("[z^%d] %s" % (k, x) for k, x in cn.items()))
        if ref:
            check("   REGRESSION: stored P_11(0, z) equals the worker's e4.log values %s" % sorted(ref), all(st[k] == x for k, x in ref.items()))
        check("   canonical: P_11(0, z) is constant", all(x == 0 for k, x in cn.items() if k > 0))
        c0 = st[0]; mult = {0: Fr(1)}
        inv = {0: 1 / st[0]}
        for k in range(5, len(ps), 5): inv[k] = -sum(st[j] * inv[k - j] for j in range(5, k + 1, 5)) / st[0]
        mult = {k: c0 * x for k, x in inv.items()}
        pn = times_c(ps, mult)
        same = all(pn[k][c].sub(pc[k][c]).trunc(min(pn[k][c].prec, pc[k][c].prec, 12)).iszero() for k in range(len(pc)) for c in range(5))
        check("   canonical (imposed in the recursion) = stored x c0/P^stored_11(0, z) (post-hoc), mod u^12", same)
        ok, c0, bad, nchk, pos = check_G(pc, Sc['b'], Sc['v'], 8, -8, 8)
        check("   (G) canonical: P Y^C_model = c0 Y^C_system on z^-8..z^8 mod u^8, c0 = %s" % c0, ok, bad[:4])
        ok2, _, bad2, _, pos2 = check_G(ps, S['b'], S['v'], 8, -8, 8)
        check("   (G) stored: P Y^C_model - c0 Y^C_system has a positive part (expected), at z^%s" % pos2, (not ok2) and 5 in pos2)
        sd = symplectic_defects(pc, 10); ss = symplectic_defects(ps, 10)
        check("   canonical gauge symplectic: P(-z)^T eta P(z) = eta through z^16 mod u^10", not sd, sd[:4])
        check("   stored gauge symplectic exactly through z^9, first defect at z^10 (c(z)c(-z) = 1 - m5^2 z^10 + ..)", ss and min(n for n, c, v in ss) == 10, ss[:3])
        log("   (%.1f s)" % (time.time() - t0))
    log("RF_G_RUN %s" % ("ALL CHECKS PASSED" if not fails else "FAILURES: %s" % fails)); LOG.close()


if __name__ == "__main__":
    main()
