# -*- coding: utf-8 -*-
"""rf_reparam.py -- B-variation is a change of coordinate on the base.  For a B-perturbed system X (variable u) let phi = v_X^-1 o v_q
(v_q = the quintic model coordinate, tau_q = 1/2F1(1,4/5;9/5;.)), u = phi(w).  Then  phi'^5 B_X(phi) = B_q  and the pulled-back system
b~ = phi' b(phi) = (f~, B_q/(f~^2 g~^2), f~, g~, g~) has the QUINTIC B, i.e. it is a univ (F, Gm) system:
    F~ = phi' f_X(phi)/t_c'(w),   Gm~ = phi' g_X(phi)/g_q(w).
This script runs that (F~, Gm~) system through the Stage-1 univ pipeline ONLY (un_gauge.pert_gauge = ur_gauge solver with the quintic
tau, un_e4_cases.normalise, un_gauge.rows_from_gauge, un_engine 'nat' base with the quintic E_M), plus the constant symplectic scalar,
and compares graph by graph with the referee's direct B-perturbed run (rf_run json):
    prediction (coordinate independence of each outer graph):  w^(2g-2) V~_i(w) = (phi(w)/w)^(2-2g) G_i(phi(w)),  V_i = u^(2-2g) G_i(u).
Usage: py -I -B rf_reparam.py X g prime [--out DIR]
"""
import os, sys, json, math, argparse, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rf_run  # noqa: sets sys.path (project, univ)
import rf_bsys as RB
from rf_ser import L, fr_mod, pmul, ppow
from rf_run import PR, SCR, in_T, rat, CG


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("X"); ap.add_argument("g", type=int); ap.add_argument("prime"); ap.add_argument("--out", default=SCR)
    a = ap.parse_args(); g = a.g; p = PR[a.prime]; K = 2 * g - 2 + 2; KZ = 3 * g - 3; t0 = time.time()
    lf = open(os.path.join(a.out, "rf_reparam_%s_g%d_%s.log" % (a.X, g, a.prime)), "w", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); lf.write(s + "\n")
    S = RB.build(a.X); Sq = RB.build("R0"); N = 30
    am = Fr(S['m']) ** 5
    vX = L(1, [S['v'].co(k) / am for k in range(1, N)], N); vq = L(1, [Sq['v'].co(k) / am for k in range(1, N)], N)
    phi = vX.revert().compose(vq); dphi = phi.der()
    # check phi'^5 B_X(phi) = B_q
    BX = S['B'].sh(1)                                              # u B_X, power series
    lhs = dphi.mul(dphi).mul(dphi).mul(dphi).mul(dphi).mul(BX.compose(phi)).mul(phi.inv())
    d = lhs.sub(Sq['B']); log("%s: phi = %r;  phi'^5 B_X(phi) = B_q: %s" % (a.X, phi, all(d.co(k) == 0 for k in range(-1, 20))))
    ft = dphi.mul(S['b'][0].compose(phi)); gt = dphi.mul(S['b'][3].compose(phi))
    Ft = ft.mul(Sq['b'][0].inv()); Gt = gt.mul(Sq['b'][3].inv())
    Ftup = tuple(Ft.co(k) for k in range(26)); Gtup = tuple(Gt.co(k) for k in range(26))
    log("   F~ = %r ;  Gm~ = %r" % (Ft, Gt))
    from un_gauge import pert_gauge, rows_from_gauge
    from un_e4_cases import normalise, p11_series
    from un_engine import get_base, get_nat, Var, collect
    from un_core import Arith, lift_rows, symplectic_check
    res = pert_gauge(12, 24, Ftup, Gtup, ("0", "0")); res, mz = normalise(res, "canon")
    q0 = [res['p'][0][c].coef(0) for c in range(5)]; C = q0[0] * q0[3]
    assert C == q0[1] * q0[2] == q0[4] ** 2
    sig = 1 / RB.qsqrt(C)
    if sig != 1: res = dict(res); res['p'] = [[s.scale(sig) for s in row] for row in res['p']]
    c0 = p11_series(res['p'])[0]; log("   univ pipeline: C = %s, sigma = %s, c0 = %s" % (C, sig, c0))
    nr = get_nat(a.prime, g, ("0", "0"), K); Qh = rows_from_gauge(res, nr['Rh'], nr['R'], KZ, K)
    ar = Arith(p, K, 1); QL = lift_rows(ar, Qh); log("   rows symplectic: %s" % (not symplectic_check(ar, QL, KZ, nr['eta_inv'])))
    v = Var(get_base(a.prime, g, ("0", "0"), K), QL, 1, None); tab = collect(v, [v.graph(t) for t in range(len(v.inner_graphs))])
    # direct run of X
    J = json.load(open(os.path.join(a.out, "rf_%s_g%d_%s_canon_none_k%d.json" % (a.X, g, a.prime, K))))
    TX = {int(i): c for i, c in J['table'].items()}
    phr = [fr_mod(phi.co(k), p) for k in range(K + 1)]; ph1 = phr[1:K + 1]       # phi = w * ph1
    pw = ppow(ph1, 2 - 2 * g, K, p)
    def comp(G):
        out = [0] * K; P_ = [1] + [0] * (K - 1); phl = phr[:K]
        for k in range(K):
            out = [(x + G[k] * y) % p for x, y in zip(out, P_)]; P_ = pmul(P_, phl, K, p)
        return out
    nbad = 0; nnz = 0
    for i in range(len(v.inner_graphs)):
        pred = pmul(pw, comp(TX.get(i, [0] * K)), K, p); got = [x % p for x in tab.get(i, [0] * K)]
        if pred != got: nbad += 1; log("      graph %d: MISMATCH pred %s got %s" % (i, pred[:3], got[:3]))
        if any(got): nnz += 1
    log("   graph by graph, w^(2g-2) V~_i(w) = (phi/w)^(2-2g) G_i(phi(w)) mod w^%d: %s (%d graphs, %d non-zero)" % (K, "ALL EQUAL" if not nbad else "%d MISMATCHES" % nbad, len(v.inner_graphs), nnz))
    sing = [i for i in range(1, len(v.inner_graphs)) if any(tab.get(i, [0] * K)[:2 * g - 2])]
    Tt = ft.integ().sc(c0); Vv = in_T(tab.get(0, [0] * K), Tt, g, K, p)
    log("   univ-pipeline run of the pulled-back system: singular edge graphs %s; vertex in T~'' %s (c_g = %s)" % (sing or "none", [rat(x, p) for x in Vv], CG[g]))
    log("RF_REPARAM %s g%d %s: %s, %.0f s" % (a.X, g, a.prime, "PASS" if not nbad else "FAIL", time.time() - t0))


if __name__ == "__main__":
    main()
