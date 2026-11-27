# -*- coding: utf-8 -*-
"""genus3_ou.py -- the formal quintic at GENUS 3: tautological factorisation M_0 = R~ mod zhat^7 and its u-corrections
M(zhat; u) = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3), tested graph by graph on the 42 stable graphs at the conifold point.

Usage:  py -I genus3_ou.py --prime 89|61 --steps 1,2,3,4,5 [--nproc 24] [--ab 0,0 | 1/3,-2/7] [--cls cov] [--placement right] [--K 6]
Modules: fq3_graphs (legged stable graphs, 42 outer graphs), fq3_mvertex (vertex CohFT M.omega^ to order 6), fq3_m (sp, classes, M_0,
u-dependent M, rows), fq3_msum (hatted theory in u, outer sum, eps-interpolation), fq3_par (multiprocessing workers), fq3_steps (steps).
Log: logs/genus3_ou_<p>[_ab].log; state: g3_state_<p>[_ab].json (residues).  Genus 3 only; exact finite computation; AI-checked only.
"""
import os, sys, time, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_steps import (frame_of, State, step1, step1_genus2, step1_sseries, step2, step3, step3_genus2, step4, pass_u, status, next_N,
                       IDX_ALL, IDX_VERTEX, IDX_BY_E, GD, total_of, fmt_u, tc_fmt)
from fq3_msum import tc_expansion3
from fq3_m import sp_mats, NSP
from fq3_par import setup
from fq_hat import tc_series


class Rec(object):
    """reconstruction with the default bound, then the 1e13 fallback; residues otherwise"""
    def __init__(self, R): self.R = R
    def fmt(self, a):
        R = self.R; a %= R.p
        r = R.ratrec(a)
        if r is None and R.p > 10 ** 20: r = R.ratrec(a, 10 ** 13)
        if r is not None: return str(r)
        for k in range(1, 5):
            r = R.ratrec(a * R.z(-k) % R.p)
            if r is not None: return "%s*zeta^%d" % (r, k)
        return "[%d mod p]" % a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", default="89"); ap.add_argument("--steps", default="1,2,3,4"); ap.add_argument("--nproc", type=int, default=24)
    ap.add_argument("--ab", default="0,0"); ap.add_argument("--K", type=int, default=6); ap.add_argument("--prec_abs", type=int, default=34)
    ap.add_argument("--cls", default="cov"); ap.add_argument("--placement", default="right"); ap.add_argument("--NN", type=int, default=21)
    ap.add_argument("--no_genus2", action="store_true"); ap.add_argument("--no_sseries", action="store_true")
    a = ap.parse_args()
    p = P2 if a.prime == "89" else P1; ab = tuple(Fr(x) for x in a.ab.split(","))
    tagab = "" if ab == (0, 0) else "_ab%s" % a.ab.replace("/", "o").replace(",", "_").replace("-", "m")
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    LOG = open(os.path.join(HERE, "logs", "genus3_ou_%s%s.log" % (a.prime, tagab)), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    steps = [s.strip() for s in a.steps.split(",")]
    frame = frame_of(p, ab, a.prec_abs, a.K); state = State(os.path.join(HERE, "g3_state_%s%s.json" % (a.prime, tagab)))
    T0 = time.time(); st = setup(frame); R = st['R']; rec = Rec(R)
    log("=" * 110); log("genus3_ou  p = 2^%s-1  (alpha, beta) = %s  prec_abs = %d  K = %d  nproc = %d  steps %s  %s" % (a.prime, ab, a.prec_abs, a.K, a.nproc, steps, time.strftime("%Y-%m-%d %H:%M")))
    log("setup (rows Nz = 6, conifold periods, SRows, symbol, HatRowsU, M_0): %.1f s" % (time.time() - T0))
    tc = tc_series(st['C'], R, 10)
    if "1" in steps:
        log("-" * 110); log("(1) hatted theory at genus 3, point frame M = 1, frozen sums, M-independence")
        d1 = step1(frame, st, state, a.nproc, log, check, rec)
        if not a.no_genus2: step1_genus2(frame, st, state, log, check, rec)
        if not a.no_sseries: step1_sseries(frame, st, state, log, check)
    tot_id = state.get('total_id')
    if "2" in steps:
        log("-" * 110); log("(2) pass A: the tautological M_0 = R~ mod zhat^7")
        step2(frame, st, state, a.nproc, log, check, rec, tot_id, tc)
    if "3" in steps:
        log("-" * 110); log("(3) frozen Jacobian J (75 directions, right placement) for the graphs with E <= 1")
        step3(frame, st, state, a.nproc, log, check, rec)
        if not a.no_genus2: step3_genus2(frame, st, state, log, check, rec)
    if "4" in steps:
        log("-" * 110); log("(4) the [u^1] systems on the classes cov / diag / circ / full, right and left placement")
        J3 = {int(i): v for i, v in state.get('J').items()}; A0 = {int(i): v for i, v in state.get('A_M0').items()}
        g1_3 = {i: A0[i][1] for i in (2, 5)}
        J2s = state.get('J2'); J2 = {int(i): v for i, v in J2s.items()} if J2s else None
        g1_2 = {2: R.of(Fr(1, 750)), 4: R.of(Fr(-13, 8000))}
        log("   g1 (genus 3, M_0): graph 2 = %s, graph 5 = %s;  genus 2: 1/750, -13/8000" % (rec.fmt(g1_3[2]), rec.fmt(g1_3[5])))
        step4(frame, st, state, log, check, rec, J3, g1_3, J2, g1_2)
    if "5" in steps:
        log("-" * 110); log("(5) u-dependent passes on class %s, placement %s: M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3)" % (a.cls, a.placement))
        J3 = {int(i): v for i, v in state.get('J').items()}; A0 = {int(i): v for i, v in state.get('A_M0').items()}
        M0 = st['M0']; Ns = []
        A = A0
        for order in (1, 2, 3):
            N, info = next_N(R, J3, A, order, a.cls, a.placement, M0, log, check, rec, "order %d" % order)
            check("order %d: equations for N_%d on class %s/%s consistent" % (order, order, a.cls, a.placement), N is not None)
            if N is None: break
            Ns.append(N)
            A = pass_u(frame, st, state, a.nproc, log, check, rec, Ns, a.placement, "%s_%s_N%d" % (a.cls, a.placement, order), tot_id, a.NN)
            stt = status(R, A)
            log("   after N_%d: non-vanishing u-orders <= 3 per edge graph: %s" % (order, {i: v for i, v in stt.items() if v}))
            check("after N_%d: 1-edge graphs vanish through u^%d" % (order, order), all(A[i][k] == 0 for i in (2, 5) for k in range(order + 1)))
            for E in (2, 3):
                bad = {i: [k for k in range(min(order + 1, 4)) if A[i][k]] for i in IDX_BY_E[E]}
                check("after N_%d: %d-edge graphs vanish through u^%d (orders <= %d are fixed by N_1..N_%d)" % (order, E, order, order, max(1, order + 1 - E)), not any(bad.values()), {i: v for i, v in bad.items() if v})
        if len(Ns) == 3:
            reg = [i for i in IDX_ALL if i != IDX_VERTEX and all(A[i][k] == 0 for k in range(4))]
            check("FINAL: all 41 edge graphs regular (O(u^4)) with M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3)", len(reg) == 41, "non-regular: %s" % [i for i in IDX_ALL if i != IDX_VERTEX and i not in reg])
            V = tc_expansion3(A[IDX_VERTEX], tc, R)
            log("   vertex graph in t_c: %s" % tc_fmt(rec, V, 5)); state.put('vertex_final_tc_%s_%s' % (a.cls, a.placement), V[:5])
            check("FINAL: vertex graph polar part = (1/25200) t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1", V[0] == R.of(Fr(1, 25200)) and V[1] == 0 and V[2] == 0 and V[3] == 0)
    log("total time %.0f s" % (time.time() - T0))
    log("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)
    LOG.close()


if __name__ == "__main__":
    main()
