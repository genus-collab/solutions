# -*- coding: utf-8 -*-
"""genus3_quad.py -- compute the O(u^2) quadrics Q_Gamma(N_1) of the five two-edge genus-3 graphs on the full sp class (75 coordinates)
by polarisation (fq3_quad), with checks against the stored full passes.  Run:  py -I genus3_quad.py --prime 89|61 [--nproc 24]
Writes g3_quad_<p>.json (raw values at the 2926 points, the coefficient data c0, b, A per graph) and appends to logs/genus3_quad_<p>.log.
The analysis (kernels, cone structure, ranks, sectors) is genus3_quadan.py.
"""
import os, sys, time, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import P1, P2
from fq3_steps import frame_of, State
from fq3_quad import job_quad, points, coefficients, IDX_2EDGE
from fq3_par import run_pool, setup
from fq3_m import sp_params, sp_mats, sp_add, sp_scale, expand, class_columns, NSP
from fq3_taut import taut_N1
from fq3_rec import Rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--nproc", type=int, default=24)
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1
    frame = frame_of(p, (0, 0), 34, 3)
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    LOG = open(os.path.join(HERE, "logs", "genus3_quad_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    st = setup(frame); R = st['R']; rec = Rec(R); T0 = time.time()
    log("=" * 110); log("genus3_quad  p = 2^%s-1  frame (0, 0)  K = 3  nproc = %d  %s" % (a.prime, a.nproc, time.strftime("%Y-%m-%d %H:%M")))
    pts = points(); log("   %d evaluation points on the 75-dimensional class, graphs %s" % (len(pts), IDX_2EDGE))
    t = time.time(); res = run_pool(job_quad, [(frame, v, IDX_2EDGE, tag) for tag, v in pts], a.nproc, log, label="quadrics")
    log("   evaluations: %.0f s" % (time.time() - t))
    vals = {i: {} for i in IDX_2EDGE}; zero01 = True
    for tag, out in res:
        for i in IDX_2EDGE:
            vals[i][tuple(tag)] = out[i][2]
            if out[i][0] != 0 or out[i][1] != 0: zero01 = False
    check("every evaluation: [u^0] = [u^1] = 0 for the two-edge graphs", zero01)
    Q = {i: coefficients(vals[i], p) for i in IDX_2EDGE}
    state = State(os.path.join(HERE, "g3_quad_%s.json" % a.prime))
    state.put('prime', p); state.put('points', [[list(tag), None] for tag, v in pts])
    state.put('Q', {str(i): dict(c0=Q[i].c0, b=Q[i].b, A=Q[i].A) for i in IDX_2EDGE})
    # checks against stored passes (g3_state): pass A (N = 0), the covariant family passes t = 0..3, the tautological control
    g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % a.prime)))
    A0 = g3['A_M0']
    check("Q_Gamma(0) = [u^2] of pass A (M_0) for the five graphs", all(Q[i].c0 == A0[str(i)][2] % p for i in IDX_2EDGE))
    M0 = st['M0']; cols, names = class_columns(R, 'cov', M0, 'right')
    sysr = g3['sys4_g3:cov:right']; x = sysr['x']
    from fq_mfit import fp_solve
    Jr = sysr['Jr']; rhs = sysr['rhs']; sol = fp_solve(Jr, rhs, p); ns = sol[1]
    N0 = sp_mats(expand(cols, x, p), p); Nf = sp_mats(expand(cols, ns[0], p), p)
    ok = True
    for tt in range(4):
        key = 'A_cov5_right_t%d' % tt
        if key not in g3: continue
        Nt = sp_add(N0, sp_scale(Nf, tt, p), p); n = sp_params(Nt, p)
        for i in IDX_2EDGE:
            if Q[i].eval(n) != g3[key][str(i)][2] % p: ok = False
    check("Q_Gamma(N_1(t)) = [u^2] of the stored full passes of step 5 (covariant family, t = 0..3, 21-node interpolation of the vertex tensors)", ok)
    N1t, chk = taut_N1(st); nt = sp_params(N1t, p); state.put('n_taut', nt)
    check("Q_Gamma(N_1^taut) = 0 for the five graphs (R' = Psi + O(u^2) there)", all(Q[i].eval(nt) == 0 for i in IDX_2EDGE), [rec.fmt(Q[i].eval(nt)) for i in IDX_2EDGE])
    check("cone structure: gradient of every Q_Gamma vanishes at N_1^taut (Q_Gamma is homogeneous quadratic in N_1 - N_1^taut)", all(not any(Q[i].gradient(nt)) for i in IDX_2EDGE))
    if 'A_taut_N1' in g3:
        check("stored tautological control pass: [u^2] = 0 for the five graphs", all(g3['A_taut_N1'][str(i)][2] == 0 for i in IDX_2EDGE))
    # a random point: evaluate directly and compare with the quadric (independent of the polarisation points)
    import random
    rng = random.Random(77); nr = [rng.randrange(p) for _ in range(NSP)]
    tag, out = job_quad((frame, nr, IDX_2EDGE, 'rand'))
    check("random N_1: direct evaluation = quadric (polarisation exact, degree <= 2 in N_1)", all(out[i][2] == Q[i].eval(nr) for i in IDX_2EDGE))
    log("   diagonal coefficients A_ii of Q_3 (first block, k = 1): %s" % [rec.fmt(Q[3].A[i][i]) for i in range(15)])
    log("total time %.0f s" % (time.time() - T0)); log("QUAD ALL CHECKS PASSED" if not FAILED else "QUAD FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
