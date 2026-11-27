# -*- coding: utf-8 -*-
"""p2_spot.py -- referee PASSES: (a) cache integrity: re-evaluate single graphs (one-edge 2, 5; optionally others) with the STORED
(n1, n2, n3) of each class through the project's job_hat (21-node interpolation) and compare all six u-coefficients with the stored
A1 / A2 / A3 tables;  (b) the N_2 / N_3 systems of genus3_u3.py tested on an independently known exact solution: the tautological
N_1^taut, N_2^taut, N_3^taut = [u^k] log(M_0^-1 R~(zhat; u)) must solve the systems built exactly as in the driver, with the right-hand
sides of the stored control pass A_taut_N1 (N_2) and of a fresh (N_1^taut, N_2^taut) evaluation of the graphs 2, 5 (N_3).
Read-only on the project.  Run: py -I p2_spot.py <61|89> [graphs, default 2,5]   -> p2_spot_<p>.json, log on stdout."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from fq3_steps import frame_of
from fq3_par import setup, job_hat
from fq3_m import sp_mats, sp_params, in_sp, class_columns, NSP, KMAX
from fq_symbol import mat_series_mul
from fq3_taut import colour_series_u


def main():
    tag = sys.argv[1]; p = P2 if tag == "89" else P1; R = FpRing(p)
    idxs = [int(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "2,5").split(",")]
    frame6 = frame_of(p, (0, 0), 34, 6); out = {}; t0 = time.time()
    def log(*x): print(" ".join(str(y) for y in x), flush=True)
    g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
    def ev(Ns, i):
        t = time.time(); idx, coefs, info = job_hat((frame6, ('M0exp', Ns, 'right'), i, 21))
        return [None if c is None else c % p for c in coefs], time.time() - t, info['maxdeg']
    # ---------------------------------------------------------------- (a) stored tables against fresh single-graph evaluations
    for cls in ("circ", "circpi", "circb"):
        d = json.load(open(os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls))))
        if 'A3' not in d: log(cls, "not finished"); continue
        N = [sp_mats(d[k], p) for k in ('n1', 'n2', 'n3')]
        for key, Ns in (('A3', N), ('A2', N[:2]), ('A1', N[:1])):
            for i in idxs:
                if key != 'A3' and i != idxs[-1]: continue          # A1, A2: one graph only (keep it light)
                c, dt, md = ev(Ns, i); st = [x % p for x in d[key][str(i)]]
                ok = (c == st); out["%s_%s_graph%d" % (cls, key, i)] = dict(equal=ok, nonzero_orders=[k for k in range(6) if c[k]], maxdeg=md, sec=round(dt))
                log("p=2^%s-1 %s %s graph %d: fresh == stored: %s; non-zero orders %s; eps-degree %s; %.0f s" % (tag, cls, key, i, ok, [k for k in range(6) if c[k]], md, dt))
    # ---------------------------------------------------------------- (b) the tautological solution through the driver's systems
    st = setup(frame6); H = st['H']; M0 = st['M0']
    C = [colour_series_u(H, k) for k in range(4)]
    X = [None] + [mat_series_mul(M0.Minv, C[k], p, KMAX) for k in (1, 2, 3)]
    mm = lambda A, B: mat_series_mul(A, B, p, KMAX)
    def lin(terms):
        o = [[[0] * 5 for _ in range(5)] for _ in range(KMAX + 1)]
        for c, A in terms:
            for k in range(KMAX + 1):
                for a in range(5):
                    for b in range(5): o[k][a][b] = (o[k][a][b] + c * A[k][a][b]) % p
        return o
    i2, i3 = pow(2, p - 2, p), pow(3, p - 2, p)
    N1t = X[1]; N2t = lin([(1, X[2]), (-i2, mm(X[1], X[1]))])
    N3t = lin([(1, X[3]), (-i2, mm(X[1], X[2])), (-i2, mm(X[2], X[1])), (i3, mm(mm(X[1], X[1]), X[1]))])
    out['taut_in_sp'] = [in_sp(N1t, p), in_sp(N2t, p), in_sp(N3t, p)]; log("tautological N_1, N_2, N_3 in sp:", out['taut_in_sp'])
    n1t, n2t, n3t = sp_params(N1t, p), sp_params(N2t, p), sp_params(N3t, p)
    nt = [x % p for x in qs['n_taut']]; out['n1t_eq_stored_n_taut'] = (n1t == nt); log("N_1^taut recomputed == stored n_taut:", n1t == nt)
    J3 = [g3['J']['2'], g3['J']['5']]; AT = g3['A_taut_N1']
    def gradQ(i, n):
        b, Am = qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']
        return [(b[a] + 2 * sum(Am[a][c] * n[c] for c in range(NSP))) % p for a in range(NSP)]
    rows = J3 + [gradQ(i, nt) for i in (3, 6, 7, 8, 17)]
    rhs = [(-AT['2'][2]) % p, (-AT['5'][2]) % p] + [(-AT[str(i)][3]) % p for i in (3, 6, 7, 8, 17)]
    out['taut_N2_system_solved_by_n2taut'] = [dot(r, n2t) == b for r, b in zip(rows, rhs)]
    out['taut_N2_rhs_nonzero'] = [bool(b) for b in rhs]
    log("N_2 system of the driver at N_1^taut (rhs from the stored pass A_taut_N1): solved by N_2^taut:", out['taut_N2_system_solved_by_n2taut'], "rhs non-zero:", out['taut_N2_rhs_nonzero'])
    cols, names = class_columns(R, 'circ')
    R15 = [[dot(r, cols[j]) for j in range(15)] for r in rows]
    from fq_mfit import fp_rref
    out['taut_N2_rank'] = len(fp_rref(R15, 15, p)[1]); log("   rank of the N_2 system at N_1^taut on the circulant class:", out['taut_N2_rank'])
    c2, dt, md = ev([N1t, N2t], 2); c5, dt5, md5 = ev([N1t, N2t], 5)
    out['taut_P2_graph2'] = [k for k in range(6) if c2[k]]; out['taut_P2_graph5'] = [k for k in range(6) if c5[k]]
    out['taut_N3_system_solved_by_n3taut'] = [dot(J3[0], n3t) == (-c2[3]) % p, dot(J3[1], n3t) == (-c5[3]) % p]
    log("fresh (N_1^taut, N_2^taut) evaluation: graph 2 non-zero orders %s, graph 5 %s (%.0f s, %.0f s); N_3 system solved by N_3^taut: %s"
        % (out['taut_P2_graph2'], out['taut_P2_graph5'], dt, dt5, out['taut_N3_system_solved_by_n3taut']))
    c2, dt, md = ev([N1t, N2t, N3t], 2); c5, dt5, md5 = ev([N1t, N2t, N3t], 5)
    out['taut_P3_graph2'] = [k for k in range(6) if c2[k]]; out['taut_P3_graph5'] = [k for k in range(6) if c5[k]]
    log("fresh (N_1^taut, N_2^taut, N_3^taut) evaluation: graph 2 non-zero orders %s, graph 5 %s" % (out['taut_P3_graph2'], out['taut_P3_graph5']))
    # ---------------------------------------------------------------- how far are the stored n2, n3 from the tautological ones?
    from fq_mfit import fp_solve
    Amat = [[cols[j][v] for j in range(15)] for v in range(NSP)]
    def coords(v):
        s = fp_solve(Amat, [x % p for x in v], p); return s[0]
    t2, t3 = coords(n2t), coords(n3t)
    out['n2taut_circulant'] = t2 is not None; out['n3taut_circulant'] = t3 is not None
    for cls in ("circ", "circpi", "circb"):
        d = json.load(open(os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls))))
        if 'n3' not in d: continue
        for nm, tt in (('n2', t2), ('n3', t3)):
            t = coords(d[nm])
            if t is None or tt is None: out["%s_%s_minus_taut" % (cls, nm)] = None; continue
            diff = [names[j] for j in range(15) if (t[j] - tt[j]) % p]
            out["%s_%s_minus_taut_nonzero_coords" % (cls, nm)] = diff
            log("%s: %s - %s^taut non-zero in the circulant coordinates %s" % (cls, nm, nm, diff))
    out['seconds'] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, "p2_spot_%s.json" % tag), "w"), indent=1)
    log("done in %d s" % out['seconds'])


if __name__ == "__main__":
    main()
