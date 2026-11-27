# -*- coding: utf-8 -*-
"""genus3_cubic.py -- the nine three-edge cubics c_Gamma(Delta) = [u^3] G_Gamma on the circulant class, with validation (genus 3, O(u^3)).

Computes (fq3_cubic.job_cubic, vertex tensors frozen at M_0, rows to O(u^3)) the homogeneous cubic forms of the nine three-edge graphs
4, 9, 10, 11, 12, 18, 19, 20, 32 in the circulant coordinates of Delta = N_1 - N_1^taut by polarisation on the 10 coordinates of zhat-orders
1..4 (220 points), and checks:
  (a) c_Gamma(0) = 0 (N_1 = N_1^taut) and [u^k] = 0 for k < 3 at every point;
  (b) the fitted forms reproduce direct evaluations at random circulant points, with and without components of orders 5, 6 (invisible);
  (c) homogeneity (c(2 Delta) = 8 c(Delta));
  (d) the frozen-vertex evaluation equals the [u^3] coefficients of the stored INDEPENDENT full passes (21-node vertex interpolation) at the
      three stored O(u^2) solutions (rational circulant, F_p circulant, full class);
  (e) full class: independence of the orders 5, 6 and affine-linearity in the order-4 coordinates at random points.
Output: g3_cubic_<p>.json (key 'C': {graph: {"i,j,k": coefficient}} on the 10 visible circulant coordinates 'vis_names'), log logs/genus3_cubic_<p>.log.
Run: py -I genus3_cubic.py --prime 89|61 [--nproc 12]
"""
import os, sys, json, argparse, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_steps import State, frame_of, GD
from fq3_par import run_pool
from fq3_m import NSP, class_columns, BLOCKS
from fq3_cubic import job_cubic, IDX_3EDGE, cubic_points, cubic_from_values, CubicForm
from fq3_rec import Rec


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--nproc", type=int, default=12)
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R)
    qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % a.prime))); nt = qs['n_taut']
    cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("=" * 110); log("genus3_cubic: the nine three-edge cubics on the circulant class, p = 2^%s-1" % a.prime)
    frame = frame_of(p, (0, 0), 34, 4)
    cols, names = class_columns(R, 'circ')
    vis = [j for j, nm in enumerate(names) if int(nm[1]) <= 4]; hid = [j for j in range(15) if j not in vis]
    assert len(vis) == 10
    def n_of(c15): return [(nt[v] + sum(cols[j][v] * c15[j] for j in range(15) if c15[j])) % p for v in range(NSP)]
    rng = random.Random(77)
    jobs = [(frame, list(nt), IDX_3EDGE, ('taut',))]
    for tag, sup in cubic_points(10):
        c = [0] * 15
        for i, val in sup.items(): c[vis[i]] = val
        jobs.append((frame, n_of(c), IDX_3EDGE, tag))
    rnd = []
    for k in range(4):
        c = [0] * 15
        for j in vis: c[j] = rng.randrange(p)
        c2 = list(c)
        for j in hid: c2[j] = rng.randrange(p)
        rnd.append((c, c2))
        jobs.append((frame, n_of(c), IDX_3EDGE, ('rv', k))); jobs.append((frame, n_of(c2), IDX_3EDGE, ('ra', k)))
    jobs.append((frame, n_of([2 * x % p for x in rnd[0][0]]), IDX_3EDGE, ('twice',)))
    for key in ('rational_solution_n', 'circ_solution_n', 'full_solution_n'): jobs.append((frame, qs[key], IDX_3EDGE, ('stored', key)))
    # full class structure: random full Delta on orders 1..4, the same plus orders 5, 6, and Delta_L + s Delta_4 (s = 0, 1, 2)
    lowf = [v for k, pairs, off in BLOCKS if k <= 3 for v in range(off, off + len(pairs))]
    o4 = [v for k, pairs, off in BLOCKS if k == 4 for v in range(off, off + len(pairs))]
    o56 = [v for k, pairs, off in BLOCKS if k >= 5 for v in range(off, off + len(pairs))]
    dL = [0] * NSP; d4 = [0] * NSP; d56 = [0] * NSP
    for v in lowf: dL[v] = rng.randrange(p)
    for v in o4: d4[v] = rng.randrange(p)
    for v in o56: d56[v] = rng.randrange(p)
    for s in range(3): jobs.append((frame, [(nt[v] + dL[v] + s * d4[v]) % p for v in range(NSP)], IDX_3EDGE, ('full4', s)))
    jobs.append((frame, [(nt[v] + dL[v] + d4[v] + d56[v]) % p for v in range(NSP)], IDX_3EDGE, ('full56',)))
    t = time.time(); res = dict(run_pool(job_cubic, jobs, a.nproc, log, label="cubic points")); log("   %d evaluations: %.0f s" % (len(jobs), time.time() - t))
    check("(a) [u^0], [u^1], [u^2] vanish for every three-edge graph at every point", all(res[tag][i][k] == 0 for tag in res for i in IDX_3EDGE for k in range(3)))
    check("(a) c_Gamma(0) = 0: the [u^3] coefficients vanish at N_1 = N_1^taut", all(res[('taut',)][i][3] == 0 for i in IDX_3EDGE))
    C = {}
    for i in IDX_3EDGE:
        C[i] = cubic_from_values({tag: res[tag][i][3] for tag, _ in cubic_points(10)}, 10, p)
    def ev(i, c15): return C[i].eval([c15[j] for j in vis])
    check("(b) fitted cubics = direct evaluation at 4 random circulant points (orders 1..4)", all(ev(i, rnd[k][0]) == res[('rv', k)][i][3] for k in range(4) for i in IDX_3EDGE))
    check("(b) the orders 5, 6 of Delta are invisible to the cubics (same values with random components of orders 5, 6)", all(res[('ra', k)][i][3] == res[('rv', k)][i][3] for k in range(4) for i in IDX_3EDGE))
    check("(c) homogeneity: c(2 Delta) = 8 c(Delta)", all(res[('twice',)][i][3] == 8 * res[('rv', 0)][i][3] % p for i in IDX_3EDGE))
    for key, vt in (('rational_solution_n', 'verify_table_rational'), ('circ_solution_n', 'verify_table_circ'), ('full_solution_n', 'verify_table_full')):
        tab = qs[vt]
        check("(d) frozen-vertex [u^3] = [u^3] of the stored independent full pass at %s, all nine graphs" % key, all(res[('stored', key)][i][3] == tab[str(i)][3] % p for i in IDX_3EDGE),
              {i: rec.fmt(tab[str(i)][3]) for i in IDX_3EDGE})
    f = {s: res[('full4', s)] for s in range(3)}
    check("(e) full class: affine-linear in the order-4 coordinates (second difference in s of c(Delta_L + s Delta_4) vanishes)", all((f[2][i][3] - 2 * f[1][i][3] + f[0][i][3]) % p == 0 for i in IDX_3EDGE))
    check("(e) full class: orders 5, 6 invisible", all(res[('full56',)][i][3] == f[1][i][3] for i in IDX_3EDGE))
    log("   order-4 dependence at the random full point (c(Delta_L + Delta_4) - c(Delta_L) non-zero for graphs): %s" % [i for i in IDX_3EDGE if (f[1][i][3] - f[0][i][3]) % p])
    vnames = [names[j] for j in vis]
    for i in IDX_3EDGE:
        sv = C[i].support_vars(); types = sorted(set(tuple(sorted(int(vnames[v][1]) for v in key)) for key in C[i].coef))
        log("   graph %2d (gv=%s, ed=%s): %3d non-zero monomials; variables %s; order types %s" % (i, GD[i]['gv'], GD[i]['ed'], len(C[i].coef), [vnames[v] for v in sv], types))
    cs.put('vis_names', vnames); cs.put('C', {str(i): {"%d,%d,%d" % k: v for k, v in C[i].coef.items()} for i in IDX_3EDGE})
    cs.put('checks_failed', FAILED)
    log("GENUS3_CUBIC ALL CHECKS PASSED" if not FAILED else "GENUS3_CUBIC FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
