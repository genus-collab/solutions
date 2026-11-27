# -*- coding: utf-8 -*-
"""ref_direct.py -- referee DIRECT: frozen-vertex graph sums (fq3_cubic.job_cubic, fq3_quad.job_quad) at the stored circulant points.
Does not use the fitted cubic forms for the evaluation itself (CircSys is used only for to_sp and for the comparison columns).
Run: py -I ref_direct.py --prime 61 --batch A|B|C [--nproc 4]
Output: ref_direct_<prime>_<batch>.json in this directory.
"""
import os, sys, json, argparse, random, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)

IDX_3 = [4, 9, 10, 11, 12, 18, 19, 20, 32]
IDX_2 = [3, 6, 7, 8, 17]
IDX_1 = [2, 5]


def worker(args):
    name, p, n75 = args
    from fq3_steps import frame_of
    from fq3_cubic import job_cubic
    from fq3_quad import job_quad
    t0 = time.time()
    frame4 = frame_of(p, (0, 0), 34, 4); frame3 = frame_of(p, (0, 0), 34, 3)
    _, c = job_cubic((frame4, list(n75), IDX_3 + IDX_2 + IDX_1, name)); t1 = time.time()
    _, q = job_quad((frame3, list(n75), IDX_2 + IDX_1, name)); t2 = time.time()
    return name, {str(i): c[i] for i in c}, {str(i): q[i] for i in q}, (t1 - t0, t2 - t1)


def build_points(S, cs, batch):
    p = S.p; rng = random.Random(20261110 + (0 if S.tag == "61" else 1))
    A = cs['Pi_A']; B = cs['Pi_B']; pts = {}
    lin = lambda x, y: [(x * A[j] + y * B[j]) % p for j in range(15)]
    m = S.m1
    # kernel of the genus-2 rows on the c3 block (columns 5, 6, 7): cross product of the two rows
    r0 = [S.L2[0][j] for j in (5, 6, 7)]; r1 = [S.L2[1][j] for j in (5, 6, 7)]
    k3 = [(r0[1] * r1[2] - r0[2] * r1[1]) % p, (r0[2] * r1[0] - r0[0] * r1[2]) % p, (r0[0] * r1[1] - r0[1] * r1[0]) % p]
    if batch == 'A':
        pts['Pi_A'] = list(A); pts['Pi_B'] = list(B); pts['Pi_A+2B'] = lin(1, 2)
        x, y = rng.randrange(p), rng.randrange(p); pts['Pi_rand'] = lin(x, y)
        pts['circ_b_sol0'] = list(cs['circ_b_solutions'][0]); pts['circ_b_extra0'] = list(cs['circ_b_extra'][0])
        pts['control_random'] = [rng.randrange(p) for _ in range(15)]
        pts['taut'] = [0] * 15
    elif batch == 'B':
        for k in (1, 2, 3): pts['circ_b_extra%d' % k] = list(cs['circ_b_extra'][k])
        pts['c6_c2line'] = [0, 0, 0, m[4] % p, (-m[3]) % p] + [rng.randrange(p) for _ in range(10)]
        pts['c6_c3line'] = [0, 0, 0, 0, 0] + k3 + [rng.randrange(p) for _ in range(7)]
        # a point of Pi shifted by a random amount in ONE low coordinate: must fail (sensitivity control)
        v = lin(3, 5); v[3] = (v[3] + 1) % p; pts['Pi_3A+5B_perturbed_c2_1'] = v
        # Pi point plus the order-6 direction c6_2 (claimed invisible)
        v = lin(1, 1); v[14] = (v[14] + rng.randrange(p)) % p; pts['Pi_A+B_plus_c6_2'] = v
        v = lin(2, 7); pts['Pi_2A+7B'] = v
    elif batch == 'C':
        pts['Pi_A+2B'] = lin(1, 2); x, y = rng.randrange(p), rng.randrange(p); pts['Pi_rand'] = lin(x, y)
        pts['circ_b_sol0'] = list(cs['circ_b_solutions'][0]); pts['circ_b_extra0'] = list(cs['circ_b_extra'][0])
    return pts


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--batch", default="A"); ap.add_argument("--nproc", type=int, default=4)
    a = ap.parse_args()
    from fq3_circsys import CircSys
    S = CircSys(a.prime); p = S.p
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime)))
    pts = build_points(S, cs, a.batch)
    out_path = os.path.join(HERE, "ref_direct_%s_%s.json" % (a.prime, a.batch)); out = {}
    jobs = [(name, p, S.to_sp(v)) for name, v in pts.items()]
    t0 = time.time()
    with mp.Pool(processes=min(a.nproc, len(jobs))) as pool:
        for name, c, q, tm in pool.imap_unordered(worker, jobs):
            v = pts[name]
            pred_c = {str(i): S.cub(i, v) for i in IDX_3}; pred_q = {str(i): S.quad(i, v) for i in IDX_2}
            out[name] = dict(v15=v, cubic_frame4=c, quad_frame3=q, pred_cub=pred_c, pred_quad=pred_q, times=tm)
            json.dump(out, open(out_path, "w"), indent=1)
            nz3 = {i: [k for k in range(4) if c[str(i)][k]] for i in IDX_3 if any(c[str(i)])}
            nz2 = {i: [k for k in range(3) if q[str(i)][k]] for i in IDX_2 if any(q[str(i)])}
            nz1 = {i: [k for k in range(2) if q[str(i)][k]] for i in IDX_1 if any(q[str(i)][:2])}
            agree3 = all(c[str(i)][3] == pred_c[str(i)] for i in IDX_3); agree2 = all(q[str(i)][2] == pred_q[str(i)] for i in IDX_2)
            cons = all(c[str(i)][k] == q[str(i)][k] for i in IDX_2 + IDX_1 for k in range(3))
            print("%-28s %.0f s (%.0f+%.0f) | 3-edge nonzero [u^0..3]: %s | 2-edge nonzero [u^0..2]: %s | 1-edge nonzero [u^0,1]: %s | direct=S.cub: %s direct=S.quad: %s | K=4 vs K=3 agree: %s"
                  % (name, time.time() - t0, tm[0], tm[1], nz3, nz2, nz1, agree3, agree2, cons), flush=True)
    print("done %.0f s -> %s" % (time.time() - t0, out_path), flush=True)


if __name__ == "__main__":
    main()
