# -*- coding: utf-8 -*-
"""ref_full.py -- referee DIRECT, bonus: INDEPENDENT full pass (fq3_par.job_hat: u-dependent vertex tensors by 21-node interpolation, NOT
frozen at M_0) with M = M_0 exp(u N_1), u-precision K = 4, at points chosen by the referee.  Graphs with E <= 3 except graph 4 (too slow).
Expected at a solution of the N_1-conditions: one-edge graphs [u^0] = [u^1] = 0, two-edge graphs [u^0..u^2] = 0, three-edge [u^0..u^3] = 0.
Run: py -I ref_full.py --prime 61 --points Pi_rand,circ_b_extra0 [--nproc 4]
"""
import os, sys, json, argparse, random, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)

IDX_3 = [9, 10, 11, 12, 18, 19, 20, 32]          # graph 4 omitted
IDX_2 = [3, 6, 7, 8, 17]
IDX_1 = [2, 5]


def worker(args):
    name, p, n75, idx = args
    from fq3_steps import frame_of
    from fq3_par import job_hat
    from fq3_m import sp_mats
    frame4 = frame_of(p, (0, 0), 34, 4)
    i, coefs, info = job_hat((frame4, ('M0exp', [sp_mats(list(n75), p)], 'right'), idx, 21))
    return name, i, coefs, info


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--points", default="Pi_rand,circ_b_extra0"); ap.add_argument("--nproc", type=int, default=4)
    a = ap.parse_args()
    from fq3_circsys import CircSys
    S = CircSys(a.prime); p = S.p
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime)))
    A = cs['Pi_A']; B = cs['Pi_B']; rng = random.Random(987654321)
    x, y = rng.randrange(p), rng.randrange(p)
    allpts = {'Pi_rand': [(x * A[j] + y * B[j]) % p for j in range(15)], 'circ_b_extra0': list(cs['circ_b_extra'][0]), 'circ_b_extra1': list(cs['circ_b_extra'][1]),
              'circ_b_extra3': list(cs['circ_b_extra'][3]), 'Pi_B': list(B), 'control_random': [rng.randrange(p) for _ in range(15)]}
    pts = {k: allpts[k] for k in a.points.split(",")}
    for k, v in pts.items(): assert not any(S.residuals(v).values()) or k == 'control_random'
    jobs = [(name, p, S.to_sp(v), i) for name, v in pts.items() for i in IDX_1 + IDX_2 + IDX_3]
    out_path = os.path.join(HERE, "ref_full_%s_%s.json" % (a.prime, "+".join(pts))); out = {k: {} for k in pts}; t0 = time.time()
    with mp.Pool(processes=min(a.nproc, len(jobs))) as pool:
        for name, i, coefs, info in pool.imap_unordered(worker, jobs):
            out[name][str(i)] = coefs; json.dump(out, open(out_path, "w"), indent=1)
            E = info['E']; need = coefs[:E + 1]
            print("%-16s graph %2d (E=%d) %.0f s (job %.0f s, eps-degree %s): [u^0..u^%d] %s | [u^0..u^3] zero pattern %s"
                  % (name, i, E, time.time() - t0, info['time'], info['maxdeg'], E, "ALL ZERO" if not any(need) else "NON-ZERO", [int(c == 0) for c in coefs]), flush=True)
    for name in pts:
        ok = all(not any(out[name][str(i)][:E + 1]) for E, idxs in ((1, IDX_1), (2, IDX_2), (3, IDX_3)) for i in idxs)
        print("%s: one-edge O(u^2), two-edge O(u^3), three-edge (without graph 4) O(u^4): %s" % (name, ok), flush=True)
        if name == 'control_random':
            v = pts[name]
            print("   control: full-pass [u^3] of three-edge graphs = S.cub: %s; [u^2] of two-edge graphs = S.quad: %s"
                  % (all(out[name][str(i)][3] == S.cub(i, v) for i in IDX_3), all(out[name][str(i)][2] == S.quad(i, v) for i in IDX_2)), flush=True)
    print("done %.0f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
