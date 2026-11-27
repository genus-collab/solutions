# -*- coding: utf-8 -*-
"""ref10_direct.py -- DIRECT graph evaluation (fq3_cubic.job_cubic, fq3_quad.job_quad: no fitted forms) at N_1 = N_1^taut + Delta for
Delta = Delta* (ref9_newsol), for a point of Pi (control: must vanish) and for a random Delta (control: must NOT vanish):
[u^0..u^3] of the nine three-edge graphs and [u^0..u^2] of the five two-edge graphs.     Run: py -I ref10_direct.py --prime 61|89"""
import sys, os, random, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sys import *


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; rng = random.Random(10); t0 = time.time()
    from fq3_steps import frame_of
    from fq3_cubic import job_cubic, IDX_3EDGE
    from fq3_quad import job_quad, IDX_2EDGE
    d = json.load(open(os.path.join(HERE, "ref9_newsol_%s.json" % a.prime))); cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime)))
    A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    pts = [("Delta*", d['delta']), ("3 Delta* (homogeneity)", [3 * x % p for x in d['delta']]), ("point A + 7 B of Pi (control)", [(A[j] + 7 * B[j]) % p for j in range(15)]),
           ("random circulant Delta on m_1 = 0, L2 (control)", S.low([1, rng.randrange(p), rng.randrange(p)], rng.randrange(p), rng.randrange(p)) + [rng.randrange(p) for _ in range(7)])]
    print("=" * 100); print("ref10_direct, p = 2^%s-1" % a.prime, flush=True)
    frame4 = frame_of(p, (0, 0), 34, 4); frame3 = frame_of(p, (0, 0), 34, 3)
    for name, v in pts:
        n75 = S.to_sp(v)
        _, c = job_cubic((frame4, n75, IDX_3EDGE, name)); print("   %s: three-edge graphs done (%.0f s)" % (name, time.time() - t0), flush=True)
        _, q = job_quad((frame3, n75, IDX_2EDGE, name))
        zc = all(c[i][k] == 0 for i in IDX_3EDGE for k in range(4)); zq = all(q[i][k] == 0 for i in IDX_2EDGE for k in range(3))
        fit = all(c[i][3] == S.cub(i, v) for i in IDX_3EDGE) and all(q[i][2] == S.quad(i, v) for i in IDX_2EDGE)
        print("   %s: [u^0..u^3] of the 9 three-edge graphs all zero: %s; [u^0..u^2] of the 5 two-edge graphs all zero: %s; direct values = fitted forms (CircSys): %s   (%.0f s)" % (name, zc, zq, fit, time.time() - t0), flush=True)
        if not (zc and zq): print("      non-zero: cubics %s, quadrics %s" % ([i for i in IDX_3EDGE if any(c[i])], [i for i in IDX_2EDGE if any(q[i])]))
    print("REF10 done (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
