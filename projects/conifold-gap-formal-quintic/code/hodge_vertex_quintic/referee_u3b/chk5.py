# -*- coding: utf-8 -*-
"""chk5.py -- referee NOTE: at the point B of Pi (c_1(0) = 0) the N_2 matrix has rank 5 (rows grad q_6 and grad q_17 vanish);
the two right-hand sides [u^3] G_6(N_1, 0), [u^3] G_17(N_1, 0) with true u-dependent vertices (one process, graphs 6 and 17 only)."""
import sys, os, json, time, argparse
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
from fq3_circsys import CircSys
from fq3_steps import frame_of
from fq3_par import job_hat
from fq3_m import sp_mats
S = CircSys(a.prime); p = S.p; cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); t0 = time.time()
frame = frame_of(p, (0, 0), 34, 4)
for name, v in (("B of Pi", cs['Pi_B']), ("extra[2] (B with another c5, c6)", cs['circ_b_extra'][2]), ("A + 2B", [(x + 2 * y) % p for x, y in zip(cs['Pi_A'], cs['Pi_B'])])):
    N1 = sp_mats(S.to_sp(v), p)
    for g in (6, 17):
        idx, coefs, info = job_hat((frame, ('M0exp', [N1], 'right'), g, 12))
        print("prime %s  %-34s graph %2d: grad q = 0: %-5s [u^0..u^2] = 0: %-5s [u^3] = 0: %-5s (eps-degree %s, %.0f s)"
              % (a.prime, name, g, not any(S.quad_grad(g, v)), all(c == 0 for c in coefs[:3]), coefs[3] == 0, info['maxdeg'], time.time() - t0), flush=True)
