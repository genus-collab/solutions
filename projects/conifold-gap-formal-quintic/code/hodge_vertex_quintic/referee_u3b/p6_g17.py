# -*- coding: utf-8 -*-
"""p6_g17.py -- referee PASSES: at N_1-solutions where rows of the N_2 system vanish identically on the circulant class (graph 17, and
graph 6 at the point B of Pi; see p5_rows.py), solvability of N_2 is the extra condition [u^3] G_Gamma(N_1, N_2 = 0) = 0 for those graphs.
Evaluate only those graphs (21-node interpolation, N_1 only) at stored solutions not covered by the driver's runs.
Run: py -I p6_g17.py <61|89> name[,name..]   (names: pi_A pi_B pi_rand0 pi_rand1 extra2 extra3 star ctrl)"""
import os, sys, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_steps import frame_of
from fq3_par import job_hat
from fq3_m import sp_mats
from fq3_circsys import CircSys
TWO = (3, 6, 7, 8, 17)
tag = sys.argv[1]; names = sys.argv[2].split(","); p = P2 if tag == "89" else P1
S = CircSys(tag); cj = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
A = [x % p for x in cj['Pi_A']]; B = [x % p for x in cj['Pi_B']]; rng = random.Random(7 + int(tag))
pts = {'pi_A': A, 'pi_B': B, 'ctrl': [(a + 2 * b) % p for a, b in zip(A, B)], 'star': cj['circ_star']}
for k in range(3):
    x, y = rng.randrange(p), rng.randrange(p); pts['pi_rand%d' % k] = [(x * a + y * b) % p for a, b in zip(A, B)]
for k, v in enumerate(cj['circ_b_extra']): pts['extra%d' % k] = v
frame6 = frame_of(p, (0, 0), 34, 6); out = {}
for nm in names:
    v = [t % p for t in pts[nm]]; bad = [k for k, r in S.residuals(v).items() if r]
    zero = [i for i in TWO if not any(S.quad_grad(i, v))]; N1 = sp_mats(S.to_sp(v), p); o = {}
    for i in zero:
        t = time.time(); idx, coefs, info = job_hat((frame6, ('M0exp', [N1], 'right'), i, 21))
        o[i] = dict(nonzero_orders=[k for k in range(6) if coefs[k] % p], u3_is_zero=(coefs[3] % p == 0), sec=round(time.time() - t))
    out[nm] = dict(N1_conditions_violated=bad, zero_rows=zero, graphs=o)
    print("p=2^%s-1 %-9s: violated %s; identically zero polar rows: graphs %s; %s" % (tag, nm, bad, zero, o), flush=True)
    json.dump(out, open(os.path.join(HERE, "p6_g17_%s_%s.json" % (tag, "_".join(names))), "w"), indent=1)
