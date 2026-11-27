# -*- coding: utf-8 -*-
"""dump M_0, M_0^-1, mu_0 and the first u-corrections of the colour form of the hatted rows (prime 61 or 89) to the scratch folder"""
import sys, os, json, time
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); sys.path.insert(0, PROJ)
HERE = os.path.dirname(os.path.abspath(__file__))
from fq_fp import P1, P2
from fq3_steps import frame_of
from fq3_par import setup
from fq3_taut import colour_series_u, taut_N1
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
p = P2 if tag == "89" else P1
t = time.time(); st = setup(frame_of(p, (0, 0), 34, 4)); print("setup %.0f s" % (time.time() - t), flush=True)
M0 = st['M0']; H = st['H']
out = dict(p=p, M=M0.M, Minv=M0.Minv, mu0=st['mu0'], zeta=st['R'].zeta)
out['C'] = [colour_series_u(H, k) for k in range(4)]
N1, chk = taut_N1(st); out['N1taut'] = N1; print(chk)
out['x'] = [H.x.coef(k) for k in range(4)]
json.dump(out, open(os.path.join(HERE, "m0_%s.json" % tag), "w"))
print("done %.0f s" % (time.time() - t))
