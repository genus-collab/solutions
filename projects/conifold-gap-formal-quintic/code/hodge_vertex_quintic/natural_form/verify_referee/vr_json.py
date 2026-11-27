# -*- coding: utf-8 -*-
"""vr_json.py -- VERIFIER: reads the referee's table G_fourier_NZ16_K12.json (already confirmed residue by residue against the project
rows at 2^89-1 by vr_proj_main.py) and prints the quoted numbers / counts."""
import json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(os.path.dirname(HERE), "referee", "G_fourier_NZ16_K12.json")))
T = {}
for key, v in d['G_fourier'].items():
    u_, z_, mm = key.split(); m, mp = (int(x) for x in mm.strip("()").split(","))
    T[(int(u_[1:]), int(z_[1:]), m, mp)] = Fr(v)
print("total entries %d, of which the identity at u^0 zhat^0: %d" % (len(T), sum(1 for k in T if k[:2] == (0, 0))))
print("nonzero Fourier components per u-order n = 1..11:", [sum(1 for k in T if k[0] == n) for n in range(1, 12)])
print("(3,3) of the zhat-linear part, u^1..u^7:", [str(T.get((n, 1, 3, 3), 0)) for n in range(1, 8)])
print("all negative for n = 2..11:", all(T[(n, 1, 3, 3)] < 0 for n in range(2, 12)))
for n in range(1, 7):
    dg = max(k[1] for k in T if k[0] == n)
    print("top coefficient of G_%d (zhat^%d), labels (a,b) = (-m-1,-m'-1): %s" % (n, dg, {((-k[2] - 1) % 5, (-k[3] - 1) % 5): str(v) for k, v in sorted(T.items()) if k[:2] == (n, dg)}))
