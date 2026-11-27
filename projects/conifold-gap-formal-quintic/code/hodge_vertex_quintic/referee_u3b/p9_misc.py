# -*- coding: utf-8 -*-
"""p9_misc.py -- referee PASSES: c_1 of the stored circpi point; pure-high kernel coordinates; counts of E per graph.  Run: py -I p9_misc.py"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_circsys import CircSys
from fq3_graphs import outer_graphs3
from p1_audit import circ_coords
G = outer_graphs3(); print("graphs by number of edges:", {E: [g['idx'] for g in G if g['E'] == E] for E in range(7)})
print("counts:", {E: sum(1 for g in G if g['E'] == E) for E in range(7)}, "gv of 17:", [g['gv'] for g in G if g['idx'] == 17], [g['ed'] for g in G if g['idx'] == 17],
      "gv of 6:", [(g['gv'], g['ed']) for g in G if g['idx'] == 6])
for tag in ("61", "89"):
    p = P2 if tag == "89" else P1; R = FpRing(p); S = CircSys(tag); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    d = json.load(open(os.path.join(PROJ, "g3_u3_%s_circpi.json" % tag)))
    dl, ok = circ_coords(S, [(d['n1'][k] - qs['n_taut'][k]) % p for k in range(75)])
    i2 = pow(2, p - 2, p); want = (-R.z(4) * i2 - 2 * R.z(2)) % p
    print("p=2^%s-1 circpi: c1 = (%s, %s, c1_2); c1_2 == -zeta^4/2 - 2 zeta^2 for the project's zeta: %s (any conjugate: %s)" % (
        tag, dl[0], dl[1], dl[2] == want, any(dl[2] == (-R.z(4 * m) * i2 - 2 * R.z(2 * m)) % p for m in range(1, 5))))
    d = json.load(open(os.path.join(PROJ, "g3_u3_%s_circ.json" % tag)))
    dl, ok = circ_coords(S, [(d['n1'][k] - qs['n_taut'][k]) % p for k in range(75)])
    sm = lambda x: x if x < p // 2 else x - p
    print("   circ: (c5_0, c5_1, c5_2, c6_1, c6_2) = %s; low coordinates all zero: %s; J3 . Delta = %s" % ([sm(x) for x in dl[10:]], not any(dl[:8]),
          [sum(S.J3[r][j] * dl[j] for j in range(15)) % p for r in range(2)]))
