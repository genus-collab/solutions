# -*- coding: utf-8 -*-
"""ref_post.py -- summary of the direct evaluations and comparison of the points used with the CURRENT stored data."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
I3 = [4, 9, 10, 11, 12, 18, 19, 20, 32]; I2 = [3, 6, 7, 8, 17]; I1 = [2, 5]
for tag, batches in (("61", "AB"), ("89", "C")):
    p = 2 ** int(tag) - 1; cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag))); A = cs['Pi_A']; B = cs['Pi_B']
    cur = {'Pi_A': A, 'Pi_B': B, 'Pi_A+2B': [(A[j] + 2 * B[j]) % p for j in range(15)], 'Pi_2A+7B': [(2 * A[j] + 7 * B[j]) % p for j in range(15)], 'circ_b_sol0': cs['circ_b_solutions'][0]}
    for k, e in enumerate(cs['circ_b_extra']): cur['circ_b_extra%d' % k] = e
    for b in batches:
        d = json.load(open(os.path.join(HERE, "ref_direct_%s_%s.json" % (tag, b))))
        for name, r in d.items():
            c = r['cubic_frame4']; q = r['quad_frame3']
            z3 = all(c[str(i)][k] == 0 for i in I3 for k in range(4)); z2 = all(q[str(i)][k] == 0 for i in I2 for k in range(3)) and all(c[str(i)][k] == 0 for i in I2 for k in range(3))
            z1 = all(q[str(i)][k] == 0 for i in I1 for k in range(2))
            nonnull = all(x is not None for i in I3 for x in c[str(i)]) and all(x is not None for i in I2 + I1 for x in q[str(i)])
            same = (cur[name] == r['v15']) if name in cur else None
            inPi = None
            if name == 'Pi_rand':
                # v = x A + y B with x = v[0], y = v[1]
                x, y = r['v15'][0], r['v15'][1]; inPi = all((x * A[j] + y * B[j] - r['v15'][j]) % p == 0 for j in range(15))
            print("p=2^%s-1 %-26s 3-edge [u^0..3] all zero: %-5s 2-edge [u^0..2] all zero: %-5s 1-edge [u^0,1] zero: %-5s (no None: %s) | point = current stored data: %s%s"
                  % (tag, name, z3, z2, z1, nonnull, same, "" if inPi is None else " in current Pi: %s" % inPi))
    print(tag, "current failure keys:", {k: cs[k] for k in cs if 'failed' in k})
