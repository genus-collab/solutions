# -*- coding: utf-8 -*-
"""C1 / C2 at random FULL-class (N_1, N_2, N_3), vertex tensors frozen at M_0 (exact for the N_2-, N_3-dependence claimed):
 three-edge [u^3] independent of N_2, N_3; two-edge [u^3] affine in N_2 with linear part grad Q(N_1).n_2, independent of N_3;
 one-edge [u^2]: + J.n_2;  one-edge [u^3]: + J.n_3."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, random, time
from evlib import *
from fq3_quad import QuadForm, IDX_2EDGE
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
p = P2 if tag == "89" else P1; rng = random.Random(4242)
qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag))); g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag)))
Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
J3 = {2: g3['J']['2'], 5: g3['J']['5']}
dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
n1 = [rng.randrange(p) for _ in range(NSP)]; n2 = [rng.randrange(p) for _ in range(NSP)]; n2b = [rng.randrange(p) for _ in range(NSP)]; n3 = [rng.randrange(p) for _ in range(NSP)]
E1 = [2, 5]; E3 = [9, 10, 11, 12, 18, 19, 20, 32] + ([4] if "--g4" in sys.argv else [])
idxs = E1 + IDX_2EDGE + E3; t = time.time()
a, _ = evaluate(tag, [n1], idxs); print("first evaluation %.0f s" % (time.time() - t), flush=True)
b, _ = evaluate(tag, [n1, n2], idxs); c, _ = evaluate(tag, [n1, n2, n3], idxs); b2, _ = evaluate(tag, [n1, n2b], idxs)
n22 = [(x + y) % p for x, y in zip(n2, n2b)]; b3, _ = evaluate(tag, [n1, n22], idxs)
print("total %.0f s" % (time.time() - t))
ok = lambda s, v: print(("PASS " if v else "FAIL ") + s)
ok("three-edge graphs %s: [u^0..u^2] = 0 and [u^3] independent of N_2, N_3 (random full N_1, N_2, N_3)" % E3, all(a[i][k] == 0 for i in E3 for k in range(3)) and all(a[i][3] == b[i][3] == c[i][3] == b2[i][3] and a[i][3] != 0 for i in E3))
ok("two-edge graphs: [u^2] = Q_Gamma(N_1) (stored quadric) and independent of N_2, N_3", all(a[i][2] == Q[i].eval(n1) == b[i][2] == c[i][2] for i in IDX_2EDGE))
ok("two-edge graphs: [u^3](N_1, N_2) - [u^3](N_1, 0) = grad Q_Gamma(N_1) . n_2 (two random N_2 and their sum: exact affine-linearity)",
   all((x[i][3] - a[i][3]) % p == dot(Q[i].gradient(n1), nn) and dot(Q[i].gradient(n1), nn) != 0 for i in IDX_2EDGE for x, nn in ((b, n2), (b2, n2b), (b3, n22))))
ok("two-edge graphs: [u^3] independent of N_3", all(c[i][3] == b[i][3] for i in IDX_2EDGE))
ok("one-edge graphs: [u^1] = g1 + J.n_1", all(a[i][1] == (g3['A_M0'][str(i)][1] + dot(J3[i], n1)) % p for i in E1))
ok("one-edge graphs: [u^2](N_1, N_2) - [u^2](N_1, 0) = J . n_2", all((b[i][2] - a[i][2]) % p == dot(J3[i], n2) for i in E1))
ok("one-edge graphs: [u^3](N_1, N_2, N_3) - [u^3](N_1, N_2, 0) = J . n_3", all((c[i][3] - b[i][3]) % p == dot(J3[i], n3) for i in E1))
ok("one-edge graphs: [u^2] independent of N_3; [u^3] NOT affine in N_2 is allowed (quadratic terms): second difference %s" % [(b3[i][3] - b[i][3] - b2[i][3] + a[i][3]) % p != 0 for i in E1], all(c[i][2] == b[i][2] for i in E1))
