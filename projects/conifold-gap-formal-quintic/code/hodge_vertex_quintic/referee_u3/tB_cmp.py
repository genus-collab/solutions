# -*- coding: utf-8 -*-
"""C1 / C2 with TRUE u-dependent vertex tensors (UVertex) at a random full-class (N_1, N_2, N_3): compare the three configurations
a = (N_1), b = (N_1, N_2), c = (N_1, N_2, N_3) and the frozen-vertex evaluation"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, random
from evlib import *
from fq3_quad import QuadForm, IDX_2EDGE
HERE = os.path.dirname(os.path.abspath(__file__)); tag = "61"; p = P1; rng = random.Random(4242)
n1 = [rng.randrange(p) for _ in range(NSP)]; n2 = [rng.randrange(p) for _ in range(NSP)]; n2b = [rng.randrange(p) for _ in range(NSP)]; n3 = [rng.randrange(p) for _ in range(NSP)]
qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag))); g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag)))
Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
J3 = {2: g3['J']['2'], 5: g3['J']['5']}; dot = lambda r, v: sum(x * y for x, y in zip(r, v)) % p
L = {c: json.load(open(os.path.join(HERE, "tB_%s.json" % c))) for c in "abc"}
a, b, c = [{int(i): v for i, v in L[x]['vals'].items()} for x in "abc"]
print("nodes %s, maximal eps-degrees found %s (bounds: 5, 10, 15 for vertices of deficit <= 5)" % ([L[x]['NN'] for x in "abc"], [L[x]['maxdeg'] for x in "abc"]))
E1 = [2, 5]; E2 = [i for i in IDX_2EDGE if i in a]; E3 = [9, 10, 11, 12, 18, 19, 20, 32]
fz, _ = evaluate(tag, [n1], E1 + E2 + E3)
ok = lambda s, v: print(("PASS " if v else "FAIL ") + s)
ok("C1 three-edge %s: [u^3] with u-dependent vertices is the same for (N_1), (N_1, N_2), (N_1, N_2, N_3) and equals the frozen-vertex value (non-zero)" % E3, all(a[i][3] == b[i][3] == c[i][3] == fz[i][3] != 0 and a[i][:3] == [0, 0, 0] for i in E3))
ok("C2 two-edge %s: [u^2] = stored quadric Q(N_1) in all three configurations" % E2, all(a[i][2] == b[i][2] == c[i][2] == Q[i].eval(n1) for i in E2))
ok("C2 two-edge %s: [u^3](N_1, N_2) - [u^3](N_1, 0) = grad Q(N_1) . n_2 with u-dependent vertices; N_3 does not enter" % E2, all((b[i][3] - a[i][3]) % p == dot(Q[i].gradient(n1), n2) and c[i][3] == b[i][3] for i in E2))
print("   the frozen-vertex [u^3] of the two-edge graphs differs from the true one (vertex u-dependence and leaves contribute, N_2-independently): %s" % [a[i][3] != fz[i][3] for i in E2])
ok("C2 one-edge: [u^2](N_1, N_2) - [u^2](N_1, 0) = J . n_2;  [u^3](N_1, N_2, N_3) - [u^3](N_1, N_2, 0) = J . n_3; [u^2] independent of N_3", all((b[i][2] - a[i][2]) % p == dot(J3[i], n2) and (c[i][3] - b[i][3]) % p == dot(J3[i], n3) and c[i][2] == b[i][2] for i in E1))
if os.path.exists(os.path.join(HERE, "tB_a3.json")) and os.path.exists(os.path.join(HERE, "tB_b3.json")):
    a3 = json.load(open(os.path.join(HERE, "tB_a3.json"))); b3 = json.load(open(os.path.join(HERE, "tB_b3.json")))
    x = a3['vals']['3']; y = b3['vals']['3']
    print("graph 3: nodes %s, eps-degrees %s (bounds 4, 8)" % ([a3['NN'], b3['NN']], [a3['maxdeg'], b3['maxdeg']]))
    ok("C2 graph 3: [u^2] = Q_3(N_1) and [u^3](N_1, N_2) - [u^3](N_1, 0) = grad Q_3(N_1) . n_2 (u-dependent vertices)", x[2] == y[2] == Q[3].eval(n1) and (y[3] - x[3]) % p == dot(Q[3].gradient(n1), n2))
