# -*- coding: utf-8 -*-
"""is Pi a stabiliser-type symmetry of the frozen vertex theory Lambda_0 = M_0.omega^ ?  Free energies and a few correlators of
exp(eps zhat D_1).Lambda_0 for D_1 in Pi (no (H, H) term) and for a control D_1 with an (H, H) term; one-point descendants of Lambda_0."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import itertools, time, random
from rlib import *
from fq3_rec import Rec
from fq3_mvertex import MVertex6, MConst
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
S, d, cs = load(tag); p = S.p; R = S.R; rec = Rec(R); z = R.z; rng = random.Random(8)
NAME = {1: "1", 2: "H", 3: "H2", 4: "H3", 0: "H4"}
def four(T, n):
    nz = {}
    for ms in itertools.product(range(5), repeat=n):
        acc = 0
        for idx in itertools.product(range(5), repeat=n):
            v = T[sum(idx[l] * 5 ** l for l in range(n))]
            if v: acc = (acc + v * z(sum(ms[l] * idx[l] for l in range(n)))) % p
        if acc: nz[ms] = acc
    return nz
def M_eps(D1, eps):
    D = zeros(); D[1] = [[x * eps % p for x in r] for r in D1]
    E = mat_series_exp(D, R, KMAX); Ei = mat_series_exp(scale(D, p - 1, p), R, KMAX)
    return MConst(R, mul(E, d['M'], p), mul(d['Minv'], Ei, p), check=True)
M0 = MConst(R, d['M'], d['Minv'], check=True); mv0 = MVertex6(R, M0)
t = time.time(); F2 = mv0.tensor(2, ())[0]; F3 = mv0.tensor(3, ())[0]
print("Lambda_0: F_2 = %s, F_3 = %s (%.0f s)" % (rec.fmt(F2), rec.fmt(F3), time.time() - t), flush=True)
# one- and two-point descendants of Lambda_0 (top degree = no M-class left is NOT required: all psi exponents with deficit >= 0)
for (g, dd) in [(1, (1,)), (2, (1,)), (2, (2,)), (2, (3,)), (2, (4,)), (1, (1, 0)), (1, (2, 0)), (1, (1, 1)), (2, (1, 0)), (2, (2, 0)), (2, (1, 1)), (3, (1,)), (3, (2,))]:
    try: nz = four(mv0.tensor(g, dd), len(dd))
    except Exception as e: print("   (%d, %s): %s" % (g, dd, repr(e)[:60])); continue
    print("   genus %d, psi exponents %s: non-zero %s" % (g, dd, {",".join(NAME[m] for m in ms): rec.fmt(v) for ms, v in sorted(nz.items())}), flush=True)
be, ga = rng.randrange(p), rng.randrange(p)
for nm, F in (("D_1 in Pi (beta, gamma random)", {(1, 2): be, (1, 4): be, (1, 0): ga, (1, 1): ga}), ("D_1 in Pi, (1, H2) propagator only", {(1, 2): 1, (1, 4): 1}), ("D_1 in Pi, (H3, H4) propagator only", {(1, 0): 1, (1, 1): 1}),
               ("control: (H, H) propagator only", {(1, 3): 1})):
    D1 = from_fourier(S, F)[1]
    for eps in (1, 2):
        mv = MVertex6(R, M_eps(D1, eps)); t = time.time()
        f2 = mv.tensor(2, ())[0]; f3 = mv.tensor(3, ())[0]
        h1 = four(mv.tensor(1, (0,)), 1); h2 = four(mv.tensor(2, (0,)), 1); h11 = four(mv.tensor(1, (0, 0)), 2)
        print("%s, eps = %d: F_2 unchanged %s, F_3 unchanged %s; <.>_1 = %s, <.>_2 = %s, <.,.>_1 = %s (%.0f s)" % (nm, eps, f2 == F2, f3 == F3,
              {NAME[m[0]]: rec.fmt(v) for m, v in h1.items()}, {NAME[m[0]]: rec.fmt(v) for m, v in h2.items()}, {",".join(NAME[x] for x in m): rec.fmt(v) for m, v in h11.items()}, time.time() - t), flush=True)
