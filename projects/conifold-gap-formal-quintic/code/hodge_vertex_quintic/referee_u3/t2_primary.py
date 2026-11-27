# -*- coding: utf-8 -*-
"""primary (no psi) correlators of the frozen vertex theory Lambda_0 = M_0.omega^ in the Fourier (flat) basis phi_m = Hhat^(m-1):
which components are non-zero?  Hypothesis (P): for stable (g, n) != (0, 3) only <H, .., H> (all m = 2) survives."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import itertools, time, json
from rlib import *
from fq3_rec import Rec
from fq3_mvertex import MVertex6, MConst
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
S, d, cs = load(tag); p = S.p; R = S.R; rec = Rec(R); z = R.z
M0 = MConst(R, d['M'], d['Minv'], check=True)
mv = MVertex6(R, M0)
NAME = {1: "1", 2: "H", 3: "H2", 4: "H3", 0: "H4"}
out = {}
for (g, n) in [(0, 3), (1, 1), (1, 2), (2, 1), (0, 4), (1, 3), (2, 2), (0, 5), (1, 4), (0, 6), (2, 3), (3, 1), (3, 2)]:
    t = time.time()
    try: T = mv.tensor(g, (0,) * n)
    except Exception as e:
        print("(%d,%d): not available: %s" % (g, n, repr(e)[:100])); continue
    nz = {}
    for ms in itertools.combinations_with_replacement(range(5), n):
        acc = 0
        for idx in itertools.product(range(5), repeat=n):
            v = T[sum(idx[l] * 5 ** l for l in range(n))]
            if v: acc = (acc + v * z(sum(ms[l] * idx[l] for l in range(n)))) % p
        if acc: nz[ms] = acc
    # symmetric tensor? (check one transposition for n >= 2)
    print("(g,n) = (%d,%d)  [%.0f s]: %d non-zero Fourier components (multisets):" % (g, n, time.time() - t, len(nz)), flush=True)
    for ms, v in sorted(nz.items()):
        ch = sum(m - 2 for m in ms) % 5
        print("     <%s>  charge %d   value %s" % (",".join(NAME[m] for m in ms), ch, rec.fmt(v)))
    out["%d,%d" % (g, n)] = {",".join(str(m) for m in ms): v for ms, v in nz.items()}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "primary_%s.json" % tag), "w"))
