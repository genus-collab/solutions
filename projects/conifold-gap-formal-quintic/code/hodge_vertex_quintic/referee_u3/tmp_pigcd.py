import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *
for tag in ("61", "89"):
    S = load(tag); p = S.p; P = load_polys(tag); rng = random.Random(1); I = IRing(p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    D = chart_data(P, p); M7 = minors(D['D7'], 4, p); ME = minors(D['E'], 3, p)
    for _ in range(3):
        s = rng.randrange(p); v = [(A[j] + s * B[j]) % p for j in range(15)]
        for nm, MM in (("D7", M7), ("E", ME)):
            tg = []
            for m in MM.values():
                if not m: continue
                ft = muni(msubs(m, {0: v[1], 1: v[2]}, p), 2, p)
                if ft: tg = ugcd(tg, ft, p) if tg else umonic(ft, p)
            print(tag, nm, "tau-gcd degree", len(tg) - 1, "roots", uroots(tg, p, rng), "mult of tau_Pi", umult(tg, v[3], p), "a4 != 0:", mevalR(D['a4'], [v[1], v[2], 0], I) != 0)
