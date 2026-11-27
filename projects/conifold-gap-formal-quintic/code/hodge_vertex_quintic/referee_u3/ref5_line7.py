# -*- coding: utf-8 -*-
"""ref5_line7.py -- task (3): an independent line test in the chart c1_0 = 1 WITHOUT the a4-scaling.
The seven conditions c_4, c_9, c_10, c_11, c_12, c_19, q_3 are affine in (c4_1, c4_2, w): rows [A1, A2, S, P] with entries in (a, b, tau);
a solution needs rank <= 3, i.e. all 35 maximal minors vanish (valid on the conic a4 = 0 as well).  Second independent family of random
lines, different random combinations, different elimination (Euclidean resultants, Newton interpolation).
Run: py -I ref5_line7.py --prime 61|89 [--lines 4]"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *
from ref4_explore import subs_line, gcdl


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--lines", type=int, default=4); a = ap.parse_args()
    S = load(a.prime); p = S.p; P = load_polys(a.prime); rng = random.Random(777); FAIL = []; I = IRing(p); z = S.R.zeta
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref5_line7, p = 2^%s-1" % a.prime)
    D = chart_data(P, p); a4 = D['a4']; M7 = minors(D['D7'], 4, p)
    print("   D7-minors (35): total degrees %s, t-degrees %s, zero minors %d" % (sorted(set(mdeg(m) for m in M7.values())), sorted(set(mdeg(m, 2) for m in M7.values())), sum(1 for m in M7.values() if not m)))
    ok = True
    for _ in range(3):
        s = rng.randrange(p); v = [(A[j] + s * B[j]) % p for j in range(15)]; ok = ok and all(mevalR(m, [v[1], v[2], v[3]], I) == 0 for m in M7.values())
    check("all 35 D7-minors vanish on Pi", ok)
    lineL = lambda aa, bb: (1 + 2 * pow(z, 3, p) * aa + 2 * z * bb) % p
    summary = []
    for line in range(a.lines):
        lam = [[rng.randrange(p) for _ in range(35)] for _ in range(5)]; DD = combos(M7, lam, p)
        a0, a1, b0, b1 = [rng.randrange(p) for _ in range(4)]
        sPi = (-lineL(a0, b0)) * pow((lineL(a1, b1) - 1) % p, p - 2, p) % p
        a4l = muni({(e[0],): c for e, c in subs_line(a4, a0, a1, b0, b1, p).items()}, 0, p)
        fl = [subs_line(f, a0, a1, b0, b1, p) for f in DD]; dt = max(mdeg(f, 2) for f in fl)
        Rk = [bires(fl[0], fl[k], 0, 2, 3, p, dt, dt) for k in (1, 2, 3, 4)]; g = gcdl(Rk, p); mPi = umult(g, sPi, p)
        rest = list(g)
        for _ in range(mPi): rest = udiv_exact(rest, [(-sPi) % p, 1], p)
        g4 = ugcd(rest, a4l, p)
        lcs = gcdl([muni({(e[0],): c for e, c in mcoeffs(f, 2)[dt].items()}, 0, p) if mdeg(f, 2) == dt else [] for f in fl], p)
        print("   line %d: tau-degree %d, degrees of Res_tau(D_1, D_k) %s; common factor of degree %d = (s - s_Pi)^%d x rest of degree %d (on a4 = 0: %d); a4 at the two conic points is a root of the common factor: %s; gcd of the leading tau-coefficients on the line: degree %d" % (line, dt, [len(r) - 1 for r in Rk], len(g) - 1, mPi, len(rest) - 1, len(g4) - 1, len(ugcd(g, a4l, p)) > 1, len(lcs) - 1))
        summary.append((len(g) - 1, mPi, len(rest) - 1))
        # at s_Pi: the tau-gcd of all 35 minors (not only the combinations) -> exactly the tau of Pi, simple or multiple?
        aP, bP = (a0 + sPi * a1) % p, (b0 + sPi * b1) % p; tg = []
        for m in M7.values():
            if not m: continue
            ft = muni(msubs(m, {0: aP, 1: bP}, p), 2, p)
            if ft: tg = ugcd(tg, ft, p) if tg else umonic(ft, p)
        vP = [(A[j] + aP * B[j]) % p for j in range(15)]
        print("      at the point of Pi on the line: gcd over all 35 minors in tau has degree %d; its root = c2_1 of the point A + a B of Pi: %s (c1_2 agrees: %s)" % (len(tg) - 1, len(tg) == 2 and (-tg[0]) % p == vP[3], vP[2] == bP))
    check("on %d random lines (new family): the common zeros of the unscaled seven-row system are the point of Pi only (double root); nothing on the conic a4 = 0, nothing else" % a.lines,
          all(s == (2, 2, 0) for s in summary), summary)
    print("REF5 failures:", FAIL)


if __name__ == "__main__":
    main()
