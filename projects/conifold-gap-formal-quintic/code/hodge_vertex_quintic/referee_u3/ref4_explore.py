# -*- coding: utf-8 -*-
"""ref4_explore.py -- task (3), first part: the chart c1_0 = 1.  Degrees of the E- and D-minors; validity on Pi; my own line tests
(E-system with my own random combinations; D-system = unscaled, valid on a4 = 0 too); degree accounting of the common factor.
Run: py -I ref4_explore.py --prime 61|89 [--lines 3]"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *


def subs_line(f, a0, a1, b0, b1, p):
    """(a, b) = (a0 + s a1, b0 + s b1): 3-variable poly (a, b, t) -> (s, -, t)"""
    out = {}; cache = {}
    def pw(base, n):
        k = (base, n)
        if k not in cache:
            r = {(0, 0, 0): 1}; lin = {(0, 0, 0): (a0 if base == 0 else b0) % p, (1, 0, 0): (a1 if base == 0 else b1) % p}
            for _ in range(n): r = mmul(r, lin, p)
            cache[k] = r
        return cache[k]
    for (i, j, k), c in f.items():
        out = madd(out, mmul(mmul(pw(0, i), pw(1, j), p), {(0, 0, k): c}, p), p)
    return out


def gcdl(fs, p):
    r = []
    for f in fs:
        f = utrim(f, p)
        if f: r = ugcd(r, f, p) if r else umonic(f, p)
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--lines", type=int, default=3); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); rng = random.Random(4242); FAIL = []; I = IRing(p); z = S.R.zeta
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref4_explore, p = 2^%s-1" % a.prime)
    D = chart_data(P, p); a4, B4 = D['a4'], D['B4']
    print("   a4: degree %d, support in (a, b) only: %s;  B4: degree %d (t-degree %d)" % (mdeg(a4), mdeg(a4, 2) <= 0, mdeg(B4), mdeg(B4, 2)))
    for k in CUB6 + ['q_3']:
        A1, A2, wc = D['rows'][k]
        print("     %-5s c4-coefficients: degrees (%d, %d), t-degrees (%d, %d); w-coefficients: degrees %s, t-degrees %s" % (k, mdeg(A1), mdeg(A2), mdeg(A1, 2), mdeg(A2, 2), [mdeg(c) for c in wc], [mdeg(c, 2) for c in wc]))
    check("the cubics 4, 9, 10, 11, 12, 19 are affine in w and their c4-coefficients depend on (a, b) only; q_3 is at most quadratic in w (in fact affine) with c4-coefficients independent of w",
          all(len(D['rows'][k][2]) <= 2 and mdeg(D['rows'][k][0], 2) <= 0 and mdeg(D['rows'][k][1], 2) <= 0 for k in CUB6) and len(D['rows']['q_3'][2]) <= 3)
    Em = minors(D['E'], 3, p); Dm = minors(D['D'], 4, p)
    print("   E-minors (20): total degrees %s, t-degrees %s, zero minors %d" % (sorted(set(mdeg(m) for m in Em.values())), sorted(set(mdeg(m, 2) for m in Em.values())), sum(1 for m in Em.values() if not m)))
    print("   D-minors (15): total degrees %s, t-degrees %s, zero minors %d" % (sorted(set(mdeg(m) for m in Dm.values())), sorted(set(mdeg(m, 2) for m in Dm.values())), sum(1 for m in Dm.values() if not m)))
    # relation between the two systems: an E-minor of three cubic rows (i, j, k) = the D-minor of (c_4, i, j, k) up to sign
    okrel = True
    for T, m in Em.items():
        if 5 in T: continue
        d = Dm[(0,) + tuple(t + 1 for t in T)]; okrel = okrel and (m == d or m == mscale(d, p - 1, p))
    check("identity: every E-minor of three cubic rows equals (up to sign) the D-minor of the same rows together with the c_4 row (the a4-scaling is exactly a Laplace expansion)", okrel)
    # Pi in this chart: the point A + s B has c1 = (1, s, c1_2), tau = c2_1
    ok = True
    for _ in range(3):
        s = rng.randrange(p); v = [(A[j] + s * B[j]) % p for j in range(15)]; pt = [v[1], v[2], v[3]]
        ok = ok and all(mevalR(m, pt, I) == 0 for m in Em.values()) and all(mevalR(m, pt, I) == 0 for m in Dm.values())
    check("all 20 E-minors and all 15 D-minors vanish at the points of Pi (necessary conditions hold on a known solution)", ok)
    lineL = lambda aa, bb: (1 + 2 * pow(z, 3, p) * aa + 2 * z * bb) % p
    lamE = [[rng.randrange(p) for _ in range(20)] for _ in range(4)]; lamD = [[rng.randrange(p) for _ in range(15)] for _ in range(4)]
    EE = combos(Em, lamE, p); DD = combos(Dm, lamD, p)
    for line in range(a.lines):
        a0, a1, b0, b1 = [rng.randrange(p) for _ in range(4)]
        # the point of Pi on this line: L(a0 + s a1, b0 + s b1) = 0
        sPi = (-lineL(a0, b0)) * pow((lineL(a1, b1) - 1) % p, p - 2, p) % p
        a4l = muni({(e[0],): c for e, c in subs_line(a4, a0, a1, b0, b1, p).items()}, 0, p)
        for nm, sysm in (("E", EE), ("D", DD)):
            fl = [subs_line(f, a0, a1, b0, b1, p) for f in sysm]; dt = max(mdeg(f, 2) for f in fl)
            Rk = [bires(fl[0], fl[k], 0, 2, 3, p, dt, dt) for k in (1, 2, 3)]; g = gcdl(Rk, p)
            mPi = umult(g, sPi, p); ga = list(g); m4 = 0
            while len(ugcd(ga, a4l, p)) > 1: h = ugcd(ga, a4l, p); ga = udiv_exact(ga, h, p); m4 += len(h) - 1
            rest = list(ga)
            for _ in range(mPi): rest = udiv_exact(rest, [(-sPi) % p, 1], p)
            # leading tau-coefficients: common zeros on the line (root at tau = infinity)
            lcs = gcdl([muni({(e[0],): c for e, c in mcoeffs(f, 2)[dt].items()}, 0, p) if mdeg(f, 2) == dt else [] for f in fl], p)
            print("   line %d, %s-system: tau-degree %d, degrees of Res_tau(f_1, f_k) %s, common factor degree %d = (s - s_Pi)^%d x (degree %d on a4 = 0) x (rest of degree %d, factors %s); gcd of the leading tau-coefficients: degree %d" % (line, nm, dt, [len(r) - 1 for r in Rk], len(g) - 1, mPi, m4, len(rest) - 1, fmt_fac(factor(rest, p, rng)) if len(rest) > 1 else [], len(lcs) - 1))
            if nm == "E": check("line %d, E-system (my combinations): common factor of degree 6 = (s - s_Pi)^2 x degree 4 on the conic a4 = 0, nothing else (script's claim)" % line, len(g) - 1 == 6 and mPi == 2 and m4 == 4 and len(rest) == 1)
            else: print("      (D-system = the six cubics only, no q_3: the extra factor of degree %d is a curve of solutions of the CUBICS ALONE; q_3 removes it)" % (len(rest) - 1))
    print("REF4 failures:", FAIL)


if __name__ == "__main__":
    main()
