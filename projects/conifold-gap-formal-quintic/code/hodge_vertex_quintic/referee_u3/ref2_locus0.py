# -*- coding: utf-8 -*-
"""ref2_locus0.py -- task (2): the locus c1_0 = 0, c1 != 0: independent eliminations.
Unscaled forms on P^2 (x : y : t) = (c1_1 : c1_2 : c2_1):  c_4 = a4 w + B4,  c_10 = a10 w + B10,  H = Hw w + H0;
   G0 = a4 B10 - a10 B4 (degree 5),   HH = a4 H0 - B4 Hw (degree 7)       [the script uses G = a4^2 G0, H = a4^4 HH after its scaling]
Eliminations:  (A) chart x = 1, eliminate tau -> U0(b) (compare with the stored U_b1);  (B) chart x = 1, eliminate b -> V(tau);
               (C) chart y = 1, eliminate tau -> U1(a);  (D) chart x = 1 after a random shear b -> b + lam tau (counts POINTS).
Run: py -I ref2_locus0.py --prime 61|89"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; P = load_polys(a.prime); F, Z = locus0_forms(P, p); rng = random.Random(77); FAIL = []; z2 = pow(S.R.zeta, 2, p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime)))
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref2_locus0, p = 2^%s-1" % a.prime)
    G0, HH, a4, B4, a10, B10, Hw, H0 = [F[k] for k in ('G0', 'HH', 'a4', 'B4', 'a10', 'B10', 'Hw', 'H0')]
    print("   degrees (total; in x, y, t): a4 %s, B4 %s, a10 %s, B10 %s, Hw %s, H0 %s, G0 %s, HH %s" % tuple((mdeg(f), [mdeg(f, i) for i in range(3)]) for f in (a4, B4, a10, B10, Hw, H0, G0, HH)))
    lin_p = {(0, 1, 0): 1, (1, 0, 0): z2}; lin_m = {(0, 1, 0): 1, (1, 0, 0): (-z2) % p}        # y + zeta^2 x ,  y - zeta^2 x
    def divides_line(f, sgn):
        """is f divisible by y - sgn zeta^2 x ?  <=> f(x, sgn zeta^2 x, t) = 0"""
        g = {}
        for (i, j, k), c in f.items():
            e = (i + j, 0, k); g[e] = (g.get(e, 0) + c * pow(sgn * z2, j, p)) % p
        return not any(g.values())
    check("G0 is divisible by (y + zeta^2 x) and HH by (y - zeta^2 x), and not the other way round", divides_line(G0, -1) and divides_line(HH, 1) and not divides_line(G0, 1) and not divides_line(HH, -1))
    # ---------- chart x = 1 : polynomials in (b, tau) = variables (1, 2)
    g = msubs(G0, {0: 1}, p); h = msubs(HH, {0: 1}, p); a4b = muni(msubs(a4, {0: 1}, p), 1, p); a10b = muni(msubs(a10, {0: 1}, p), 1, p)
    dgt, dht = mdeg(g, 2), mdeg(h, 2); dgb, dhb = mdeg(g, 1), mdeg(h, 1)
    print("   chart x = 1: G0 degrees (b %d, tau %d), HH degrees (b %d, tau %d); a4(b) degree %d with roots in F_p: %d; a10(b) degree %d; gcd(a4, a10) degree %d" % (dgb, dgt, dhb, dht, len(a4b) - 1, len(uroots(a4b, p, rng)), len(a10b) - 1, len(ugcd(a4b, a10b, p)) - 1))
    gc = [muni(c, 1, p) for c in mcoeffs(g, 2)]; hc = [muni(c, 1, p) for c in mcoeffs(h, 2)]
    def gcdl(fs):
        r = []
        for f in fs: r = ugcd(r, f, p) if r else umonic(f, p)
        return r
    contG = gcdl(gc); contH = gcdl(hc); lead = ugcd(gc[-1], hc[-1], p)
    print("   tau-contents: G0: %s (roots %s), HH: %s (roots %s); gcd of the leading tau-coefficients: degree %d, factors %s" % (len(contG) - 1, ["-zeta^2" if r == (-z2) % p else "zeta^2" if r == z2 else r for r in uroots(contG, p, rng)], len(contH) - 1, ["-zeta^2" if r == (-z2) % p else "zeta^2" if r == z2 else r for r in uroots(contH, p, rng)], len(lead) - 1, fmt_fac(factor(lead, p, rng))))
    check("content of G0 = (b + zeta^2), content of HH = (b - zeta^2) exactly", contG == [z2, 1] and contH == [(-z2) % p, 1])
    U0 = bires(g, h, 1, 2, 3, p); V = bires(g, h, 2, 1, 3, p)
    print("   (A) U0(b) = Res_tau(G0, HH): degree %d;   (B) V(tau) = Res_b(G0, HH): degree %d" % (len(U0) - 1, len(V) - 1))
    # compare with the script's U_b1 = Res_tau(a4^2 G0, a4^4 HH) = a4^(2*3 + 4*2) U0 up to a constant
    Ust = utrim(cs['U_b1'], p); t = list(U0)
    for _ in range(2 * dht + 4 * dgt): t = umul(t, a4b, p)
    check("(A) the stored U_b1 (degree %d) is proportional to a4^%d U0: the script's scaling only adds powers of a4" % (len(Ust) - 1, 2 * dht + 4 * dgt), umonic(Ust, p) == umonic(t, p))
    # genuine part: remove from the squarefree part the roots of a4 (w undetermined there) and of the common leading coefficient (tau = infinity,
    # i.e. the point x = y = 0 of P^2 where c1 = 0)
    sq = usqfree(U0, p); fs = factor(sq, p, rng)
    print("   U0: squarefree part of degree %d, irreducible factors of degrees %s; multiplicities of the F_p-roots in U0: %s" % (len(sq) - 1, fmt_fac(fs), {("zeta^2" if r == z2 else "-zeta^2" if r == (-z2) % p else r): umult(U0, r, p) for r in uroots(U0, p, rng)}))
    Ug = uremove(uremove(uremove(uremove(sq, a4b, p), lead, p), contG, p), contH, p); fg = factor(Ug, p, rng); rg = uroots(Ug, p, rng)
    print("   GENUINE part of U0 (a4, contents and common leading coefficient removed): degree %d, factor degrees %s, F_p-roots %d" % (len(Ug) - 1, fmt_fac(fg), len(rg)))
    check("(A) genuine degree 13 with exactly one F_p-root, equal to the stored solution's c1_2", len(Ug) - 1 == 13 and len(rg) == 1 and rg[0] == cs['circ_b_solutions'][0][2] and cs['U_b1_genuine_degree'] == 13)
    # is anything genuine lost by removing 'lead'?  at a root b of lead the two curves share the tangent line y = b x at (0:0:1); a finite common
    # root tau may exist in addition: test gcd of the tau-polynomials at each root of lead (over the root field)
    for f in factor(lead, p, rng):
        K = KRing(f, p); b = K.gen()
        gt = [mevalR({(e[1],): c for e, c in mcoeffs(g, 2)[k].items()}, [b], K) for k in range(dgt + 1)]
        ht = [mevalR({(e[1],): c for e, c in mcoeffs(h, 2)[k].items()}, [b], K) for k in range(dht + 1)]
        while gt and K.iszero(gt[-1]): gt.pop()
        while ht and K.iszero(ht[-1]): ht.pop()
        # resultant over K at the true degrees through the Sylvester determinant over K (fraction-free: use rank_field on the Sylvester matrix)
        m, n = len(gt) - 1, len(ht) - 1; M = []
        for i in range(n): M.append([K.zero] * i + gt[::-1] + [K.zero] * (n - 1 - i))
        for i in range(m): M.append([K.zero] * i + ht[::-1] + [K.zero] * (m - 1 - i))
        rk = rank_field(K, M) if M else 0
        a4v = mevalR({(i,): c for i, c in enumerate(a4b)}, [b], K)
        print("      root field of degree %d of the common leading coefficient: true tau-degrees (%d, %d), Sylvester rank %d of %d -> %s; a4 %s" % (K.d, m, n, rk, m + n, "NO finite common root" if rk == m + n else "COMMON ROOT", "= 0" if K.iszero(a4v) else "!= 0"))
        check("no genuine solution hidden at the root(s) (degree %d) of the common leading tau-coefficient" % K.d, rk == m + n)
    # ---------- (B) eliminate b: tau-values
    sqV = usqfree(V, p); print("   (B) V(tau): squarefree part degree %d, factor degrees %s" % (len(sqV) - 1, fmt_fac(factor(sqV, p, rng))))
    # tau-values of the genuine points: for each irreducible factor of Ug compute tau in the root field and its minimal polynomial degree
    pts = []; tau_ok = True; deg_sum = 0
    for f in fg:
        K = KRing(f, p); b = K.gen()
        gt = [mevalR({(e[1],): c for e, c in mcoeffs(g, 2)[k].items()}, [b], K) for k in range(dgt + 1)]
        ht = [mevalR({(e[1],): c for e, c in mcoeffs(h, 2)[k].items()}, [b], K) for k in range(dht + 1)]
        # Euclid over K: gcd of the two tau-polynomials
        A_, B_ = ht, gt
        def ktrim(v):
            v = list(v)
            while v and K.iszero(v[-1]): v.pop()
            return v
        A_, B_ = ktrim(A_), ktrim(B_)
        while B_:
            inv = K.inv(B_[-1]); r = list(A_)
            while len(r) >= len(B_):
                fct = K.mul(r[-1], inv); s = len(r) - len(B_)
                for i, x in enumerate(B_): r[s + i] = K.sub(r[s + i], K.mul(fct, x))
                r = ktrim(r)
            A_, B_ = B_, r
        dg_ = len(A_) - 1; tau = K.mul(K.sub(K.zero, A_[0]), K.inv(A_[1])) if dg_ == 1 else None
        vv = mevalR({(i,): c for i, c in enumerate(V)}, [tau], K) if tau is not None else None
        print("      genuine factor of degree %d: gcd_tau(G0, HH) has degree %d; V(tau) = 0 there: %s" % (K.d, dg_, vv is not None and K.iszero(vv)))
        tau_ok = tau_ok and dg_ == 1 and K.iszero(vv); deg_sum += K.d * max(dg_, 0); pts.append((f, tau))
    check("every genuine root b carries exactly ONE tau (gcd degree 1 over the root field) and that tau is a root of V: 13 genuine points", tau_ok and deg_sum == 13, deg_sum)
    json.dump({'Ug': Ug, 'factors': [[f, t] for f, t in pts], 'U0': U0, 'V': V}, open(os.path.join(HERE, "ref2_genuine_%s.json" % a.prime), "w"))
    # ---------- (C) chart y = 1
    g1 = msubs(G0, {1: 1}, p); h1 = msubs(HH, {1: 1}, p); a4a = muni(msubs(a4, {1: 1}, p), 0, p)
    U1 = bires(g1, h1, 0, 2, 3, p); g1c = [muni(c, 0, p) for c in mcoeffs(g1, 2)]; h1c = [muni(c, 0, p) for c in mcoeffs(h1, 2)]
    U1g = usqfree(U1, p)
    for f in (a4a, ugcd(g1c[-1], h1c[-1], p), gcdl(g1c), gcdl(h1c)): U1g = uremove(U1g, f, p)
    # reciprocal of Ug
    rec_ = umonic(list(reversed(Ug)), p)
    check("(C) chart y = 1 (c1 = (0, a, 1)): degree of U1 %d, genuine part of degree %d = the reciprocal polynomial of the genuine part of chart x = 1 (same 13 points, a = 1/b)" % (len(U1) - 1, len(U1g) - 1), U1g == rec_)
    print("      U1(0) %s 0, a4(a = 0) %s 0: the point c1 = (0, 0, 1) %s" % ("=" if ueval(U1, 0, p) == 0 else "!=", "=" if ueval(a4a, 0, p) == 0 else "!=", "carries NO solution" if ueval(U1, 0, p) and ueval(a4a, 0, p) else "NEEDS EXAMINATION"))
    check("c1 = (0, 0, 1): no solution (U1(0) != 0 and a4 != 0 there)", ueval(U1, 0, p) != 0 and ueval(a4a, 0, p) != 0)
    # ---------- (D) shear: count points
    lam = rng.randrange(p)
    def shear(f):
        out = {}
        for (i, j, k), c in f.items():                                  # b -> b + lam tau
            term = {(i, 0, k): c}
            for _ in range(j): term = mmul(term, {(0, 1, 0): 1, (0, 0, 1): lam}, p)
            out = madd(out, term, p)
        return out
    gs = shear(g); hs = shear(h); Us = bires(gs, hs, 1, 2, 3, p); sqs = usqfree(Us, p)
    print("   (D) sheared chart (b -> b + lam tau): tau-degrees (%d, %d), Res degree %d, squarefree degree %d, factor degrees %s" % (mdeg(gs, 2), mdeg(hs, 2), len(Us) - 1, len(sqs) - 1, fmt_fac(factor(sqs, p, rng))))
    print("REF2 failures:", FAIL)


if __name__ == "__main__":
    main()
