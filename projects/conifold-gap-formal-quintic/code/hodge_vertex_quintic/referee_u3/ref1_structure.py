# -*- coding: utf-8 -*-
"""ref1_structure.py -- task (1): the whole circulant system on {m_1 = 0, genus-2 rows} as explicit polynomials in the 12 variables
   z = c1_0, x = c1_1, y = c1_2, t = c2_1, w = c3_2, c4_1, c4_2, c5_0, c5_1, c5_2, c6_1, c6_2       (indices 0..11)
then the structure on the locus z = 0 (exact, symbolic; no random sampling).  Writes ref_polys_<prime>.json for the later scripts.
Run: py -I ref1_structure.py --prime 61|89"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sys import *

VN = ['z', 'x', 'y', 't', 'w', 'c4_1', 'c4_2', 'c5_0', 'c5_1', 'c5_2', 'c6_1', 'c6_2']; NV = 12


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); R = MRing(NV, p); rng = random.Random(5); FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref1_structure, p = 2^%s-1" % a.prime)
    # homogeneity of the stored quadrics (CircSys keeps only the matrix A of the restricted QuadForm)
    from fq3_quad import QuadForm
    qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % a.prime)))
    hom = True
    for i in QUA:
        q = QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']).restrict(S.nt, S.cols)
        hom = hom and q.c0 == 0 and not any(q.b) and all(q.A[r][c] == q.A[c][r] for r in range(15) for c in range(15))
    check("the five two-edge quadrics restricted to the circulant class at N_1^taut have no constant and no linear part (homogeneous in Delta), symmetric matrices", hom)
    v = Y.low(R, R.var(0), R.var(1), R.var(2), R.var(3), R.var(4)) + [R.var(k) for k in range(5, 12)]
    P = Y.residuals(R, v)
    check("m_1 and the two genus-2 rows vanish identically on my parametrisation", not P['m1'] and not P['L2_0'] and not P['L2_1'])
    # cross-check with the project evaluator at random points
    ok = True
    for _ in range(5):
        vals = [rng.randrange(p) for _ in range(NV)]; num = S.low(vals[:3], vals[3], vals[4]) + vals[5:]; res = S.residuals(num)
        mine = {k: (list(msubs(f, dict(enumerate(vals)), p).values()) or [0])[0] for k, f in P.items()}
        ok = ok and all(mine[k] == res[k] for k in res)
        num2 = [ (list(msubs(f, dict(enumerate(vals)), p).values()) or [0])[0] for f in v]; ok = ok and num2 == [x % p for x in num]
    check("my symbolic polynomials agree with CircSys.residuals / CircSys.low at 5 random points (all 18 conditions)", ok)
    check("c_32 = kappa m_1^3: c_32 vanishes identically on m_1 = 0; so do c_18, c_20 (C3, on m_1 = 0 AND the genus-2 rows)", not P['c_32'] and not P['c_18'] and not P['c_20'],
          {k: len(P[k]) for k in ('c_32', 'c_18', 'c_20')})
    zvar = R.var(0)
    def divisible_by_z(f): return all(k[0] >= 1 for k in f)
    check("c_11 and c_19 are divisible by c1_0 on {m_1 = 0, L2}", divisible_by_z(P['c_11']) and divisible_by_z(P['c_19']))
    print("   general chart (z free): per condition: #terms, total degree, support, degree per variable")
    for k in sorted(P):
        f = P[k]
        if k in ('m1', 'L2_0', 'L2_1'): continue
        print("     %-5s terms %4d deg %2d  %s" % (k, len(f), mdeg(f), {VN[i]: mdeg(f, i) for i in msupport(f)}))
    # ---------------- the locus z = 0
    Z = {k: msubs(f, {0: 0}, p) for k, f in P.items()}
    print("   locus z = c1_0 = 0:")
    for k in sorted(Z):
        f = Z[k]
        if k in ('m1', 'L2_0', 'L2_1'): continue
        print("     %-5s terms %4d deg %2d  %s" % (k, len(f), mdeg(f), {VN[i]: mdeg(f, i) for i in msupport(f)}))
    zero = sorted(k for k, f in Z.items() if not f and k not in ('m1', 'L2_0', 'L2_1'))
    check("on z = 0 exactly the cubics 11, 12, 18, 19, 20, 32 and the quadrics 6, 17 vanish identically", zero == sorted(['c_11', 'c_12', 'c_18', 'c_19', 'c_20', 'c_32', 'q_6', 'q_17']), zero)
    low = {1, 2, 3, 4}; c4 = {5, 6}; hi = {7, 8, 9, 10, 11}
    check("c_4, c_10 involve (x, y, t, w) only", set(msupport(Z['c_4'])) <= low and set(msupport(Z['c_10'])) <= low)
    check("c_9, q_3, q_8 involve (x, y, t, w, c4) only and are affine in c4", all(set(msupport(Z[k])) <= low | c4 and all(e[5] + e[6] <= 1 for e in Z[k]) for k in ('c_9', 'q_3', 'q_8')))
    check("c_4, c_10 are affine in w; the w-coefficients involve (x, y) only", all(mdeg(Z[k], 4) == 1 and set(msupport(mcoeffs(Z[k], 4)[1])) <= {1, 2} for k in ('c_4', 'c_10')))
    used = {k: sorted(set(msupport(Z[k])) & hi) for k in Z if Z[k]}
    print("   high variables (c5, c6) present: %s" % {k: [VN[i] for i in u] for k, u in used.items()})
    check("(c5, c6) enter only q_7 and the two genus-3 rows", all(not u for k, u in used.items() if k not in ('q_7', 'J3_0', 'J3_1')))
    # (c5, c6) system: rows J3_0, J3_1, q_7 ; coefficient matrix 3 x 5 with polynomial entries
    def coef(f, i):
        cs = mcoeffs(f, i); assert len(cs) <= 2; return cs[1] if len(cs) == 2 else {}
    rows = [[coef(Z[k], i) for i in (7, 8, 9, 10, 11)] for k in ('J3_0', 'J3_1', 'q_7')]
    print("   (c5, c6)-coefficients: J3 rows constant: %s; q_7 row: supports %s, degrees %s" % (all(mdeg(e) <= 0 for r in rows[:2] for e in r), [[VN[i] for i in msupport(e)] for e in rows[2]], [mdeg(e) for e in rows[2]]))
    import itertools
    minors = [mdet([[rows[r][c] for c in T] for r in range(3)], p, NV) for T in itertools.combinations(range(5), 3)]
    nz = [m for m in minors if m]
    print("   3 x 3 minors of the (c5, c6) coefficient matrix: %d non-zero of 10; supports %s" % (len(nz), sorted(set(VN[i] for m in nz for i in msupport(m)))))
    # the minors are linear forms in (x, y): they vanish simultaneously only where all are proportional and zero
    lin = [[m.get(tuple(1 if j == i else 0 for j in range(NV)), 0) for i in (1, 2)] for m in nz]
    rk2 = any((l1[0] * l2[1] - l1[1] * l2[0]) % p for l1 in lin for l2 in lin)
    allhom1 = all(mdeg(m) == 1 and set(msupport(m)) <= {1, 2} for m in nz)
    check("the (c5, c6) system (J3_0, J3_1, q_7) has rank 3 at EVERY point of the locus with c1 != 0 (minors are linear forms in (x, y) spanning both)", allhom1 and rk2,
          "minors linear in (x,y): %s; two independent: %s" % (allhom1, rk2))
    if allhom1 and not rk2:
        l = lin[0]; print("      all minors proportional to %d x + %d y: rank drops on that line of c1" % (l[0], l[1]))
    # c4 system: rows c_9, q_3, q_8 -> H = det [A | value]
    A = [[coef(Z[k], 5), coef(Z[k], 6), msubs(Z[k], {5: 0, 6: 0}, p)] for k in ('c_9', 'q_3', 'q_8')]
    print("   c4-rows: coefficient degrees %s, value degrees %s; supports of the coefficients %s" % ([[mdeg(r[0]), mdeg(r[1])] for r in A], [mdeg(r[2]) for r in A], [[ [VN[i] for i in msupport(r[c])] for c in (0, 1)] for r in A]))
    H = mdet(A, p, NV)
    print("   H = det[A | value]: terms %d, degree %d, per variable %s" % (len(H), mdeg(H), {VN[i]: mdeg(H, i) for i in msupport(H)}))
    # does J3 / q_7 constrain c4?  (they are solved by (c5, c6) whatever c4 is, by the rank statement)
    out = {'VN': VN, 'P': {k: [[list(e), c] for e, c in f.items()] for k, f in P.items()}, 'H0': [[list(e), c] for e, c in H.items()]}
    json.dump(out, open(os.path.join(HERE, "ref_polys_%s.json" % a.prime), "w"))
    print("REF1 failures:", FAIL)


if __name__ == "__main__":
    main()
