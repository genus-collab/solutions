# -*- coding: utf-8 -*-
"""ref3_lift.py -- task (2), second part: every candidate point of the locus c1_0 = 0, c1 != 0 lifted to a full 15-vector over its ROOT FIELD
(F_p or an extension F_p[x]/(f)) and verified against every condition; Jacobian ranks (isolatedness); the special fibres
b = +-zeta^2, a4 = 0; comparison with the stored circ_b_solutions, circ_b_extra, Pi_A, Pi_B.
Run: py -I ref3_lift.py --prime 61|89   (after ref1_structure.py and ref2_locus0.py)"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); F, Z = locus0_forms(P, p); rng = random.Random(99); FAIL = []; z2 = pow(S.R.zeta, 2, p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); gen = json.load(open(os.path.join(HERE, "ref2_genuine_%s.json" % a.prime)))
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref3_lift, p = 2^%s-1" % a.prime)
    def lift(K, b, tau, label):
        """point (x : y : t) = (1 : b : tau) with a4 != 0"""
        xyz = [K.one, b, tau]; a4v = mevalR(F['a4'], xyz, K); B4v = mevalR(F['B4'], xyz, K)
        if K.iszero(a4v): print("      %s: a4 = 0" % label); return None
        w = K.sub(K.zero, K.mul(B4v, K.inv(a4v))); low = [K.zero, K.one, b, tau, w]
        G0v = mevalR(F['G0'], xyz, K); HHv = mevalR(F['HH'], xyz, K)
        rkA = rank_field(K, [[mevalR(e, xyz, K) for e in r] for r in F['A']])
        v, rk, nfree, bad = full_lift(Y, P, low, K, p)
        jr = jac_rank(P, low + v[8:], K, p) if v is not None else None
        print("      %s: G0 = 0: %s, HH = 0: %s; rank of the c4-block %d; high system rank %d, free directions %d; non-zero residuals %s; Jacobian rank %s (12 variables)" % (label, K.iszero(G0v), K.iszero(HHv), rkA, rk, nfree, bad, jr))
        return dict(v=v, bad=bad, jr=jr, rkA=rkA, nfree=nfree)
    # ---------- (i) the 13 genuine points
    allok = True; njr = []; sol_fp = None
    for f, tau in gen['factors']:
        K = KRing(f, p); r = lift(K, K.gen(), utrim(tau, p), "genuine factor of degree %d" % K.d)
        allok = allok and r is not None and r['v'] is not None and not r['bad'] and r['rkA'] == 2; njr.append(r['jr'] if r else None)
        if K.d == 1: sol_fp = [(x[0] if x else 0) for x in r['v']]
    check("ALL 13 genuine points (root fields of degrees %s) lift to full solutions: every condition (m_1, 4 linear rows, 5 quadrics, 9 cubics) vanishes over the root field" % [len(f) - 1 for f, _ in gen['factors']], allok)
    check("Jacobian rank 9 = 12 - 3 at each genuine point (cone direction + 2 free high directions only): reduced and ISOLATED in the low part", all(j == 9 for j in njr), njr)
    st = cs['circ_b_solutions'][0]
    check("the F_p-rational genuine point reproduces the stored circ_b_solutions[0] on the coordinates 0..9 (c5, c6 differ at most by the free fibre)", sol_fp is not None and sol_fp[:10] == [x % p for x in st[:10]])
    I = IRing(p)
    check("the stored circ_b_solutions[0] satisfies all conditions with my evaluator", not any(Y.residuals(I, [x % p for x in st]).values()))
    # ---------- (ii) the special lines b = +-zeta^2
    g = msubs(F['G0'], {0: 1}, p); h = msubs(F['HH'], {0: 1}, p); found = []
    for nm, b0, pol in (("b = zeta^2 (HH = 0 identically; roots of G0)", z2, g), ("b = -zeta^2 (G0 = 0 identically; roots of HH)", (-z2) % p, h)):
        ft = muni(msubs(pol, {1: b0}, p), 2, p); fs = factor(ft, p, rng)
        print("   %s: tau-polynomial of degree %d, squarefree factors %s, root multiplicities %s" % (nm, len(ft) - 1, fmt_fac(fs), [umult(ft, (-f_[0]) % p, p) for f_ in fs if len(f_) == 2]))
        for f_ in fs:
            K = KRing(f_, p); r = lift(K, K.const(b0), K.gen(), "tau in a field of degree %d" % K.d)
            if r and r['v'] is not None and not r['bad']: found.append((b0, K.d, r))
    ex = [[x % p for x in v] for v in cs['circ_b_extra']]
    mine = sorted([(x[0] if x else 0) for x in r['v'][:10]] for b0, d, r in found if d == 1)
    def norm(v):
        inv = pow(v[1], p - 2, p); return [x * inv % p for x in v[:10]]
    check("special lines: %d solution points over the closure (%d over F_p), the F_p ones = the 4 stored circ_b_extra (coordinates 0..9, normalised c1_1 = 1)" % (sum(d for _, d, _ in found), len(mine)),
          mine == sorted(norm(v) for v in ex) and len(mine) == 4)
    B = [x % p for x in cs['Pi_B']]; A = [x % p for x in cs['Pi_A']]
    isB = [r for b0, d, r in found if d == 1 and [(x[0] if x else 0) for x in r['v'][:10]] == B[:10]]
    check("exactly one of the special points is the point B of Pi (at b = -zeta^2)", len(isB) == 1)
    print("      Jacobian ranks at the special points: %s  (B of Pi: %s)" % ([r['jr'] for _, _, r in found], isB[0]['jr'] if isB else None))
    # ---------- (iii) the fibres a4 = 0 (w free): need B4 = 0 and, for some w, a10 w + B10 = 0 and Hw w + H0 = 0
    a4b = muni(msubs(F['a4'], {0: 1}, p), 1, p); okA = True
    for f_ in factor(a4b, p, rng):
        K = KRing(f_, p); b = K.gen()
        def tpoly(form):
            cf = mcoeffs(msubs(form, {0: 1}, p), 2); return [mevalR({(e[1],): c for e, c in c_.items()}, [b], K) for c_ in cf]
        if K.d != 1: print("      a4 has an irreducible factor of degree %d: fibre not examined here" % K.d); okA = False; continue
        tp = lambda form: utrim([(x[0] if x else 0) for x in tpoly(form)], p)
        B4t, B10t, Hwt, H0t = tp(F['B4']), tp(F['B10']), tp(F['Hw']), tp(F['H0']); a10v = ueval(muni(msubs(F['a10'], {0: 1}, p), 1, p), (-f_[0]) % p, p)
        Rt = usub(uscale(H0t, a10v, p), umul(B10t, Hwt, p), p); gd = ugcd(B4t, Rt, p)
        print("      a4 = 0 at b = %d: B4(tau) degree %d, a10 %s 0, residual condition a10 H0 - B10 Hw of degree %d, gcd degree %d" % ((-f_[0]) % p, len(B4t) - 1, "!=" if a10v else "=", len(Rt) - 1, len(gd) - 1))
        okA = okA and len(gd) == 1 and a10v != 0
    check("no solution on the two fibres a4 = 0 (B4 and the eliminated (c_10, H) condition have no common root; a10 != 0 there)", okA)
    # ---------- (iv) generic points of Pi: Jacobian rank (is Pi a reduced component of the expected dimension?)
    def to12(v):                                                    # 15-vector -> (z, x, y, t, w, high7)
        return [v[0], v[1], v[2], v[3], v[7]] + list(v[8:])
    jr = []; okPi = True
    for _ in range(4):
        s, t = rng.randrange(p), rng.randrange(p); v = [(s * A[j] + t * B[j]) % p for j in range(15)]
        okPi = okPi and not any(Y.residuals(I, v).values()); jr.append(jac_rank(P, to12(v), I, p))
    jA = jac_rank(P, to12(A), I, p); jB = jac_rank(P, to12(B), I, p)
    check("C4 (my evaluator): all conditions incl. m_1 vanish at 4 random points of Pi = span(A, B)", okPi)
    print("   Jacobian rank at 4 random points of Pi: %s; at A: %d; at B: %d   (rank 9 <=> local dimension 3 = dim Pi + 1 invisible direction, reduced)" % (jr, jA, jB))
    check("Pi is a reduced component: Jacobian rank 9 at random points of Pi (tangent space = Pi + the invisible direction)", all(j == 9 for j in jr), jr)
    print("REF3 failures:", FAIL)


if __name__ == "__main__":
    main()
