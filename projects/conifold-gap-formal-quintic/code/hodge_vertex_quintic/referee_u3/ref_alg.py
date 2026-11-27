# -*- coding: utf-8 -*-
"""ref_alg.py -- referee DIRECT, algebraic part (own polynomial code; no graph sums):
 (1) test of fq3_cubic.cubic_points / cubic_from_values / CubicForm on a random cubic;
 (3) C3 from the stored forms: c_32 = kappa m_1^3 (exact), m_1 | c_18, c_20; c1_0 | c_11, c_19; further divisibility facts;
 (4) C6 from the CircSys data;
 (L) the four linear rows with own dot products at the stored points.
Run: py -I ref_alg.py --prime 61|89
"""
import os, sys, json, argparse, random, itertools
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
P = {"61": 2 ** 61 - 1, "89": 2 ** 89 - 1}
IDX_3 = [4, 9, 10, 11, 12, 18, 19, 20, 32]
IDX_2 = [3, 6, 7, 8, 17]


# ---------------------------------------------------------------- own sparse polynomials: {sorted index tuple: coef}
def padd(a, b, p, s=1):
    out = dict(a)
    for k, v in b.items():
        x = (out.get(k, 0) + s * v) % p
        if x: out[k] = x
        else: out.pop(k, None)
    return out


def pmul(a, b, p):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(sorted(k1 + k2)); x = (out.get(k, 0) + v1 * v2) % p
            if x: out[k] = x
            else: out.pop(k, None)
    return out


def lin(coefs, p):
    return {(j,): c % p for j, c in enumerate(coefs) if c % p}


def psubst(f, subs, p):
    """substitute variable j -> polynomial subs[j] (for j in subs), other variables unchanged"""
    out = {}
    for k, v in f.items():
        term = {(): v}
        for j in k:
            term = pmul(term, subs[j] if j in subs else {(j,): 1}, p)
        out = padd(out, term, p)
    return out


def peval(f, x, p):
    s = 0
    for k, v in f.items():
        t = v
        for j in k: t = t * x[j] % p
        s = (s + t) % p
    return s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args(); tag = a.prime; p = P[tag]
    rng = random.Random(4242); inv = lambda x: pow(x % p, p - 2, p)
    res = {}
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True); res[name] = bool(ok)
    # ------------------------------------------------------------ (1) polarisation test
    from fq3_cubic import cubic_points, cubic_from_values, CubicForm
    for nvar in (4, 10):
        keys = [k for k in itertools.combinations_with_replacement(range(nvar), 3)]
        f = {k: rng.randrange(p) for k in keys}
        pts = cubic_points(nvar)
        vals = {}
        for t, sup in pts:
            x = [0] * nvar
            for i, val in sup.items(): x[i] = val
            vals[t] = peval(f, x, p)
        C = cubic_from_values(vals, nvar, p)
        check("(1) cubic_from_values recovers a random homogeneous cubic in %d variables exactly (%d points, %d monomials)" % (nvar, len(pts), len(keys)),
              C.coef == {k: v for k, v in f.items() if v})
        x = [rng.randrange(p) for _ in range(nvar)]
        check("(1) CubicForm.eval = own evaluation at a random point (n = %d)" % nvar, C.eval(x) == peval(f, x, p))
        # gradient by own differentiation
        g = []
        for i in range(nvar):
            s = 0
            for k, v in f.items():
                c = k.count(i)
                if c:
                    rest = list(k); rest.remove(i); t = v * c
                    for j in rest: t = t * x[j] % p
                    s = (s + t) % p
            g.append(s)
        check("(1) CubicForm.gradient = own derivative (n = %d)" % nvar, C.gradient(x) == g)
        B = [[rng.randrange(p) for _ in range(nvar)] for _ in range(3)]; s3 = [rng.randrange(p) for _ in range(3)]
        y = [sum(s3[m] * B[m][i] for m in range(3)) % p for i in range(nvar)]
        check("(1) CubicForm.restrict consistent with eval (n = %d)" % nvar, C.restrict(B).eval(s3) == peval(f, y, p))
        # a NON-homogeneous cubic is NOT detected by the fit itself (needs the separate homogeneity / random-point checks)
    check("(1) number of polarisation points for 10 variables is 220 = number of cubic monomials", len(cubic_points(10)) == 220 and len(list(itertools.combinations_with_replacement(range(10), 3))) == 220)
    # ------------------------------------------------------------ data
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag)))
    names = cs['vis_names']
    C = {i: {tuple(sorted(int(t) for t in k.split(","))): v % p for k, v in cs['C'][str(i)].items() if v % p} for i in IDX_3}
    from fq3_circsys import CircSys
    S = CircSys(tag); m = [x % p for x in S.m1]
    order = lambda j: int(S.names[j][1])
    print("names:", S.names)
    print("m_1 support:", [S.names[j] for j in range(15) if m[j]])
    # weight bound from the stored supports
    for i in IDX_3:
        types = sorted(set(tuple(sorted(order(j) for j in k)) for k in C[i]))
        wmax = max(sum(order(j) - 1 for j in k) for k in C[i])
        print("   c_%d: %d monomials, max weight %d, order types %s" % (i, len(C[i]), wmax, types))
    check("C1-support: every stored monomial has weight Sum (k - 1) <= 3", all(sum(order(j) - 1 for j in k) <= 3 for i in IDX_3 for k in C[i]))
    # ------------------------------------------------------------ (3) C3
    m1 = lin(m[:10], p); m1c = pmul(pmul(m1, m1, p), m1, p)
    k0 = next(iter(m1c)); kappa = C[32].get(k0, 0) * inv(m1c[k0]) % p
    diff = padd(C[32], {k: kappa * v % p for k, v in m1c.items()}, p, -1)
    check("C3: c_32 = kappa m_1^3 as an exact polynomial identity (kappa != 0)", not diff and kappa != 0, "kappa = %d; c_32 has %d monomials, m_1^3 has %d" % (kappa, len(C[32]), len(m1c)))
    try:
        from fq3_rec import Rec
        rec = Rec(S.R); print("   kappa =", rec.fmt(kappa), "; m_1 =", ["%s: %s" % (S.names[j], rec.fmt(m[j])) for j in range(15) if m[j]])
    except Exception as e: print("   (Rec unavailable: %s)" % e)
    # substitution c2_2 -> -(m_1 - m4 c2_2)/m4
    sub = {4: {(j,): (-m[j] * inv(m[4])) % p for j in (1, 2, 3) if m[j]}}
    div_m1 = {}
    for i in IDX_3:
        r = psubst(C[i], sub, p); div_m1[i] = not r
    print("   cubics divisible by m_1:", [i for i in IDX_3 if div_m1[i]])
    check("C3: c_18 and c_20 vanish identically on m_1 = 0 (m_1 divides them)", div_m1[18] and div_m1[20])
    div_c10 = {i: all(0 in k for k in C[i]) for i in IDX_3}
    print("   cubics divisible by c1_0:", [i for i in IDX_3 if div_c10[i]])
    check("C3: c_11 and c_19 are divisible by c1_0", div_c10[11] and div_c10[19])
    # cubicb claim: on {c1_0 = 0, m_1 = 0, L2 rows} the cubics 11, 12, 18, 19, 20, 32 vanish identically
    sub0 = dict(sub); sub0[0] = {}
    on_locus = {i: not psubst(psubst(C[i], {0: {}}, p), sub, p) for i in IDX_3}
    print("   cubics vanishing identically on {c1_0 = 0, m_1 = 0} (without the genus-2 rows):", [i for i in IDX_3 if on_locus[i]])
    # with the genus-2 rows: eliminate c3_0, c3_1 by S.rref (pivots on columns 5, 6; free c3_2 and coords 0..4)
    rr = S.rref
    sub2 = {0: {}, 4: sub[4]}
    for r, jv in ((0, 5), (1, 6)):
        e = {}
        if rr[r][2] % p: e[(7,)] = (-rr[r][2]) % p
        for j in range(5):
            if rr[r][3 + j] % p: e[(j,)] = (-rr[r][3 + j]) % p
        sub2[jv] = e
    def on_full_locus(f):
        r = psubst(f, {5: sub2[5], 6: sub2[6]}, p)        # c3_0, c3_1 -> linear in (c1, c2, c3_2)
        r = psubst(r, {0: {}}, p); r = psubst(r, {4: sub[4]}, p)
        return not r
    onl = {i: on_full_locus(C[i]) for i in IDX_3}
    print("   cubics vanishing identically on {c1_0 = 0, m_1 = 0, genus-2 rows}:", [i for i in IDX_3 if onl[i]])
    check("cubicb structure: cubics 11, 12, 18, 19, 20, 32 vanish identically on {c1_0 = 0, m_1 = 0, genus-2 rows}", all(onl[i] for i in (11, 12, 18, 19, 20, 32)))
    res['_div_m1'] = [i for i in IDX_3 if div_m1[i]]; res['_div_c10'] = [i for i in IDX_3 if div_c10[i]]; res['_on_locus_noL2'] = [i for i in IDX_3 if on_locus[i]]
    res['_on_locus_L2'] = [i for i in IDX_3 if onl[i]]
    # ------------------------------------------------------------ (4) C6: the stratum c1 = 0
    ok_sup = True; bin_cub = {}
    for i in IDX_3:
        r = {k: v for k, v in C[i].items() if not (set(k) & {0, 1, 2})}
        ok_sup = ok_sup and all(set(k) <= {3, 4} for k in r); bin_cub[i] = r
    check("C6(a): at c1 = 0 every cubic depends on (c2_1, c2_2) only (exact support of the stored forms)", ok_sup, {i: len(bin_cub[i]) for i in IDX_3})
    v2 = [m[4], (-m[3]) % p]                             # the m_1 = 0 direction of c2 at c1 = 0
    x = [0] * 10; x[3], x[4] = v2
    vals = {i: peval(C[i], x, p) for i in IDX_3}
    check("C6(b): on c1 = 0, m_1 = 0 (c2 = s v2) some cubic is non-zero at s = 1, so c2 = 0", any(vals.values()), "non-zero for graphs %s" % [i for i in IDX_3 if vals[i]])
    check("C6(b'): c_32 itself vanishes at c1 = 0, c2 = v2 (consistency of m_1)", vals[32] == 0)
    # independent of m_1: common zeros of the nine binary cubics in (c2_1 : c2_2)
    def gcd_uni(f, g):
        def trim(a):
            a = [x % p for x in a]
            while a and a[-1] == 0: a.pop()
            return a
        f = trim(f); g = trim(g)
        while g:
            ig = inv(g[-1]); f = list(f)
            while len(f) >= len(g):
                c = f[-1] * ig % p; s = len(f) - len(g)
                for t, y in enumerate(g): f[s + t] = (f[s + t] - c * y) % p
                f = trim(f)
                if not f: break
            f, g = g, f
        return f
    unis = []
    for i in IDX_3:
        if bin_cub[i]:
            u = [0] * 4
            for k, v in bin_cub[i].items(): u[k.count(4)] = v       # c2_1 = 1, c2_2 = t
            unis.append((i, u))
    g = None
    for i, u in unis: g = list(u) if g is None else gcd_uni(g, u)
    all_inf = all(u[3] % p == 0 for i, u in unis)                    # common zero at (0 : 1)?
    check("C6(b''): the nine cubics restricted to c1 = 0 (binary cubics in c2) have NO common zero in P^1 at all", g is not None and len(g) == 1 and not all_inf,
          "non-zero binary cubics for graphs %s; gcd degree %s; common zero at (0:1): %s" % ([i for i, u in unis], (len(g) - 1) if g else None, all_inf))
    # c1 = c2 = 0: genus-2 rows on c3
    L2 = S.L2
    check("C6(c0): the genus-2 rows involve the coordinates of orders <= 3 only", all(L2[r][j] == 0 for r in range(2) for j in range(8, 15)))
    r0 = [L2[0][j] for j in (5, 6, 7)]; r1 = [L2[1][j] for j in (5, 6, 7)]
    k3 = [(r0[1] * r1[2] - r0[2] * r1[1]) % p, (r0[2] * r1[0] - r0[0] * r1[2]) % p, (r0[0] * r1[1] - r0[1] * r1[0]) % p]
    check("C6(c): the genus-2 rows have rank 2 on the c3 block: at c1 = c2 = 0 they leave a line of c3", any(k3))
    hi = list(range(8, 15)); c3 = [5, 6, 7]
    pure = all(S.Q[i][a_][b_] == 0 for i in IDX_2 for a_ in range(5, 15) for b_ in hi)
    check("C6(d): at c1 = c2 = 0 every two-edge quadric depends on c3 only (matrix entries (c3|high) and (high|high) vanish exactly)", pure)
    sym = all(S.Q[i][a_][b_] == S.Q[i][b_][a_] for i in IDX_2 for a_ in range(15) for b_ in range(15)); check("C6(d'): stored quadric matrices symmetric", sym)
    qv = {i: sum(k3[a_] * S.Q[i][5 + a_][5 + b_] * k3[b_] for a_ in range(3) for b_ in range(3)) % p for i in IDX_2}
    check("C6(e): some two-edge quadric is non-zero on the c3-line, so c3 = 0", any(qv.values()), "non-zero for graphs %s" % [i for i in IDX_2 if qv[i]])
    # ------------------------------------------------------------ (L) the four linear rows at the stored points
    R_of = lambda fr: fr.numerator * inv(fr.denominator) % p
    rows = [g3['J2']['2'], g3['J2']['4'], g3['J']['2'], g3['J']['5']]
    rhs4 = [(-R_of(Fr(1, 750))) % p, (-R_of(Fr(-13, 8000))) % p, (-g3['A_M0']['2'][1]) % p, (-g3['A_M0']['5'][1]) % p]
    dot = lambda r, v: sum(x * y for x, y in zip(r, v)) % p
    nt = qs['n_taut']
    check("(L) N_1^taut satisfies the four O(u) rows with the right-hand sides of genus3_quadsolve3.py", [dot(r, nt) for r in rows] == rhs4)
    check("(L) rhs of the genus-3 rows = -(1463/600000), -(-232351/120000000) (note s.4)", rhs4[2] == R_of(Fr(1463, 600000)) and rhs4[3] == R_of(Fr(-232351, 120000000)),
          "(sign convention: J.n = -g1)")
    A = cs['Pi_A']; B = cs['Pi_B']; x, y = rng.randrange(p), rng.randrange(p)
    pts = {'Pi_A': A, 'Pi_B': B, 'Pi_A+2B': [(A[j] + 2 * B[j]) % p for j in range(15)], 'Pi_rand': [(x * A[j] + y * B[j]) % p for j in range(15)],
           'circ_b_sol0': cs['circ_b_solutions'][0]}
    for k, e in enumerate(cs['circ_b_extra']): pts['circ_b_extra%d' % k] = e
    for nm, v in pts.items():
        n = S.to_sp(v)
        own = [(nt[k] + sum(S.cols[j][k] * v[j] for j in range(15))) % p for k in range(75)]
        ok = [dot(r, n) for r in rows] == rhs4 and own == n
        rs = S.residuals(v); bad = [k for k, val in rs.items() if val]
        m1v = sum(m[j] * v[j] for j in range(15)) % p
        check("(L) %s: the four rows hold at N_1 = N_1^taut + Delta (own dot products); stored-form residuals all zero; m_1 = 0" % nm, ok and not bad and m1v == 0, "bad=%s m1=%s c1=%s" % (bad, m1v != 0, v[:3]))
    z = S.R.zeta; i2 = inv(2)
    check("C4: c1(A) = (1, 0, -zeta^4/2), c1(B) = (0, 1, -zeta^2)", A[:3] == [1, 0, (-i2 * pow(z, 4, p)) % p] and B[:3] == [0, 1, (-pow(z, 2, p)) % p])
    check("C4: low parts (orders 1..3) of A and B linearly independent", any((A[i] * B[j] - A[j] * B[i]) % p for i in range(8) for j in range(8)))
    sol = cs['circ_b_solutions'][0]
    check("C5: stored solution has c1 = (0, 1, b), b = stored genuine root", sol[0] == 0 and sol[1] == 1 and [sol[2]] == cs['U_b1_genuine_roots'])
    z2 = pow(z, 2, p)
    print("   circ_b_extra: c1_2 / c1_1 =", ["+zeta^2" if e[2] * inv(e[1]) % p == z2 else ("-zeta^2" if e[2] * inv(e[1]) % p == (-z2) % p else "other") for e in cs['circ_b_extra']])
    check("C5: the four extra points have c1_0 = 0 and c1_2/c1_1 = +-zeta^2", all(e[0] == 0 and e[2] * inv(e[1]) % p in (z2, (-z2) % p) for e in cs['circ_b_extra']))
    # are the extra points related to Pi (B has c1 = (0, 1, -zeta^2))?
    for k, e in enumerate(cs['circ_b_extra']):
        sc = e[1]; d = [(e[j] - sc * B[j]) % p for j in range(15)]
        print("   extra%d - c1_1 * Pi_B: zero coordinates %s (low part equal: %s)" % (k, [S.names[j] for j in range(15) if d[j] == 0], not any(d[:8])))
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ref_alg_%s.json" % tag), "w"), indent=1)
    print("FAILED:", [k for k, v in res.items() if v is False])


if __name__ == "__main__":
    main()
