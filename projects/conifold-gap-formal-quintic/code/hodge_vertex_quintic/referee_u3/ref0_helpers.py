# -*- coding: utf-8 -*-
"""ref0_helpers.py -- task (4): the helpers of fq3_circsys.py (interp_simple, sylvester_res with formal degrees, pgcd) and of fq3_quadlib
(det_fp, roots_fp, pdiv, pmod) against independent implementations (ref_lib) and against first principles."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_sys import *
from fq3_circsys import interp_simple, sylvester_res, pgcd, ptrim, peval
import fq3_quadlib as L
from fq_fp import P1, P2

rng = random.Random(20261110); FAIL = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
    if not ok: FAIL.append(name)


def rpoly(d, p): return [rng.randrange(p) for _ in range(d)] + [rng.randrange(1, p)]


def from_roots(rs, p, lc=1):
    f = [lc]
    for r in rs: f = umul(f, [(-r) % p, 1], p)
    return f


for p in (P1, P2, 101):
    tag = "p=%d bits" % p.bit_length()
    # ---- interp_simple
    ok = True
    for trial in range(60):
        n = rng.randrange(1, 40); d = rng.randrange(0, n); f = rpoly(d, p) if rng.random() < 0.9 else []
        xs = sorted(set(rng.randrange(p) for _ in range(n))); n = len(xs)
        if len(f) > n: f = f[:n]
        ys = [ueval(f, x, p) for x in xs]; g = ptrim(interp_simple(xs, ys, p), p)
        ok = ok and g == utrim(f, p) and g == newton_interp(xs, ys, p)
    check("%s interp_simple reproduces random polynomials of degree < n (60 trials, arbitrary nodes) and agrees with Newton interpolation" % tag, ok)
    # interpolation when the true degree exceeds the node count: aliasing (documented pitfall, not a bug)
    f = rpoly(7, p); xs = list(range(1, 6)); g = ptrim(interp_simple(xs, [ueval(f, x, p) for x in xs], p), p)
    check("%s interp_simple with too few nodes returns the remainder mod the master polynomial (silent aliasing)" % tag, g == umod(f, from_roots(xs, p), p))
    # ---- sylvester_res
    ok = ok2 = ok3 = ok4 = True
    for trial in range(80):
        df = rng.randrange(1, 7); dg = rng.randrange(1, 7); mf = rng.randrange(0, df + 1); mg = rng.randrange(0, dg + 1)
        f = rpoly(mf, p); g = rpoly(mg, p)
        if rng.random() < 0.3:
            c = rpoly(rng.randrange(1, 3), p); f2 = umul(f, c, p); g2 = umul(g, c, p)
            if len(f2) - 1 <= df and len(g2) - 1 <= dg: f, g = f2, g2
        r1 = sylvester_res(f, g, df, dg, p); r2 = ures(f, g, df, dg, p)
        # first principles: my own Sylvester matrix (columns instead of rows) and my own determinant
        n = df + dg; M = [[0] * n for _ in range(n)]; fp_ = (f + [0] * (df + 1))[:df + 1]; gp_ = (g + [0] * (dg + 1))[:dg + 1]
        for i in range(dg):
            for j in range(df + 1): M[i][i + j] = fp_[df - j]
        for i in range(df):
            for j in range(dg + 1): M[dg + i][i + j] = gp_[dg - j]
        r3 = det_modp(M, p)
        ok = ok and r1 == r2 == r3
        common = len(ugcd(f, g, p)) > 1; bothdrop = (len(f) - 1 < df and len(g) - 1 < dg)
        ok2 = ok2 and ((r1 == 0) == (common or bothdrop))
    check("%s sylvester_res(f, g, df, dg) = Euclidean resultant with formal degrees = det of an independently built Sylvester matrix (80 trials, random degree drops)" % tag, ok)
    check("%s sylvester_res vanishes iff (common root) or (both leading coefficients vanish)" % tag, ok2)
    for trial in range(20):
        ra = [rng.randrange(p) for _ in range(rng.randrange(1, 5))]; rb = [rng.randrange(p) for _ in range(rng.randrange(1, 5))]
        la = rng.randrange(1, p); lb = rng.randrange(1, p); f = from_roots(ra, p, la); g = from_roots(rb, p, lb)
        prod = pow(la, len(rb), p) * pow(lb, len(ra), p) % p
        for x in ra:
            for y in rb: prod = prod * (x - y) % p
        ok3 = ok3 and sylvester_res(f, g, len(ra), len(rb), p) == prod
    check("%s sylvester_res = lc(f)^dg lc(g)^df prod (alpha_i - beta_j) on split polynomials" % tag, ok3)
    f = rpoly(5, p); g = rpoly(2, p)
    check("%s sylvester_res silently TRUNCATES an input longer than the formal degree (pitfall: caller must guarantee deg <= formal degree)" % tag,
          sylvester_res(f, g, 3, 2, p) == sylvester_res(f[:4], g, 3, 2, p))
    # ---- pgcd
    ok = True
    for trial in range(80):
        c = from_roots([rng.randrange(p) for _ in range(rng.randrange(0, 4))], p); a = umul(rpoly(rng.randrange(0, 6), p), c, p); b = umul(rpoly(rng.randrange(0, 6), p), c, p)
        if rng.random() < 0.15: a = []
        if rng.random() < 0.15: b = [0, 0]
        g1 = pgcd(a, b, p); g2 = ugcd(a, b, p)
        ok = ok and list(g1) == list(g2)
        if a and utrim(b, p): ok = ok and len(g1) >= len(c)
    check("%s pgcd = monic gcd (80 trials incl. zero inputs, planted common factors)" % tag, ok)
    # ---- det_fp, roots_fp, pdiv, pmod
    ok = True
    for trial in range(30):
        n = rng.randrange(1, 9); M = [[rng.randrange(p) for _ in range(n)] for _ in range(n)]
        if rng.random() < 0.3 and n > 1: M[-1] = [(3 * x + 5 * y) % p for x, y in zip(M[0], M[1 % n])]
        ok = ok and L.det_fp(M, p) == det_modp(M, p)
    check("%s det_fp agrees with an independent determinant (30 trials, some singular)" % tag, ok)
    if p > 200:
        ok = True
        for trial in range(6):
            rs = sorted(set(rng.randrange(p) for _ in range(rng.randrange(1, 5)))); irr = rpoly(rng.randrange(0, 6), p)
            f = umul(umul(from_roots(rs, p), from_roots(rs[:1], p), p), irr, p)             # a double root included
            r1 = L.roots_fp(f, p, rng); r2 = uroots(f, p, rng)
            ok = ok and r1 == r2 and set(rs) <= set(r1) and all(ueval(f, x, p) == 0 for x in r1)
        check("%s roots_fp = independent root finder (distinct roots, multiplicities ignored), planted roots found" % tag, ok)
        a = rpoly(9, p); b = rpoly(4, p); q, r = udivmod(a, b, p)
        check("%s pdiv / pmod = quotient / remainder" % tag, list(L.pdiv(umul(a, b, p), b, p)) == utrim(a, p) and [x for x in L.pmod(a, b, p)] in (r, r or [0]))
        fs = factor(umul(umul(from_roots([5, 7], p), [3, 0, 1], p), [2, 1, 0, 0, 1], p), p, rng)
        check("%s ref_lib.factor: product of the factors = squarefree part" % tag, sum(len(f) - 1 for f in fs) == 8, fmt_fac(fs))
print("REF0 failures:", FAIL)
