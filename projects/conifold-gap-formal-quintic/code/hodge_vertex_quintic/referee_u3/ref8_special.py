# -*- coding: utf-8 -*-
"""ref8_special.py -- task (3): the special points of the chart c1_0 = 1 where the line L of Pi meets the conic a4 = 0 (the root of
multiplicity 21 of g in ref6_plane.py), and, more generally, a solver that treats w as an unknown:
all conditions except q_7 are jointly affine in u = (w, c4_1, c4_2, c5_0, c5_1, c5_2, c6_1, c6_2) with coefficients in (a, b, tau); q_7 is
quadratic in w.  At a point (a, b, tau) over a field K: solve the 11 x 8 linear system, then restrict q_7 to the solution space.
Run: py -I ref8_special.py --prime 61|89"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *
from ref_kpoly import *

UV = [4, 5, 6, 7, 8, 9, 10, 11]
LINC = ['J3_0', 'J3_1', 'q_3', 'q_6', 'q_8', 'c_4', 'c_9', 'c_10', 'c_11', 'c_12', 'c_19']


def affine_tables(P, p):
    """for the 11 jointly affine conditions: {name: (list of 8 coefficient polys in (a, b, t), constant poly)} on z = 1"""
    out = {}
    for k in LINC:
        f = msubs(P[k], {0: 1}, p); assert all(sum(e[i] for i in UV) <= 1 for e in f), "not jointly affine: " + k
        t3 = lambda g: {(e[1], e[2], e[3]): c for e, c in g.items()}
        out[k] = ([t3({e: c for e, c in f.items() if e[i]}) for i in UV], t3({e: c for e, c in f.items() if not any(e[i] for i in UV)}))
    return out


def ev3K(K, c, pt):
    s = K.zero
    for (i, j, k), v in c.items(): s = K.add(s, K.scal(v, K.mul(K.mul(kpow(K, pt[0], i), kpow(K, pt[1], j)), kpow(K, pt[2], k))))
    return s


def solve_point(K, Y, P, T, pt, p):
    """returns a description string and the list of full solutions found (15-vectors over K)"""
    rows = [[ev3K(K, c, pt) for c in T[k][0]] for k in LINC]; rhs = [K.sub(K.zero, ev3K(K, T[k][1], pt)) for k in LINC]
    sol, rk, free = solve_field(K, rows, rhs)
    if sol is None: return "linear system (11 conditions, 8 unknowns) INCONSISTENT (rank of the coefficient matrix %d)" % rk, []
    # nullspace directions
    ns = []
    for fcol in free:
        s2, _, _ = solve_field(K, rows + [[K.one if j == c else K.zero for j in range(8)] for c in free], [K.zero] * len(rows) + [K.one if c == fcol else K.zero for c in free])
        ns.append(s2)
    wdirs = [n for n in ns if not K.iszero(n[0])]
    def full(u): return Y.low(K, K.one, pt[0], pt[1], pt[2], u[0]) + list(u[1:])
    res = Y.residuals(K, full(sol)); bad = sorted(k for k, x in res.items() if not K.iszero(x))
    msg = "linear system consistent: rank %d, %d free directions (%d of them move w); residuals at the particular solution: %s" % (rk, len(free), len(wdirs), bad)
    sols = []
    if not bad: sols.append(full(sol))
    elif bad == ['q_7'] and wdirs:
        # q_7 along sol + s n (n moving w): quadratic in s; sample s = 0, 1, 2 -> coefficients; roots need a square root: report the discriminant
        n = wdirs[0]; q = [Y.residuals(K, full([K.add(x, K.scal(s_, y)) for x, y in zip(sol, n)]))['q_7'] for s_ in (0, 1, 2)]
        msg += "; q_7 along the w-direction is a non-trivial polynomial (values %s zero at s = 0, 1, 2): finitely many further points on this line" % [K.iszero(x) for x in q]
    return msg, sols


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); rng = random.Random(8); z = S.R.zeta; I = IRing(p)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    print("=" * 100); print("ref8_special, p = 2^%s-1" % a.prime)
    D = chart_data(P, p); M7 = minors(D['D7'], 4, p); T = affine_tables(P, p)
    # sanity of the solver on Pi and on a random point
    v = [(A[j] + 5 * B[j]) % p for j in range(15)]; msg, sols = solve_point(I, Y, P, T, [v[1], v[2], v[3]], p)
    print("   solver at the point A + 5 B of Pi: %s -> %d solution(s); reproduces the point up to the free direction: %s" % (msg, len(sols), bool(sols) and sols[0][:10] == v[:10]))
    msg, sols = solve_point(I, Y, P, T, [rng.randrange(p) for _ in range(3)], p); print("   solver at a random (a, b, tau): %s" % msg)
    # L and the conic: b = -(1 + 2 z^3 a) / (2 z)
    i2z = pow(2 * z, p - 2, p); bl = [(-i2z) % p, (-2 * pow(z, 3, p) * i2z) % p]                      # b(a) = bl[0] + bl[1] a
    qa = []
    for (i, j, _), c in D['a4'].items():
        t = [c]
        for _ in range(i): t = umul(t, [0, 1], p)
        for _ in range(j): t = umul(t, bl, p)
        qa = uadd(qa, t, p)
    fs = factor(qa, p, rng); disc = (qa[1] * qa[1] - 4 * qa[0] * qa[2]) % p if len(qa) == 3 else None
    print("   a4 restricted to L: degree %d in a, discriminant %s 0, irreducible factors of degrees %s, root multiplicities %s" % (len(qa) - 1, "=" if disc == 0 else "!=", fmt_fac(fs), [umult(qa, (-f[0]) % p, p) for f in fs if len(f) == 2]))
    total = 0
    for f in fs:
        if len(f) != 2: print("   intersection point in a field of degree %d: NOT examined" % (len(f) - 1)); continue
        a_s = (-f[0]) % p; b_s = (bl[0] + bl[1] * a_s) % p; vP = [(A[j] + a_s * B[j]) % p for j in range(15)]
        print("   special point (a, b) on L with a4 = 0: the point of Pi there has tau = c2_1 = %d; B4 at the Pi point: %s 0" % (vP[3], "=" if mevalR(D['B4'], [a_s, b_s, vP[3]], I) == 0 else "!="))
        Tg = []; nz = 0
        for m in M7.values():
            if not m: continue
            ft = muni(msubs(m, {0: a_s, 1: b_s}, p), 2, p)
            if ft: Tg = ugcd(Tg, ft, p) if Tg else umonic(ft, p); nz += 1
        ft_ = factor(Tg, p, rng)
        print("      tau-gcd of the 35 minors: degree %d (%d non-zero minors); irreducible factors of degrees %s; multiplicity of tau_Pi: %d" % (len(Tg) - 1, nz, fmt_fac(ft_), umult(Tg, vP[3], p)))
        B4t = muni(msubs(D['B4'], {0: a_s, 1: b_s}, p), 2, p); print("      B4(tau) at this (a, b): degree %d, gcd with the tau-gcd of degree %d" % (len(B4t) - 1, len(ugcd(B4t, Tg, p)) - 1))
        for h in ft_:
            K = KRing(h, p); tau = K.gen(); msg, sols = solve_point(K, Y, P, T, [K.const(a_s), K.const(b_s), tau], p)
            isPi = (K.d == 1 and (-h[0]) % p == vP[3])
            print("      tau in a field of degree %d (multiplicity %d in the gcd)%s: %s -> %d solution(s)%s" % (K.d, 0 if K.d > 1 else umult(Tg, (-h[0]) % p, p), " [= tau of Pi]" if isPi else "", msg, len(sols), "" if not sols else ("; it is the point of Pi (coordinates 0..9): %s" % ([(x[0] if x else 0) for x in sols[0][:10]] == vP[:10]))))
            total += len(sols) * K.d if not isPi else 0
    print("REF8: solutions over the special point(s) of L other than Pi: %d" % total)


if __name__ == "__main__":
    main()
