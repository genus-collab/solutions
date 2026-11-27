# -*- coding: utf-8 -*-
"""ref9_newsol.py -- the additional solution Delta* found by ref8_special.py over the point where the line L of Pi touches the conic a4 = 0:
reconstruction, verification with the PROJECT evaluator CircSys.residuals, position relative to Pi, Jacobian rank, coordinates; rank of the
linear system along the curve of Pi (no hidden fibres); writes ref9_newsol_<prime>.json.     Run: py -I ref9_newsol.py --prime 61|89"""
import sys, os, random, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref8_special import *
from ref3b_pi import nullspace


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); rng = random.Random(9); z = S.R.zeta; I = IRing(p); FAIL = []
    from fq3_rec import Rec
    rec = Rec(S.R)
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime))); A = [x % p for x in cs['Pi_A']]; B = [x % p for x in cs['Pi_B']]
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref9_newsol, p = 2^%s-1" % a.prime)
    D = chart_data(P, p); M7 = minors(D['D7'], 4, p); T = affine_tables(P, p)
    i2z = pow(2 * z, p - 2, p); bl = [(-i2z) % p, (-2 * pow(z, 3, p) * i2z) % p]; qa = []
    for (i, j, _), c in D['a4'].items():
        t = [c]
        for _ in range(i): t = umul(t, [0, 1], p)
        for _ in range(j): t = umul(t, bl, p)
        qa = uadd(qa, t, p)
    check("the line L of Pi (c1_0 + 2 zeta^3 c1_1 + 2 zeta c1_2 = 0) is TANGENT to the conic a4 = 0 (double root)", len(qa) == 3 and (qa[1] * qa[1] - 4 * qa[0] * qa[2]) % p == 0)
    a_s = (-qa[1]) * pow(2 * qa[2], p - 2, p) % p; b_s = (bl[0] + bl[1] * a_s) % p; vP = [(A[j] + a_s * B[j]) % p for j in range(15)]
    Tg = []
    for m in M7.values():
        if not m: continue
        ft = muni(msubs(m, {0: a_s, 1: b_s}, p), 2, p)
        if ft: Tg = ugcd(Tg, ft, p) if Tg else umonic(ft, p)
    roots = uroots(Tg, p, rng); new = [r for r in roots if r != vP[3]]
    check("tau-gcd of the 35 minors at the tangency point: degree 4 = (tau - tau_Pi)^3 (tau - tau*)", len(Tg) == 5 and umult(Tg, vP[3], p) == 3 and len(new) == 1)
    msg, sols = solve_point(I, Y, P, T, [a_s, b_s, new[0]], p); v = sols[0]
    res = S.residuals(v); m1v = sum(S.m1[j] * v[j] for j in range(15)) % p
    check("Delta*: ALL conditions vanish with the PROJECT evaluator CircSys.residuals (2 genus-2 rows, 2 genus-3 rows, 5 quadrics, 9 cubics) and m_1 = 0", not any(res.values()) and m1v == 0, [k for k, x in res.items() if x])
    rk = rank_field(I, [A, B, v]); hk = nullspace([[S.J3[r][j] if j >= 8 else 0 for j in range(15)] for r in range(2)] + [[1 if j == i else 0 for j in range(15)] for i in range(8)], p)
    rk2 = rank_field(I, [A[:10], B[:10], v[:10]])
    check("Delta* does not lie in Pi, not even modulo the orders 5, 6 (rank of (A, B, Delta*) on the coordinates 0..9 is 3)", rk == 3 and rk2 == 3)
    to12 = lambda x: [x[0], x[1], x[2], x[3], x[7]] + list(x[8:]); jr = jac_rank(P, to12(v), I, p)
    hr, hh, names = high_system(P, to12(v)[:5], I, p); sol_, hrk, free = solve_field(I, hr, hh)
    print("   Delta*: Jacobian rank %d of 12 variables (10 = isolated and reduced up to scaling and the invisible direction); high system rank %d, free high directions %d" % (jr, hrk, len(free)))
    check("Delta* is an isolated reduced solution: Jacobian rank 10 = 12 - 1 (scaling) - 1 (invisible c6 direction)", jr == 10 and len(free) == 1)
    print("   Delta* (c1_0 = 1, invisible direction fixed by the solver):")
    for j in range(15): print("      %-5s %s      [Pi point over the same c1: %s]" % (S.names[j], rec.fmt(v[j]), rec.fmt(vP[j])))
    print("   tangency point: c1 = (1, %s, %s);  tau_Pi = %s, tau* = %s" % (rec.fmt(a_s), rec.fmt(b_s), rec.fmt(vP[3]), rec.fmt(new[0])))
    dv = [(v[j] - vP[j]) % p for j in range(15)]; print("   Delta* - (point of Pi over the same c1): non-zero coordinates %s" % [S.names[j] for j in range(15) if dv[j]])
    # a4, B4 at this c1
    print("   at this c1: a4 = %d, B4(tau) identically zero: %s (so c_4 = 0 for every tau and w)" % (mevalR(D['a4'], [a_s, b_s, 0], I), not muni(msubs(D['B4'], {0: a_s, 1: b_s}, p), 2, p)))
    json.dump(dict(delta=v, pi_point=vP, a=a_s, b=b_s, tau=new[0]), open(os.path.join(HERE, "ref9_newsol_%s.json" % a.prime), "w"))
    # ---------- rank of the 11 x 8 linear system along the curve of Pi: where does it drop below 7 ?
    xs = list(range(1, 40)); dets = [[], [], []]
    XY = [([[rng.randrange(p) for _ in range(11)] for _ in range(7)], [[rng.randrange(p) for _ in range(7)] for _ in range(8)]) for _ in range(3)]
    for s in xs:
        w_ = [(A[j] + s * B[j]) % p for j in range(15)]; pt = [w_[1], w_[2], w_[3]]
        M = [[ev3K(I, c, pt) for c in T[k][0]] for k in LINC]
        for q, (X, Yv) in enumerate(XY):
            XM = [[sum(X[i][r] * M[r][c] for r in range(11)) % p for c in range(8)] for i in range(7)]
            XMY = [[sum(XM[i][c] * Yv[c][j] for c in range(8)) % p for j in range(7)] for i in range(7)]; dets[q].append(det_modp(XMY, p))
    dp = [newton_interp(xs, d, p) for d in dets]; g = []
    for d in dp: g = ugcd(g, d, p) if g else umonic(d, p)
    print("   along Pi (point A + s B): three random 7 x 7 compressions of the 11 x 8 coefficient matrix: degrees %s in s (nodes %d); gcd degree %d" % ([len(d) - 1 for d in dp], len(xs), len(g) - 1))
    check("the linear system for (w, c4, c5, c6) has rank 7 (unique solution up to the invisible direction) at EVERY point of the chart part of Pi", len(g) == 1 and all(len(d) - 1 < len(xs) - 3 for d in dp))
    print("REF9 failures:", FAIL)


if __name__ == "__main__":
    main()
