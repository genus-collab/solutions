# -*- coding: utf-8 -*-
"""ref6_plane.py -- task (3): a TWO-DIMENSIONAL elimination in the chart c1_0 = 1 (isolated solutions).
Unscaled seven-row system (c_4, c_9, c_10, c_11, c_12, c_19, q_3 affine in (c4_1, c4_2, w)): D_1..D_4 = random combinations of the 35 maximal
minors, polynomials in (a, b, tau).  Shear a = alpha + mu b.  R_k(alpha, b) = Res_tau(D_1, D_k), k = 2, 3, 4 (plane curves of degree 42, all
divisible by L^2, L = the line of Pi);  R'_k = R_k / L^2;   U_23(alpha) = Res_b(R'_2, R'_3),  U_24(alpha) = Res_b(R'_2, R'_4);
g(alpha) = gcd(U_23, U_24): its roots are the alpha-coordinates of ALL common zeros of the seven-row system in the chart (off the generic
point of Pi), over the algebraic closure, plus possibly spurious ones (examined by ref7_cands.py).
Run: py -I ref6_plane.py --prime 61|89 [--seed 1]"""
import sys, os, random, json, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_chart import *


def shear_a(f, mu, p):
    """a = alpha + mu b in a 3-variable polynomial (a, b, t)"""
    out = {}; lin = {(1, 0, 0): 1, (0, 1, 0): mu % p}; pw = [{(0, 0, 0): 1}]
    for _ in range(mdeg(f, 0)): pw.append(mmul(pw[-1], lin, p))
    for (i, j, k), c in f.items(): out = madd(out, mmul(pw[i], {(0, j, k): c}, p), p)
    return out


def newton_consecutive(x0, ys, p):
    """Newton interpolation on the consecutive nodes x0, x0 + 1, ...; returns coefficients low -> high"""
    n = len(ys); c = [y % p for y in ys]; inv = [0] + [pow(j, p - 2, p) for j in range(1, n)]
    for j in range(1, n):
        ij = inv[j]
        for i in range(n - 1, j - 1, -1): c[i] = (c[i] - c[i - 1]) * ij % p
    out = [c[n - 1]]
    for i in range(n - 2, -1, -1):
        xi = x0 + i; new = [0] * (len(out) + 1)
        for k, v in enumerate(out): new[k + 1] = v
        for k, v in enumerate(out): new[k] = (new[k] - xi * v) % p
        new[0] = (new[0] + c[i]) % p; out = new
    return utrim(out, p)


def setup(tag, seed):
    S = load(tag); p = S.p; P = load_polys(tag); rng = random.Random(1000 + seed); z = S.R.zeta
    D = chart_data(P, p); M7 = minors(D['D7'], 4, p); lam = [[rng.randrange(p) for _ in range(35)] for _ in range(4)]; mu = rng.randrange(p)
    DD = [shear_a(f, mu, p) for f in combos(M7, lam, p)]; dt = max(mdeg(f, 2) for f in DD)
    cf = []
    for f in DD:
        c = mcoeffs(f, 2); cf += c + [{}] * (dt + 1 - len(c))
    L = [(1, 2 * pow(z, 3, p) % p), ((2 * pow(z, 3, p) * mu + 2 * z) % p,)]          # L = (1 + 2 z^3 alpha) + (2 z^3 mu + 2 z) b
    return dict(S=S, p=p, P=P, D=D, M7=M7, lam=lam, mu=mu, DD=DD, dt=dt, fast=Fast2(cf, p), Lc=(1, 2 * pow(z, 3, p) % p, (2 * pow(z, 3, p) * mu + 2 * z) % p), rng=rng)


def Rprime_at(st, alpha, nb, dR, mL=2):
    """the three polynomials R'_k(alpha, b) (k = 2, 3, 4), degree <= dR - mL, by interpolation on b = 1..nb"""
    p = st['p']; dt = st['dt']; us = st['fast'].at_a(alpha); vals = [[], [], []]; l0 = (st['Lc'][0] + st['Lc'][1] * alpha) % p; l1 = st['Lc'][2]
    for b in range(1, nb + 1):
        co = [ueval(u, b, p) for u in us]; f = [utrim(co[k * (dt + 1):(k + 1) * (dt + 1)], p) for k in range(4)]
        Lv = (l0 + l1 * b) % p; iL = pow(pow(Lv, mL, p), p - 2, p)
        for k in range(3): vals[k].append(ures(f[0], f[k + 1], dt, dt, p) * iL % p)
    return [newton_consecutive(1, v, p) for v in vals]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--nodes", type=int, default=0)
    a = ap.parse_args(); st = setup(a.prime, a.seed); p = st['p']; t0 = time.time()
    print("=" * 100); print("ref6_plane, p = 2^%s-1, seed %d" % (a.prime, a.seed))
    dt = st['dt']; dtot = max(mdeg(f) for f in st['DD']); print("   D_k: total degree %d, tau-degree %d" % (dtot, dt))
    # degree of R_k in b at generic alpha, multiplicity of L, common factor
    al = 12345; nb0 = 2 * dtot * dt + 6; R0 = Rprime_at(st, al, nb0, None, mL=0); dR = max(len(r) - 1 for r in R0)
    g0 = []
    for r in R0: g0 = ugcd(g0, r, p) if g0 else umonic(r, p)
    bL = (-(st['Lc'][0] + st['Lc'][1] * al)) * pow(st['Lc'][2], p - 2, p) % p; mL = umult(g0, bL, p)
    print("   at a generic alpha: b-degrees of R_2, R_3, R_4: %s (node count %d); common factor of degree %d = (b - b_L)^%d" % ([len(r) - 1 for r in R0], nb0, len(g0) - 1, mL))
    assert len(g0) - 1 == mL == 2 and dR < nb0 - 3, "unexpected common factor"
    dRp = dR - mL; nb = dRp + 1 + 2
    chk = Rprime_at(st, al, nb, dR); lin = [(-bL) % p, 1]
    assert all(len(r) - 1 <= dRp for r in chk) and all(uscale(umul(umul(chk[k], lin, p), lin, p), pow(st['Lc'][2], 2, p), p) == R0[k] for k in range(3)), "division by L^2 failed"
    print("   R'_k = R_k / L^2: b-degree %d (verified with 2 spare nodes; R'_k L^2 = R_k at the test alpha)" % dRp)
    N = a.nodes or (dRp * dRp + 1 + 8); U23 = []; U24 = []
    for i in range(N):
        Rp = Rprime_at(st, 1 + i, nb, dR)
        U23.append(ures(Rp[0], Rp[1], dRp, dRp, p)); U24.append(ures(Rp[0], Rp[2], dRp, dRp, p))
        if i % 200 == 0: print("      node %d / %d  (%.0f s)" % (i, N, time.time() - t0), flush=True)
    P23 = newton_consecutive(1, U23, p); P24 = newton_consecutive(1, U24, p)
    print("   U_23, U_24: degrees %d, %d from %d nodes (bound %d): resolved: %s   (%.0f s)" % (len(P23) - 1, len(P24) - 1, N, dRp * dRp, len(P23) - 1 <= N - 5 and len(P24) - 1 <= N - 5, time.time() - t0))
    g = ugcd(P23, P24, p); rngf = random.Random(5)
    print("   g = gcd(U_23, U_24): degree %d  (%.0f s)" % (len(g) - 1, time.time() - t0))
    if len(g) > 1:
        sq = usqfree(g, p); fs = factor(sq, p, rngf)
        print("   squarefree part of g: degree %d, irreducible factors of degrees %s; multiplicities of the F_p-roots: %s" % (len(sq) - 1, fmt_fac(fs), [umult(g, (-f[0]) % p, p) for f in fs if len(f) == 2]))
    else: fs = []
    json.dump(dict(seed=a.seed, mu=st['mu'], lam=st['lam'], g=g, factors=fs, dR=dR, dRp=dRp, degU=[len(P23) - 1, len(P24) - 1], N=N), open(os.path.join(HERE, "ref6_plane_%s_seed%d.json" % (a.prime, a.seed)), "w"))
    print("REF6 done (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
