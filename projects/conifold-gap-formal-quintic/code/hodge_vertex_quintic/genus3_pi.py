# -*- coding: utf-8 -*-
"""genus3_pi.py -- the plane Pi: a 2-dimensional LINEAR space of circulant Delta = N_1 - N_1^taut, with non-zero low-order part, on which every
N_1-condition of genus 2 and 3 vanishes identically (4 linear rows, 5 two-edge quadrics, 9 three-edge cubics).

Found by the line test of genus3_cubicc.py in the chart c_1(0) = 1: with w eliminated by c_4 = 0 (scaling by a_4), the 3 x 3 minors of the six
c4-rows (cubics 9, 10, 11, 12, 19 and q_3) are polynomials E_k(a, b; tau), c1 = (1, a, b); on a random line of the (a, b)-plane the
tau-resultants Res(E_1, E_k) have a common factor of degree 6: a factor of degree 4 supported on the conic a_4 = 0 (spurious) and ONE more
point (a double root), F_p-rational on every line, at which the whole system is solvable.  These points span a 2-dimensional linear space.
This script: collects such points on random lines, checks that they span a plane, normalises a basis
   A: (c1_0, c1_1) = (1, 0),   B: (c1_0, c1_1) = (0, 1)      (high part: the particular solution of the linear system, free directions zero),
checks that all quadrics and cubics vanish IDENTICALLY on x A + y B (5 points of P^1 for forms of degree <= 3), describes the line of c1
(a Fourier mode of c_1 vanishes), accounts for the degree-6 common factor, and stores Pi (keys 'Pi_A', 'Pi_B' of g3_cubic_<p>.json).
Run: py -I genus3_pi.py --prime 61|89
"""
import os, sys, json, argparse, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_solve
from fq3_steps import State
from fq3_quad import rank_fp, IDX_2EDGE
from fq3_cubic import IDX_3EDGE
from fq3_rec import Rec
import fq3_quadlib as L
from fq3_circsys import CircSys, ptrim, peval, pgcd, interp_simple, sylvester_res
from genus3_cubicb import full_point


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61")
    a = ap.parse_args(); S = CircSys(a.prime); p = S.p; rec = Rec(S.R); cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime)); names = S.names
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("-" * 110); log("genus3_pi: the plane Pi of circulant solutions, p = 2^%s-1" % a.prime); rng = random.Random(23)
    lam = [[rng.randrange(p) for _ in range(20)] for _ in range(4)]; triples = list(itertools.combinations(range(6), 3))
    def rows_at(c1, tau):
        l0 = S.low(c1, tau, 0); l1 = S.low(c1, tau, 1); dl = [(l1[j] - l0[j]) % p for j in range(8)]
        B4 = S.cub(4, l0 + [0] * 7); a4 = (S.cub(4, l1 + [0] * 7) - B4) % p; rows = []
        for g_ in (9, 10, 11, 12, 19):
            P = S.cub(g_, l0 + [0, 0]); Sg = (S.cub(g_, l1 + [0, 0]) - P) % p
            rows.append([(S.cub(g_, l0 + [1, 0]) - P) % p, (S.cub(g_, l0 + [0, 1]) - P) % p, (a4 * P - B4 * Sg) % p])
        v0 = l0 + [0] * 7; g3 = S.quad_grad(3, v0); q0 = S.quad(3, v0); q2 = S.quad(3, dl + [0] * 7); q1 = (S.quad(3, l1 + [0] * 7) - q0 - q2) % p
        rows.append([a4 * g3[8] % p, a4 * g3[9] % p, (a4 * a4 * q0 - a4 * B4 * q1 + B4 * B4 * q2) % p])
        return rows, a4, B4, l0, dl
    def epolys(c1, n=6):
        xt = list(range(1, n + 1)); vals = []
        for t in xt:
            rows = rows_at(c1, t)[0]; ms = [L.det_fp([rows[x] for x in T], p) for T in triples]
            vals.append([sum(lam[k][i] * ms[i] for i in range(20)) % p for k in range(4)])
        return [ptrim(interp_simple(xt, [x[k] for x in vals], p), p) for k in range(4)]
    def gcdl(fs):
        g = []
        for f in fs:
            if f: g = pgcd(g, f, p) if g else list(f)
        return g
    sols = []; accounted = True
    for line in range(4):
        a0, a1, b0, b1 = [rng.randrange(p) for _ in range(4)]; ns = 60; xs = list(range(1, ns + 1)); vals = [[], [], []]; a4s = []; lcs = [[], [], [], []]
        for s_ in xs:
            c1 = [1, (a0 + s_ * a1) % p, (b0 + s_ * b1) % p]; E = epolys(c1); a4s.append(rows_at(c1, 1)[1])
            for k in range(3): vals[k].append(sylvester_res(E[0], E[k + 1], 3, 3, p))
            for k in range(4): lcs[k].append(E[k][3] if len(E[k]) > 3 else 0)
        Us = [ptrim(interp_simple(xs, v, p), p) for v in vals]; g = gcdl(Us); a4l = ptrim(interp_simple(xs, a4s, p), p)
        lcg = gcdl([ptrim(interp_simple(xs, v, p), p) for v in lcs])            # common factor of the tau^3-coefficients (common root at tau = infinity)
        good = []
        for s_ in L.roots_fp(g, p, rng):
            c1 = [1, (a0 + s_ * a1) % p, (b0 + s_ * b1) % p]
            if rows_at(c1, 1)[1] == 0: continue
            for tau in L.roots_fp(gcdl(epolys(c1)), p, rng):
                rows, a4, B4, l0, dl = rows_at(c1, tau); w = (-B4) * pow(a4, p - 2, p) % p
                v, nfree = full_point(S, [(l0[j] + w * dl[j]) % p for j in range(8)])
                if v is not None and not any(S.residuals(v).values()): sols.append(v); good.append(s_)
        # the common factor: (s - s0) times a factor supported on the roots of a_4 restricted to the line
        rest = list(g)
        mult = 0
        for s_ in good:
            while len(rest) > 1 and not any(L.pmod(rest, [(-s_) % p, 1], p)): rest = ptrim(L.pdiv(rest, [(-s_) % p, 1], p), p); mult += 1
        def pmul_(x_, y_):
            out = [0] * (len(x_) + len(y_) - 1)
            for i, x in enumerate(x_):
                for j, y in enumerate(y_): out[i + j] = (out[i + j] + x * y) % p
            return out
        spur = [1]
        for _ in range(4): spur = pmul_(pmul_(spur, a4l), lcg)
        acc = len(pgcd(rest, spur, p)) == len(rest)
        log("   line %d: resultant degrees %s, common factor of degree %d: %d solution point(s) in F_p (multiplicity %d as a root); a_4|line has degree %d, the common factor of the tau^3-coefficients degree %d; the remaining factor of degree %d is supported on a_4 = 0: %s" % (line, [len(u) - 1 for u in Us], len(g) - 1, len(good), mult, len(a4l) - 1, len(lcg) - 1, len(rest) - 1, acc))
        accounted = accounted and acc and len(good) == 1
    check("on every random line of the chart c_1(0) = 1: exactly one solution point (F_p-rational, a double root of the common factor); the rest of the common factor is supported on the conic a_4 = 0", accounted)
    check("the solution points span a 2-dimensional linear space (rank of the 15-vectors, of their low parts, of their c1 parts)", rank_fp(sols, p) == 2 and rank_fp([v[:8] for v in sols], p) == 2 and rank_fp([v[:3] for v in sols], p) == 2,
          (rank_fp(sols, p), rank_fp([v[:8] for v in sols], p), rank_fp([v[:3] for v in sols], p)))
    def comb(target):
        co = fp_solve([[v[0] for v in sols[:2]], [v[1] for v in sols[:2]]], target, p)[0]
        return [sum(co[k] * sols[k][j] for k in range(2)) % p for j in range(15)]
    A = comb([1, 0]); B = comb([0, 1])
    # (the high parts of A and B are the combinations of the high parts of generic solution points, unique up to the invisible direction,
    #  normalised by c6_2 = 0; re-solving the fibre system at B, where the fibre is larger, would pick a different point of the fibre)
    ok = True
    for (x, y) in ((1, 0), (0, 1), (1, 1), (1, 2), (3, p - 1), (5, 7), (rng.randrange(p), rng.randrange(p))):
        v = [(x * A[j] + y * B[j]) % p for j in range(15)]; ok = ok and not any(S.residuals(v).values())
    check("ALL N_1-conditions vanish identically on Pi = span(A, B): 7 points of P^1 (forms of degree <= 3 in (x : y))", ok)
    log("   A: %s" % ["%s: %s" % (names[j], rec.fmt(A[j])) for j in range(15)])
    log("   B: %s" % ["%s: %s" % (names[j], rec.fmt(B[j])) for j in range(15)])
    z = S.R.zeta; inv2 = pow(2, p - 2, p)
    check("the c1-part of Pi is the line c1_0 + 2 zeta^3 c1_1 + 2 zeta c1_2 = 0 (the Fourier mode m = 3 of c_1 vanishes): A has c1 = (1, 0, -zeta^4/2), B has c1 = (0, 1, -zeta^2)",
          A[:3] == [1, 0, (-inv2 * pow(z, 4, p)) % p] and B[:3] == [0, 1, (-pow(z, 2, p)) % p])
    m = S.m1; check("Pi lies in m_1 = 0", all(sum(m[j] * v[j] for j in range(15)) % p == 0 for v in (A, B)))
    # is Pi + (pure-high kernel of the genus-3 rows) also made of solutions?  (it is not expected to be)
    hk = fp_solve([[S.J3[r][j] for j in range(8, 15)] for r in range(2)], [0, 0], p)[1]
    vv = [(A[j] + 3 * B[j]) % p for j in range(15)]; cnt = 0
    for h in hk:
        w_ = list(vv)
        for j in range(7): w_[8 + j] = (w_[8 + j] + h[j]) % p
        cnt += 0 if any(S.residuals(w_).values()) else 1
    log("   adding pure-high directions (kernel of the genus-3 rows, %d of them) to a point of Pi: %d keep all conditions (none is expected: the only free direction over a generic point of Pi is ONE combination of order-6 coordinates in the kernel of the genus-3 rows, not a basis vector)" % (len(hk), cnt))
    cs.put('Pi_A', A); cs.put('Pi_B', B); cs.put('pi_failed', FAILED)
    log("GENUS3_PI ALL CHECKS PASSED" if not FAILED else "GENUS3_PI FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
