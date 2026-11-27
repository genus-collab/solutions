# -*- coding: utf-8 -*-
"""genus3_cubicc.py -- circulant class at O(u^3): the strata around the isolated solution of genus3_cubicb.py (genus 3).

  NOTE (referee, g3_referee_u3.json): of the four solutions found at b = +-zeta^2 one (at b = -zeta^2) is the point of the plane Pi with
  c_1(0) = 0; the removed factors are examined at their F_p-roots only (they all split at both primes); the 12 non-rational roots of the
  genuine part were lifted and verified by the referee, not here.
  (1) chart c1 = (0, 1, b): the genuine part of U(b) = Res_tau(G, H): remove the roots of a_4(b), of the contents of G and H (b with G(b, .) = 0
      or H(b, .) = 0) and of the common leading coefficients; degree and number of F_p-roots; the special values of b examined directly.
  (2) the point c1 = (0, 0, 1).
  (3) the stratum c1 = 0: m_1 = 0 fixes the direction of c2; the (2,2,2) parts of the cubics; then the genus-2 rows and the quadrics on c3.
  (4) the chart c_1(0) = 1: at random c1 = (1, a, b) random combinations of the 3 x 3 minors of the six c4-rows (q_3, c_9, c_10, c_11, c_12,
      c_19), polynomials in tau after the c_4-scaling, have no common root: for GENERIC c1 there is no solution.  The solutions of this chart
      (the plane Pi) are found and tested in genus3_pi.py and genus3_cubicline.py.
Run: py -I genus3_cubicc.py --prime 61|89
"""
import os, sys, json, argparse, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_solve
from fq3_steps import State
from fq3_quad import rank_fp, IDX_2EDGE
from fq3_cubic import IDX_3EDGE
from fq3_rec import Rec
import fq3_quadlib as L
from fq3_circsys import CircSys, ptrim, peval, pgcd, interp_simple, sylvester_res
from genus3_cubicb import point_data, full_point


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61")
    a = ap.parse_args(); S = CircSys(a.prime); p = S.p; rec = Rec(S.R); cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("-" * 110); log("genus3_cubicc: circulant class, strata, p = 2^%s-1" % a.prime); rng = random.Random(23)
    def tau_poly(c1, n, which, raw=False):
        xs = list(range(1, n + 1)); ys = [point_data(S, c1, t)[which] for t in xs]; c = interp_simple(xs, ys, p)
        return c if raw else ptrim(c, p)
    def try_point(c1, tau, label):
        G1, H1, aux = point_data(S, c1, tau)
        if not any(aux['lowt']): log("      %s: scaled low part is zero (degenerate)" % label); return None
        v, nfree = full_point(S, aux['lowt'])
        if v is None: log("      %s: high system inconsistent" % label); return None
        bad = [k for k, x in S.residuals(v).items() if x]
        log("      %s: point reconstructed, non-zero residuals %s" % (label, bad))
        return v if not bad else None
    # ---------------- (1) genuine part of U on the chart c1 = (0, 1, b)
    nb = 70; xs = list(range(2, nb + 2)); Gc = [[], [], []]; Hc = [[], [], [], []]; a4v = []
    for b in xs:
        g = tau_poly([0, 1, b], 3, 0, raw=True); h = tau_poly([0, 1, b], 4, 1, raw=True)
        for k in range(3): Gc[k].append(g[k])
        for k in range(4): Hc[k].append(h[k])
        a4v.append(point_data(S, [0, 1, b], 1)[2]['a4'])
    Gp = [ptrim(interp_simple(xs, ys, p), p) for ys in Gc]; Hp = [ptrim(interp_simple(xs, ys, p), p) for ys in Hc]; a4p = ptrim(interp_simple(xs, a4v, p), p)
    check("b-degrees of the tau-coefficients of G, H and of a_4 resolved with %d nodes" % nb, all(len(f) < nb - 3 for f in Gp + Hp + [a4p]), ([len(f) - 1 for f in Gp], [len(f) - 1 for f in Hp], len(a4p) - 1))
    def gcdl(fs):
        g = []
        for f in fs:
            if f: g = pgcd(g, f, p) if g else list(f)
        return g
    contG = gcdl(Gp); contH = gcdl(Hp); lead = pgcd(Gp[2], Hp[3], p)
    U = cs.d['U_b1']; dU = [(k * U[k]) % p for k in range(1, len(U))]; sq = ptrim(L.pdiv(U, pgcd(U, dU, p), p), p)
    Ug = list(sq)
    for f in (a4p, contG, contH, lead):
        while True:
            g = pgcd(Ug, f, p)
            if len(g) <= 1: break
            Ug = ptrim(L.pdiv(Ug, g, p), p)
    rg = L.roots_fp(Ug, p, rng)
    log("   degrees: U %d, squarefree part %d; a_4 %d, content of G %d, content of H %d, gcd of the leading coefficients %d; GENUINE part %d with %d root(s) in F_p" % (len(U) - 1, len(sq) - 1, len(a4p) - 1, len(contG) - 1, len(contH) - 1, len(lead) - 1, len(Ug) - 1, len(rg)))
    cs.put('U_b1_genuine_degree', len(Ug) - 1); cs.put('U_b1_genuine_roots', rg)
    sol = cs.d.get('circ_b_solutions', [])
    check("the stored solution sits at a root of the genuine part", bool(sol) and sol[0][2] in rg)
    if len(Ug) - 1 > len(rg):
        # degrees of the irreducible factors of the non-rational part
        rest = list(Ug)
        for r in rg: rest = ptrim(L.pdiv(rest, [(-r) % p, 1], p), p)
        degs = []; xq = [0, 1]; d = 0
        while len(rest) > 1 and d < 40:
            d += 1; xq = L.ppowmod(xq, p, rest, p); g = pgcd(rest, ptrim([(xq[k] if k < len(xq) else 0) - (1 if k == 1 else 0) for k in range(max(len(xq), 2))], p), p)
            if len(g) > 1: degs.append((d, (len(g) - 1) // d)); rest = ptrim(L.pdiv(rest, g, p), p); xq = L.pmod(xq, rest, p) if len(rest) > 1 else xq
        log("   non-rational genuine roots: irreducible factors (degree, multiplicity count) %s" % degs); cs.put('U_b1_genuine_factors', degs)
    extra = []
    for nm, f, which, n in (("content of H (H(b, .) = 0)", contH, 0, 3), ("content of G (G(b, .) = 0)", contG, 1, 4)):
        for b in L.roots_fp(f, p, rng):
            pol = tau_poly([0, 1, b], n, which); rts = L.roots_fp(pol, p, rng) if len(pol) > 1 else []
            log("   special b = %s (%s): the other polynomial has degree %d, %d root(s) in F_p" % (rec.fmt(b), nm, len(pol) - 1, len(rts)))
            for t in rts:
                v = try_point([0, 1, b], t, "tau = %s" % rec.fmt(t))
                if v is not None: extra.append(v)
    for b in L.roots_fp(a4p, p, rng):
        # a_4 = 0: c_4 = B_4(tau) (no w); c_10 affine in w: scale by its w-coefficient; then H
        def lowt2(t):
            l0 = S.low([0, 1, b], t, 0); l1 = S.low([0, 1, b], t, 1); dl = [(l1[j] - l0[j]) % p for j in range(8)]
            B10 = S.cub(10, l0 + [0] * 7); a10 = (S.cub(10, l1 + [0] * 7) - B10) % p
            return [(a10 * l0[j] - B10 * dl[j]) % p for j in range(8)], S.cub(4, l0 + [0] * 7)
        def h2(t):
            lt, B4 = lowt2(t); v = lt + [0] * 7; b9 = S.cub(9, v); g3 = S.quad_grad(3, v); g8 = S.quad_grad(8, v)
            return L.det_fp([[(S.cub(9, lt + [1, 0]) - b9) % p, (S.cub(9, lt + [0, 1]) - b9) % p, b9], [g3[8], g3[9], S.quad(3, v)], [g8[8], g8[9], S.quad(8, v)]], p), B4
        xt = list(range(1, 25)); vals = [h2(t) for t in xt]
        Hh = ptrim(interp_simple(xt, [x[0] for x in vals], p), p); B4p = ptrim(interp_simple(xt, [x[1] for x in vals], p), p); g = pgcd(Hh, B4p, p) if (Hh and B4p) else (B4p or Hh)
        log("   a_4(b) = 0 at b = %s: B_4(tau) has degree %d, the residual condition degree %d, gcd degree %d" % (rec.fmt(b), len(B4p) - 1, len(Hh) - 1, len(g) - 1))
        check("no solution on the fibre a_4 = 0, b = %s (B_4 and the residual condition have no common root)" % rec.fmt(b), len(g) <= 1)
    cs.put('circ_b_extra', extra)
    log("   additional solutions at special b over F_p: %d" % len(extra))
    # ---------------- (2) c1 = (0, 0, 1)
    g = tau_poly([0, 0, 1], 6, 0); h = tau_poly([0, 0, 1], 8, 1); d = pgcd(g, h, p) if (g and h) else []
    log("   c1 = (0, 0, 1): a_4 %s; deg G %d, deg H %d, gcd degree %s" % ("= 0" if point_data(S, [0, 0, 1], 1)[2]['a4'] == 0 else "!= 0", len(g) - 1, len(h) - 1, (len(d) - 1) if d else None))
    check("c1 = (0, 0, 1): no solution (G and H non-zero with trivial gcd, a_4 != 0)", bool(g) and bool(h) and len(d) == 1 and point_data(S, [0, 0, 1], 1)[2]['a4'] != 0)
    # ---------------- (3) c1 = 0
    m = S.m1; v2 = [m[4] % p, (-m[3]) % p]                                   # c2 direction with m_1 = 0
    def lowc(s, w): return S.low([0, 0, 0], s * v2[0] % p, w)
    l = lowc(1, 0); check("c1 = 0: m_1 = 0 fixes c2 up to scale (c2_2 from m_1 agrees)", l[4] == v2[1] % p)
    c222 = {g_: S.cub(g_, lowc(1, 0) + [0] * 7) for g_ in IDX_3EDGE}
    log("   c1 = 0, c2 = s v2: cubics at s = 1, w = 0: non-zero for graphs %s" % [g_ for g_ in IDX_3EDGE if c222[g_]])
    indep_w = all(S.cub(g_, lowc(1, 5) + [0] * 7) == c222[g_] for g_ in IDX_3EDGE)
    check("c1 = 0: the cubics do not depend on c3 and some cubic is non-zero on the c2-line, so c2 = 0", indep_w and any(c222.values()))
    k3 = lowc(0, 1) + [0] * 7; qv = {i: S.quad(i, k3) for i in IDX_2EDGE}
    cross = all(S.quad(i, lowc(0, 1) + [rng.randrange(p) for _ in range(7)]) == qv[i] for i in IDX_2EDGE)
    log("   c1 = c2 = 0: the genus-2 rows leave the c3-line; quadrics on it: %s (independent of the orders >= 4: %s)" % ({i: rec.fmt(qv[i]) for i in IDX_2EDGE}, cross))
    check("c1 = c2 = 0: some two-edge quadric is non-zero on the c3-line and independent of the high part, so c3 = 0: the stratum c1 = 0 contains only the pure-high solutions", cross and any(qv.values()))
    # ---------------- (4) chart c1_0 = 1
    import itertools
    lam = [[rng.randrange(p) for _ in range(20)] for _ in range(4)]; triples = list(itertools.combinations(range(6), 3))
    def minors(c1, tau):
        """four fixed random combinations of the 20 maximal minors of the six c4-rows (cubics 9, 10, 11, 12, 19 and q_3), without denominators:
        unknown a_4 c4; cubic rows [A(c1) | a_4 P - B_4 S], q_3 row [a_4 grad | a_4^2 q0 - a_4 B_4 q1 + B_4^2 q2]"""
        l0 = S.low(c1, tau, 0); l1 = S.low(c1, tau, 1); dl = [(l1[j] - l0[j]) % p for j in range(8)]
        B4 = S.cub(4, l0 + [0] * 7); a4 = (S.cub(4, l1 + [0] * 7) - B4) % p; rows = []
        for g_ in (9, 10, 11, 12, 19):
            P = S.cub(g_, l0 + [0, 0]); Sg = (S.cub(g_, l1 + [0, 0]) - P) % p
            rows.append([(S.cub(g_, l0 + [1, 0]) - P) % p, (S.cub(g_, l0 + [0, 1]) - P) % p, (a4 * P - B4 * Sg) % p])
        v0 = l0 + [0] * 7; g3 = S.quad_grad(3, v0); q0 = S.quad(3, v0); q2 = S.quad(3, dl + [0] * 7); q1 = (S.quad(3, l1 + [0] * 7) - q0 - q2) % p
        rows.append([a4 * g3[8] % p, a4 * g3[9] % p, (a4 * a4 * q0 - a4 * B4 * q1 + B4 * B4 * q2) % p])
        ms = [L.det_fp([rows[x] for x in T], p) for T in triples]
        return [sum(lam[k][t] * ms[t] for t in range(20)) % p for k in range(4)]
    def epolys(c1, n=12):
        xt = list(range(1, n + 1)); vals = [minors(c1, t) for t in xt]
        return [ptrim(interp_simple(xt, [x[k] for x in vals], p), p) for k in range(4)]
    degs = []; gds = []
    for _ in range(6):
        E = epolys([1, rng.randrange(p), rng.randrange(p)]); degs.append([len(e) - 1 for e in E]); g = gcdl(E); gds.append(len(g) - 1)
    log("   chart c_1(0) = 1: tau-degrees of the four combinations of minors at random c1: %s; gcd degrees %s" % (degs[0], gds))
    check("chart c_1(0) = 1: at 6 random c1 = (1, a, b) the combinations of minors have no common root in tau (the cubic rows and q_3 are inconsistent)", all(d < 10 for d in degs[0]) and all(g == 0 for g in gds))
    log("   (generic c1 only: the positive-dimensional part of the solution set in this chart is the plane Pi, see genus3_pi.py and genus3_cubicline.py)")
    cs.put('cubicc_failed', FAILED)
    log("GENUS3_CUBICC ALL CHECKS PASSED" if not FAILED else "GENUS3_CUBICC FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
