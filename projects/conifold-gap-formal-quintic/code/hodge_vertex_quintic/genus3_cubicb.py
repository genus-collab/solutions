# -*- coding: utf-8 -*-
"""genus3_cubicb.py -- circulant class at O(u^3): the locus c_1(0) = 0, m_1 = 0, c_1 != 0 (genus 3).

c_32 = kappa m_1^3 forces m_1 = 0.  On {c_1(0) = 0, m_1 = 0, genus-2 rows} (coordinates (c1_1 : c1_2 : tau = c2_1 : w = c3_2), a P^3):
  * the cubics 11, 12, 18, 19, 20, 32 and the quadrics 6, 17 vanish identically (in all 15 coordinates);
  * c_4 and c_10 do not involve the orders >= 4: two conditions on the low part;
  * c4 (2 unknowns) enters c_9, q_3, q_8: one consistency condition  det [A | value]_{c_9, q_3, q_8} = 0;
  * (c5, c6) enter q_7 and the two genus-3 rows only: always solvable (2-dimensional fibre, one direction invisible).
So the solutions with c_1(0) = 0, c_1 != 0 are the common zeros of THREE forms on P^3 (degrees 3, 3, 5): finitely many points expected.
Elimination without denominators (homogeneity): c_4 is affine in w with coefficient a_4(c1); scale the low part by a_4 so that c_4 = 0:
  G(b, tau) = c_10,   H(b, tau) = det [A | value]   at c1 = (0, 1, b);   U(b) = Res_tau(G, H).
F_p-roots of U -> candidate points -> the full Delta is reconstructed and EVERY condition (4 linear rows, 5 quadrics, 9 cubics) is verified.
Run: py -I genus3_cubicb.py --prime 61|89
"""
import os, sys, json, argparse, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_solve
from fq3_steps import State
from fq3_quad import rank_fp
from fq3_cubic import IDX_3EDGE
from fq3_quad import IDX_2EDGE
from fq3_rec import Rec
import fq3_quadlib as L
from fq3_circsys import CircSys, ptrim, peval, pgcd, interp_simple, sylvester_res

ROWS_C = (9, 10, 12)


def point_data(S, c1, tau):
    """returns (G, H, aux) at c1 (3 coordinates with c1[0] = 0) and tau"""
    p = S.p
    l0 = S.low(c1, tau, 0); l1 = S.low(c1, tau, 1); dl = [(l1[j] - l0[j]) % p for j in range(8)]
    B4 = S.cub(4, l0 + [0] * 7); a4 = (S.cub(4, l1 + [0] * 7) - B4) % p
    lowt = [(a4 * l0[j] - B4 * dl[j]) % p for j in range(8)]; v = lowt + [0] * 7
    G = S.cub(10, v)
    b9 = S.cub(9, v); A9 = [(S.cub(9, lowt + [1, 0]) - b9) % p, (S.cub(9, lowt + [0, 1]) - b9) % p]
    g3 = S.quad_grad(3, v); g8 = S.quad_grad(8, v)
    M = [A9 + [b9], [g3[8], g3[9], S.quad(3, v)], [g8[8], g8[9], S.quad(8, v)]]
    H = L.det_fp(M, p)
    return G, H, dict(a4=a4, B4=B4, lowt=lowt, A9=A9, b9=b9)


def full_point(S, lowt):
    """solve the linear system for y = (c4, c5, c6) (2 genus-3 rows, 5 quadrics, the c4-rows of the cubics 9, 10, 11, 12, 19); returns the
    15-vector or None"""
    p = S.p; v = list(lowt) + [0] * 7; rows = []; rhs = []
    for r in range(2): rows.append([S.J3[r][j] for j in range(8, 15)]); rhs.append((-sum(S.J3[r][j] * v[j] for j in range(8))) % p)
    for i in IDX_2EDGE:
        g = S.quad_grad(i, v); rows.append([g[j] for j in range(8, 15)]); rhs.append((-S.quad(i, v)) % p)
    for g in (9, 10, 11, 12, 19):
        b0 = S.cub(g, v); rows.append([(S.cub(g, lowt + [1, 0]) - b0) % p, (S.cub(g, lowt + [0, 1]) - b0) % p] + [0] * 5); rhs.append((-b0) % p)
    sol = fp_solve(rows, rhs, p)
    if sol[0] is None: return None, None
    return list(lowt) + sol[0], len(sol[1])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--nb", type=int, default=0)
    a = ap.parse_args(); S = CircSys(a.prime); p = S.p; rec = Rec(S.R); cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("-" * 110); log("genus3_cubicb: circulant class, the locus c_1(0) = 0, m_1 = 0, p = 2^%s-1" % a.prime)
    rng = random.Random(11)
    # structure checks at random points of the locus (all 15 coordinates, with the linear rows imposed on the low part only)
    v = S.low([0, rng.randrange(p), rng.randrange(p)], rng.randrange(p), rng.randrange(p)) + [rng.randrange(p) for _ in range(7)]
    res = S.residuals(v)
    check("on {c_1(0) = 0, m_1 = 0, L2}: the cubics 11, 12, 18, 19, 20, 32 and the quadrics 6, 17 vanish identically (all 15 coordinates random)",
          all(res['c_%d' % i] == 0 for i in (11, 12, 18, 19, 20, 32)) and res['q_6'] == 0 and res['q_17'] == 0)
    v2 = v[:8] + [rng.randrange(p) for _ in range(7)]; res2 = S.residuals(v2)
    check("on the locus c_4 and c_10 do not depend on the orders >= 4", res2['c_4'] == res['c_4'] and res2['c_10'] == res['c_10'])
    v3 = v[:10] + [rng.randrange(p) for _ in range(5)]; res3 = S.residuals(v3)
    check("on the locus c_9, q_3, q_8 depend on the orders <= 4 only", all(res3[k] == res[k] for k in ('c_9', 'q_3', 'q_8')))
    # degrees in tau
    def tau_poly(b, n, which):
        xs = list(range(1, n + 1)); ys = [point_data(S, [0, 1, b], t)[which] for t in xs]
        return ptrim(interp_simple(xs, ys, p), p)
    dG = dH = 0
    for _ in range(3):
        b = rng.randrange(p); dG = max(dG, len(tau_poly(b, 14, 0)) - 1); dH = max(dH, len(tau_poly(b, 30, 1)) - 1)
    log("   tau-degrees: deg G = %d, deg H = %d" % (dG, dH))
    check("tau-degrees resolved by the node counts (14, 30)", dG < 12 and dH < 28)
    def U_at(b):
        g = tau_poly(b, dG + 1, 0); h = tau_poly(b, dH + 1, 1)
        return sylvester_res(g, h, dG, dH, p), g, h
    nb = a.nb or 900
    t = time.time(); xs = list(range(2, nb + 2)); ys = [U_at(b)[0] for b in xs]; U = ptrim(interp_simple(xs, ys, p), p)
    log("   U(b) = Res_tau(G, H): interpolated on %d nodes in %.0f s: degree %d" % (nb, time.time() - t, len(U) - 1))
    check("b-degree of U resolved (degree < nodes - 5)", len(U) - 1 < nb - 5)
    cs.put('U_b1', U)
    roots = L.roots_fp(U, p, rng); log("   U has %d distinct roots in F_p" % len(roots))
    # distinct-degree information: gcd(U, x^p - x) already used; squarefree part and its degree
    dU = [(k * U[k]) % p for k in range(1, len(U))]; sq = L.pdiv(U, pgcd(U, dU, p), p); log("   squarefree part of U has degree %d" % (len(ptrim(sq, p)) - 1))
    sols = []
    for b in roots:
        Ub, g, h = U_at(b); d = pgcd(g, h, p) if (g and h) else []
        G0, H0, aux = point_data(S, [0, 1, b], 1)
        info = "a4 %s, G(b,.) %s, H(b,.) %s, gcd degree %s" % ("=0" if aux['a4'] == 0 else "!=0", "zero" if not g else "deg %d" % (len(g) - 1), "zero" if not h else "deg %d" % (len(h) - 1), (len(d) - 1) if d else None)
        log("   root b = %s: %s" % (rec.fmt(b), info))
        if not d or len(d) < 2: continue
        for tau in L.roots_fp(d, p, rng):
            G1, H1, aux = point_data(S, [0, 1, b], tau)
            if aux['a4'] == 0: log("      tau = %s: degenerate (a4 = 0)" % rec.fmt(tau)); continue
            v, nfree = full_point(S, aux['lowt'])
            if v is None: log("      tau = %s: high system inconsistent" % rec.fmt(tau)); continue
            if v[1] == 0: log("      tau = %s: c1_1 = 0 after scaling (degenerate)" % rec.fmt(tau)); continue
            inv = pow(v[1], p - 2, p); v = [x * inv % p for x in v]          # normalise c1_1 = 1
            res = S.residuals(v); bad = [k for k, x in res.items() if x]
            log("      tau = %s: full point reconstructed (%d free high directions); non-zero residuals: %s" % (rec.fmt(tau), nfree, bad))
            if not bad: sols.append(v)
    check("genuine solutions with c_1(0) = 0, c_1 != 0 found over F_p (all 4 linear rows, 5 quadrics, 9 cubics vanish)", len(sols) > 0, len(sols))
    for k, v in enumerate(sols):
        log("   SOLUTION %d (normalised c1_1 = 1): %s" % (k, ["%s: %s" % (S.names[j], rec.fmt(v[j])) for j in range(15)]))
    cs.put('circ_b_solutions', sols); cs.put('circ_b_roots', roots)
    log("GENUS3_CUBICB done; failures: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
