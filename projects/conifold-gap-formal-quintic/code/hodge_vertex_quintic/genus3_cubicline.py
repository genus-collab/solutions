# -*- coding: utf-8 -*-
"""genus3_cubicline.py -- chart c_1(0) = 1 of the circulant class at O(u^3): what else lies on a random line of the (a, b)-plane?  (genus 3)

On a random line c1 = (1, a0 + s a1, b0 + s b1) the tau-resultants of the minors E_k (consistency of the six c4-rows: cubics 9, 10, 11, 12,
19 and q_3, after eliminating w by c_4 = 0) have a common factor g(s) of degree 6:
    g = (s - s_Pi)^2 . (factor of degree 4 supported on the conic a_4 = 0),
i.e. off a_4 = 0 the only common zero is the point of the plane Pi (genus3_pi.py), a double root.  The remaining condition is tested through
F(tau), the determinant expressing the consistency of the (c5, c6) system (2 genus-3 rows, q_6, q_7, q_8) at the Cramer solution c4 of the
rows 9, 10 (everything scaled, no denominators): Res_tau(E_k, F) vanishes at s_Pi for all k.  Any further factor c(s) of g would be tested by
divisibility c | Res_tau(E_k, F).  A line test sees positive-dimensional families only.  The isolated solutions of this chart were determined by the referee's two-dimensional
elimination (referee_u3/ref6_plane.py ..): exactly one, Delta*, over the point where the c1-line of Pi touches the conic a_4 = 0.
Run: py -I genus3_cubicline.py --prime 61|89
"""
import os, sys, json, argparse, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq3_steps import State
import fq3_quadlib as L
from fq3_circsys import CircSys, ptrim, pgcd, interp_simple, sylvester_res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--lines", type=int, default=2)
    a = ap.parse_args(); S = CircSys(a.prime); p = S.p; cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("-" * 110); log("genus3_cubicline: chart c_1(0) = 1, random lines, p = 2^%s-1" % a.prime); rng = random.Random(23)
    lam = [[rng.randrange(p) for _ in range(20)] for _ in range(4)]; triples = list(itertools.combinations(range(6), 3))
    def data(c1, tau):
        l0 = S.low(c1, tau, 0); l1 = S.low(c1, tau, 1); dl = [(l1[j] - l0[j]) % p for j in range(8)]
        B4 = S.cub(4, l0 + [0] * 7); a4 = (S.cub(4, l1 + [0] * 7) - B4) % p; rows = []
        for g_ in (9, 10, 11, 12, 19):
            P = S.cub(g_, l0 + [0, 0]); Sg = (S.cub(g_, l1 + [0, 0]) - P) % p
            rows.append([(S.cub(g_, l0 + [1, 0]) - P) % p, (S.cub(g_, l0 + [0, 1]) - P) % p, (a4 * P - B4 * Sg) % p])
        v0 = l0 + [0] * 7; g3 = S.quad_grad(3, v0); q0 = S.quad(3, v0); q2 = S.quad(3, dl + [0] * 7); q1 = (S.quad(3, l1 + [0] * 7) - q0 - q2) % p
        rows.append([a4 * g3[8] % p, a4 * g3[9] % p, (a4 * a4 * q0 - a4 * B4 * q1 + B4 * B4 * q2) % p])
        ms = [L.det_fp([rows[x] for x in T], p) for T in triples]; E = [sum(lam[k][i] * ms[i] for i in range(20)) % p for k in range(4)]
        # Cramer solution of the rows 9, 10 and the (c5, c6) determinant, everything scaled by delta * a_4
        r9, r10 = rows[0], rows[1]; delta = (r9[0] * r10[1] - r9[1] * r10[0]) % p
        c4 = [(-(r9[2] * r10[1] - r10[2] * r9[1])) % p, (-(r9[0] * r10[2] - r10[0] * r9[2])) % p]
        v10 = [delta * ((a4 * l0[j] - B4 * dl[j]) % p) % p for j in range(8)] + c4
        A5, rhs = S.high_system(v10, rows_q=(6, 7, 8))
        F = L.det_fp([A5[r][:4] + [rhs[r]] for r in range(5)], p)
        return E, F, a4
    def polys(c1, nE=6, nF=0):
        n = max(nE, nF); xt = list(range(1, n + 1)); vals = [data(c1, t) for t in xt]
        E = [ptrim(interp_simple(xt[:nE], [x[0][k] for x in vals[:nE]], p), p) for k in range(4)]
        F = ptrim(interp_simple(xt, [x[1] for x in vals], p), p) if nF else None
        return E, F
    dF = len(polys([1, rng.randrange(p), rng.randrange(p)], 6, 40)[1]) - 1
    log("   tau-degree of the (c5, c6) determinant F: %d" % dF); check("tau-degree of F resolved", 0 < dF < 36)
    def gcdl(fs):
        g = []
        for f in fs:
            if f: g = pgcd(g, f, p) if g else list(f)
        return g
    verdicts = []
    for line in range(a.lines):
        a0, a1, b0, b1 = [rng.randrange(p) for _ in range(4)]; ns = 6 * dF + 3 * 40 + 20; xs = list(range(1, ns + 1)); UE = [[], [], []]; UF = [[], [], [], []]; a4s = []
        for s_ in xs:
            c1 = [1, (a0 + s_ * a1) % p, (b0 + s_ * b1) % p]; E, F = polys(c1, 6, dF + 1); a4s.append(data(c1, 1)[2])
            for k in range(3): UE[k].append(sylvester_res(E[0], E[k + 1], 3, 3, p))
            for k in range(4): UF[k].append(sylvester_res(E[k], F, 3, dF, p))
        UEp = [ptrim(interp_simple(xs, v, p), p) for v in UE]; UFp = [ptrim(interp_simple(xs, v, p), p) for v in UF]; a4l = ptrim(interp_simple(xs, a4s, p), p)
        ok_deg = all(len(u) - 1 < ns - 3 for u in UEp + UFp)
        g = gcdl(UEp); a4pow = [1]
        for _ in range(5):
            new = [0] * (len(a4pow) + len(a4l) - 1)
            for i, x in enumerate(a4pow):
                for j, y in enumerate(a4l): new[i + j] = (new[i + j] + x * y) % p
            a4pow = new
        rest = ptrim(L.pdiv(g, pgcd(g, a4pow, p), p), p)                         # roots of g off the conic a_4 = 0
        rts = L.roots_fp(rest, p, rng); c = list(rest)
        mult = []
        for r in rts:
            m_ = 0
            while len(c) > 1 and not any(L.pmod(c, [(-r) % p, 1], p)): c = ptrim(L.pdiv(c, [(-r) % p, 1], p), p); m_ += 1
            mult.append(m_)
        gF = [len(pgcd(c, u, p)) - 1 if len(c) > 1 else 0 for u in UFp]; gPi = [all(L.pmod(u, [(-r) % p, 1], p) in ([], [0]) or not any(L.pmod(u, [(-r) % p, 1], p)) for u in UFp) for r in rts]
        log("   line %d: nodes %d, degrees of Res(E_1, E_k) %s and of Res(E_k, F) %s (resolved: %s); common factor g of degree %d; off a_4 = 0: degree %d with %d F_p-root(s) of multiplicities %s; remaining (non-rational) factor c of degree %d" % (line, ns, [len(u) - 1 for u in UEp], [len(u) - 1 for u in UFp], ok_deg, len(g) - 1, len(rest) - 1, len(rts), mult, len(c) - 1))
        log("      at the F_p-root(s) (the point of Pi): every Res(E_k, F) vanishes: %s;  gcd(c, Res(E_k, F)) degrees: %s" % (gPi, gF))
        verdicts.append(dict(ok_deg=ok_deg, deg_c=len(c) - 1, gcd_c_F=gF, pi_ok=gPi, nroots=len(rts)))
    cs.put('chart_a_lines', verdicts)
    check("degrees resolved; on each line exactly one F_p-root off a_4 = 0, at which the (c5, c6) condition holds (the point of Pi)", all(v['ok_deg'] and v['nroots'] == 1 and all(v['pi_ok']) for v in verdicts))
    nomore = all(v['deg_c'] == 0 for v in verdicts)
    check("on each random line the common zeros of the minors off the conic a_4 = 0 are exactly the point of Pi (a double root): no other positive-dimensional family of solutions in the chart c_1(0) = 1", nomore)
    full = nomore; none = False
    cs.put('chart_a_only_Pi_on_lines', nomore); cs.put('cubicline_failed', FAILED)
    log("GENUS3_CUBICLINE done; failures: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
