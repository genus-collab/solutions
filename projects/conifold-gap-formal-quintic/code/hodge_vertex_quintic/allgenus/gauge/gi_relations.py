# -*- coding: utf-8 -*-
"""gi_relations.py -- the u-expansion of the gauge identity (G), its symplectic dual (G*), and the symplecticity of the
canonical gauge (GAUGE_IDENTITY.md s.3, s.4).  Exact over Q, stdlib only.  (Second version, 2026-11-13: see DIAGNOSIS below.)

Write P = Sum_n u^n P^(n)(z) (canonical gauge), Y^C(v(u), z) = Sum_n u^n eta_n(z), Y^C_Q = Sum_n u^n theta_n(z), eta_0 = theta_0 = -e_1.
  (G)_n :   P^(n) e_1 = Sum_{j<n} P^(j) eta_{n-j} - theta_n            (column 1 = the H column, from the other columns of lower order)
  (G*)  :   Y^C_Q(-z)^T eta P(z) = Y^C(-z)^T eta,  eta_ab = 5 delta_{a+b=3 mod 5}
            [(G*) = (G) + symplecticity: Y^C_Q(-z) = P(-z) Y^C(-z) and P(-z)^T eta P(z) = eta]
Checks: (S)  the canonical gauge is symplectic, P(z) eta^-1 P(-z)^T = eta^-1;
        (S') the stored gauge P_s = P_can c(z), c(z) = P^s_11(0, z), has the EXACT defect P_s(z) eta^-1 P_s(-z)^T - eta^-1 =
             (c(z) c(-z) - 1) eta^-1, which is even in z, hence starts at z^10, not at z^5 (the first version expected z^5: wrong);
        (G)_n, (G*)_n for n = 1..NMAX on every complete z-power;  (R1)-(R3);  the explicit row-2 relations at u^3, u^4.
DIAGNOSIS of the (G*) failure of the first version: its covector map used (-1) ** m with m < 0 (powers of 1/z), which is a Python
FLOAT (-1.0); every product with it became a float and a few slots failed by rounding (exact comparison float vs Fraction).  The
statement was right.  Here the sign is an integer and every summand is asserted to be exact; the float version is re-run as a
diagnostic and its failing slots are printed with the type and size of the defect.
Run: py -I -B gi_relations.py [K] [NU] [alpha] [beta]      (defaults 20, 66, 0, 0)      log: gi_relations_<frame>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import Q, Ser, gauge_frame, model_vanishing, quintic_vanishing, p11_at_zero, inv_poly, normalise
PAIR = {0: 3, 1: 2, 2: 1, 3: 0, 4: 4}


def sgn(m):
    return 1 if m % 2 == 0 else -1                       # (-1)^m for every integer m, as an int


def zpoly_add(A, B, c=1):
    out = dict(A)
    for e, v in B.items():
        out[e] = out.get(e, Fr(0)) + c * v
    return {e: v for e, v in out.items() if v != 0}


def sympl_defect(g, K, NS):
    """{(k, a1, a2): Ser} of [z^k] (P(z) eta^-1 P(-z)^T - eta^-1), k <= K, entries mod u^NS"""
    out = {}
    for k in range(K + 1):
        acc = {}
        for k1 in range(k + 1):
            k2 = k - k1
            for b in range(5):
                if not g[k1][b].c or not g[k2][PAIR[b]].c: continue
                a1 = (b + k1) % 5; a2 = (PAIR[b] + k2) % 5
                t = g[k1][b].mul(g[k2][PAIR[b]]).truncate(NS).scale(Fr(sgn(k2), 5))
                acc[(a1, a2)] = t if (a1, a2) not in acc else acc[(a1, a2)].add(t)
        if k == 0:
            for a in range(5):
                key = (a, PAIR[a]); t = Ser.const(Q, Fr(-1, 5), NS)
                acc[key] = t if key not in acc else acc[key].add(t)
        for key, s in acc.items():
            s = s.truncate(NS)
            if s.c: out[(k,) + key] = s
    return out


def gstar(P_, theta, eta, K, NMAX, sign):
    """(G*)_n for n = 0..NMAX; returns (number of tested identities, list of failing (n, b, e, lhs, rhs))"""
    star = lambda vec: [{m: 5 * c * sign(m) for m, c in vec[PAIR[a]].items()} for a in range(5)]
    tot = 0; bad = []
    for n in range(NMAX + 1):
        acc = [dict() for _ in range(5)]
        for j in range(n + 1):
            ts = star(theta[n - j])
            for (a, b), pz in P_[j].items():
                for k, c in pz.items():
                    for m, y in ts[a].items(): acc[b][k + m] = acc[b].get(k + m, Fr(0)) + c * y
        es = star(eta[n]); mneg = (5 * n) // 4                  # deepest power of 1/z in theta_{n-j}(-z): floor(5n/4)
        for b in range(5):
            for e in set(acc[b]) | set(es[b]):
                if e + mneg > K: continue
                tot += 1
                l = acc[b].get(e, 0); r = es[b].get(e, 0)
                if l != r: bad.append((n, b, e, l, r))
    return tot, bad


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    NU = int(sys.argv[2]) if len(sys.argv) > 2 else 66
    al = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); be = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "gi_relations_%s.log" % tag), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(name)

    t0 = time.time(); NMAX = 8
    res = gauge_frame(K, NU, al, be); p = res['p']; v = res['v']; bq = res['bq']
    c11 = p11_at_zero(p); pc = normalise(p, inv_poly(c11, K))
    log("gi_relations (v2): frame (%s, %s), canonical gauge to z^%d; v = %s" % (al, be, K, v.fmt(5)))
    NS = min(min(x.prec for x in pk) for pk in p) - 1
    # ---- (S) and (S')
    dc = sympl_defect(pc, K, NS)
    check("(S) the canonical gauge is symplectic: P(z) eta^-1 P(-z)^T = eta^-1, all entries, z^0..z^%d, mod u^%d" % (K, NS), not dc, sorted(dc)[:4])
    ds = sympl_defect(p, K, NS)
    cc = {}
    for k1, x1 in c11.items():
        for k2, x2 in c11.items():
            if k1 + k2 <= K: cc[k1 + k2] = cc.get(k1 + k2, Fr(0)) + x1 * x2 * sgn(k2)
    cc[0] -= 1; cc = {k: x for k, x in cc.items() if x != 0}
    pred = {}
    for k, x in cc.items():
        for a in range(5): pred[(k, a, PAIR[a])] = x / 5
    ok = set(ds) == set(pred) and all(ds[key].sub(Ser.const(Q, pred[key], NS)).truncate(NS).iszero() for key in pred)
    log("   stored gauge: c(z) = P_11(0, z) = %s" % " + ".join("(%s) z^%d" % (c11[k], k) for k in sorted(c11)))
    log("   c(z) c(-z) - 1 = %s   (even in z: no z^5, z^15 terms)" % " + ".join("(%s) z^%d" % (cc[k], k) for k in sorted(cc)))
    check("(S') the stored gauge has exactly the defect (c(z) c(-z) - 1) eta^-1: first at z^%s, constant in u, all entries mod u^%d"
          % (min(k for (k, a, b) in ds) if ds else None, NS), ok and min(k for (k, a, b) in ds) == 10)
    # ---- coefficient matrices in u:  P_[n][(a, b)] = {z-power: Fraction}
    def Pn(n):
        out = {}
        for k in range(K + 1):
            for b in range(5):
                c = pc[k][b].coef(n)
                if c != 0: out.setdefault(((b + k) % 5, b), {})[k] = c
        return out
    P_ = [Pn(n) for n in range(NMAX + 1)]
    YC = model_vanishing(v, NMAX + 1, 2 * NMAX + 4); YQ = quintic_vanishing(bq, NMAX + 1, 2 * NMAX + 4)
    def ucoef(Y, n):
        out = [dict() for _ in range(5)]
        for m, vec in Y.items():
            for a in range(5):
                c = vec[a].coef(n) if n < vec[a].prec else None
                if c: out[a][m] = c
        return out
    eta = [ucoef(YC, n) for n in range(NMAX + 1)]; theta = [ucoef(YQ, n) for n in range(NMAX + 1)]
    fz = lambda d: " + ".join("(%s) z^%d" % (d[e], e) for e in sorted(d, reverse=True)) if d else "0"
    log("   eta_n(z)  (u^n of Y^C(v(u), z)) and theta_n(z) (u^n of Y^C_Q), components a = 0..4:")
    for n in range(0, 6):
        log("      eta_%d   = %s" % (n, [fz(eta[n][a]) for a in range(5)])); log("      theta_%d = %s" % (n, [fz(theta[n][a]) for a in range(5)]))
    def mulvec(Pm, vec):
        out = [dict() for _ in range(5)]
        for (a, b), pz in Pm.items():
            for k, c in pz.items():
                for m, y in vec[b].items(): out[a][k + m] = out[a].get(k + m, Fr(0)) + c * y
        return [{e: x for e, x in d.items() if x != 0} for d in out]
    # ---- (G)_n
    tot = 0; bad = []
    for n in range(1, NMAX + 1):
        rhs = [dict() for _ in range(5)]
        for j in range(n):
            t = mulvec(P_[j], eta[n - j]); rhs = [zpoly_add(rhs[a], t[a]) for a in range(5)]
        rhs = [zpoly_add(rhs[a], theta[n][a], -1) for a in range(5)]
        mneg = (5 * n) // 4
        for a in range(5):
            lhs = P_[n].get((a, 1), {})
            for e in set(lhs) | set(rhs[a]):
                if e + mneg > K: continue
                tot += 1
                if lhs.get(e, 0) != rhs[a].get(e, 0): bad.append((n, a, e))
    check("(G)_n: P^(n) e_1 = Sum_{j<n} P^(j) eta_{n-j} - theta_n for n = 1..%d, all rows, all complete z-powers (%d non-trivial coefficient identities)" % (NMAX, tot), not bad, bad[:4])
    # ---- (G*)_n, exact
    tot, bad = gstar(P_, theta, eta, K, NMAX, sgn)
    exact = all(isinstance(x, (int, Fr)) for (n, b, e, l, r) in bad for x in (l, r))
    check("(G*)_n: Y^C_Q(-z)^T eta P(z) = Y^C(-z)^T eta at u^n, n = 0..%d, all columns, complete z-powers (%d coefficient identities, exact)" % (NMAX, tot),
          not bad, [(n, b, e, str(l), str(r)) for (n, b, e, l, r) in bad[:4]])
    # ---- diagnostic: the first version's sign (-1) ** m (a float for m < 0)
    tot2, bad2 = gstar(P_, theta, eta, K, NMAX, lambda m: (-1) ** m)
    log("   DIAGNOSTIC: the same identities with the first version's sign (-1) ** m: %d of %d fail; failing slots (n, column, z-power, type, |lhs - rhs|): %s"
        % (len(bad2), tot2, [(n, b, e, type(l).__name__, "%.1e" % abs(float(l) - float(r))) for (n, b, e, l, r) in bad2[:6]]))
    check("DIAGNOSIS: every failure of the float version is a float rounding residue below 1e-12 in a slot where the exact identity holds",
          bool(bad2) and all(isinstance(l, float) or isinstance(r, float) for (n, b, e, l, r) in bad2)
          and all(abs(float(l) - float(r)) < 1e-12 for (n, b, e, l, r) in bad2) and not bad)
    # ---- (R1)-(R3), and the explicit row-2 relations at u^3, u^4
    E = lambda n, a, b: P_[n].get((a, b), {})
    log("   (R1) column 1 of P^(0): %s ;  row 2 of P^(0): %s" % ([fz(E(0, a, 1)) for a in range(5)], [fz(E(0, 2, b)) for b in range(5)]))
    check("(R1) = (G)_0 and (G*)_0: P(0, z) e_1 = e_1, e_2^T P(0, z) = e_2^T", all(E(0, a, 1) == ({0: 1} if a == 1 else {}) for a in range(5)) and all(E(0, 2, b) == ({0: 1} if b == 2 else {}) for b in range(5)))
    pi = E(0, 3, 0)
    sh = lambda d, s, c=1: {e + s: c * x for e, x in d.items()}
    lim = lambda d, top: {e: x for e, x in d.items() if e <= top}
    log("   pi(z) = P_30(0, z) = %s" % fz(pi))
    check("(R2) = row 2 of (G)_1: P^(1)_21 = -P^(0)_20/z = 0", E(1, 2, 1) == {})
    check("(R3b) = row 3 of (G)_1: P^(1)_31 = -pi/z", lim(E(1, 3, 1), K - 2) == lim(sh(pi, -1, -1), K - 2))
    check("(R3a) = column 0 of (G*)_1: P^(1)_20 = +pi/z", lim(E(1, 2, 0), K - 2) == lim(sh(pi, -1), K - 2))
    check("(R3c) = row 2 of (G)_2: P^(2)_21 = -P^(1)_20/z = -pi/z^2", lim(E(2, 2, 1), K - 3) == lim(sh(pi, -2, -1), K - 3))
    v2, v3 = v.coef(2), v.coef(3)
    r3 = zpoly_add(zpoly_add(sh(E(2, 2, 0), -1, -1), sh(E(1, 2, 0), -1, -v2)), sh(E(1, 2, 4), -2, Fr(1, 2)))
    log("   u^3, row 2:  P^(3)_21 = -P^(2)_20/z - v_2 P^(1)_20/z + P^(1)_24/(2 z^2)")
    log("        P^(3)_21 = %s" % fz(lim(E(3, 2, 1), 8))); log("        rhs      = %s" % fz(lim(r3, 8)))
    check("row 2 of (G)_3", lim(E(3, 2, 1), K - 4) == lim(r3, K - 4))
    r4 = zpoly_add(sh(E(3, 2, 0), -1, -1), sh(E(2, 2, 0), -1, -v2))
    r4 = zpoly_add(r4, sh(E(1, 2, 0), -1, -v3)); r4 = zpoly_add(r4, sh(E(2, 2, 4), -2, Fr(1, 2))); r4 = zpoly_add(r4, sh(E(1, 2, 4), -2, v2))
    r4 = zpoly_add(r4, sh(E(1, 2, 3), -3, Fr(-1, 6))); r4 = zpoly_add(r4, sh(zpoly_add({0: Fr(1)}, E(0, 2, 2), -1), -4, Fr(1, 24)))
    log("   u^4, row 2:  P^(4)_21 = -P^(3)_20/z - v_2 P^(2)_20/z - v_3 P^(1)_20/z + P^(2)_24/(2 z^2) + v_2 P^(1)_24/z^2 - P^(1)_23/(6 z^3) + (1 - P^(0)_22)/(24 z^4)")
    log("        P^(4)_21 = %s" % fz(lim(E(4, 2, 1), 8))); log("        rhs      = %s" % fz(lim(r4, 8)))
    check("row 2 of (G)_4 (the term (1 - P_22(0, z))/(24 z^4) vanishes: canonical normalisation)", lim(E(4, 2, 1), K - 6) == lim(r4, K - 6))
    # the same u^4 relation in the STORED normalisation fails exactly by the normalisation term
    log("   in the stored normalisation (1 - P_22(0, z))/(24 z^4) = %s" % fz(lim(sh({k: -x for k, x in c11.items() if k}, -4, Fr(1, 24)), 6)))
    log("   u^5, row 2 (general formula): P^(5)_21 = %s" % fz(lim(E(5, 2, 1), 9)))
    log("GI_RELATIONS %s: %s   (%.0f s)" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
