# -*- coding: utf-8 -*-
"""v2_gr_and_frame.py -- ADVERSARIAL checks of the proof ingredients of LG_MODEL.md s.4 and of the canonical-frame statements.
 (G) the grading deg*(z^m L^e D^r) = e + 5r - 4m on the conjugated Picard-Fuchs operator Pcal of fq_rows.WKB (the operator the project
     really uses): max deg* = 0, and the deg* = 0 part equals  u (z Dbar + l)^5 - 1 - 3 z (z Dbar + l)^4,  u = L^-5, Dbar L^k = (k/5) L^{k+5},
     built here with an independent small operator algebra.  Next non-zero degree is -5.
 (C) conifold-frame leading terms Cf(C0) = varpi0 = 1 + O(u), Cf(C1) = D t_c = -1 + O(u) for six frames (alpha, beta), using only
     fq_conifold.Conifold (.varpi0, .tc) and D = -(1-u) d/du.
 (E) the canonical-frame equation (5/4)[C_{n+1}, Z] = (V - n) C_n with the PROJECT's own M_0 (fq3_par.setup -> M0.M[k], k <= 6), not the
     report's dump; and the alternative sign/transposition variants as controls (they must fail).
 (U) uniqueness: a random circulant-covariant perturbation is not needed; instead the recursion is SOLVED here from (U, V) alone
     (off-diagonal from order n, diagonal from order n+1) and compared with M0.M[k], k <= 6.
 (F) genus one: C_1[i][i] from the project's M0 and dF_1 = (1/24 - (5/2) C_1[0][0]) du/u.
Run:  py -I v2_gr_and_frame.py
"""
import sys
sys.dont_write_bytecode = True
import os, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from fq_rows import WKB
from fq_conifold import Conifold

MU = (Fr(-1, 4), Fr(-1, 2), Fr(1, 2), Fr(1, 4), Fr(0))


# ---------------------------------------------------------------- graded operator algebra: {(m, r): {e: c}} = Sum z^m c L^e Dbar^r
def g_add(A, B, p, c=1):
    out = {k: dict(v) for k, v in A.items()}
    for k, v in B.items():
        d = out.setdefault(k, {})
        for e, x in v.items():
            y = (d.get(e, 0) + c * x) % p
            if y: d[e] = y
            elif e in d: del d[e]
        if not d: del out[k]
    return out


def g_mul(A, B, p, inv5):
    """composition with Dbar L^k = (k/5) L^{k+5} + L^k Dbar"""
    import math
    out = {}
    for (m, r), a in A.items():
        for (n, s), b in B.items():
            Dt = dict(b)
            for t in range(r + 1):
                key = (m + n, r - t + s); d = out.setdefault(key, {})
                for e1, x in a.items():
                    for e2, y in Dt.items():
                        v = (d.get(e1 + e2, 0) + math.comb(r, t) * x * y) % p
                        if v: d[e1 + e2] = v
                        elif e1 + e2 in d: del d[e1 + e2]
                if not d: del out[key]
                Dt = {e + 5: y * e % p * inv5 % p for e, y in Dt.items() if e % p}
    return out


def main():
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    for tag, p in (("89", P2), ("61", P1)):
        R = FpRing(p); inv5 = R.inv(5)
        for i in (0, 1, 3):
            w = WKB(R, i, 1); P = w.P; l = R.z(i)
            degs = {}
            for (m, r), a in P.items():
                for e, c in a.items(): degs.setdefault(e + 5 * r - 4 * m, []).append(((m, r), e, c))
            gr0 = {}
            for ((m, r), e, c) in degs.get(0, []): gr0.setdefault((m, r), {})[e] = c
            Ml = {(1, 1): {0: 1}, (0, 0): {1: l}}
            pw = [{(0, 0): {0: 1}}]
            for _ in range(5): pw.append(g_mul(pw[-1], Ml, p, inv5))
            uop = {(0, 0): {-5: 1}}
            pred = g_mul(uop, pw[5], p, inv5)
            pred = g_add(pred, {(0, 0): {0: 1}}, p, -1)
            pred = g_add(pred, g_mul({(1, 0): {0: 3}}, pw[4], p, inv5), p, -1)
            check("(G) p = 2^%s-1, colour %d: max deg*(Pcal) = 0, next degree -5; gr_0(Pcal) = u (z Dbar + l)^5 - 1 - 3 z (z Dbar + l)^4"
                  % (tag, i), max(degs) == 0 and sorted(degs, reverse=True)[1] == -5 and gr0 == pred,
                  "degrees present: %s" % sorted(degs, reverse=True)[:4])
            # control: the coefficient 3 is essential
            pred2 = g_add(pred, g_mul({(1, 0): {0: 1}}, pw[4], p, inv5), p, 1)
            check("(G) control: with 2 z (..)^4 in place of 3 z (..)^4 the comparison fails", gr0 != pred2)
    # ---------------------------------------------------------------- (C)
    R = FpRing(P2)
    for (al, be) in ((0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(5), Fr(-11, 3)), (Fr(-7, 2), Fr(9)), (Fr(1), Fr(1)), (Fr(-13, 5), Fr(0))):
        C = Conifold(R, Nu=10, alpha=Fr(al), beta=Fr(be))
        v0 = C.varpi0; tc = C.tc
        dtc = [-(k + 1) * tc.coef(k + 1) + k * tc.coef(k) for k in range(3)]          # D t_c = -(1-u) t_c'
        check("(C) frame (alpha, beta) = (%s, %s): varpi0 = 1 + alpha u + O(u^2), t_c = u + O(u^2), D t_c = -1 + O(u)" % (al, be),
              v0.coef(0) == 1 and v0.coef(1) == Fr(al) and tc.coef(0) == 0 and tc.coef(1) == 1 and dtc[0] == -1,
              "varpi0 = %s, %s; D t_c = %s" % (v0.coef(0), v0.coef(1), [str(x) for x in dtc]))
    # ---------------------------------------------------------------- (E), (U), (F): the project's M_0
    from fq3_par import setup
    for tag, p in (("89", P2), ("61", P1)):
        st = setup((p, "0", "0", 34, 6)); R = st['R']; M = st['M0'].M; z = R.zeta
        f = R.of; zp = [pow(z, k, p) for k in range(5)]; inv5 = R.inv(5)
        V = [[sum(zp[((i - j) * (a + 1)) % 5] * f(MU[a]) for a in range(5)) * inv5 % p for i in range(5)] for j in range(5)]
        def test(sign_n, transpose, coef):
            bad = 0
            for n in range(6):
                A = M[n]; B = M[n + 1]
                for j in range(5):
                    for i in range(5):
                        lhs = f(coef) * B[j][i] * (zp[i] - zp[j]) % p
                        if transpose: rhs = (sum(A[j][k] * V[k][i] for k in range(5)) + sign_n * n * A[j][i]) % p
                        else: rhs = (sum(V[j][k] * A[k][i] for k in range(5)) + sign_n * n * A[j][i]) % p
                        if lhs != rhs: bad += 1
            return bad
        check("(E) p = 2^%s-1: project M_0 satisfies (5/4)[C_{n+1}, Z] = (V - n) C_n, n <= 5" % tag, test(-1, False, Fr(5, 4)) == 0)
        check("(E) controls fail: (V + n), right action C_n V, coefficient -5/4",
              test(1, False, Fr(5, 4)) > 0 and test(-1, True, Fr(5, 4)) > 0 and test(-1, False, Fr(-5, 4)) > 0)
        # (U) solve the recursion from (Z, V) alone
        C = [[[1 if i == j else 0 for i in range(5)] for j in range(5)]]
        c54 = f(Fr(5, 4))
        for n in range(6):
            A = C[n]; B = [[0] * 5 for _ in range(5)]
            for j in range(5):
                for i in range(5):
                    if i == j: continue
                    rhs = (sum(V[j][k] * A[k][i] for k in range(5)) - n * A[j][i]) % p
                    B[j][i] = rhs * R.inv(c54 * (zp[i] - zp[j]) % p) % p
            for i in range(5):                                   # (n+1) B_ii = Sum_{k != i} V_ik B_ki   (V_ii = 0)
                B[i][i] = sum(V[i][k] * B[k][i] for k in range(5) if k != i) * R.inv(n + 1) % p
            C.append(B)
        check("(U) p = 2^%s-1: the recursion solved from (Z, V, mu) alone reproduces the project's M_0 through zhat^6 (150 entries); V_ii = 0"
              % tag, all(C[n] == [[x % p for x in r] for r in M[n]] for n in range(7)) and all(V[i][i] == 0 for i in range(5)))
        c11 = M[1][0][0]
        check("(F) p = 2^%s-1: C_1[i][i] = zeta^{-i}/20 and 1/24 - (5/2) C_1[0][0] = -1/12" % tag,
              all(M[1][i][i] == f(Fr(1, 20)) * pow(z, (-i) % 5, p) % p for i in range(5)) and (f(Fr(1, 24)) - f(Fr(5, 2)) * c11) % p == f(Fr(-1, 12)))
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
