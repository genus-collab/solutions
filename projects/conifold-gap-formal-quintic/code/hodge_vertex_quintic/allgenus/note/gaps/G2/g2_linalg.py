# -*- coding: utf-8 -*-
"""g2_linalg.py -- evidence for lem:gau-comm (G1, G2): the linear algebra of a weighted cyclic shift At = Sum_a bt_a E_(a,a+1)
over a field, tested exactly over Q at random points bt_a in Q^x (the statements are identities of rational functions in the
bt_a, used in the note over K = Q((u))).  Checks: characteristic polynomial x^5 - B (B = prod bt_a); rank ad At = 20; the
commutant is span(1, At, .., At^4); the trace form on Q[At] (Gram determinant 5^5 B^4); the projection formula
pr(X) = (1/5) [Sum_a xi_a / calB_a^(m)] At^m on the shift-m part, with X - pr(X) in im(ad At); its y-form; the lifting
formula of lem:gau-lift.  Controls: the factor 1/4 instead of 1/5; a vector not in the image.
Usage: py -I -B g2_linalg.py [number of random points].  Stdlib only.
"""
import sys, random
from fractions import Fraction as Fr

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def zero():
    return [[Fr(0)] * 5 for _ in range(5)]


def mm(X, Y):
    return [[sum(X[i][k] * Y[k][j] for k in range(5)) for j in range(5)] for i in range(5)]


def lin(*terms):
    out = zero()
    for c, X in terms:
        for i in range(5):
            for j in range(5):
                out[i][j] += c * X[i][j]
    return out


def rank(rows):
    rows = [list(r) for r in rows]
    rk = 0; ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def flat(X):
    return [X[i][j] for i in range(5) for j in range(5)]


def det(M):
    M = [list(r) for r in M]; n = len(M); d = Fr(1)
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]; d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def main():
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    rnd = random.Random(20261115)
    print("g2_linalg.py: %d random points, exact over Q" % npts)
    for t in range(npts):
        bt = [Fr(rnd.choice([-1, 1]) * rnd.randint(1, 9), rnd.randint(1, 9)) for _ in range(5)]
        B = bt[0] * bt[1] * bt[2] * bt[3] * bt[4]
        At = zero()
        for a in range(5):
            At[a][(a + 1) % 5] = bt[a]
        I = [[Fr(int(i == j)) for j in range(5)] for i in range(5)]
        pw = [I]
        for m in range(1, 6):
            pw.append(mm(pw[-1], At))
        calB = lambda a, m: eval_prod(bt, a, m)
        ok_cp = all(det(lin((Fr(x), I), (Fr(-1), At))) == Fr(x) ** 5 - B for x in (-2, 0, 1, 3, Fr(1, 2), 7))
        check("pt %d: det(x - At) = x^5 - B at six x" % t, ok_cp)
        ad = []
        for i in range(5):
            for j in range(5):
                E = zero(); E[i][j] = Fr(1)
                ad.append(flat(lin((Fr(1), mm(At, E)), (Fr(-1), mm(E, At)))))
        adT = [list(col) for col in zip(*ad)]       # 25 x 25: columns = images of the basis
        check("pt %d: rank ad At = 20" % t, rank(ad) == 20)
        check("pt %d: At^m commute with At, span of dim 5; At^5 = B" % t,
              all(all(x == 0 for x in flat(lin((Fr(1), mm(At, pw[m])), (Fr(-1), mm(pw[m], At))))) for m in range(5))
              and rank([flat(pw[m]) for m in range(5)]) == 5 and pw[5] == lin((B, I)))
        G = [[sum(mm(pw[i], pw[j])[k][k] for k in range(5)) for j in range(5)] for i in range(5)]
        check("pt %d: trace form on Q[At]: Gram det = 5^5 B^4" % t, det(G) == 5 ** 5 * B ** 4)
        X = [[Fr(rnd.randint(-9, 9), rnd.randint(1, 5)) for _ in range(5)] for _ in range(5)]

        def pr(X, fac=Fr(1, 5)):
            out = zero()
            for m in range(5):
                lam = fac * sum(X[a][(a + m) % 5] / calB(a, m) for a in range(5))
                out = lin((Fr(1), out), (lam, pw[m]))
            return out
        R = lin((Fr(1), X), (Fr(-1), pr(X)))
        check("pt %d: X - pr(X) lies in im(ad At) (rank stays 20)" % t, rank(ad + [flat(R)]) == 20)
        R4 = lin((Fr(1), X), (Fr(-1), pr(X, Fr(1, 4))))
        check("pt %d: control: with 1/4 in place of 1/5, X - pr(X) is NOT in the image" % t, rank(ad + [flat(R4)]) == 21)
        check("pt %d: control: a non-zero element of Q[At] is not in the image" % t, rank(ad + [flat(pw[2])]) == 21)
        Y = [[Fr(0)] * 5 for _ in range(5)]
        for a in range(5):
            Y[(a + 1) % 5][a] = 1 / bt[a]
        ypw = [I]
        for j in range(1, 5):
            ypw.append(mm(ypw[-1], Y))
        ok_y = True
        for m in range(1, 5):
            lhs = lin(((Fr(1, 5) * sum(X[a][(a + m) % 5] * B / calB(a, m) for a in range(5))), ypw[5 - m]))
            rhs = lin(((Fr(1, 5) * sum(X[a][(a + m) % 5] / calB(a, m) for a in range(5))), pw[m]))
            ok_y = ok_y and lhs == rhs
        check("pt %d: y-form (1/5)[Sum xi_a B/calB_a^(m)] y^(5-m) = the At^m form; supports of y^j disjoint" % t,
              ok_y and all(sum(1 for j in range(5) if ypw[j][r][c] != 0) == 1 for r in range(5) for c in range(5)))
        ok_l = True
        for d in range(5):                          # lifting: shift d+1 part of R (pr(R) = 0) from Pi of shift d
            xi = [R[a][(a + d + 1) % 5] for a in range(5)]
            sg = [xi[a] / calB(a, d + 1) for a in range(5)]
            ok_l = ok_l and sum(sg) == 0
            rho = {1: Fr(0)}
            for a in (1, 2, 3, 4):
                rho[(a + 1) % 5] = rho[a] + sg[a]
            Pi = zero()
            for a in range(5):
                Pi[a][(a + d) % 5] = rho[a] * calB(a, d)
            C = lin((Fr(1), mm(At, Pi)), (Fr(-1), mm(Pi, At)))
            ok_l = ok_l and all(C[a][(a + d + 1) % 5] == xi[a] for a in range(5))
        check("pt %d: lifting formula pi_a = rho_a calB_a^(d), rho_(a+1) - rho_a = xi_a/calB_a^(d+1) solves [At, Pi] = X - pr X" % t, ok_l)
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


def eval_prod(bt, a, m):
    out = Fr(1)
    for i in range(m):
        out *= bt[(a + i) % 5]
    return out


if __name__ == "__main__":
    main()
