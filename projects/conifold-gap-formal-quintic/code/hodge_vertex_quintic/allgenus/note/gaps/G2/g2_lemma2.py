# -*- coding: utf-8 -*-
"""g2_lemma2.py -- evidence for lem:gau-scal (G3): the two projection facts
   pr((Y^j)') = -(j/5) (log calB)' Y^j,     pr(Theta^-1 Theta' Y^j) = (tr(Theta^-1 Theta')/5) Y^j     (Y = At^-1, j = 0..4),
for the pulled-back model At and the quintic (Theta = Theta_U = diag(d_a), d_0 = 1, d_(a+1) = d_a bt_a/b_a), exact over Q in
truncated Laurent series.  pr is computed here from the trace form, pr(X) = Sum_m [tr(X At^(5-m)) / tr(At^m At^(5-m))] At^m
(matrix products and traces, not the closed formula), and compared with the right-hand sides.  Also: the residue j/5 of
-(j/5)(log calB)', and pr(Z) of a holomorphic Z has at most a simple pole in the Y^j-components (j >= 1).
Controls: the factor j/4; replacing Theta^-1 Theta' by a non-diagonal holomorphic matrix (the trace rule fails).
Usage: py -I -B g2_lemma2.py [N].  Stdlib only; imports g2_core.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2_core import Ser, systems, prod

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def mat_zero(N):
    return [[Ser(N, [], N) for _ in range(5)] for _ in range(5)]


def mm(X, Y):
    return [[sum((X[i][k] * Y[k][j] for k in range(5)), Ser(10 ** 6, [], 10 ** 6)) for j in range(5)] for i in range(5)]


def tr(X):
    return sum((X[i][i] for i in range(5)), Ser(10 ** 6, [], 10 ** 6))


def add(X, Y, c=1):
    return [[X[i][j] + Y[i][j].scale(Fr(c)) for j in range(5)] for i in range(5)]


def scal(s, X):
    return [[s * X[i][j] for j in range(5)] for i in range(5)]


def deriv(X):
    return [[X[i][j].d() for j in range(5)] for i in range(5)]


def iszero(X):
    return all(X[i][j].iszero() for i in range(5) for j in range(5))


def shift_part(X, m):
    out = mat_zero(10 ** 6)
    for a in range(5):
        out[a][(a + m) % 5] = X[a][(a + m) % 5]
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print("g2_lemma2.py: exact over Q, u-precision %d, frame (0,0)" % N)
    S = systems(N, 0, 0)
    b, bt = S["b"], S["bt"]
    At = mat_zero(N); Y = mat_zero(N)
    for a in range(5):
        At[a][(a + 1) % 5] = bt[a]
        Y[(a + 1) % 5][a] = bt[a].inv()
    one = [[Ser.const(int(i == j), N) for j in range(5)] for i in range(5)]
    Apw = [one]
    for m in range(1, 6):
        Apw.append(mm(Apw[-1], At))
    Ypw = [one]
    for j in range(1, 5):
        Ypw.append(mm(Ypw[-1], Y))
    B = prod(bt)
    check("At Y = 1 and Y holomorphic", iszero(add(mm(At, Y), one, -1)) and all(Y[i][j].v >= 0 for i in range(5) for j in range(5)))
    norms = [tr(mm(Apw[m], Apw[(5 - m) % 5])) for m in range(5)]       # tr(At^m At^(5-m)) = 5 B (m != 0), 5 (m = 0)

    def pr(X):
        out = mat_zero(N)
        for m in range(5):
            lam = tr(mm(shift_part(X, m), Apw[(5 - m) % 5])) * norms[m].inv()
            out = add(out, scal(lam, Apw[m]))
        return out
    d = [Ser.const(1, N)]
    for a in range(4):
        d.append(d[a] * bt[a] * b[a].inv())
    Dl = mat_zero(N)
    for a in range(5):
        Dl[a][a] = d[a].d() * d[a].inv()
    trDl = tr(Dl)
    dlogB = B.d() * B.inv()
    for j in range(5):
        lhs = pr(deriv(Ypw[j]))
        rhs = scal(dlogB.scale(Fr(-j, 5)), Ypw[j])
        check("pr((Y^%d)') = -(%d/5)(log calB)' Y^%d" % (j, j, j), iszero(add(lhs, rhs, -1)))
        if j > 0:
            ctl = scal(dlogB.scale(Fr(-j, 4)), Ypw[j])
            check("control: -(%d/4)(log calB)' Y^%d differs" % (j, j), not iszero(add(lhs, ctl, -1)))
        lhs2 = pr(mm(Dl, Ypw[j]))
        rhs2 = scal(trDl.scale(Fr(1, 5)), Ypw[j])
        check("pr(Theta^-1 Theta' Y^%d) = (tr/5) Y^%d" % (j, j), iszero(add(lhs2, rhs2, -1)))
    kap = dlogB.scale(Fr(-1, 5))
    check("residue of -(1/5)(log calB)' is 1/5 (so -(j/5)(log calB)' = j/(5u) + holomorphic)", kap.resid() == Fr(1, 5) and kap.v == -1)
    check("tr(Theta^-1 Theta') holomorphic", trDl.v >= 0)
    Z = mat_zero(N)
    for i in range(5):
        for k in range(5):
            Z[i][k] = Ser(0, [Fr(i + 2 * k + 1, 3), Fr(k - i, 2), Fr(1)], N)
    PZ = pr(Z)
    poles = [min(PZ[i][k].v for i in range(5) for k in range(5) if (i - k) % 5 == j) for j in range(5)]
    ok = True
    for j in range(1, 5):            # coefficient of Y^j = (Y^j-shaped entry) * calB_a^(j): at most a simple pole
        a = 2                        # window {2..2+j-1} avoids the pole index 1, so calB_2^(j) is a unit
        cj = PZ[(a + j) % 5][a] * prod([bt[(a + i) % 5] for i in range(j)])
        ok = ok and cj.v >= -1
    c0 = PZ[0][0]
    check("pr(Z), Z holomorphic: Y^j-coefficients have at most a simple pole (j >= 1), Y^0-coefficient holomorphic",
          ok and c0.v >= 0, "(entry valuations by shift %s)" % poles)
    Wnd = mat_zero(N)
    Wnd[0][1] = Ser.const(1, N)
    lhs3 = pr(mm(add(Dl, Wnd), Ypw[1]))
    rhs3 = scal(tr(add(Dl, Wnd)).scale(Fr(1, 5)), Ypw[1])
    check("control: for a non-diagonal matrix in place of Theta^-1 Theta' the trace rule fails", not iszero(add(lhs3, rhs3, -1)))
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
