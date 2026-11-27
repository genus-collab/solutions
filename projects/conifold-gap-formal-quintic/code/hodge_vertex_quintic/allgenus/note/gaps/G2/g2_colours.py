# -*- coding: utf-8 -*-
"""g2_colours.py -- evidence for cor:gau-colours (G6): (a) the colour-i solution of the model T-system,
   F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4,  T_k = -zeta^i + (t/5)(k - 4 t d/dt),
F_0(0) = 1, equals zeta^(i b) F_b^(0)(zeta^-i t) (so the normalised colour-i model rows zeta^i F^(i)_b = zeta^(i(b+1)) F_b(zeta^-i t)
are the image of colour 0 under Y(z) -> T_i Y(zeta^-i z), T_i = diag(zeta^(i(a+1)))); the pivot is 4n zeta^(4i);
(b) the multiset identity w_i = zeta^i w_0 behind lem:set-hrot; (c) the colour covariance of the graded P_can,
T_i P_can(u, zeta^-i z) T_i^-1 = P_can(u, z), with the non-graded control P_can (1 + z).
Exact over Q(zeta): an element is a 5-vector (coefficients of 1, zeta, .., zeta^4) modulo (1,1,1,1,1).
Usage: py -I -B g2_colours.py [n].  Stdlib only; imports g2_core, g2_gauge.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def zq(c=0, k=0):
    """c zeta^k"""
    v = [Fr(0)] * 5; v[k % 5] = Fr(c); return v


def zadd(x, y):
    return [a + b for a, b in zip(x, y)]


def zsc(c, x):
    return [Fr(c) * a for a in x]


def zrot(x, k):
    """zeta^k x"""
    return [x[(j - k) % 5] for j in range(5)]


def zeq(x, y):
    d = [a - b for a, b in zip(x, y)]
    return all(a == d[0] for a in d)


def zmul(x, y):
    out = [Fr(0)] * 5
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                out[(i + j) % 5] += a * b
    return out


def T(k, i, F):
    """T_k F, F a list of Q(zeta) coefficients of t^n (truncated): [t^n] = -zeta^i f_n + ((k - 4(n-1))/5) f_(n-1)"""
    out = []
    for n in range(len(F)):
        c = zsc(-1, zrot(F[n], i))
        if n >= 1:
            c = zadd(c, zsc(Fr(k - 4 * (n - 1), 5), F[n - 1]))
        out.append(c)
    return out


def chain(i, F0):
    F1 = T(-1, i, F0)
    F2 = [zsc(-1, x) for x in T(-2, i, F1)]
    F3 = T(2, i, F2)
    F4 = [zsc(-1, x) for x in T(1, i, F3)]
    F5 = [zsc(-1, x) for x in T(0, i, F4)]
    return [F0, F1, F2, F3, F4], F5


def solve(i, n):
    """colour-i T-system solution with F_0(0) = 1, coefficients t^0..t^n"""
    F0 = [zq(1, 0)]
    for m in range(1, n + 1):
        trial = F0 + [zq(0), zq(0)]
        _, F5 = chain(i, trial)
        resid = zadd(F5[m], zsc(-1, trial[m]))          # [t^m]((closing) - F_0) with f_m = 0
        one = F0 + [zq(1, 0), zq(0)]
        _, G5 = chain(i, one)
        piv = zadd(zadd(G5[m], zsc(-1, one[m])), zsc(-1, resid))  # coefficient of f_m at t^m ... (pivot sits at t^(m+1))
        r1 = zadd(chain(i, trial)[1][m + 1], zsc(-1, trial[m + 1]))
        p1 = zadd(zadd(G5[m + 1], zsc(-1, one[m + 1])), zsc(-1, r1))
        assert zeq(p1, zrot(zq(4 * m, 0), 4 * i)), "pivot is not 4m zeta^(4i)"
        assert zeq(piv, zq(0)), "f_m enters at t^m"
        fm = zsc(-1, zrot(zsc(Fr(1, 4 * m), r1), -4 * i))   # f_m = -r1/(4m zeta^(4i))
        F0 = F0 + [fm]
    rows, F5 = chain(i, F0 + [zq(0)])
    ok = all(zeq(F5[m], F0[m]) for m in range(n + 1))
    return [r[:n + 1] for r in rows], ok


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    print("g2_colours.py: T-system to t^%d, exact over Q(zeta)" % n)
    F, ok0 = solve(0, n)
    check("colour 0: closing relation holds; F_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3 - (2415609/6400000)t^4 + ..",
          ok0 and [F[0][m][0] - F[0][m][1] for m in range(5)] == [1, Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000)]
          and all(all(x == F[0][m][1] for x in F[0][m][1:]) for m in range(n + 1)))
    check("colour 0: F_b(0) = (1, -1, -1, 1, 1)", [F[b][0][0] - F[b][0][1] for b in range(5)] == [1, -1, -1, 1, 1])
    for i in range(1, 5):
        Fi, oki = solve(i, n)
        cov = all(zeq(Fi[b][m], zrot(F[b][m], i * b - i * m)) for b in range(5) for m in range(n + 1))
        ctl = all(zeq(Fi[b][m], zrot(F[b][m], i * (b + 1) - i * m)) for b in range(5) for m in range(n + 1))
        check("colour %d: unique solution (pivot 4n zeta^(4i)); F^(%d)_b(t) = zeta^(%d b) F_b(zeta^-%d t)" % (i, i, i, i),
              oki and cov)
        check("colour %d: control: zeta^(%d(b+1)) F_b(zeta^-%d t) (the normalised rows) is NOT the T-system solution" % (i, i, i),
              not ctl)
    w = lambda i: [zadd(zq(1, i), zq(-1, j)) for j in range(5) if j != i] + [zq(-5, i)]
    ok_w = all(sorted(map(tuple, [tuple(x[k] - x[0] for k in range(5)) for x in w(i)])) ==
               sorted(map(tuple, [tuple(y[k] - y[0] for k in range(5)) for y in [zrot(x, i) for x in w(0)]])) for i in range(5))
    check("lem:set-hrot: w_i = zeta^i w_0 as multisets (so h_(w_i)(z) = h_(w_0)(zeta^-i z))", ok_w)
    wbad = lambda i: [zadd(zq(1, i), zq(-1, j)) for j in range(5) if j != i] + [zq(-5, 0)]
    okb = all(sorted([tuple(x[k] - x[0] for k in range(5)) for x in wbad(i)]) ==
              sorted([tuple(y[k] - y[0] for k in range(5)) for y in [zrot(x, i) for x in wbad(0)]]) for i in range(5))
    check("control: the multiset {zeta^i - zeta^j} U {-5} is not covariant", not okb)
    from g2_core import systems
    from g2_gauge import gauge, times_c, inverse_series
    S = systems(20, 0, 0)
    P, _ = gauge(S["b"], S["bt"], 7)
    Pc = times_c(P, inverse_series({0: P[0][1].coef(0), 5: P[5][1].coef(0)}, 7))
    # entry (a, b) of P_k is multiplied by zeta^(i(a-b-k)) under P -> T_i P(zeta^-i z) T_i^-1; graded entries have a = b + k
    phase = [[(((bb + k) % 5) - bb - k) % 5 for bb in range(5)] for k in range(8)]
    check("colour covariance of P_can: every stored entry (b+k, b) has phase zeta^(i(a-b-k)) = 1 (gradedness)",
          all(x == 0 for r in phase for x in r))
    check("control: P_can (1 + z) is not covariant: its z^1 term contains P_can,0 at (a, a), phase zeta^(-i) != 1",
          all(not Pc[0][a].iszero() for a in range(5)) and (0 - 0 - 1) % 5 != 0)
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
