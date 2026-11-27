# -*- coding: utf-8 -*-
"""g2_G.py -- evidence for thm:gau-G Cor. (a) (G9) and the identification of the constructed gauge with the canonical one.
(1) the gauge identity P_can Y^C(v(u), z) = Y^C_Q in Q((z))[[u]] (z-powers e <= 3, mod u^NC), for the gauge built by
g2_gauge.gauge (normalised by P_11(0, z) = 1); (2) control: P_can (1 + c z^5) gives (1 + c z^5) Y^C_Q, which differs from Y^C_Q
already at u^0 (the uniqueness argument of G9); (3) the values pi(z) = P_30(0, z) and p_(1,b)(0) quoted in the sources.
Usage: py -I -B g2_G.py [K] [N] [NC].  Stdlib only; imports g2_core, g2_gauge.
"""
import os, sys
from fractions import Fraction as Fr
from math import factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2_core import Ser, systems
from g2_gauge import gauge, times_c, inverse_series

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)
FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def YC_model(v, NC):
    """Y^C(v(u), z) as {n: [5 Ser]} = coefficient of z^(-n); Y^C_b = eps_b I^C_(nu_b), I^C_m = Sum_i (-1/4z)^i/i! (-v/z)^j/j!"""
    out = {}
    vpow = [Ser.const(1, NC)]
    for j in range(1, NC):
        vpow.append(vpow[-1] * v)
    for bb in range(5):
        m = NU[bb]
        i = 0
        while True:
            j = m + 4 * i + 1
            if j >= NC:
                break
            if j >= 0:
                coef = Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** j, factorial(j)) * EPS[bb]
                n = i + j
                out.setdefault(n, [Ser(NC, [], NC) for _ in range(5)])
                out[n][bb] = out[n][bb] + vpow[j].scale(coef)
            i += 1
    return out


def integrate0(f):
    assert f.v >= 0 or f.iszero(), "integrand has a pole"
    return Ser(f.v + 1 if not f.iszero() else f.N + 1,
               [x / (f.v + 1 + k) for k, x in enumerate(f.c)], f.N + 1)


def YC_quintic(b, NC, nmax):
    """Y^C_Q = Sum_n z^(-n) y^(n), y^(0) = -e_1, y^(n+1)_a = Int_0^u b_a y^(n)_(a+1)"""
    y = [Ser(NC, [], NC) for _ in range(5)]
    y[1] = Ser.const(-1, NC)
    out = {0: y}
    for n in range(nmax):
        y = [integrate0(b[a] * y[(a + 1) % 5]) for a in range(5)]
        y = [Ser(x.v, x.c, min(x.N, NC)) for x in y]
        out[n + 1] = y
    return out


def apply_gauge(P, Y, e, NC):
    """[z^e] (P Y): Sum_k P_k [z^(e-k)] Y, P graded (entry (b+k, b) = P[k][b]); Y = {n: vector of z^(-n)}"""
    res = [Ser(NC, [], NC) for _ in range(5)]
    for k in range(len(P)):
        n = k - e
        if n not in Y:
            continue
        for bb in range(5):
            res[(bb + k) % 5] = res[(bb + k) % 5] + P[k][bb] * Y[n][bb]
    return [Ser(x.v, x.c, min(x.N, NC)) for x in res]


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 32
    NC = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    print("g2_G.py: K = %d, N = %d, check mod u^%d (exact over Q)" % (K, N, NC))
    emax = K + 1 - (5 * NC + 3) // 4
    print("  z-powers tested: e <= %d (complete: omitted z-orders k > K enter only at u^>= NC)" % emax)
    for (al, be) in [(0, 0), (Fr(1, 3), Fr(-2, 7))]:
        S = systems(N, al, be)
        b, bt, v = S["b"], S["bt"], S["v"]
        P, _ = gauge(b, bt, K)
        c11 = {5 * q: P[5 * q][1].coef(0) for q in range(K // 5 + 1)}
        Pc = times_c(P, inverse_series(c11, K))
        print(" frame (%s, %s): pi(z) = P_30(0,z) = %s z^3 + %s z^8 + ..; p_(1,b)(0) = %s" % (
            al, be, Pc[3][0].coef(0), Pc[8][0].coef(0), [str(Pc[1][x].coef(0)) for x in range(5)]))
        YC = YC_model(v, NC)
        YQ = YC_quintic(b, NC, (5 * NC) // 4 + 2)
        check("Y^C(0, z) = -e_1 and Y^C_Q(0) = -e_1", [x.coef(0) for x in YC[0]] == [0, -1, 0, 0, 0]
              and [x.coef(0) for x in YQ[0]] == [0, -1, 0, 0, 0])
        check("[z^-1] Y^C_Q = -t_c e_0", (YQ[1][0] + S["tc"]).iszero() and all(YQ[1][a].iszero() for a in range(1, 5)))
        bad = 0; tested = 0
        for e in range(-((5 * NC) // 4) - 1, emax + 1):
            lhs = apply_gauge(Pc, YC, e, NC)
            rhs = YQ.get(-e, [Ser(NC, [], NC)] * 5) if e <= 0 else [Ser(NC, [], NC)] * 5
            for a in range(5):
                tested += 1
                if not (lhs[a] - rhs[a]).iszero():
                    bad += 1
        check("(G) for P_can: [z^e] P_can Y^C = [z^e] Y^C_Q (all e tested, all slots, mod u^%d)" % NC, bad == 0,
              "%d slots" % tested)
        cc = Fr(1, 7)
        Pd = times_c(Pc, {0: Fr(1), 5: cc})
        lhs5 = apply_gauge(Pd, YC, 5, NC)      # its u^0 coefficient is complete for K >= 5 (Y^C(0) = -e_1)
        check("control: P_can (1 + z^5/7) Y^C is not Y^C_Q: [z^5 u^0] slot 1 = -1/7 (Y^C_Q has no z^5)",
              lhs5[1].coef(0) == Fr(-1, 7) and all(lhs5[a].coef(0) == 0 for a in (0, 2, 3, 4)))
        bad5 = 0
        for e in range(-((5 * NC) // 4) - 1, emax + 1):
            l = apply_gauge(Pd, YC, e, NC)
            r1 = YQ.get(-e, [Ser(NC, [], NC)] * 5) if e <= 0 else [Ser(NC, [], NC)] * 5
            r2 = YQ.get(5 - e, [Ser(NC, [], NC)] * 5) if e - 5 <= 0 else [Ser(NC, [], NC)] * 5
            for a in range(5):
                if not (l[a] - r1[a] - r2[a].scale(cc)).iszero():
                    bad5 += 1
        check("(G) for P_can c(z): P Y^C = c(z) Y^C_Q, c = 1 + z^5/7 (z-powers <= %d)" % emax, bad5 == 0)
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
