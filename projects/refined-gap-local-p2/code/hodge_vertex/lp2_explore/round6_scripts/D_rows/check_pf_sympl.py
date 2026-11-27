"""Theorem A, Prop 2.10: exact checks with controls (evidence only).
1. PF stripped form on a = vS: symbolic in v, multiply by (Ev-3)^3 to clear y. Through z^N.
2. Pi-tilde_ij at rational t in Q(omega), through z^N (all i,j), value -E delta_ij.
3. Full rows of Def 1.3 with arbitrary (C1i, XX) as free rationals: pairing = delta_ij (frame B).
4. Controls: bad r_l; bad denominators; different c1 at two colours; wrong eta; d2 not of the form E(C1i^2-1)/(2 Delta_c).
"""
import sys, time
from fractions import Fraction as Fr
from rows_core import *

def pf_check(E, N, r=None, sign=1):
    S = series_S(E, N, r, sign)
    a = [lp_mul(x, {1: Fr(1)}) for x in S]
    T1 = Th_plain(a, E, N)
    T2 = Th_plain(T1, E, N)
    T3 = Th_plain(T2, E, N)
    # Delta_c Th^3 + (E+27 y z) Th^2 + (E^2/3 + 6 y z^2) Th + E^3/27
    # y = 1/(Ev-3)^3 ; multiply through by C=(Ev-3)^3
    C = lp_mul(lp_mul({1: E, 0: Fr(-3)}, {1: E, 0: Fr(-3)}), {1: E, 0: Fr(-3)})
    Dc = lp_add(C, {0: Fr(27)})  # C*Delta_c = C+27
    res = []
    for n in range(N + 1):
        tot = lp_mul(Dc, T3[n])
        tot = lp_add(tot, lp_mul(C, lp_scale(T2[n], E)))
        if n >= 1:
            tot = lp_add(tot, lp_scale(T2[n - 1], Fr(27)))
        tot = lp_add(tot, lp_mul(C, lp_scale(T1[n], E * E / 3)))
        if n >= 2:
            tot = lp_add(tot, lp_scale(T1[n - 2], Fr(6)))
        tot = lp_add(tot, lp_mul(C, lp_scale(a[n], E ** 3 / 27)))
        res.append(tot)
    bad = [n for n, x in enumerate(res) if x]
    return bad

def evalrow(row, v):
    return [lp_eval(x, v) for x in row]

def colour_data(t, i, E):
    ti = (OMEGA ** i) * t
    v = (1 + 3 * ti) / (E * ti)
    Y = (-E / 3) + 1 / v
    return ti, v, Y

def pair_tilde(E, t, N, i, j, r=None, sign=1):
    """Pi-tilde_ij through z^N at rational t."""
    y = t ** 3
    Dc = 1 + 27 * y
    a, b, Thb = abc_vs(E, N, r, sign)
    a2, b2, Thb2 = abc_flip(E, N, r, sign)
    _, vi, _ = colour_data(t, i, E)
    _, vj, _ = colour_data(t, j, E)
    au = evalrow(a, vi); bu = evalrow(b, vi); cu = [Dc * x for x in evalrow(Thb, vi)]
    aw = evalrow(a2, vj); bw = evalrow(b2, vj); cw = [Dc * x for x in evalrow(Thb2, vj)]
    out = []
    for n in range(N + 1):
        s = Qw(0)
        for k in range(n + 1):
            s = s + (-(E ** 3) / 27) * au[k] * aw[n - k] + E * bu[k] * bw[n - k] + bu[k] * cw[n - k] + cu[k] * bw[n - k]
        out.append(s)
    return out

def full_rows(E, N, v, C1i, XX, Dc, r=None, flip=False):
    """Def 1.3 rows with sqrt(c0) stripped: P0=a, P1=C1i b, P2=(Dc/C1i^2)[Th P1 - d2 P1]
    Th P1 = C1i(Th b) - z C1i XX b  (at +z); with z -> -z handled by flip."""
    d2 = E * (C1i * C1i - 1) / (2 * Dc)
    if not flip:
        a, b, Thb = abc_vs(E, N, r)
        sg = 1
    else:
        a, b, Thb = abc_flip(E, N, r)
        sg = -1
    av = evalrow(a, v); bv = evalrow(b, v); Tv = evalrow(Thb, v)
    P0 = av
    P1 = [C1i * x for x in bv]
    P2 = []
    for n in range(N + 1):
        s = C1i * Tv[n] - d2 * C1i * bv[n]
        if n >= 1:
            s = s - sg * XX * C1i * bv[n - 1]
        P2.append(s * Dc / (C1i * C1i))
    return P0, P1, P2

def pair_full(E, e1, e2, t, N, i, j, ci, xi, cj=None, xj=None, eta_err=False, r=None):
    if cj is None:
        cj, xj = ci, xi
    y = t ** 3; Dc = 1 + 27 * y
    _, vi, _ = colour_data(t, i, E)
    _, vj, _ = colour_data(t, j, E)
    A = full_rows(E, N, vi, ci, xi, Dc, r)
    B = full_rows(E, N, vj, cj, xj, Dc, r, flip=True)
    c0 = 3 * e1 * e2 * E
    n00 = e1 * e2 * E ** 3 / 9
    n11 = -3 * e1 * e2 * E * (2 if eta_err else 1)
    n12 = -3 * e1 * e2
    out = []
    for n in range(N + 1):
        s = Qw(0)
        for k in range(n + 1):
            s = s + n00 * A[0][k] * B[0][n - k] + n11 * A[1][k] * B[1][n - k] + n12 * (A[1][k] * B[2][n - k] + A[2][k] * B[1][n - k])
        out.append(s * (1 / c0))
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    weights = [(Fr(2, 3), Fr(5, 7)), (Fr(-11, 5), Fr(3, 2)), (Fr(1, 1), Fr(1, 13))]
    t0 = time.time()
    for e1, e2 in weights:
        E = e1 + e2
        print("weights", e1, e2, "E", E)
        print(" PF stripped form nonzero z-orders (want []):", pf_check(E, N))
    E = Fr(2, 3) + Fr(5, 7)
    print(" control PF, bad r_3:", pf_check(E, 12, rec_r(14, (3, Fr(1, 1000))))[:6])
    print(" control PF, bad denominators:", pf_check(E, 12, None, -1)[:6])
    print("time", time.time() - t0)

    ts = [Fr(-1, 7), Fr(3, 11), Fr(-5, 2)]
    for e1, e2 in weights[:2]:
        E = e1 + e2
        for t in ts:
            for i in range(3):
                for j in range(i, 3):
                    out = pair_tilde(E, t, N, i, j)
                    nz = [n for n, x in enumerate(out) if not Qw._c(x).iszero()]
                    if i == j:
                        ok = nz == [0] and out[0] == -E
                    else:
                        ok = nz == []
                    print(f" Pitilde t={t} ij={i}{j} ok={ok}", "" if ok else (nz[:6], out[0]))
    print("time", time.time() - t0)
    # generic c1 data
    e1, e2 = weights[0]; E = e1 + e2
    t = Fr(-1, 7)
    for (ci, xi) in [(Fr(3, 5), Fr(7, 9)), (Fr(-1, 1), Fr(-2, 9)), (Fr(101, 17), Fr(0))]:
        for i in range(3):
            for j in range(3):
                out = pair_full(E, e1, e2, t, min(N, 30), i, j, ci, xi)
                nz = [n for n, x in enumerate(out) if not Qw._c(x).iszero()]
                ok = (nz == [0] and out[0] == 1) if i == j else nz == []
                print(f" full rows c1i={ci} XX={xi} ij={i}{j} ok={ok}", "" if ok else (nz[:6], out[0]))
    print(" controls:")
    o = pair_full(E, e1, e2, t, 6, 0, 0, Fr(3, 5), Fr(7, 9), Fr(4, 5), Fr(7, 9))
    print("  different C1i at the two colours: first nonzero orders", [n for n, x in enumerate(o) if not Qw._c(x).iszero()][:5], "z^0 =", o[0])
    o = pair_full(E, e1, e2, t, 6, 0, 0, Fr(3, 5), Fr(7, 9), Fr(3, 5), Fr(1, 9))
    print("  different XX at the two colours: orders", [n for n, x in enumerate(o) if not Qw._c(x).iszero()][:5], "z^0 =", o[0])
    o = pair_full(E, e1, e2, t, 6, 0, 1, Fr(3, 5), Fr(7, 9), eta_err=True)
    print("  wrong eta^11, i!=j: orders", [n for n, x in enumerate(o) if not Qw._c(x).iszero()][:5])
    o = pair_tilde(E, t, 8, 0, 0, rec_r(12, (2, Fr(1, 1000))))
    print("  bad r_2, Pitilde_00: orders", [n for n, x in enumerate(o) if not Qw._c(x).iszero()][:6])
    o = pair_tilde(E, t, 8, 0, 1, rec_r(12, (2, Fr(1, 1000))))
    print("  bad r_2, Pitilde_01: orders", [n for n, x in enumerate(o) if not Qw._c(x).iszero()][:6])
    print("time", time.time() - t0)
