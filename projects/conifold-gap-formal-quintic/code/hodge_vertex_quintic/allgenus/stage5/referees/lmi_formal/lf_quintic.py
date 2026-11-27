"""lf_quintic N N0: the two-variable chart of REFEREE s.3.3 for the quintic, exact over Q.

Q0  power sums: P_r = p_r + (1 - x)^-r; P_r(0) = 5 [5 | r]; p_r(x) against explicit Q(zeta_5) sums (Cyc5).
Q1  h = h_(w_0) from its DEFINITION in Q(zeta_5): exp(-Sum_k b_k z^(2k-1) Sum_(a in w_0) a^(1-2k)).
Q2  chart at s* = 1 (alpha0 = -1, M0 = 0) = one-variable c(z; L) at six rational L, through z^N;
    total Hessian determinant = -g''(L) a0(L) - alpha0^2 and the Schur relation det = -a0(L) S_0''(L).
Q3  chart at L = 1 (Q = 0) = h through z^N0.
Q4  chart at s* = 0 (alpha0 = 5, M0 = 1) = c(z; L) at four rational L (z^N); chart at L = 0 = Ginf (z^N0).
Failing controls: Q e^(2 mu) for Q e^mu; the singular Stirling series dropped in c(z; L); twist weight +5 in h;
M0 = 0 at the orbifold point; Bernoulli polynomials at 0, 4/5, -1/5 instead of 1/5 in Ginf.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from lf_core import (bk, bern_poly, ser_exp, Cyc5, zeta_pow, power_sums_all, power_sums_reg)
from lf_wick import wick_1, wick_2sep, cov_2
from lf_models import q_onevar, q_chart1, q_chart0
from lf_log import Log, first_diff

Z, O = Fr(0), Fr(1)


def h_definition(N, twist=-5):
    """h_(w_0), w_0 = {1 - zeta^j, j = 1..4} U {twist}, computed in Q(zeta_5)."""
    s = [Fr(0)] * (N + 1)
    inv = [(Cyc5([1]) + zeta_pow(j) * (-1)).inv() for j in range(1, 5)]     # (1 - zeta^j)^-1
    k = 1
    while 2 * k - 1 <= N:
        e = 2 * k - 1
        tot = Cyc5([Fr(twist) ** (1 - 2 * k)])
        for u in inv:
            pw = Cyc5([1])
            for _ in range(e):
                pw = pw * u
            tot = tot + pw
        s[e] -= bk(k) * tot.rational()
        k += 1
    return ser_exp(s, N + 1)


def g_shift(N, h):
    """exp(Sum_(5 | n) 5 (-1)^(n+1) B_(n+1)(h) z^n/(n(n+1))); h = 1/5 is Ginf."""
    s = [Fr(0)] * (N + 1)
    for n in range(5, N + 1, 5):
        s[n] = 5 * Fr((-1) ** (n + 1)) * bern_poly(n + 1, h) / (n * (n + 1))
    return ser_exp(s, N + 1)


def ginf(N):
    s = [Fr(0)] * (N + 1)
    for l in range(1, N // 5 + 1):
        s[5 * l] = Fr((-1) ** (l + 1)) * bern_poly(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1))
    return ser_exp(s, N + 1)


def chart(data, N):
    phs, las, phm, lam, H = data
    Css, Csm, Cmm, det = cov_2(*H)
    return wick_2sep(phs, las, phm, lam, (Css, Csm, Cmm), N, Z, O), det


def main():
    N, N0 = int(sys.argv[1]), int(sys.argv[2])
    log = Log("lf_quintic_%d_%d" % (N, N0))
    # Q0
    ok = True
    for x in [Fr(2), Fr(-2, 3), Fr(1, 2), Fr(3)]:
        P, p = power_sums_all(x, 40), power_sums_reg(x, 40)
        ok &= all(P[r] == p[r] + (1 - x) ** (-r) for r in range(1, 41))
    log.check("Q0 P_r(x) = p_r(x) + (1 - x)^-r, x = 2, -2/3, 1/2, 3, r <= 40", ok)
    P0 = power_sums_all(Fr(0), 40)
    log.check("Q0 P_r(0) = 5 if 5 | r else 0, r <= 40", all(P0[r] == (5 if r % 5 == 0 else 0) for r in range(1, 41)))
    ok = True
    for x in [Fr(1), Fr(2)]:
        p = power_sums_reg(x, 25)
        invs = [(zeta_pow(j) + Cyc5([-x])).inv() for j in range(1, 5)]
        for r in range(1, 26):
            tot = Cyc5([0])
            for u in invs:
                pw = Cyc5([1])
                for _ in range(r):
                    pw = pw * u
                tot = tot + pw
            ok &= (tot.rational() == p[r])
    log.check("Q0 p_r(x) = explicit Sum_(j=1..4) (zeta^j - x)^-r in Q(zeta_5), x = 1, 2, r <= 25", ok)
    # Q1
    h = h_definition(max(N0, N))
    log.say("Q1 h_1..h_5 = %s" % ", ".join(str(x) for x in h[1:6]))
    log.check("Q1 h_1 = -3/20, h_2 = 9/800 (values quoted in CHART_PROOF s.3.6 and REFEREE_2 s.4.2)",
              h[1] == Fr(-3, 20) and h[2] == Fr(9, 800))
    # Q2
    for L in [Fr(2), Fr(-2, 3), Fr(1, 2), Fr(3), Fr(-1), Fr(5, 4)]:
        ph3, logA, C = q_onevar(L, N)
        c = wick_1(ph3, logA, C, N, Z, O)
        data = q_chart1(L, N)
        ch, det = chart(data, N)
        gpp, a0 = data[4][0], 1 - L
        log.check("Q2 chart(s*=1) = c(z; L) through z^%d at L = %s" % (N, L), ch == c,
                  "c_1 = %s" % c[1] if L == 2 else "")
        if L == 2 and N >= 5:
            stored = [Fr(-12, 5), Fr(72, 25), Fr(354278, 625), Fr(-73040834, 3125), Fr(11495493197, 31500)]
            log.check("Q2 c_1..c_5(2) = the values quoted in REFEREE_2 s.4.2 (E1) and CHART_PROOF s.1.7", c[1:6] == stored)
        S0pp = 5 / L + power_sums_all(L, 2)[1]
        log.check("Q2 det H = -g''(L) a0(L) - 1 = -a0(L) S_0''(L) at L = %s" % L,
                  det == -gpp * a0 - 1 == -a0 * S0pp)
        if L in (Fr(2), Fr(-2, 3)):
            ch2, _ = chart(q_chart1(L, N, mu_scale=2), N)
            log.control("Q2 chart with Q e^(2 mu) at L = %s" % L, ch2 == c, "first diff z^%s" % first_diff(ch2, c))
            ph3x, logAx, Cx = q_onevar(L, N, sing_stirling=False)
            cx = wick_1(ph3x, logAx, Cx, N, Z, O)
            log.control("Q2 c(z; L) without the singular Stirling series vs chart, L = %s" % L, cx == ch,
                        "first diff z^%s" % first_diff(cx, ch))
    # Q3
    ch1, det1 = chart(q_chart1(Fr(1), N0), N0)
    log.check("Q3 det H at L = 1 equals -alpha0^2 = -1", det1 == -1)
    log.check("Q3 chart at L = 1 (Q = 0) = h_(w_0) from its definition, through z^%d" % N0, ch1 == h[:N0 + 1])
    hp = h_definition(N0, twist=5)
    log.control("Q3 twist weight +5 in h", ch1 == hp, "first diff z^%s" % first_diff(ch1, hp))
    # Q4
    for L in [Fr(1, 2), Fr(-1, 3), Fr(2), Fr(-2)]:
        ph3, logA, C = q_onevar(L, N)
        c = wick_1(ph3, logA, C, N, Z, O)
        ch, det = chart(q_chart0(L, N), N)
        log.check("Q4 chart(s*=0) = c(z; L) through z^%d at L = %s" % (N, L), ch == c)
    ch0, det0 = chart(q_chart0(Fr(0), N0), N0)
    G = ginf(N0)
    log.check("Q4 det H at L = 0 equals -alpha0^2 = -25", det0 == -25)
    log.check("Q4 chart at L = 0 = Ginf through z^%d" % N0, ch0 == G)
    log.check("Q4 Ginf = g_shift(1/5) (two codings of the closed form)", G == g_shift(N0, Fr(1, 5)))
    log.say("Q4 [z^5] Ginf = %s, [z^10] Ginf = %s" % (G[5], G[10] if N0 >= 10 else "-"))
    log.check("Q4 [z^5] Ginf = B_6(1/5)/6 = 4537/3937500", G[5] == Fr(4537, 3937500))
    d = q_chart0(Fr(0), N0)
    chM, _ = chart((d[0], d[1], d[2], {}, d[4]), N0)
    log.control("Q4 orbifold chart with M0 = 0 (no e^mu) vs Ginf", chM == G, "first diff z^%s" % first_diff(chM, G))
    for hh in [Fr(0), Fr(4, 5), Fr(-1, 5)]:
        Gx = g_shift(N0, hh)
        log.control("Q4 Bernoulli polynomials at %s instead of 1/5" % hh, Gx == ch0,
                    "first diff z^%s" % first_diff(Gx, ch0))
    log.close()


if __name__ == "__main__":
    main()
