# -*- coding: utf-8 -*-
"""ref_g1_loc.py -- referee evidence (not proof) for gap cluster G1, items 1 and 2.

Independent of the prover's code.  Exact arithmetic (fractions.Fraction; Q(zeta_5) as Q[x]/(x^4+x^3+x^2+x+1)).
Python stdlib only.  Run:  py -I -B ref_g1_loc.py      (writes ref_g1_loc.log next to itself; prints relative names only)

Conventions (those of GAPS_G1.md s.1.1): H|_(p_i) = lambda_i, T_(p_i)P^4 has weights lambda_i - lambda_k, O(5)|_(p_i) has weight 5 lambda_i
(lambda' = 0).  Degree-1 genus-0 moduli: Mbar_(0,1)(P^4,1) = P(S) over G(2,5), Mbar_(0,2)(P^4,1) = P(S) x_G P(S) (S tautological);
torus-fixed points (line L_ij, marked points at p_i or p_j).

Check A (theorem (M) at order q, all orders in 1/z).  Localisation gives
    J_1|_(p_i) = Sum_(j != i) prod_(m=1..5)((5-m) lam_i + m lam_j) / ( prod_(k not in {i,j})(lam_j - lam_k) (z + lam_i - lam_j) ),
the q^1 coefficient of e^(-TH/z) J(TH,z) (J = e^(TH/z)(z + Sum_d e^(dT) ev_*(e(E'_d)/(z - psi)))).  Givental Thm 11.6 (CY case,
r = 1, l = 5, lambda' = 0), rewritten for s_i, predicts at order q:
    J_1|_(p_i) = z A_1 - 120 z - 770 lam_i - 120 e_1(lam),   A_1 = prod_(m=1..5)(5 lam_i + m z)/prod_j (lam_i - lam_j + z)
(120 = [q]C_0, 770 = [q] I~_1, 120 = [q] g_1).  Tested at random rational weights and at lam_i = zeta^i, as rational functions of z
(exact values at 25 rational points z).
Controls that MUST fail: (A1) Cor. 11.7 read literally ("Sum_alpha l_alpha" = l_1 = 5 in the exponent): extra term + 5*120 = +600;
(A2) the e_1-term omitted (generic weights only; at the zeta point e_1 = 0, so it must PASS there); (A3) E_d in place of E'_d
(no division by e(O(5)) in the eta-dual push-forward).

Check B (prop:set-cyclicLR at order q).  <H^a, H^b>_(0,2,1) = Sum_(i<j) e(Sym^5 S^vee) h_(a-1)(lam_i,lam_j) h_(b-1)(lam_i,lam_j)/e(T G(2,5))
at lam_i = zeta^i, compared with the q^1 coefficients of 5 kappa_1 and 5 kappa_3 (<<H,H,H>> = 5 kappa_1, <<H,H^3,H^4>> = 5 kappa_3,
divisor equation, e^T = q + O(q^2)) from the closed forms kappa_1 = L^5/(C_0^2 C_1^3), kappa_3 = C_0/C_1.  Z/5 grading: <H^2,H^3>_(0,2,1)
and <H^4,H^4>_(0,2,1) vanish at zeta (and not at generic weights: control).  Control (B1): the guess kappa_3 = C_1/C_0 must fail.
"""
import os, random
from fractions import Fraction as Fr
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
LOGNAME = "ref_g1_loc.log"


class Z5:
    """element of Q(zeta), zeta = e^(2 pi i/5); stored mod x^5 - 1 as 5 coefficients, compared mod Phi_5."""
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = tuple(Fr(x) for x in c)

    @staticmethod
    def const(a):
        return Z5((a, 0, 0, 0, 0))

    @staticmethod
    def zpow(k):
        c = [0] * 5
        c[k % 5] = 1
        return Z5(c)

    def __add__(self, o):
        o = o if isinstance(o, Z5) else Z5.const(o)
        return Z5([a + b for a, b in zip(self.c, o.c)])
    __radd__ = __add__

    def __neg__(self):
        return Z5([-a for a in self.c])

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Z5) else Z5.const(o)))

    def __rsub__(self, o):
        return Z5.const(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Z5) else Z5.const(o)
        r = [Fr(0)] * 5
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    if b:
                        r[(i + j) % 5] += a * b
        return Z5(r)
    __rmul__ = __mul__

    def red(self):
        c = list(self.c)
        return tuple(x - c[4] for x in c[:4])

    def gal(self, k):
        r = [Fr(0)] * 5
        for m, a in enumerate(self.c):
            r[(k * m) % 5] += a
        return Z5(r)

    def inv(self):
        num = self.gal(2) * self.gal(3) * self.gal(4)
        nrm = (self * num).red()
        assert all(x == 0 for x in nrm[1:]) and nrm[0] != 0, "not invertible"
        return num * (1 / nrm[0])

    def __truediv__(self, o):
        o = o if isinstance(o, Z5) else Z5.const(o)
        return self * o.inv()

    def __eq__(self, o):
        o = o if isinstance(o, Z5) else Z5.const(o)
        return self.red() == o.red()

    def rational(self):
        r = self.red()
        return r[0] if all(x == 0 for x in r[1:]) else None

    def __pow__(self, n):
        r = Z5.const(1)
        for _ in range(n):
            r = r * self
        return r


def prod(it, one):
    p = one
    for x in it:
        p = p * x
    return p


def J1_loc(lam, i, z, full_E=False):
    """localisation value of J_1|_(p_i) at the number z (full_E: control A3, E_d instead of E'_d)."""
    one = lam[0] * 0 + 1
    tot = lam[0] * 0
    for j in range(5):
        if j == i:
            continue
        num = prod((((5 - m) * lam[i] + m * lam[j]) for m in range(1, 6)), one)
        if full_E:
            num = num * (5 * lam[i])
        den = prod(((lam[j] - lam[k]) for k in range(5) if k not in (i, j)), one) * (z + lam[i] - lam[j])
        tot = tot + num / den
    return tot


def J1_mirror(lam, i, z, extra=0, with_e1=True):
    one = lam[0] * 0 + 1
    e1 = sum(lam[1:], lam[0])
    A1z = prod(((5 * lam[i] + m * z) for m in range(1, 6)), one) / prod(((lam[i] - lam[j] + z) for j in range(5) if j != i), one)
    r = A1z - 120 * z - 770 * lam[i] + extra
    if with_e1:
        r = r - 120 * e1
    return r


def hcomp(k, a, b):
    """complete homogeneous symmetric polynomial h_k(a, b); h_(-1) = 0."""
    if k < 0:
        return a * 0
    return sum((a ** m * b ** (k - m) for m in range(k + 1)), a * 0)


def two_point(lam, a, b):
    one = lam[0] * 0 + 1
    tot = lam[0] * 0
    for i in range(5):
        for j in range(i + 1, 5):
            eE = prod((((5 - m) * lam[i] + m * lam[j]) for m in range(6)), one)
            eT = prod(((lam[l] - lam[k]) for l in (i, j) for k in range(5) if k not in (i, j)), one)
            tot = tot + eE * hcomp(a - 1, lam[i], lam[j]) * hcomp(b - 1, lam[i], lam[j]) / eT
    return tot


def qseries(N):
    """C_0, I~_1, C_1, L^5 to O(q^(N+1)) from the closed formulas (5d)!/(d!)^5 and harmonic numbers (not from Fhat)."""
    from math import factorial
    a = [Fr(factorial(5 * d), factorial(d) ** 5) for d in range(N + 1)]
    It = [5 * a[d] * sum(Fr(1, m) for m in range(d + 1, 5 * d + 1)) for d in range(N + 1)]
    def mul(x, y):
        return [sum(x[k] * y[n - k] for k in range(n + 1)) for n in range(N + 1)]
    def inv(x):
        r = [Fr(0)] * (N + 1)
        r[0] = 1 / x[0]
        for n in range(1, N + 1):
            r[n] = -sum(x[k] * r[n - k] for k in range(1, n + 1)) / x[0]
        return r
    T1 = mul(It, inv(a))                                  # T - log q
    C1 = [Fr(1)] + [n * T1[n] for n in range(1, N + 1)]
    L5 = [Fr(3125) ** n for n in range(N + 1)]
    k1 = mul(L5, inv(mul(mul(a, a), mul(C1, mul(C1, C1)))))
    k3 = mul(a, inv(C1))
    eT = [Fr(0)] + [Fr(0)] * N                            # e^T = q exp(T1): coefficient of q^1 and q^2 suffice
    eT[1] = Fr(1)
    if N >= 2:
        eT[2] = T1[1]
    return a, It, C1, k1, k3, eT, mul, inv


def main():
    lines = []
    def log(s):
        print(s); lines.append(s)
    bad = [0]
    def check(name, ok, must_fail=False):
        good = (not ok) if must_fail else ok
        if not good:
            bad[0] += 1
        tag = ("PASS (control fails as required)" if good else "FAIL (control did not fail)") if must_fail else ("PASS" if ok else "FAIL")
        log("%s  %s" % (tag, name))
    log("ref_g1_loc: degree-1 localisation checks of (M) and of prop:set-cyclicLR (referee, independent code)")
    rng = random.Random(20261115)
    zs = [Fr(rng.randint(-400, 400), rng.randint(1, 37)) + Fr(1, 1000003) for _ in range(25)]
    # A: generic rational weights (three draws) and the zeta point
    for trial in range(3):
        lam = [Fr(rng.randint(-50, 50), rng.randint(1, 9)) for _ in range(5)]
        while len(set(lam)) < 5 or sum(lam) == 0:
            lam = [Fr(rng.randint(-50, 50), rng.randint(1, 9)) for _ in range(5)]
        okA = all(J1_loc(lam, i, z) == J1_mirror(lam, i, z) for i in range(5) for z in zs)
        check("A  generic weights #%d: J_1 (localisation) = z A_1 - 120 z - 770 lam_i - 120 e_1" % trial, okA)
        check("A1 generic #%d: Cor. 11.7 as printed (extra +600)" % trial,
              all(J1_loc(lam, i, z) == J1_mirror(lam, i, z, extra=600) for i in range(5) for z in zs), must_fail=True)
        check("A2 generic #%d: e_1 term omitted" % trial,
              all(J1_loc(lam, i, z) == J1_mirror(lam, i, z, with_e1=False) for i in range(5) for z in zs), must_fail=True)
        check("A3 generic #%d: E_d instead of E'_d" % trial,
              all(J1_loc(lam, i, z, full_E=True) == J1_mirror(lam, i, z) for i in range(5) for z in zs), must_fail=True)
    lamz = [Z5.zpow(i) for i in range(5)]
    zsz = [Z5.const(z) for z in zs[:12]]
    check("A  zeta point: J_1 (localisation) = z A_1 - 120 z - 770 zeta^i",
          all(J1_loc(lamz, i, z) == J1_mirror(lamz, i, z) for i in range(5) for z in zsz))
    check("A1 zeta point: Cor. 11.7 as printed (extra +600)",
          all(J1_loc(lamz, i, z) == J1_mirror(lamz, i, z, extra=600) for i in range(5) for z in zsz), must_fail=True)
    check("A2' zeta point: e_1 term omitted (e_1 = 0 here, so this must PASS)",
          all(J1_loc(lamz, i, z) == J1_mirror(lamz, i, z, with_e1=False) for i in range(5) for z in zsz))
    # B: two-point invariants at the zeta point
    a, It, C1, k1, k3, eT, mul, inv = qseries(3)
    I11 = two_point(lamz, 1, 1).rational()
    I34 = two_point(lamz, 3, 4).rational()
    I23 = two_point(lamz, 2, 3).rational()
    I44 = two_point(lamz, 4, 4).rational()
    log("   <H,H>_(0,2,1) = %s, <H^3,H^4>_(0,2,1) = %s, <H^2,H^3>_(0,2,1) = %s, <H^4,H^4>_(0,2,1) = %s  (zeta point)" % (I11, I34, I23, I44))
    log("   [q]C_0 = %s, [q]I~_1 = %s, [q]C_1 = %s, [q]kappa_1 = %s, [q]kappa_3 = %s" % (a[1], It[1], C1[1], k1[1], k3[1]))
    check("B  <H,H>_(0,2,1) = 2875 = 5 [q] kappa_1", I11 == 2875 and I11 == 5 * k1[1])
    check("B  <H^3,H^4>_(0,2,1) = 5 [q] kappa_3  (kappa_3 = C_0/C_1)", I34 == 5 * k3[1])
    check("B1 <H^3,H^4>_(0,2,1) = 5 [q] (C_1/C_0)  (wrong guess)", I34 == 5 * mul(C1, inv(a))[1], must_fail=True)
    check("B  Z/5 grading: <H^2,H^3>_(0,2,1) = <H^4,H^4>_(0,2,1) = 0 at zeta", I23 == 0 and I44 == 0)
    lamg = [Fr(x) for x in (0, 1, 3, 7, 12)]
    check("B2 grading control: <H^2,H^3>_(0,2,1) = 0 at generic weights", two_point(lamg, 2, 3) == 0, must_fail=True)
    log("TOTAL: %d unexpected result(s)" % bad[0])
    with open(os.path.join(HERE, LOGNAME), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("log written to " + LOGNAME)


if __name__ == "__main__":
    main()
