# -*- coding: utf-8 -*-
"""rc_t2_tc.py -- part (2) of the re-check: the change of variable u -> t_c, re-derived independently,
the gap, the triangularity/perturbation analysis, and the closed form G = c_5 (u/t_c)^8.
"""
import os
import sys
import json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rc_pf
import rc_rec

RES = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..",
                                   "cluster_genus5_quintic", "results"))
G = 5
M = 2 * (G - 1)          # 8


def load(tag):
    with open(os.path.join(RES, tag + ".json")) as f:
        return json.load(f)


def fmt(x):
    x = Fr(x)
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def main():
    tot = [Fr(s) for s in json.load(open(os.path.join(HERE, "rc_total_rational.json")))['total']]
    ver = [Fr(s) for s in json.load(open(os.path.join(HERE, "rc_total_rational.json")))['vertex']]
    c5 = abs(rc_rec.bern(10)) / (10 * 8 * 5 ** 4)
    print("B_10 =", fmt(rc_rec.bern(10)), " c_5 = |B_10|/(10*8*5^4) =", fmt(c5), " 1/660000 =", fmt(Fr(1, 660000)))
    print("c_5 == 1/660000:", c5 == Fr(1, 660000))
    print()
    tcc, _ = rc_pf.tc_coeffs(0, 0, kmax=10)
    print("t_c = u + ... (own Frobenius solution), u^1..u^10:", [fmt(x) for x in tcc])
    print()

    for label, ser in (("total", tot), ("vertex graph row", ver)):
        T = rc_pf.polar_in_tc(ser, 0, 0, G)
        print("--- %s in t_c  (coefficients of t_c^-8 .. t_c^-1)" % label)
        for k, x in enumerate(T):
            print("    t_c^%-3d %s" % (k - M, fmt(x)))
        print("    all seven subleading coefficients vanish:", not any(T[1:]))
        print("    leading coefficient equals c_5:", T[0] == c5)
        print()

    print("--- closed form: is the u-series exactly c_5 (u/t_c)^8 mod u^8 ?")
    cf = rc_pf.u_over_tc_pow(M, 0, 0, kmax=M)
    pred = [c5 * x for x in cf]
    for k in range(M):
        print("    u^%d  data %-34s  c_5 (u/t_c)^8 %-34s  %s"
              % (k, fmt(tot[k]), fmt(pred[k]), "EQUAL" if tot[k] == pred[k] else "DIFFERENT"))
    print("    identical in all 8 coefficients:", pred == tot)
    print()

    print("--- Jacobian dT_k / da_j of the change of variable (a = u-coefficients, T = t_c-coefficients)")
    rows = []
    for j in range(M):
        e = [Fr(0)] * M
        e[j] = Fr(1)
        rows.append(rc_pf.polar_in_tc(e, 0, 0, G))
    print("    (row j = perturbation of a_j by +1; column k = t_c^(k-8))")
    hdr = "     a_j \\ T_k " + "".join("%14s" % ("t_c^%d" % (k - M)) for k in range(M))
    print(hdr)
    for j in range(M):
        print("     a_%-9d" % j + "".join("%14s" % fmt(rows[j][k]) for k in range(M)))
    print("    lower triangular with 1 on the diagonal:",
          all(rows[j][k] == 0 for j in range(M) for k in range(j)) and all(rows[j][j] == 1 for j in range(M)))
    print()

    print("--- perturbation test: change ONE u-coefficient of the reconstructed total by +1 and see which t_c coefficients move")
    for j in range(M):
        s = list(tot)
        s[j] = s[j] + 1
        T = rc_pf.polar_in_tc(s, 0, 0, G)
        moved = [k - M for k in range(M) if T[k] != 0]
        print("    perturb a_%d: non-zero t_c exponents become %s" % (j, moved))
    print()

    print("--- the unique u-series with a_0 = c_5 and vanishing subleading polar part")
    sol = [c5] + [Fr(0)] * (M - 1)
    for k in range(1, M):
        T = rc_pf.polar_in_tc(sol, 0, 0, G)
        sol[k] = sol[k] - T[k]          # diagonal of the Jacobian is 1
    T = rc_pf.polar_in_tc(sol, 0, 0, G)
    print("    forced series:", [fmt(x) for x in sol])
    print("    check: its polar part is", [fmt(x) for x in T])
    print("    equals the reconstructed total:", sol == tot)
    print()

    print("--- second conifold frame (alpha,beta) = (1/3,-2/7), p89 only: the same test mod p89")
    d2 = load("g5_89_nat_ab1o3_m2o7")
    p = rc_rec.P89
    rows2 = []
    for j in range(M):
        e = [Fr(0)] * M
        e[j] = Fr(1)
        rows2.append(rc_pf.polar_in_tc(e, Fr(1, 3), Fr(-2, 7), G))
    for label, res in (("total", d2['info']['total']), ("vertex row", d2['table']['0'])):
        Tp = [sum(rc_rec.tofp(rows2[j][k], p) * res[j] for j in range(M)) % p for k in range(M)]
        print("    %s: t_c residues mod p89 =" % label, Tp[:1], "then", Tp[1:])
        print("       subleading all zero:", not any(Tp[1:]),
              "  leading == c_5 mod p89:", Tp[0] == rc_rec.tofp(c5, p))
    print()
    print("--- cross-check: same mod-p route applied to the main p89 run (frame (0,0))")
    d1 = load("g5_89_nat")
    Tp = [sum(rc_rec.tofp(rows[j][k], p) * d1['info']['total'][j] for j in range(M)) % p for k in range(M)]
    print("    subleading all zero:", not any(Tp[1:]), "  leading == c_5 mod p89:", Tp[0] == rc_rec.tofp(c5, p))
    print("    matches info['total_tc'] stored in the JSON:", Tp == d1['info']['total_tc'])
    print("    matches info['vertex'] stored in the JSON:", Tp == d1['info']['vertex'])
    d61 = load("g5_61_nat")
    rows61 = rows
    q = rc_rec.P61
    Tq = [sum(rc_rec.tofp(rows61[j][k], q) * d61['info']['total'][j] for j in range(M)) % q for k in range(M)]
    print("    p61 run: subleading all zero:", not any(Tq[1:]),
          "  leading == c_5 mod p61:", Tq[0] == rc_rec.tofp(c5, q),
          "  matches stored total_tc:", Tq == d61['info']['total_tc'])


if __name__ == "__main__":
    main()
