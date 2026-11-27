# -*- coding: utf-8 -*-
"""xy_legs.py -- pairings of the TR differentials of x = z^4 + v/z, y = z (v = 4) with leg functions
   <omega_{g,n}, Phi_1 .. Phi_n> = Sum_{a_1..a_n} Res ... Res prod_j Phi_j(z_j) omega_{g,n}(z_1..z_n)       (direct TR, xy_tr),
   Phi = A(z) + log(z) f(x(z)),  A a Laurent polynomial, f a polynomial (log a_i dropped: Res_a f(x) omega = 0).
Homogeneity: z -> c z, v -> c^5 v; a leg of degree d gives <..> = coeff * v^((5(2-2g-n) + Sum d_j)/5).
Class L (claimed regular): Phi in Q[v][z] + Q[v][x] (single valued, no log, regular at z = 0 modulo functions of x).
Usage: py -I -B xy_legs.py [fp89|fp61|qz]
"""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb
from xy_field import QZ, make_fp, P89, P61
from xy_tr import TR
from xy_jets import binom_gen

V = 4
LEGS = [  # name, A {exp: coeff}, f [coeffs in x], degree, class
    ("unit", {1: Fr(-1)}, [], 1, "L"), ("H4", {2: Fr(-1, 2)}, [], 2, "L"), ("H3", {3: Fr(-1, 3)}, [], 3, "L"),
    ("H2=v/z", {-1: Fr(V)}, [], 4, "L"), ("eul=-4z^5/5", {5: Fr(-4, 5)}, [], 5, "L"), ("z^9", {9: Fr(1)}, [], 9, "L"),
    ("z^14", {14: Fr(1)}, [], 14, "L"), ("x^2+3z^8", {8: Fr(4), 3: Fr(2 * V), -2: Fr(V * V)}, [], 8, "L"),
    ("z^10", {10: Fr(1)}, [], 10, "L"), ("z^15", {15: Fr(1)}, [], 15, "L"), ("z^19", {19: Fr(1)}, [], 19, "L"),
    ("z^20", {20: Fr(1)}, [], 20, "L"),
    ("H=-log z", {}, [Fr(-1)], 0, "bad: log"), ("dil=-4z^5/5+v log z", {5: Fr(-4, 5)}, [Fr(V)], 5, "bad: log"),
    ("tau1H=x log z-z^4/4+v/z", {4: Fr(-1, 4), -1: Fr(V)}, [Fr(0), Fr(1)], 4, "bad: log"),
    ("z^-5", {-5: Fr(1)}, [], -5, "bad: 1/v"), ("z^-1", {-1: Fr(1)}, [], -1, "bad: 1/v"),
]


def leg_series(F, A, f, i, K):
    """[t^k] Phi(a_i + t), k = 0..K, a_i = zeta^i (log a_i dropped)"""
    z = lambda e: F.zeta((i * e) % 5)
    out = [F.zero() for _ in range(K + 1)]
    for m, c in A.items():
        for k in range(K + 1): out[k] = out[k] + z(m - k).scale(c * binom_gen(m, k))
    if f:
        xs = [F.zero() for _ in range(K + 1)]          # x(a+t) = (a+t)^4 + V/(a+t)
        for k in range(K + 1):
            xs[k] = z(4 - k).scale(binom_gen(4, k)) + z(-1 - k).scale(Fr(V) * binom_gen(-1, k))
        fx = [F.zero() for _ in range(K + 1)]; pw = [F.one()] + [F.zero()] * K
        for c in f:
            for k in range(K + 1): fx[k] = fx[k] + pw[k].scale(c)
            pw = [sum((pw[j] * xs[k - j] for j in range(k + 1)), F.zero()) for k in range(K + 1)]
        lg = [F.zero()] + [z(-k).scale(Fr((-1) ** (k + 1), k)) for k in range(1, K + 1)]
        for k in range(K + 1):
            out[k] = out[k] + sum((lg[j] * fx[k - j] for j in range(k + 1)), F.zero())
    return out


def pairing(tr, g, n, legs):
    F = tr.F; K = tr.kmax + 2
    ser = [{i: leg_series(F, A, f, i, K) for i in range(5)} for (_, A, f, _, _) in legs]
    tot = F.zero()
    for key, cf in tr.om[(g, n)].items():
        t = cf
        for j, (i, k) in enumerate(key): t = t * ser[j][i][k]
        tot = tot + t
    return tot


def show(F, val, e):
    """value = coeff * V^e; return coeff as a rational (reconstructed for FP)"""
    c = val.scale(Fr(1) / Fr(V) ** e)
    if F is QZ: return c.rational()
    return c.rational()


def main():
    fld = sys.argv[1] if len(sys.argv) > 1 else "fp89"
    F = QZ if fld == "qz" else make_fp(P89 if fld == "fp89" else P61)
    tr = TR(F, 2, 3)
    for (g, n) in [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2)]: tr.step(g, n)
    print("field", F.name, " v =", V)
    for (g, n) in [(1, 1), (2, 1)]:
        for leg in LEGS:
            d = leg[3]; e = Fr(5 * (2 - 2 * g - n) + d, 5)
            val = pairing(tr, g, n, [leg])
            if val.is_zero():
                if e.denominator == 1 and e < 0 and leg[4] == "L":
                    print("  <omega_{%d,%d}, %-26s> = 0  (homogeneity allows v^(%s): nontrivial zero)" % (g, n, leg[0], e))
                continue
            c = show(F, val, e)
            reg = (e >= 0)
            print("  <omega_{%d,%d}, %-26s> = %s * v^(%s)   [%s]  %s" % (g, n, leg[0], c, e, leg[4],
                  "regular" if reg else "SINGULAR"))
    nz = {}
    for (g, n) in [(1, 2), (2, 2)]:
        for a, b in itertools.combinations_with_replacement(range(len(LEGS)), 2):
            la, lb = LEGS[a], LEGS[b]; e = Fr(5 * (2 - 2 * g - n) + la[3] + lb[3], 5)
            if e.denominator != 1: continue
            val = pairing(tr, g, n, [la, lb])
            cls = "L,L" if la[4] == "L" and lb[4] == "L" else "mixed/bad"
            if val.is_zero():
                if e < 0 and cls == "L,L": nz[(g, n)] = nz.get((g, n), 0) + 1
                continue
            c = show(F, val, e)
            print("  <omega_{%d,%d}, %s | %s> = %s * v^(%s)   [%s]  %s" % (g, n, la[0], lb[0], c, e, cls,
                  "regular" if e >= 0 else "SINGULAR"))
    print("  (L,L) pairs with integer v-exponent < 0 that vanish (nontrivial zeros): %s" % nz)
    # vanishing periods: Int_inf^0 omega_{g,1} and the double period of omega_{g,2}
    for g in (1, 2):
        per = F.zero()
        for key, cf in tr.om[(g, 1)].items():
            (i, k), = key
            per = per + cf * F.zeta((-i * k) % 5).scale(Fr(-(-1) ** k, k))
        print("  Int_inf^0 omega_{%d,1} = %s * v^(%d)" % (g, show(F, per, 1 - 2 * g), 1 - 2 * g))


if __name__ == "__main__":
    main()
