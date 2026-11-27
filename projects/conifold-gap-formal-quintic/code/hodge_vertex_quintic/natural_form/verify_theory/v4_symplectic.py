# -*- coding: utf-8 -*-
"""v4_symplectic.py -- a property NOT used in UNIFORM_REDUCTION.md: the regular gauge is SYMPLECTIC,
        P(u, z) eta^-1 P(u, -z)^T = eta^-1 x s(z^10),   eta^{b b'} = 1/5 for b + b' = 3 mod 5,
with s = 1 for z-orders 1..9 (in any normalisation: the constant c_5 drops out of the odd order 5), and at order z^10 a constant
scalar defect delta that depends on the normalisation constants: for P -> P (1 + c_5 z^5 + c_10 z^10) the defect becomes
delta + 2 c_10 - c_5^2.  Consequence: for the gauge of the actual (symplectic) R-matrix the EVEN constants c_10, c_20, .. are fixed by
symplecticity; only c_5, c_15, .. are genuinely free.  With (R1) (P_11(0,z) = P_22(0,z) =: phi) the (1,2) entry at u = 0 reads
phi(z) phi(-z) = 1, so 'P_11(0,z) = 1 + O(z^10)' automatically improves to 'O(z^15)'.
VERIFIED here exactly over Q for k <= 10 (frame (0,0)) and k <= 8 (frame (1/3,-2/7)); not proved.
Run: py -I -B v4_symplectic.py
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from v2_gauge import gauge
OUT = open(os.path.join(HERE, "v4_symplectic.log"), "w", encoding="utf-8")


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def form(p, n):
    """entries of [z^n] P(z) eta^-1 P(-z)^T times 5: {(a, a'): series}"""
    acc = {}
    for k in range(n + 1):
        kp = n - k
        for b in range(5):
            bp = (3 - b) % 5
            key = ((b + k) % 5, (bp + kp) % 5)
            t = (p[k][b] * p[kp][bp]).scale(Fr((-1) ** kp))
            acc[key] = t if key not in acc else acc[key] + t
    return acc


def main():
    ok_all = True
    for name, al, be, K, N in (("frame (0,0)", 0, 0, 10, 32), ("frame (1/3,-2/7)", Fr(1, 3), Fr(-2, 7), 8, 24)):
        p = gauge(K, N, al, be)['p']
        f0 = form(p, 0)
        ok = all(v.coef(0) == 1 and all(v.coef(j) == 0 for j in range(1, v.prec)) for v in f0.values()) and len(f0) == 5
        log(("PASS " if ok else "FAIL ") + "%s: order z^0: P_0 eta^-1 P_0^T = eta^-1" % name); ok_all = ok_all and ok
        for n in range(1, min(K, 9) + 1):
            f = form(p, n); ok = all(not v.c for v in f.values())
            ok_all = ok_all and ok
            log(("PASS " if ok else "FAIL ") + "%s: order z^%d of P(z) eta^-1 P(-z)^T vanishes identically in u (precision >= O(u^%d))" % (name, n, min(v.prec for v in f.values())))
        if K >= 10:
            f = form(p, 10); c5 = -p[5][1].coef(0); p101 = p[10][1].coef(0)
            vals = {key: v for key, v in f.items() if v.c}
            const = all(all(v.coef(j) == 0 for j in range(1, v.prec)) for v in vals.values())
            dl = set(v.coef(0) for v in vals.values())
            ok = const and len(dl) == 1 and sorted(vals) == [(0, 3), (1, 2), (2, 1), (3, 0), (4, 4)] and dl == {2 * p101 - c5 * c5}
            log(("PASS " if ok else "FAIL ") + "%s: order z^10 (normalisation p_{5,0}(0) = p_{10,0}(0) = 0): defect = delta eta^-1, constant in u, delta = 2 p_{10,1}(0) - c_5^2 = %s"
                % (name, 2 * p101 - c5 * c5)); ok_all = ok_all and ok
            log("   => symplectic normalisation c_10 = (c_5^2 - delta)/2 = %s, for which P_11(0,z) = P_22(0,z) = 1 + O(z^15) if c_5 = %s" % ((c5 * c5 - (2 * p101 - c5 * c5)) / 2, c5))
    log("V4_SYMPLECTIC " + ("ALL PASS" if ok_all else "FAILURES")); OUT.close()


if __name__ == "__main__":
    main()
