# -*- coding: utf-8 -*-
"""s15_rational.py -- display only: rational reconstruction (CRT of the residues mod 2^89-1 and 2^61-1, Wang bound sqrt(p1 p2 / 2))
of the genus-3 flat graph values, frame (0, 0), and where the stored-normalisation control (variant 'stored') differs.
Reconstructed rationals are NOT independent evidence (they are read off the same two residue tables); the evidence is s15_compare.
Run: py -I -B s15_rational.py      Log: logs/s15_rational.log (overwritten).
"""
import os, sys, json, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from s15_flat import SCRATCH
P1, P2 = (1 << 61) - 1, (1 << 89) - 1


def crt(a1, a2):
    m = P1 * P2; x = (a1 * P2 * pow(P2, -1, P1) + a2 * P1 * pow(P1, -1, P2)) % m
    return x, m


def ratrec(a, m):
    B = math.isqrt(m // 2); r0, r1, s0, s1 = m, a % m, 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    return Fr(r1, s1) if s1 and abs(s1) <= B else None


def load(fn): return {int(k): v for k, v in json.load(open(os.path.join(SCRATCH, fn)))['A'].items()}


def main():
    LOG = open(os.path.join(HERE, "logs", "s15_rational.log"), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n")
    A61, A89 = load("flat_g3_61_can.json"), load("flat_g3_89_can.json")
    log("genus-3 flat graph values, frame (0,0): coefficients of u^-4 .. u^3 (CRT of the 2^61-1 and 2^89-1 residues, rational reconstruction)")
    for i in sorted(A89):
        vals = []
        for k in range(8):
            x, m = crt(A61[i][k], A89[i][k]); r = ratrec(x, m); vals.append(str(r) if r is not None else "?")
        if any(v != "0" for v in vals): log("  graph %2d: %s" % (i, vals))
        else: log("  graph %2d: identically 0 through u^3" % i)
    S89, S61 = load("flat_g3_89_stored.json"), load("flat_g3_61_stored.json")
    log("stored-normalisation control: (graph, k) with flat(P_stored) != flat(P_can) [k = index of u^(k-4)], difference reconstructed")
    for i in sorted(A89):
        for k in range(8):
            if (S89[i][k] - A89[i][k]) % P2:
                x, m = crt((S61[i][k] - A61[i][k]) % P1, (S89[i][k] - A89[i][k]) % P2); log("  graph %2d, u^%d: difference %s" % (i, k - 4, ratrec(x, m)))
    LOG.close()


if __name__ == "__main__":
    main()
