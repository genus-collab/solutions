# -*- coding: utf-8 -*-
"""c8_thirdprime.py -- a THIRD prime (2^521-1, also = 1 mod 5, the one an earlier referee used) on a sample of
genus-5 inner graphs, point frame (variant id).  Two purposes:
  (1) arithmetic independence: the p89/p61 residues of the shipped id table must be the reductions of ONE rational
      that the three primes together reconstruct.  Two Mersenne primes run through the same code are not an
      independent check of the code; a 521-bit prime at least rules out modulus-specific accidents.
  (2) height: 2^89-1 and 2^61-1 give only 150 bits, and most genus-5 per-graph coefficients do not reconstruct
      from them.  With 671 bits the true rationals become visible, so the margin of the two-prime ZERO test
      (a non-zero rational can vanish mod both primes only if its numerator exceeds p89*p61 ~ 2^150) can be judged.
Run:  py -I -B c8_thirdprime.py [n]
"""
import json, os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
P89 = (1 << 89) - 1
P61 = (1 << 61) - 1
P521 = (1 << 521) - 1


def crt_list(pairs):
    M = 1
    x = 0
    for a, p in pairs:
        # x = x mod M, a mod p
        t = (a - x) * pow(M % p, p - 2, p) % p
        x = x + M * t
        M *= p
    return x % M, M


def ratrec(a, M, B):
    a %= M
    if a == 0:
        return Fr(0)
    r0, r1, s0, s1 = M, a, 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B:
        return None
    return Fr(r1, s1)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    from g5_engine import make_split, graph_id
    cost = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))
    a89 = json.load(open(os.path.join(PKG, "results", "g5_89_id.json")))["table"]
    a61 = json.load(open(os.path.join(PKG, "results", "g5_61_id.json")))["table"]
    order = sorted(range(4555), key=lambda t: cost[str(t)]["est_split"])
    sample = order[:n] + [1958]
    t0 = time.time()
    sp = make_split(P521, 5, ("0", "0"), 8, "id")
    print("Split at 2^521-1 built in %.0f s" % (time.time() - t0), flush=True)
    M2 = P89 * P61
    bad = []
    hs = []
    nrec = 0
    nz = 0
    t0 = time.time()
    for t in sample:
        row = graph_id(sp, t)
        for k in range(8):
            x, M = crt_list([(a61[str(t)][k], P61), (a89[str(t)][k], P89), (row[k], P521)])
            r = ratrec(x, M, 1 << 320)
            if r is None:
                bad.append((t, k, "no reconstruction even with 671 bits"))
                continue
            nrec += 1
            if r:
                nz += 1
            h = max(abs(r.numerator), r.denominator)
            hs.append(h)
            # consistency: the rational must reduce to each stored residue
            for a, p in ((a89[str(t)][k], P89), (a61[str(t)][k], P61), (row[k], P521)):
                if r.numerator % p * pow(r.denominator % p, p - 2, p) % p != a % p:
                    bad.append((t, k, "reduction mismatch"))
    print("sample of %d inner graphs, %d coefficients, %.0f s" % (len(sample), 8 * len(sample), time.time() - t0))
    print("  reconstructed from the three primes: %d of %d (%d non-zero)" % (nrec, 8 * len(sample), nz))
    print("  inconsistencies: %d %s" % (len(bad), bad[:5]))
    if hs:
        hs.sort()
        print("  height: median %.2e, max %.3e = 2^%.1f" % (hs[len(hs) // 2], hs[-1], hs[-1].bit_length()))
        print("  reconstructible from p89*p61 alone (height <= 2^74): %d of %d"
              % (sum(1 for h in hs if h.bit_length() <= 74), len(hs)))
        print("  a false zero at BOTH p89 and p61 needs a numerator >= %.2e = 2^%d; largest numerator in the sample"
              " is 2^%d" % (M2, M2.bit_length(), max(h.bit_length() for h in hs)))


main()
