# -*- coding: utf-8 -*-
"""c1_json_audit.py -- audit of the copied genus-5 result JSONs ALONE (no project code, no engine).
Checks: table sums vs info['total']; nat total vs id total; CRT of the two primes; the leading coefficient
against 1/660000; the by-E accounting of what the regularity claim actually tested.
Run:  py -I -B c1_json_audit.py
"""
import json, os
from fractions import Fraction as Fr

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic", "results")
P89 = (1 << 89) - 1
P61 = (1 << 61) - 1


def load(tag):
    return json.load(open(os.path.join(D, "%s.json" % tag)))


def ofrac(x, p):
    return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def crt(a, b):
    """x mod P89*P61"""
    M = P89 * P61
    return (a * P61 * pow(P61, P89 - 2, P89) + b * P89 * pow(P89, P61 - 2, P61)) % M


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


def rowsum(d, p):
    K = 8
    s = [0] * K
    for v in d["table"].values():
        for k in range(K):
            s[k] = (s[k] + v[k]) % p
    return s


def main():
    print("=" * 78)
    print("1. table sums vs info['total'] (is info['total'] really the sum of the stored table?)")
    for tag, p in [("g5_89_nat", P89), ("g5_61_nat", P61), ("g5_89_id", P89), ("g5_61_id", P61),
                   ("g5_89_nat_ab1o3_m2o7", P89)]:
        d = load(tag)
        s = rowsum(d, p)
        print("   %-24s sum==info['total']: %s" % (tag, s == d["info"]["total"]))

    print("2. natural total == point-frame total, prime by prime")
    for a, b in [("g5_89_nat", "g5_89_id"), ("g5_61_nat", "g5_61_id")]:
        print("   %-24s vs %-12s %s" % (a, b, load(a)["info"]["total"] == load(b)["info"]["total"]))
    print("   second frame vs frame (0,0) at p89: u-series equal?",
          load("g5_89_nat_ab1o3_m2o7")["info"]["total"] == load("g5_89_nat")["info"]["total"])
    print("   second frame vs frame (0,0) at p89: t_c polar equal?",
          load("g5_89_nat_ab1o3_m2o7")["info"]["total_tc"] == load("g5_89_nat")["info"]["total_tc"])

    print("3. CRT over the two primes of the natural total (u-series of u^8 * term)")
    M = P89 * P61
    B = 10 ** 24
    a = load("g5_89_nat")["info"]["total"]
    b = load("g5_61_nat")["info"]["total"]
    for k in range(8):
        x = crt(a[k], b[k])
        print("      [u^%d] %s" % (k, ratrec(x, M, B)))
    print("   leading coefficient == 1/660000 mod p89:", a[0] == ofrac(Fr(1, 660000), P89))
    print("   leading coefficient == 1/660000 mod p61:", b[0] == ofrac(Fr(1, 660000), P61))
    print("   1/660000 == |B_10|/(10*8*5^4):", Fr(1, 660000) == Fr(5, 66) / (10 * 8 * 5 ** 4))

    print("4. what the regularity claim tested, by number of edges of the outer graph")
    d = load("g5_89_nat")
    out = d["outer"]
    tab = d["table"]
    byE = {}
    for i, G in enumerate(out):
        byE.setdefault(G["E"], [0, 0])[0] += 1
    for i in tab:
        byE[out[int(i)]["E"]][1] += 1
    print("      E   graphs  with a stored row (i.e. some non-zero piece)")
    for E in sorted(byE):
        print("     %2d   %5d   %5d" % (E, byE[E][0], byE[E][1]))
    free = sum(n for E, (n, s) in byE.items() if E >= 8)
    canc = sum(s for E, (n, s) in byE.items() if 1 <= E <= 7)
    sel = sum(n for E, (n, s) in byE.items() if 1 <= E <= 7) - canc
    print("     4554 edge graphs = %d free (E >= 8: no outer-edge set of that size is generated, cap |S| <= K-1 = 7)"
          % free)
    print("                      + %d with only identically vanishing pieces (no row stored at all)" % sel)
    print("                      + %d that vanish by cancellation between stored pieces" % canc)
    print("     polar coefficients genuinely at stake = sum over E<=7 of (8-E) =",
          sum(n * (8 - E) for E, (n, s) in byE.items() if 1 <= E <= 7), "(headline would be 4554*8 =", 4554 * 8, ")")

    print("5. controls: does the leading coefficient c_5 = 1/660000 discriminate at all?")
    for tag in ("g5_89_nat", "g5_89_nat_ab1o3_m2o7", "g5_89_M0_out2", "g5_89_nat_t7_out2"):
        d = load(tag)
        i = d["info"]
        print("   %-24s failed=%d  c_vertex_equals_pattern=%s  vertex t_c polar nonzero at %s"
              % (tag, len(d["failed"]), i.get("c_vertex_equals_pattern"),
                 [k for k, x in enumerate(i["vertex"]) if x] if "vertex" in i else "-"))

    print("6. is additivity recorded anywhere in the JSON?")
    for tag in ("g5_89_nat", "g5_61_nat", "g5_89_nat_ab1o3_m2o7"):
        print("   %-24s info keys: %s" % (tag, sorted(load(tag)["info"].keys())))


main()
