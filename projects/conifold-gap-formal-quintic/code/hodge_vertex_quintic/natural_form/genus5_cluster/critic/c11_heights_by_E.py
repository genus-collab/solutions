# -*- coding: utf-8 -*-
"""c11_heights_by_E.py -- height of the per-graph point-frame coefficients as a function of the number of edges,
with the three primes 2^61-1, 2^89-1, 2^521-1.  The point: the regularity claim is "this residue is 0 mod p89 and
mod p61"; a non-zero rational can do that only if its numerator exceeds p89*p61 ~ 2^150, so the margin of the
test is (150 - max height in bits).  At genus 4 the per-graph values were small; at genus 5 they are not.
Run:  py -I -B c11_heights_by_E.py [k per edge count]
"""
import json, os, sys, time
from fractions import Fraction as Fr
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
P89 = (1 << 89) - 1
P61 = (1 << 61) - 1
P521 = (1 << 521) - 1


def crt_list(pairs):
    M = 1
    x = 0
    for a, p in pairs:
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
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    from g5_engine import make_split, graph_id
    cost = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))
    a89 = json.load(open(os.path.join(PKG, "results", "g5_89_id.json")))["table"]
    a61 = json.load(open(os.path.join(PKG, "results", "g5_61_id.json")))["table"]
    sp = make_split(P521, 5, ("0", "0"), 8, "id")
    byE = {}
    for t, G in enumerate(sp.inner_graphs):
        byE.setdefault(G["E"], []).append(t)
    print("%3s %6s  %-10s %-10s %s" % ("E", "graph", "min bits", "median", "MAX bits (margin to 150)"), flush=True)
    worst = 0
    for E in sorted(byE):
        sel = sorted(byE[E], key=lambda t: cost[str(t)]["est_split"])[:k]
        for t in sel:
            t0 = time.time()
            row = graph_id(sp, t)
            hs = []
            nofit = 0
            for j in range(8):
                x, M = crt_list([(a61[str(t)][j], P61), (a89[str(t)][j], P89), (row[j], P521)])
                r = ratrec(x, M, 1 << 320)
                if r is None:
                    nofit += 1
                    continue
                hs.append(max(abs(r.numerator), r.denominator).bit_length())
            if not hs:
                print("%3d %6d  no reconstruction" % (E, t), flush=True)
                continue
            hs.sort()
            worst = max(worst, hs[-1])
            print("%3d %6d  %-10d %-10d %d  (margin %d bits)%s"
                  % (E, t, hs[0], hs[len(hs) // 2], hs[-1], 150 - hs[-1],
                     "  [%d coefficients beyond 2^320]" % nofit if nofit else ""), flush=True)
    print("largest height over the sample: 2^%d; a false zero at both primes needs > 2^150" % worst)


main()
