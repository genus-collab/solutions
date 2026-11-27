# -*- coding: utf-8 -*-
"""c7_height.py -- how big are the per-graph rationals? A non-zero rational a/b can be 0 mod BOTH primes only if
p89*p61 | a, i.e. |a| >= 2^150 ~ 1.4e45.  This script measures the actual height of the per-graph coefficients
(point-frame terms, CRT of the two id runs) to see how much margin the two-prime zero test really has.
Run:  py -I -B c7_height.py
"""
import json, os
from fractions import Fraction as Fr
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic", "results")
P89 = (1 << 89) - 1
P61 = (1 << 61) - 1
M = P89 * P61


def crt(a, b):
    return (a * P61 * pow(P61, P89 - 2, P89) + b * P89 * pow(P89, P61 - 2, P61)) % M


def ratrec(a, B):
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


a = json.load(open(os.path.join(D, "g5_89_id.json")))["table"]
b = json.load(open(os.path.join(D, "g5_61_id.json")))["table"]
B = 1 << 70
ok = 0
fail = 0
hmax = 0
worst = None
for t in a:
    for k in range(8):
        x = crt(a[t][k], b[t][k])
        r = ratrec(x, B)
        if r is None:
            fail += 1
            continue
        ok += 1
        h = max(abs(r.numerator), r.denominator)
        if h > hmax:
            hmax = h
            worst = (t, k, r)
print("per-graph point-frame coefficients: %d of %d reconstruct with height <= 2^70" % (ok, ok + fail))
print("largest height seen: %.3e  (graph %s, [u^%d] = %s)" % (hmax, worst[0], worst[1], worst[2]))
print("a non-zero rational vanishes mod both primes only if its numerator is >= p89*p61 = %.3e" % M)
print("margin: %.1f orders of magnitude" % (len(str(M)) - len(str(hmax))))
num = 0
for t in a:
    for k in range(8):
        r = ratrec(crt(a[t][k], b[t][k]), B)
        if r is not None and r.numerator:
            num = max(num, abs(r.numerator))
print("largest NUMERATOR seen: %.3e (what would have to exceed %.3e for a false zero)" % (num, M))
