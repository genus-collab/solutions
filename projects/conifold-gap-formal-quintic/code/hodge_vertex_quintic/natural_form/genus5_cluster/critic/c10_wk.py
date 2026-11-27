# -*- coding: utf-8 -*-
"""c10_wk.py -- independent checks of the Witten-Kontsevich numbers fq_hodge.wk used for the (one-colour, trivial
theory) vertices.  These are RATIONAL inputs: they are identical at every prime, so no number of primes tests them.
fq_hodge._wk is the DVV/Virasoro recursion with the string equation built in, so a DVV-based check would be
circular.  Checks used here instead:
  (1) <tau_{3g-2}>_g = 1/(24^g g!)  for g = 1..5 (Witten);
  (2) genus 0: <tau_{d_1}..tau_{d_n}>_0 = (n-3)! / prod d_i!;
  (3) the dilaton equation <tau_1 X>_g = (2g-2+n) <X>_g, on all partitions needed at genus <= 5, n <= 5;
  (4) a few published genus-2 and genus-3 values.
Run:  py -I -B c10_wk.py
"""
import os
import sys, math
from fractions import Fraction as Fr
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
from fq_hodge import wk

FAIL = []


def ck(name, ok, d=""):
    print(("PASS " if ok else "FAIL ") + name + ((": " + str(d)) if d else ""))
    if not ok:
        FAIL.append(name)


# (1)
bad = []
for g in range(1, 6):
    v = wk(g, (3 * g - 2,))
    e = Fr(1, 24 ** g * math.factorial(g))
    if v != e:
        bad.append((g, v, e))
ck("<tau_{3g-2}>_g = 1/(24^g g!) for g = 1..5", not bad, bad)

# (2)
def parts(n, s):
    if n == 1:
        yield (s,)
        return
    for x in range(s + 1):
        for r in parts(n - 1, s - x):
            yield (x,) + r
bad = []
for n in range(3, 8):
    for d in parts(n, n - 3):
        v = wk(0, d)
        e = Fr(math.factorial(n - 3), math.prod(math.factorial(x) for x in d))
        if v != e:
            bad.append((d, v, e))
ck("genus 0, n <= 7: <prod tau_d>_0 = (n-3)!/prod d!", not bad, bad[:5])

# (3) dilaton
bad = []
cnt = 0
for g in range(1, 6):
    for n in range(1, 5):
        tot = 3 * g - 3 + n
        if tot < 0:
            continue
        for d in parts(n, tot):
            if 2 * g - 2 + n <= 0:
                continue
            cnt += 1
            if wk(g, (1,) + d) != (2 * g - 2 + n) * wk(g, d):
                bad.append((g, d))
ck("dilaton equation on %d correlators with g <= 5, n <= 4" % cnt, not bad, bad[:5])

# (4) published values
known = {(2, (4,)): Fr(1, 1152), (2, (0, 5)): Fr(1, 1152), (2, (1, 4)): Fr(1, 384),
         (2, (2, 3)): Fr(29, 5760), (2, (2, 2, 2)): Fr(7, 240), (3, (7,)): Fr(1, 82944),
         (1, (1,)): Fr(1, 24), (3, (0, 8)): Fr(1, 82944), (1, (0, 2)): Fr(1, 24),
         (0, (0, 0, 0)): Fr(1, 1)}
bad = [(g, d, wk(g, d), v) for (g, d), v in known.items() if wk(g, d) != v]
ck("published low-genus values", not bad, bad)
print("WK CHECKS:", "ALL PASSED" if not FAIL else FAIL)
