"""Referee: reproduce the exact counts quoted in ../DICT_PROOF.md s.8 and the JSON with rf_core (own code).
Usage: py -I -B rf_counts.py
  T2: monomials v^m p^j, w <= 80, 0 <= m <= w//5 + w//4 + 3: 861 in L, 859 not in L (worker's range).
  T3: basis monomials of M_w (p^a q^c and log p v^m x^n), w <= 80: 1042.
  X3: monomials p^a q^c with a + 4c <= 40: 231.
  K_0 totals Sum_w d(w): 2640 (w <= 160), 10400 (w <= 320), 16200 (w <= 400).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import Ar, in_L, in_Qpq  # noqa: E402

ar = Ar(None)


def d(w):
    return w // 5 + 1 - (1 if w % 5 == 0 else 0)


if __name__ == "__main__":
    nin = nout = agree = 0
    for w in range(81):
        for m in range(0, w // 5 + w // 4 + 4):
            j = w - 5 * m
            ok = in_L(ar, {j: ar.q(1)}, w)
            nin += ok
            nout += not ok
            agree += ok == in_Qpq({j: 1}, w)
    nb = sum(sum(1 for c in range(w // 4 + 1)) + sum(1 for n in range(w // 4 + 1) if (w - 4 * n) % 5 == 0) for w in range(81))
    nx3 = sum(1 for a in range(41) for c in range(11) if a + 4 * c <= 40)
    tots = [sum(d(w) for w in range(W + 1)) for W in (160, 320, 400)]
    print("T2 (worker's range): %d in L, %d not in L (worker: 861, 859); support rule agrees on %d/%d" % (nin, nout, agree, nin + nout))
    print("T3 basis monomials of M_w, w <= 80: %d (worker: 1042)" % nb)
    print("X3 monomials p^a q^c, a + 4c <= 40: %d (worker: 231)" % nx3)
    print("Sum_w d(w) for w <= 160, 320, 400: %s (worker/critic: 2640, 10400, 16200)" % tots)
    ok = (nin, nout, agree, nb, nx3, tots) == (861, 859, 1720, 1042, 231, [2640, 10400, 16200])
    print("RESULT:", "ALL COUNTS REPRODUCED" if ok else "MISMATCH")
