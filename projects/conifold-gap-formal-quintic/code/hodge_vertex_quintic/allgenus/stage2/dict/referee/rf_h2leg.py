"""Referee: sensitivity of Theorem D(iii)(d) to the primary H^2 leg (the only primary leg not fixed by its weight).
Usage: py -I -B rf_h2leg.py
In weight 4, M_4 = <p^4, q> (+) log p <x>.  Lemma 3 predicts: with Phi_0^2 := a p^4 + b q (q = v/p), the leg map
S_w -> N_w is bijective for all w iff a != b/4, i.e. iff Phi_0^2 is not a multiple of x = p^4/4 + q modulo nothing else.
Checked by generic elimination for w <= 40, exact over Q.
"""
import os
import sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import Ar, legs_iter, nf_vx, x_coords, triples, rank, wt  # noqa: E402

ar = Ar(None)

if __name__ == "__main__":
    agree = True
    for a, b in [(0, 1), (1, 0), (1, 4), (Fr(1, 4), 1), (2, 3), (Fr(1, 2), 2), (-3, 5)]:
        phi = {0: ({1: 1}, {}), 1: ({}, {0: 1}), 2: (ar.lift({4: a, -1: b}), {}), 3: ({3: 1}, {}), 4: ({2: 1}, {})}
        T = legs_iter(ar, 12, phi)
        bad = []
        for w in range(41):
            S = triples(w)
            rows = [dict([(("p", j), c) for j, c in nf_vx(ar, T[(bb, k)][0], wt(bb, k))[0].items()] +
                         [(("L", n), c) for n, c in x_coords(ar, T[(bb, k)][1], wt(bb, k))[0].items()]) for (bb, k, m) in S]
            if rank(ar, rows) < len(S):
                bad.append(w)
        pred = Fr(a) != Fr(b, 4)
        agree = agree and pred == (not bad)
        print("Phi_0^2 = (%s) p^4 + (%s) q: predicted bijective %s; failing weights (w <= 40): %s" % (a, b, pred, bad[:6]))
    print("RESULT:", "PREDICTION CONFIRMED" if agree else "PREDICTION FAILED")
