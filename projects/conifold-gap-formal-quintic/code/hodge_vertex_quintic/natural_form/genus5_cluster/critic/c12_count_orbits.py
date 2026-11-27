# -*- coding: utf-8 -*-
"""c12_count_orbits.py -- independent COUNT of the leg-free stable graphs of genus g, with no canonical form at all.
The set of labelled structures (gv as a vector on [V], symmetric multiplicity matrix M) carries an S_V action whose
orbits are the isomorphism classes, so the number of classes with V vertices is
      (1 / V!) * sum over labelled structures of |Aut_vertices(gv, M)|
(orbit counting / Burnside).  |Aut_vertices| is computed by backtracking (c2_graphlib.vertex_auts).  This is a route
independent of g5_graphs.canon_count, of fq3_graphs.degenerations and of the critic's own canon_brute.
Also accumulates sum 1/|Aut| per edge count, a second invariant to compare with the shipped list.
Run:  py -I -B c12_count_orbits.py [g]
"""
import sys, os, math, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c2_graphlib as L
from c3_enumerate import degseqs, fill, genvecs


def count(g):
    byE = {}
    autsum = {}
    for V in range(1, max(1, 2 * g - 2) + 1):
        slots = [(i, j) for i in range(V) for j in range(i, V)]
        fac = math.factorial(V)
        for s in range(0, g + 1):
            E = g - 1 + V - s
            if E < 0 or (V > 1 and E < V - 1):
                continue
            for gv in genvecs(V, s):
                for d in degseqs(V, gv, E):
                    for M in fill(V, list(d), slots, 0, {}):
                        if L.nedges(M) != E or not L.connected(V, M) or not L.stable(V, gv, M):
                            continue
                        if L.genus(V, gv, M) != g:
                            continue
                        a = L.vertex_auts(V, gv, M)
                        byE[E] = byE.get(E, Fr(0)) + Fr(a, fac)
                        autsum[E] = autsum.get(E, Fr(0)) + Fr(a, fac) / L.aut_total(V, gv, M)
    return byE, autsum


if __name__ == "__main__":
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    t0 = time.time()
    byE, autsum = count(g)
    tot = sum(byE.values())
    print("genus %d: %s classes in %.0f s" % (g, tot, time.time() - t0))
    print("  by number of edges:", [(E, byE[E]) for E in sorted(byE)])
    print("  sum 1/|Aut| by number of edges:", [(E, str(autsum[E])) for E in sorted(autsum)])
