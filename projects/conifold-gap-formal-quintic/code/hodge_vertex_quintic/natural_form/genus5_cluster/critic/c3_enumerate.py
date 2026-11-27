# -*- coding: utf-8 -*-
"""c3_enumerate.py -- independent enumeration of the leg-free stable graphs of genus g (critic's own code).
Route: for each V, each genus vector gv (sum s), E = g - 1 + V - s edges; each admissible degree sequence
(d_v >= 3 if gv_v = 0, >= 1 otherwise, sum d = 2E); then every symmetric multiplicity matrix realising d;
keep connected ones; canonicalise with c2_graphlib.canon_brute.  Output: dict key -> |Aut|.
Run:  py -I -B c3_enumerate.py [g]
"""
import sys, os
from itertools import combinations_with_replacement
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c2_graphlib as L


def degseqs(V, gv, E):
    lo = [3 if gv[v] == 0 else (1 if gv[v] == 1 else 0) for v in range(V)]
    tot = 2 * E
    base = sum(lo)
    if base > tot:
        return
    extra = tot - base
    # distribute `extra` among V vertices, each d_v <= tot
    def rec(v, left, cur):
        if v == V:
            if left == 0:
                yield tuple(cur)
            return
        for x in range(left + 1):
            cur.append(lo[v] + x)
            yield from rec(v + 1, left - x, cur)
            cur.pop()
    yield from rec(0, extra, [])


def fill(V, d, slots, k, M):
    """assign multiplicities to the slots (i<=j) in order, respecting residual degrees d"""
    if k == len(slots):
        if all(x == 0 for x in d):
            yield dict(M)
        return
    i, j = slots[k]
    cap = d[i] // 2 if i == j else min(d[i], d[j])
    # prune: the remaining slots must be able to absorb the residual degree
    for m in range(cap, -1, -1):
        if i == j:
            d[i] -= 2 * m
        else:
            d[i] -= m
            d[j] -= m
        if m:
            M[(i, j)] = m
        rest = sum(d[x] for x in range(V) if all(x not in s for s in slots[k + 1:]))
        if rest == 0:
            yield from fill(V, d, slots, k + 1, M)
        if m:
            del M[(i, j)]
        if i == j:
            d[i] += 2 * m
        else:
            d[i] += m
            d[j] += m


def genvecs(V, s):
    def rec(v, left, cur):
        if v == V:
            if left == 0:
                yield tuple(cur)
            return
        for x in range(left + 1):
            cur.append(x)
            yield from rec(v + 1, left - x, cur)
            cur.pop()
    yield from rec(0, s, [])


def enumerate_graphs(g):
    out = {}
    Vmax = max(1, 2 * g - 2)
    for V in range(1, Vmax + 1):
        for s in range(0, g + 1):
            E = g - 1 + V - s
            if E < 0 or (V > 1 and E < V - 1):
                continue
            for gv in genvecs(V, s):
                slots = [(i, j) for i in range(V) for j in range(i, V)]
                for d in degseqs(V, gv, E):
                    for M in fill(V, list(d), slots, 0, {}):
                        if L.nedges(M) != E:
                            continue
                        if not L.connected(V, M):
                            continue
                        if not L.stable(V, gv, M):
                            continue
                        if L.genus(V, gv, M) != g:
                            continue
                        key = L.canon_brute(V, gv, M)
                        if key not in out:
                            out[key] = L.aut_total(V, gv, M)
    return out


if __name__ == "__main__":
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    import time
    t0 = time.time()
    d = enumerate_graphs(g)
    byE = {}
    for k in d:
        byE[len(k[1])] = byE.get(len(k[1]), 0) + 1
    print("genus %d: %d leg-free stable graphs (%.0f s)" % (g, len(d), time.time() - t0))
    print("  by number of edges:", sorted(byE.items()))
