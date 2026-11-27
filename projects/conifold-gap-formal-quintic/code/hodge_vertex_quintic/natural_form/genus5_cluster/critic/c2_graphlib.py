# -*- coding: utf-8 -*-
"""c2_graphlib.py -- independent library for genus-g leg-free stable graphs (critic's own code, shares nothing
with g5_graphs.py / fq3_graphs.py). Representation: gv = tuple of vertex genera, M = symmetric multiplicity
dict {(i,j): m} with i <= j (loops i == j).  |Aut| of the stable graph = (vertex automorphisms) *
prod_loops 2^l l! * prod_{i<j} m!.  Automorphisms found by backtracking, NOT by refinement + individualisation.
"""
import math
from itertools import permutations


def mat_of(V, ed):
    M = {}
    for u, w in ed:
        k = (u, w) if u <= w else (w, u)
        M[k] = M.get(k, 0) + 1
    return M


def deg(V, M):
    d = [0] * V
    for (i, j), m in M.items():
        if i == j:
            d[i] += 2 * m
        else:
            d[i] += m
            d[j] += m
    return d


def loops(V, M):
    return [M.get((v, v), 0) for v in range(V)]


def connected(V, M):
    if V == 1:
        return True
    adj = [[] for _ in range(V)]
    for (i, j) in M:
        if i != j:
            adj[i].append(j)
            adj[j].append(i)
    seen = {0}
    st = [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return len(seen) == V


def nedges(M):
    return sum(M.values())


def genus(V, gv, M):
    return sum(gv) + nedges(M) - V + 1


def stable(V, gv, M):
    d = deg(V, M)
    return all(d[v] >= (3 if gv[v] == 0 else (1 if gv[v] == 1 else 0)) for v in range(V))


def _apply(perm, M):
    out = {}
    for (i, j), m in M.items():
        a, b = perm[i], perm[j]
        out[(a, b) if a <= b else (b, a)] = m
    return out


def vertex_auts(V, gv, M):
    """number of vertex permutations preserving gv and M, by backtracking on vertex images."""
    d = deg(V, M)
    lo = loops(V, M)
    sig = [(gv[v], d[v], lo[v]) for v in range(V)]
    nbr = [{} for _ in range(V)]
    for (i, j), m in M.items():
        if i != j:
            nbr[i][j] = m
            nbr[j][i] = m
    cnt = [0]
    img = [-1] * V
    used = [False] * V

    def rec(v):
        if v == V:
            cnt[0] += 1
            return
        for w in range(V):
            if used[w] or sig[w] != sig[v]:
                continue
            ok = True
            for u in range(v):
                if nbr[v].get(u, 0) != nbr[w].get(img[u], 0):
                    ok = False
                    break
            if ok:
                img[v] = w
                used[w] = True
                rec(v + 1)
                used[w] = False
                img[v] = -1
    rec(0)
    return cnt[0]


def aut_total(V, gv, M):
    r = vertex_auts(V, gv, M)
    for (i, j), m in M.items():
        r *= math.factorial(m)
        if i == j:
            r *= 2 ** m
    return r


def canon_brute(V, gv, M):
    """canonical key by minimising over the vertex permutations compatible with the local signature
    (brute force over the compatible set, grouped by signature; independent of refinement)."""
    d = deg(V, M)
    lo = loops(V, M)
    sig = [(gv[v], d[v], lo[v]) for v in range(V)]
    order = sorted(range(V), key=lambda v: sig[v])
    groups = []
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and sig[order[j]] == sig[order[i]]:
            j += 1
        groups.append(order[i:j])
        i = j
    best = None
    slots = []
    pos = 0
    for gr in groups:
        slots.append((gr, list(range(pos, pos + len(gr)))))
        pos += len(gr)

    def build(perm):
        gv2 = [0] * V
        for v in range(V):
            gv2[perm[v]] = gv[v]
        e2 = sorted(tuple(sorted((perm[i], perm[j]))) for (i, j), m in M.items() for _ in range(m))
        return (tuple(gv2), tuple(e2))

    def rec(k, perm):
        nonlocal best
        if k == len(slots):
            key = build(perm)
            if best is None or key < best:
                best = key
            return
        gr, tgt = slots[k]
        for p in permutations(tgt):
            for v, t in zip(gr, p):
                perm[v] = t
            rec(k + 1, perm)
    rec(0, [0] * V)
    return best


def iso(V, gv1, M1, gv2, M2):
    return canon_brute(V, gv1, M1) == canon_brute(V, gv2, M2)
