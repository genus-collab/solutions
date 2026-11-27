# Stable graphs of genus g without legs, by iterated degeneration from the single vertex. Evidence only.
import sys, os, json
from itertools import permutations, product


def canon(gv, A):
    V = len(gv)
    # colour refinement
    col = [(gv[i], A[i][i], sum(A[i])) for i in range(V)]
    for _ in range(3):
        col = [(col[i], tuple(sorted((A[i][j], col[j]) for j in range(V) if j != i))) for i in range(V)]
        vals = sorted(set(col))
        idx = {c: k for k, c in enumerate(vals)}
        col = [idx[c] for c in col]
    classes = {}
    for i, c in enumerate(col):
        classes.setdefault(c, []).append(i)
    keys = sorted(classes)
    best = None
    for ps in product(*[permutations(classes[k]) for k in keys]):
        order = [x for p in ps for x in p]
        t = (tuple(gv[i] for i in order), tuple(tuple(A[i][j] for j in order) for i in order))
        if best is None or t < best:
            best = t
    return best


def degenerations(gv, A):
    V = len(gv)
    out = []
    for v in range(V):
        if gv[v] >= 1:
            g2 = list(gv)
            g2[v] -= 1
            A2 = [list(r) for r in A]
            A2[v][v] += 1
            out.append((g2, A2))
    for v in range(V):
        nbrs = [w for w in range(V) if w != v and A[v][w] > 0]
        l = A[v][v]
        for g1 in range(gv[v] + 1):
            g2v = gv[v] - g1
            for split in product(*[range(A[v][w] + 1) for w in nbrs]):
                for l1 in range(l + 1):
                    for l2 in range(l - l1 + 1):
                        l12 = l - l1 - l2
                        n1 = sum(split) + 2 * l1 + l12 + 1
                        n2 = sum(A[v][w] for w in nbrs) - sum(split) + 2 * l2 + l12 + 1
                        if 2 * g1 - 2 + n1 <= 0 or 2 * g2v - 2 + n2 <= 0:
                            continue
                        # new vertex index V goes with part 2
                        B = [list(r) + [0] for r in A] + [[0] * (V + 1)]
                        for w, s in zip(nbrs, split):
                            B[v][w] = B[w][v] = s
                            B[V][w] = B[w][V] = A[v][w] - s
                        B[v][v] = l1
                        B[V][V] = l2
                        B[v][V] = B[V][v] = l12 + 1
                        g3 = list(gv) + [g2v]
                        g3[v] = g1
                        out.append((g3, B))
    return out


def stable_graphs(g):
    level = {canon([g], [[0]])}
    allg = set(level)
    for e in range(3 * g - 3):
        nxt = set()
        for gv, A in level:
            for g2, A2 in degenerations(list(gv), [list(r) for r in A]):
                nxt.add(canon(g2, A2))
        nxt -= allg
        allg |= nxt
        level = nxt
    return sorted(allg, key=lambda t: (len(t[0]), t))


if __name__ == '__main__':
    g = int(sys.argv[1])
    gs = stable_graphs(g)
    print(g, len(gs))
    if len(sys.argv) > 2:
        json.dump([[list(a), [list(r) for r in b]] for a, b in gs], open(sys.argv[2], 'w'))
