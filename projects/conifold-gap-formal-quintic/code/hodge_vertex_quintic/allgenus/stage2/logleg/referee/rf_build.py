"""Referee: build omega_{g,n} as symmetric tensors (dict sorted-index-tuple -> value mod P), with a symmetry check."""
import time
from collections import defaultdict, Counter
from math import comb

def distinct(S):
    out, last = [], None
    for q, x in enumerate(S):
        if x != last:
            out.append((x, S[:q] + S[q + 1:])); last = x
    return out

def merge(A, B):
    return tuple(sorted(A + B))

def stable(g, n):
    return 2 * g - 2 + n > 0

def build(cv, order, log=print):
    P = cv.P
    W = {}
    for (g, n1) in order:          # compute omega_{g,n1}, n1 = n+1
        t0 = time.time()
        n = n1 - 1
        out = defaultdict(int)
        if (g, n1) == (0, 3):
            for b in range(5):
                for k0 in range(cv.KM + 1):
                    for k1 in range(cv.KM + 1):
                        for k2 in range(k1, cv.KM + 1):   # one labelled assignment per key
                            f = [0] * (cv.E - cv.V0 + 1); f[0 - cv.V0] = 1
                            r = (-2 * (2 * k1 + 1) * (2 * k2 + 1)) * cv._resk(b, k0, f, 2 * k1 + 2 * k2) % P
                            if r:
                                out[(cv.idx(b, k0), (cv.idx(b, min(k1, k2)), cv.idx(b, max(k1, k2))))] += r
        elif (g, n1) == (1, 1):
            for b in range(5):
                t = cv.tser[b]; N = cv.N
                td = [(i * t[i]) % P for i in range(1, N)] + [0]
                tdm = [(-x) % P if i % 2 else x for i, x in enumerate(td)]
                D = [(2 * t[i + 1]) % P if i % 2 == 0 else 0 for i in range(N - 1)] + [0]
                from rf_tr import smul, sinv
                Di = sinv(D, N, P); G = smul(smul(td, tdm, N, P), smul(Di, Di, N, P), N, P)
                # B(q, sigma q) = -zeta^-2 G ; array over exponents
                f = [0] * (cv.E - cv.V0 + 1)
                for i in range(N):
                    ex = i - 2
                    if cv.V0 <= ex <= cv.E: f[ex - cv.V0] = (-G[i]) % P
                for k0 in range(cv.KM + 1):
                    r = cv._resk(b, k0, f, 0)
                    if r: out[(cv.idx(b, k0), ())] += r
        else:
            if g >= 1 and stable(g - 1, n + 2):
                for S, w in W[(g - 1, n + 2)].items():
                    for i, S1 in distinct(S):
                        for j, J in distinct(S1):
                            for a0, val in cv.Clist(i, j):
                                out[(a0, J)] += val * w
            for g1 in range(g + 1):
                g2 = g - g1
                for m1 in range(n + 1):
                    m2 = n - m1
                    if not (stable(g1, m1 + 1) and stable(g2, m2 + 1)): continue
                    L1 = [(i, I, w) for S, w in W[(g1, m1 + 1)].items() for i, I in distinct(S)]
                    L2 = [(j, I, w) for S, w in W[(g2, m2 + 1)].items() for j, I in distinct(S)]
                    for i, I1, w1 in L1:
                        c1 = Counter(I1)
                        for j, I2, w2 in L2:
                            cl = cv.Clist(i, j)
                            if not cl: continue
                            A = merge(I1, I2)
                            cA = Counter(A); wt = 1
                            for x, mu in c1.items(): wt *= comb(cA[x], mu)
                            f = w1 * w2 * wt % P
                            for a0, val in cl:
                                out[(a0, A)] += val * f
            if n >= 1 and stable(g, n):
                for S, w in W[(g, n)].items():
                    for j, R in distinct(S):
                        for a0, u, val in cv.Blist(j):
                            A = tuple(sorted(R + (u,)))
                            out[(a0, A)] += val * w * A.count(u)
        # symmetry check and symmetric storage
        full = {}
        for (a0, A), val in out.items():
            S = tuple(sorted(A + (a0,)))
            full.setdefault(S, None)
        Wn, bad = {}, 0
        for S in full:
            vals = set(out.get((a0, R), 0) % P for a0, R in distinct(S))
            if len(vals) != 1: bad += 1
            v = vals.pop()
            if v: Wn[S] = v
        maxk = max((max(cv.bk(i)[1] for i in S) for S in Wn), default=-1)
        sumk = max((sum(cv.bk(i)[1] for i in S) for S in Wn), default=-1)
        log(f"  omega_({g},{n1}): {len(Wn)} symmetric dxi-coefficients, symmetry violations {bad}, max k {maxk}, "
            f"max Sum k {sumk} (bound 3g-3+n = {3*g-3+n1}), {time.time()-t0:.1f} s")
        assert bad == 0 and maxk < cv.KM
        W[(g, n1)] = Wn
    return W

def closure(targets):
    need, stack = set(), list(targets)
    while stack:
        g, n = stack.pop()
        if (g, n) in need or not stable(g, n) or n < 1: continue
        need.add((g, n))
        stack.append((g - 1, n + 1)) if g >= 1 else None
        stack.append((g, n - 1))
        for g1 in range(g + 1):
            for m1 in range(n):
                stack += [(g1, m1 + 1), (g - g1, n - 1 - m1 + 1)]
    return sorted(need, key=lambda t: (2 * t[0] - 2 + t[1], t[0]))
