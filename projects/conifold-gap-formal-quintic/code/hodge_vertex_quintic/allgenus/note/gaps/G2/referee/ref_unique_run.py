# -*- coding: utf-8 -*-
"""ref_unique_run.py -- driver of ref_unique_core.py.  For random one-pole systems with (H), (EP) (two seeds):
(1) per Z/5-degree delta, the null space of the truncated gauge equation projected onto the window k <= K-1,
    m <= M - k - margin has dimension #{k <= K-1 : k = delta mod 5}  (Theorem U (b)+(c): V(A, At) = P^nat k[[z]]);
(2) lem:gau-lead: in degree delta the coordinates with k < delta vanish on the window and those with k = delta span the line
    of frak d Theta_U (computed independently from the weights);
(3) the pair (At, At): the window projection is spanned by z^r 1 (consequence 2 of s.5.4: V(At, At) = k[[z]] 1);
(4) unsplit cross-check on a small box: the null space of the full system has the sum of the per-degree dimensions.
Failing controls: no pole (b_1, bt_1 units, (EP) kept): K per degree, i.e. rank 5 instead of 1; unequal products: 0.
Usage: py -I -B ref_unique_run.py [K] [M] [margin].  Stdlib only; imports ref_unique_core only.
"""
import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_unique_core import PR, systems, equations, nullspace, rank_of, theta_frakd

FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def window(K, M, margin, slope=1):
    return [(k, bb, m) for k in range(K) for bb in range(5) for m in range(M + 1) if m <= M - slope * k - margin]


def solve_degree(b, bt, K, M, delta):
    cols = [(k, bb, m) for k in range(K + 1) for bb in range(5) for m in range(M + 1)]
    return nullspace(equations(b, bt, K, M, delta), cols)


def full_equations(b, bt, K, M):
    rows = []
    for k in range(K + 1):
        for a in range(5):
            for bb in range(5):
                for j in range(-1, M):
                    r = {}
                    if k >= 1 and 0 < j + 1 <= M:
                        r[(k - 1, a, bb, j + 1)] = (j + 1) % PR
                    for e, c in b[a].items():
                        if 0 <= j - e <= M:
                            key = (k, (a + 1) % 5, bb, j - e); r[key] = (r.get(key, 0) - c) % PR
                    for e, c in bt[(bb - 1) % 5].items():
                        if 0 <= j - e <= M:
                            key = (k, a, (bb - 1) % 5, j - e); r[key] = (r.get(key, 0) + c) % PR
                    r = {c: v for c, v in r.items() if v}
                    if r:
                        rows.append(r)
    return rows


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    margin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    N = M + 3
    print("ref_unique_run.py: random one-pole systems mod 2^61-1; unknowns k <= %d, m <= %d; window k <= %d, m <= M - k - %d"
          % (K, M, K - 1, margin))
    W = window(K, M, margin)
    for seed in (20261115, 7):
        rng = random.Random(seed)
        b, bt = systems(rng, N)
        fdT = theta_frakd(b, bt, N)
        print(" seed %d: (H), (EP) random system" % seed)
        for delta in range(5):
            basis = solve_degree(b, bt, K, M, delta)
            pred = len([k for k in range(K) if k % 5 == delta])
            dw = rank_of(basis, W)
            check("degree %d: window dimension %d = #{k <= %d : k = %d mod 5} = %d" % (delta, dw, K - 1, delta, pred),
                  dw == pred, "(truncated null space: %d)" % len(basis))
            low = [c for c in W if c[0] < delta]
            lead = [c for c in W if c[0] == delta and c[2] <= M - delta - margin]
            r0 = rank_of(basis, low) if low else 0
            r1 = rank_of(basis, lead)
            fdv = {(delta, bb, m): fdT[bb].get(m, 0) for bb in range(5) for m in range(M + 1)}
            r1b = rank_of(basis + [fdv], lead)
            check("degree %d: lem:gau-lead: P_k = 0 for k < %d, P_%d spans the line of frak d Theta_U" % (delta, delta, delta),
                  r0 == 0 and r1 == 1 and r1b == 1, "(ranks %d, %d, %d)" % (r0, r1, r1b))
        rng = random.Random(seed + 1)
        b2, _ = systems(rng, N, same=True)
        ok = True; info = []
        for delta in range(5):
            basis = solve_degree(b2, b2, K, M, delta)
            pred = [r for r in range(K) if r % 5 == delta]
            exp = [{(r, bb, 0): 1 for bb in range(5)} for r in pred]   # z^r 1: degree r, entries (bb, bb), u^0
            dw = rank_of(basis, W); dwe = rank_of(basis + exp, W)
            ok = ok and dw == len(pred) and dwe == len(pred); info.append((dw, dwe, len(pred)))
        check("pair (At, At): window spanned by z^r 1, every degree (V(At, At) = k[[z]] 1)", ok, str(info))
        rng = random.Random(seed + 2)
        b3, bt3 = systems(rng, N, pole=False)
        ok = True; info = []
        for delta in range(5):
            dw = rank_of(solve_degree(b3, bt3, K, M, delta), W)
            ok = ok and dw == K; info.append(dw)
        check("control: no pole ((H) violated, (EP) kept): every degree has window dimension K = %d (rank 5, not 1)" % K, ok, str(info))
        rng = random.Random(seed + 3)
        b4, bt4 = systems(rng, N, equal=False)
        mb = M - 12                      # unequal products: P_k = 0 is forced only ~5 u-orders per z-order below the truncation
        box = [(k, bb, m) for k in range(2) for bb in range(5) for m in range(mb + 1)]
        ok = True; info = []
        for delta in range(5):
            dw = rank_of(solve_degree(b4, bt4, K, M, delta), box)
            de = rank_of(solve_degree(b, bt, K, M, delta), box)
            pe = len([k for k in range(2) if k % 5 == delta])
            ok = ok and dw == 0 and de == pe; info.append((dw, de))
        check("control: unequal products: dimension 0 on the box k <= 1, m <= %d (equal products, same box: 1, 1, 0, 0, 0)" % mb,
              ok, "(unequal, equal) = %s" % info)
        if seed == 7:
            for mg in (0, 1):
                dw = rank_of(solve_degree(b, bt, K, M, 1), window(K, M, mg))
                print("  [info] margin %d (window touching the truncation), degree 1: dimension %d (prediction %d)"
                      % (mg, dw, len([k for k in range(K) if k % 5 == 1])))
    Ks, Ms = 4, 8
    rng = random.Random(99)
    b, bt = systems(rng, Ms + 3)
    cols = [(k, a, bb, m) for k in range(Ks + 1) for a in range(5) for bb in range(5) for m in range(Ms + 1)]
    full = len(nullspace(full_equations(b, bt, Ks, Ms), cols))
    split = sum(len(solve_degree(b, bt, Ks, Ms, d)) for d in range(5))
    check("unsplit system (k <= %d, m <= %d): null space %d = sum of the five degree parts %d (lem:gau-degree)" % (Ks, Ms, full, split),
          full == split)
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
