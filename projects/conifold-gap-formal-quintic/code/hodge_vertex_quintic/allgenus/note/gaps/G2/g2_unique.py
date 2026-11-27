# -*- coding: utf-8 -*-
"""g2_unique.py -- evidence for the uniqueness clause of Theorem U (G4) INCLUDING non-graded gauges, by a route independent
of the construction: the gauge equation z P' = A P - P At written as a linear system for the unknown coefficients
x(k,a,b,m) = [u^m](P_k)_(ab), k <= K, m <= M, with every equation that involves only these unknowns, solved mod p = 2^61-1.
The system splits by the Z/5-degree delta = k + b - a (mod 5) of the unknowns.  Prediction of Theorem U (every gauge is
P_can c(z)): the solutions of degree delta are P_can z^r Q[[z^5]], r = delta mod 5, so the projection of the solution space
onto the coordinates k <= K-1, m <= M - k - margin has dimension #{k <= K-1 : k = delta mod 5}.  Controls: unequal products (prediction 0),
the residue equations (u^-1) dropped, and windows reaching the truncation (artefacts expected).
Usage: py -I -B g2_unique.py [K] [M] [margin]; projection window k <= K-1, m <= M-k-margin.  Stdlib only; imports g2_core, g2_gauge.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g2_core import Ser, systems
from g2_gauge import gauge, times_c, inverse_series

PR = (1 << 61) - 1
FAIL = []


def check(name, ok, detail=""):
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)


def modp(x):
    return x.numerator % PR * pow(x.denominator % PR, PR - 2, PR) % PR


def coeffs(s, lo, hi):
    """{j: [u^j] s mod p} for lo <= j <= hi"""
    return {j: modp(s.coef(j)) for j in range(lo, hi + 1) if s.coef(j) != 0}


def equations(b, bt, K, M, delta, drop_residue=False):
    """rows (sparse dicts {column: value}) of the degree-delta part; columns x(k, b, m) with a = k + b - delta"""
    bc = [coeffs(x, -1, M + 1) for x in b]
    btc = [coeffs(x, -1, M + 1) for x in bt]
    col = lambda k, bb, m: (k, bb % 5, m)
    rows = []
    for k in range(K + 1):
        for bb in range(5):
            a = (k - 1 + bb - delta) % 5          # equation entry (a, bb) at z^k
            for j in range(0 if drop_residue else -1, M):
                r = {}
                if k >= 1 and j + 1 <= M:          # (P_(k-1)')_(a,bb) at u^j = (j+1) x(k-1, a, bb, j+1); entry (a,bb) has degree delta at k-1
                    r[col(k - 1, bb, j + 1)] = (j + 1) % PR
                for m in range(0, M + 1):          # - b_a (P_k)_(a+1, bb)
                    c = bc[a].get(j - m)
                    if c:
                        key = col(k, bb, m); r[key] = (r.get(key, 0) - c) % PR
                for m in range(0, M + 1):          # + (P_k)_(a, bb-1) bt_(bb-1)
                    c = btc[(bb - 1) % 5].get(j - m)
                    if c:
                        key = col(k, bb - 1, m); r[key] = (r.get(key, 0) + c) % PR
                r = {c: v for c, v in r.items() if v}
                if r:
                    rows.append(r)
    return rows


def ranks(rows, columns_first, columns_last):
    """forward elimination with the columns in the order columns_first + columns_last; returns (rank of the first block, total rank)"""
    cols = columns_first + columns_last
    order = {c: i for i, c in enumerate(cols)}
    nfirst = len(columns_first)
    rws = [dict(r) for r in rows]
    piv_first = piv_tot = 0
    byc = {}
    for i, r in enumerate(rws):
        lead = min(r, key=lambda c: order[c])
        byc.setdefault(order[lead], []).append(i)
    for ci in range(len(cols)):
        c0 = cols[ci]
        lst = [i for i in byc.get(ci, []) if rws[i] and min(rws[i], key=lambda c: order[c]) == c0]
        if not lst:
            continue
        p = lst[0]; pr_ = rws[p]
        inv = pow(pr_[c0], PR - 2, PR)
        piv_tot += 1
        if ci < nfirst:
            piv_first += 1
        for i in lst[1:]:
            r = rws[i]; f = r[c0] * inv % PR
            for c, v in pr_.items():
                w = (r.get(c, 0) - f * v) % PR
                if w:
                    r[c] = w
                else:
                    r.pop(c, None)
            if r:
                lead = min(r, key=lambda c: order[c])
                byc.setdefault(order[lead], []).append(i)
    return piv_first, piv_tot


def projected_dim(rows, K, M, Kp, margin, box=None):
    """dimension of the projection of the truncated solution space onto the triangle k <= Kp, m <= M - k - margin
    (or onto the box k <= box[0], m <= box[1])"""
    allc = [(k, bb, m) for k in range(K + 1) for bb in range(5) for m in range(M + 1)]
    if box is None:
        S = [c for c in allc if c[0] <= Kp and c[2] <= M - c[0] - margin]
    else:
        S = [c for c in allc if c[0] <= box[0] and c[2] <= box[1]]
    sset = set(S)
    Sc = [c for c in allc if c not in sset]
    rSc, rE = ranks(rows, Sc, S)
    return len(S) - rE + rSc, len(allc) - rE


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    margin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    Kp = K - 1
    print("g2_unique.py: unknowns k <= %d, m <= %d; projection onto the triangle k <= %d, m <= M - k - %d; mod 2^61-1" % (
        K, M, Kp, margin))
    for (al, be) in [(0, 0), (Fr(1, 3), Fr(-2, 7))]:
        S = systems(M + 8, al, be)
        b, bt = S["b"], S["bt"]
        print(" frame (%s, %s):" % (al, be))
        S2 = systems(M + K + 12, al, be)
        Pg, _ = gauge(S2["b"], S2["bt"], K)
        c11 = {5 * q: Pg[5 * q][1].coef(0) for q in range(K // 5 + 1)}
        Pc = times_c(Pg, inverse_series(c11, K))
        for delta in range(5):
            rows = equations(b, bt, K, M, delta)
            dproj, dW = projected_dim(rows, K, M, Kp, margin)
            pred = len([k for k in range(Kp + 1) if k % 5 == delta])
            check("degree %d: projected dimension %d = #{k <= %d : k = %d mod 5} = %d" % (delta, dproj, Kp, delta, pred),
                  dproj == pred, "(truncated solution space: dim %d)" % dW)
            vec = {(k, bb, m): (modp(Pc[k - delta][bb].coef(m)) if k >= delta else 0)
                   for k in range(K + 1) for bb in range(5) for m in range(M + 1)}
            bad = sum(1 for r in rows if sum(v * vec[c] for c, v in r.items()) % PR)
            check("degree %d: z^%d P_can solves every truncated equation (so it spans the predicted directions)" % (
                delta, delta), bad == 0)
        if al != 0:
            continue
        print(" controls (frame (0,0)):")
        btu = [bt[0], bt[1], bt[2], bt[3] * Ser(0, [1, 1], M + 8), bt[4]]
        box = (5, M - 14)    # for unequal products the inverse of X -> A X - X At has poles in u: truncation artefacts
        for delta in range(5):  # travel about 5 u-orders per z-order, so the control uses a box far from the truncation
            rows = equations(b, btu, K, M, delta)
            dproj, _ = projected_dim(rows, K, M, Kp, margin, box)
            rows_eq = equations(b, bt, K, M, delta)
            deq, _ = projected_dim(rows_eq, K, M, Kp, margin, box)
            pred = len([k for k in range(box[0] + 1) if k % 5 == delta])
            check("unequal products (bt_3 -> (1+u) bt_3), degree %d, box k <= %d, m <= %d: dimension 0 (equal products: %d = %d)" % (
                delta, box[0], box[1], deq, pred), dproj == 0 and deq == pred, "(found %d)" % dproj)
        for delta in (0, 1):
            rows = equations(b, bt, K, M, delta, drop_residue=True)
            dproj, _ = projected_dim(rows, K, M, Kp, margin)
            pred = len([k for k in range(Kp + 1) if k % 5 == delta])
            check("residue equations (u^-1) dropped, degree %d: the count changes (%d != %d)" % (delta, dproj, pred),
                  dproj != pred)
        for mg in (0, 1):
            rows = equations(b, bt, K, M, 1)
            dproj, _ = projected_dim(rows, K, M, Kp, mg)
            pred = len([k for k in range(Kp + 1) if k % 5 == 1])
            print("  [info] margin %d (window touching the truncation), degree 1: projected dimension %d (prediction %d)" % (
                mg, dproj, pred))
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
