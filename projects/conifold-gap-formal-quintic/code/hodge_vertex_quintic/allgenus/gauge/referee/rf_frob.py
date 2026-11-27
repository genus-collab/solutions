# -*- coding: utf-8 -*-
"""rf_frob.py -- REFEREE: the regular gauge by a route DIFFERENT from Theorem U's z-recursion (ur_gauge) and from gi_*:
the u-recursion of the gauge equation with an UNKNOWN initial value C(z) = P(0, z), and linear algebra over Q.

z P' = A P - P At,  A = N/u + Sum_j A_j u^j,  At = N~/u + ...,  N = r E_12, N~ = r~ E_12.  Order u^(n-1):
    (n z - L) P^(n) = Sum_{j<n} (A_j P^(n-1-j) - P^(n-1-j) At_j),   L X = N X - X N~,  L^3 = 0,
so P^(n) = Sum_{k=0..2} L^k(rhs) / (n z)^(k+1) is determined by C (Lemma-1 type uniqueness).  C must commute with the residues
(order u^-1: (R1)).  REGULARITY (no negative power of z in any P^(n)) is a linear system for the coefficients of C.  Theorem U predicts
that its solution space is exactly C -> C c(z^5).  Unknowns: graded slots (a, b, k), a = b + k mod 5, k <= KC, with (R1) imposed.
A coefficient [z^e] P^(n)_ab is used only if COMPLETE: e + reach_n(a, b) <= KC, reach = largest downward z-shift of any path from C.
"""
from fractions import Fraction as Fr


def slots(KC):
    idx = {}; lst = []
    for k in range(KC + 1):
        for b in range(5):
            a = (b + k) % 5
            if b == 1 and a != 1: continue
            if a == 2 and b != 2: continue
            if (a, b) == (2, 2): continue
            idx[(a, b, k)] = len(lst); lst.append((a, b, k))
    for k in range(0, KC + 1, 5): idx[(2, 2, k)] = idx[(1, 1, k)]
    return idx, lst


def fadd(acc, form, c):
    for i, x in form.items():
        y = acc.get(i, 0) + c * x
        if y: acc[i] = y
        else: acc.pop(i, None)


def matrices(b, NU):
    """A_j (j >= 0) as {(a, a+1): coefficient} and the residue r, from b = [b0, u b1, b2, b3, b4]"""
    A = []
    for j in range(NU - 1):
        A.append({(a, (a + 1) % 5): (b[a][j] if a != 1 else b[1][j + 1]) for a in range(5)})
    return A, b[1][0]


def recursion(b, bt, KC, NMAX):
    idx, lst = slots(KC)
    A, r = matrices(b, NMAX + 2); At, rt = matrices(bt, NMAX + 2)
    P = [{}]
    for (a, bb, k), i in idx.items(): P[0][(a, bb, k)] = {i: Fr(1)}
    reach = [{(a, c): 0 for a in range(5) for c in range(5)}]
    for n in range(1, NMAX + 1):
        R = {}; rR = {(a, c): -10 ** 6 for a in range(5) for c in range(5)}
        for j in range(n):
            X = P[n - 1 - j]; rX = reach[n - 1 - j]
            for (a, c, e), form in X.items():
                ap = (a - 1) % 5                                    # (A_j X)[ap][c] = A_j[ap][a] X[a][c]
                fadd(R.setdefault((ap, c, e), {}), form, A[j][(ap, a)])
                cp = (c + 1) % 5                                    # (X At_j)[a][cp] = X[a][c] At_j[c][cp]
                fadd(R.setdefault((a, cp, e), {}), form, -At[j][(c, cp)])
            for a in range(5):
                for c in range(5):
                    rR[(a, c)] = max(rR[(a, c)], rX[((a + 1) % 5, c)], rX[(a, (c - 1) % 5)])
        new = {}; rnew = {(a, c): -10 ** 6 for a in range(5) for c in range(5)}
        Lk = R; rL = dict(rR)
        for k in range(3):
            f = Fr(1, n ** (k + 1))
            for (a, c, e), form in Lk.items():
                if form: fadd(new.setdefault((a, c, e - k - 1), {}), form, f)
            for key in rL: rnew[key] = max(rnew[key], rL[key] + k + 1)
            # L X = r E12 X - rt X E12
            nxt = {}; rn = {(a, c): -10 ** 6 for a in range(5) for c in range(5)}
            for (a, c, e), form in Lk.items():
                if a == 2: fadd(nxt.setdefault((1, c, e), {}), form, r)
                if c == 1: fadd(nxt.setdefault((a, 2, e), {}), form, -rt)
            for c in range(5): rn[(1, c)] = max(rn[(1, c)], rL[(2, c)])
            for a in range(5): rn[(a, 2)] = max(rn[(a, 2)], rL[(a, 1)])
            Lk = nxt; rL = rn
        assert not any(Lk.values()), "L^3 != 0"
        P.append({key: f for key, f in new.items() if f}); reach.append(rnew)
    return idx, lst, P, reach


def solve(rows, ncols, pivot_order):
    """rows: list of (form dict, rhs).  Echelon form with pivots chosen by pivot_order (low k first); returns (pivots, inconsistent)"""
    piv = {}; rank_order = []; bad = 0; born = {}
    pos = {c: i for i, c in enumerate(pivot_order)}
    for form, rhs in rows:
        f = dict(form); s = Fr(rhs)
        while True:
            cols = [c for c in f if c in piv]
            if not cols: break
            c = min(cols, key=lambda x: born[x])                 # earliest pivot first: terminates
            pf, ps = piv[c]; m = f[c]
            for i, x in pf.items():
                y = f.get(i, 0) - m * x
                if y: f[i] = y
                else: f.pop(i, None)
            s -= m * ps
        if not f:
            if s != 0: bad += 1
            continue
        c = min(f, key=lambda x: pos[x]); m = f[c]
        piv[c] = ({i: x / m for i, x in f.items()}, s / m); born[c] = len(rank_order); rank_order.append(c)
    # back substitution: express each pivot through free columns
    val = {}
    for c in reversed(rank_order):
        pf, ps = piv[c]; v_const = ps; v_free = {}
        for i, x in pf.items():
            if i == c: continue
            if i in val:
                vc, vf = val[i]; v_const -= x * vc
                for j, y in vf.items(): v_free[j] = v_free.get(j, 0) - x * y
            else:
                v_free[i] = v_free.get(i, 0) - x
        val[c] = (v_const, {j: y for j, y in v_free.items() if y})
    return val, bad
