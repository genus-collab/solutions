"""Referee re-derivation of the VERIFIED claims of ../DICT_PROOF.md with independent code (rf_core only).
Usage: py -I -B rf_check.py FIELD WMAX KLEG WRANK KCLOSED      FIELD in {Q, 61, 89} (Q exact, or mod 2^61-1 / 2^89-1)
  R0 counts (Lemma 0, Lemma 1(a),(d))                     w <= 4000 (FIELD Q only)
  R1 Lemma 1(c): L-membership (elimination test) == support rule j >= -m; monomials and random elements, w <= 100
  R2 Lemma 2 (i), (ii) closed formulas and I(M) in M on the monomial basis of M_w, w <= 100
  R3 legs_iter in M (k <= KLEG); legs_closed == legs_iter modulo Q[v][x] (k <= KCLOSED)
  R4 Lemma 3: mod-v images (closed forms) and the 2x2 determinants, k <= KLEG
  R5 the full leg map S_w -> N_w has rank |S_w| = dim N_w (generic elimination), w <= WRANK
  R6 K_0 := ker(log map) (generic elimination): dim d(w); legs log-free, in Q[p,q], in L; nf_vx linear; NF rank d(w), w <= WMAX
  R7 Cor. 4 (b=1 log map triangular, diagonal (-1)^k/k!) and Cor. 5 (dim K cap S_w = d(w) + [5|w]), w <= WMAX
"""
import os
import random
import sys
import time
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import (Ar, P61, P89, xpow, lin, mul, integ, I_step, legs_iter, legs_closed, wt, in_Qpq,  # noqa: E402
                     x_coords, nf_vx, in_L, in_M, triples, rank, left_kernel)

FIELD = sys.argv[1] if len(sys.argv) > 1 else "Q"
WMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 80
KLEG = int(sys.argv[3]) if len(sys.argv) > 3 else 40
WRANK = int(sys.argv[4]) if len(sys.argv) > 4 else 80
KCL = int(sys.argv[5]) if len(sys.argv) > 5 else 30
ar = Ar({"Q": None, "61": P61, "89": P89}[FIELD])
FAILS = []
T0 = time.time()


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg, flush=True)


def d_(w):
    return w // 5 + 1 - (1 if w % 5 == 0 else 0)


def e_(w):
    return sum(1 for n in range(w // 4 + 1) if (w - 4 * n) % 5 == 0) if w >= 0 else 0


def r0():
    for w in range(4001):
        S = triples(w) if w <= 1000 else None
        if S is not None:
            check(sum(1 for t in S if t[0] != 1) == d_(w) and sum(1 for t in S if t[0] == 1) == e_(w), "R0 S_w w=%d" % w)
        check(w // 4 + 1 == d_(w) + e_(w), "R0 floor(w/4)+1 = d+e, w=%d" % w)
        check(w // 4 + 1 == sum(e_(w - i) for i in range(5)), "R0 Q[p,q] free on 1..p^4 count, w=%d" % w)
        check(sum(e_(w - i) for i in range(1, 5)) == d_(w), "R0 Q[p,q]/Q[v,x] count = d(w), w=%d" % w)
        check(sum(1 for j in range(1, w + 1) if (w - j) % 5 == 0) == d_(w), "R0 Q[v][p]/Q[v] count, w=%d" % w)
    print("R0 counts w <= 4000 (S_w for w <= 1000): done", flush=True)


def r1():
    rng = random.Random(20261114)
    nin = nout = nrand = 0
    for w in range(101):
        for m in range(0, w // 5 + w // 4 + 5):
            j = w - 5 * m
            ok = in_L(ar, {j: ar.q(1)}, w)
            check(ok == (j >= -m) == in_Qpq({j: 1}, w), "R1 monomial v^%d p^%d" % (m, j))
            nin += ok
            nout += not ok
        keys = [j for j in range(-w - 10, w + 30) if (j - w) % 5 == 0]
        good = [j for j in keys if in_Qpq({j: 1}, w)]
        for trial in range(6):
            A = {j: ar.q(rng.randint(-9, 9)) for j in rng.sample(good, min(len(good), 4))}
            if trial % 2:
                A[rng.choice([j for j in keys if j not in good])] = ar.q(rng.choice([-2, -1, 1, 3]))
            A = ar.clean(A)
            check(in_L(ar, A, w) == in_Qpq(A, w), "R1 random element w=%d %s" % (w, sorted(A)))
            nrand += 1
    print("R1 Lemma 1(c): monomials %d in L (j >= -m), %d not in L (j < -m); %d random elements agree" % (nin, nout, nrand),
          flush=True)


def closed_I(a, c):
    out0, out1 = {}, {}
    for e, cf, key_c in ((a - c + 4, -1, 0), (a - c - 1, 1, 1)):
        if e == 0:
            out1[0] = out1.get(0, 0) + ar.q(cf)
        else:
            out0[e] = out0.get(e, 0) + ar.q(cf, e)
    return ar.clean(out0), ar.clean(out1)


def r2():
    nb = 0
    for w in range(101):
        for c in range(w // 4 + 1):
            a = w - 4 * c
            r = I_step(ar, {a - c: ar.q(1)}, {})
            check(r == closed_I(a, c), "R2 (i) I(p^%d q^%d)" % (a, c))
            check(in_M(ar, r[0], r[1], w + 4), "R2 I(p^%d q^%d) not in M" % (a, c))
            nb += 1
        for n in range(w // 4 + 1):
            if (w - 4 * n) % 5:
                continue
            r = I_step(ar, {}, xpow(ar, n))
            G = lin(ar, (ar.q(1, n + 1), xpow(ar, n + 1)))
            b0, b1 = integ(ar, mul(ar, G, {-1: ar.q(1)}), {})
            want = (b0, lin(ar, (-1, G), (1, b1)))
            check(r == want, "R2 (ii) I(log p v^m x^%d), w=%d" % (n, w))
            check(in_M(ar, r[0], r[1], w + 4), "R2 I(log p x^%d) not in M" % n)
            nb += 1
    print("R2 Lemma 2 (i), (ii) on the monomial basis of M_w, w <= 100: %d basis elements" % nb, flush=True)


def r3(T):
    for (b, k), (A0, A1) in T.items():
        w = wt(b, k)
        check(all((w - j) % 5 == 0 and j <= w for j in list(A0) + list(A1)), "R3 leg b=%d k=%d weight" % (b, k))
        check(in_M(ar, A0, A1, w), "R3 leg b=%d k=%d not in M" % (b, k))
    TC = legs_closed(ar, KCL)
    nd = 0
    for (b, k), (C0, C1) in TC.items():
        D0, D1 = lin(ar, (1, T[(b, k)][0]), (-1, C0)), lin(ar, (1, T[(b, k)][1]), (-1, C1))
        xc, ok = x_coords(ar, D0, wt(b, k))
        check(not D1 and ok, "R3 closed vs iter b=%d k=%d: difference not in Q[v][x]" % (b, k))
        nd += bool(D0)
    print("R3 legs_iter in M for k <= %d; legs_closed = legs_iter mod Q[v][x] for k <= %d (%d of %d differ by a "
          "non-zero element of Q[v][x])" % (KLEG, KCL, nd, len(TC)), flush=True)


def r4(T, NFP, LX):
    n4 = 0
    for k in range(KLEG + 1):
        for b in range(5):
            w = wt(b, k)
            rem = NFP[(b, k)]
            lv0 = LX[(b, k)].get(w // 4, 0) if w % 4 == 0 else 0
            if b in (0, 4, 3):
                j = {0: 1, 4: 2, 3: 3}[b]
                den = 1
                for i in range(1, k + 1):
                    den *= j + 4 * i
                check(rem.get(w, 0) == ar.q((-1) ** k, den) and lv0 == 0, "R4 b=%d k=%d" % (b, k))
            elif b == 2:
                den = 4
                for i in range(1, k + 1):
                    den *= 4 + 4 * i
                check(rem.get(w, 0) == ar.q((-1) ** (k + 1), den) and lv0 == 0, "R4 b=2 k=%d" % k)
                det = rem.get(w, 0) * LX[(1, k + 1)].get(k + 1, 0) if k + 1 <= KLEG else None
                check(det is None or ar.nz(det), "R4 2x2 determinant weight %d" % w)
            else:
                check(lv0 == ar.q((-1) ** k, factorial(k)), "R4 b=1 k=%d log part mod v" % k)
            n4 += 1
    print("R4 Lemma 3 mod-v images and 2x2 determinants: %d checks (k <= %d)" % (n4, KLEG), flush=True)


def rows_N(S, NFP, LX):
    return [dict([(("p", j), c) for j, c in NFP[(b, k)].items()] + [(("L", n), c) for n, c in LX[(b, k)].items()])
            for (b, k, m) in S]


def r5_r7(T, NFP, LX):
    totK, bad = 0, []
    for w in range(WMAX + 1):
        S = sorted(triples(w), key=lambda t: (t[0] != 1, t))
        if w <= WRANK:
            rk = rank(ar, rows_N(S, NFP, LX))
            check(rk == len(S) == d_(w) + e_(w), "R5 w=%d rank %d |S| %d" % (w, rk, len(S)))
        kern = left_kernel(ar, [LX[(b, k)] for (b, k, m) in S])
        ok = len(kern) == d_(w)
        nfrows = []
        for y in kern:
            A0 = lin(ar, *[(c, T[S[i][:2]][0]) for i, c in y.items()])
            A1 = lin(ar, *[(c, T[S[i][:2]][1]) for i, c in y.items()])
            rem, negv = nf_vx(ar, A0, w)
            lin_nf = lin(ar, *[(c, NFP[S[i][:2]]) for i, c in y.items()])
            ok = ok and not A1 and in_Qpq(A0, w) and in_L(ar, A0, w) and not negv and rem == lin_nf
            nfrows.append(rem)
        ok = ok and rank(ar, nfrows) == d_(w)
        b1 = [(k, LX[(b, k)]) for (b, k, m) in S if b == 1]
        ok = ok and all(max(r) == k and r[k] == ar.q((-1) ** k, factorial(k)) for k, r in b1)
        ok = ok and rank(ar, [r for k, r in b1]) == e_(w)
        rk5 = rank(ar, [{n: c for n, c in LX[(b, k)].items() if n >= 1} for (b, k, m) in S])
        ok = ok and len(S) - rk5 == d_(w) + (1 if w % 5 == 0 else 0)
        totK += len(kern)
        if not ok:
            bad.append(w)
        if w <= 12 or w % 100 == 0 or not ok:
            print(" w=%4d |S|=%3d dim K_0=%3d (d(w)=%3d) ok=%s  [%.0fs]" % (w, len(S), len(kern), d_(w), ok,
                                                                           time.time() - T0), flush=True)
    check(not bad, "R6/R7 failing weights %s" % bad[:10])
    print("R5 rank = |S_w| for w <= %d; R6/R7 w <= %d: %d K_0 basis vectors, failing weights %s" % (WRANK, WMAX, totK, bad),
          flush=True)


if __name__ == "__main__":
    print("rf_check.py FIELD=%s WMAX=%d KLEG=%d WRANK=%d KCLOSED=%d" % (FIELD, WMAX, KLEG, WRANK, KCL), flush=True)
    if FIELD == "Q":
        r0()
    r1()
    r2()
    K = max(KLEG, WMAX // 4 + 1)
    T = legs_iter(ar, K)
    print("legs computed k <= %d [%.0fs]" % (K, time.time() - T0), flush=True)
    NFP, LX = {}, {}
    for (b, k), (A0, A1) in T.items():
        rem, negv = nf_vx(ar, A0, wt(b, k))
        check(not negv, "nf_vx needed Q(v)[x] for leg b=%d k=%d" % (b, k))
        NFP[(b, k)] = rem
        LX[(b, k)] = x_coords(ar, A1, wt(b, k))[0]
    r3({key: val for key, val in T.items() if key[1] <= KLEG})
    r4(T, NFP, LX)
    r5_r7(T, NFP, LX)
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS), "[%.0fs]" % (time.time() - T0))
