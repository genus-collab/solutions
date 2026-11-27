"""Stage 2, D(iii): cross-checks for DICT_PROOF.md.
Usage: py -I -B dict_xcheck.py [KMAX]
 X1  my legs (dict_core.legs_A) = the critic's legs (critic/cr_dictionary.py leg recursion via critic/cr_laurent.py,
     imported READ-ONLY), k <= KMAX, all b: identical exact dictionaries.
 X2  the closed mod-v images used in Lemma 3 of DICT_PROOF.md, k <= KMAX:
       b = 0, 4, 3 (Phi_0 = p^j, j = 1, 2, 3):  lambda(tau_k) = (-1)^k p^(j+4k) / prod_{i=1..k} (j + 4i)       mod (v, Q[v,x])
       b = 2 (Phi_0 = q = v/p):                  lambda(tau_k) = (-1)^(k+1) p^(4k+4) / (4 prod_{i=1..k} (4+4i))
       b = 1 (Phi_0 = log p):                    log-part of lambda(tau_k) = (-1)^k x^k / k!  (v = 0)
 X4  log part of the leg of tau_k(H) has x-degree exactly k, leading coefficient (-1)^k/k! (Corollary 4), k <= KMAX.
 X3  Q[p,q] = Q[v][p] + Q[v][x] (Lemma 1(a)): each monomial p^a q^c, a + 4c <= 40, reduced to Q[v][p] by subtracting Q[v][x].
 C4  control: over Q(v) the statement fails: tau_0(H^2)/v has leg 1/p, not in L.
"""
import os
import sys
from fractions import Fraction as Fr
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CRITIC = os.path.normpath(os.path.join(HERE, "..", "..", "critic"))
sys.path.insert(1, CRITIC)
from dict_core import legs_A, split, xpow, mul, add  # noqa: E402
from dict_nf import NF, L_test  # noqa: E402
import cr_laurent as CL  # noqa: E402  (read-only import of Stage-1 critic code)

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg)


def critic_legs(kmax):
    X = {(4, 0, 0): Fr(1, 4), (-1, 0, 1): Fr(1)}
    XP = CL.deriv(X)
    LEG0 = {0: {(1, 0, 0): Fr(1)}, 1: {(0, 1, 0): Fr(1)}, 2: {(-1, 0, 1): Fr(1)}, 3: {(3, 0, 0): Fr(1)},
            4: {(2, 0, 0): Fr(1)}}
    T = {}
    for b in range(5):
        cur = LEG0[b]
        for k in range(kmax + 1):
            T[(b, k)] = {(j, m, l): c for (j, l, m), c in cur.items()}   # critic keys (j, log, v) -> mine (j, v, log)
            cur = CL.scal(CL.integ(CL.mul(cur, XP)), -1)
    return T


def at_v0(nf):
    return {k: c for k, c in nf.items() if (k[2] if k[0] == 'p' else k[1]) == 0}


def main():
    print("dict_xcheck.py KMAX=%d" % KMAX)
    T = legs_A(KMAX)
    TC = critic_legs(KMAX)
    nd = sum(1 for key in T if T[key] != TC[key])
    check(nd == 0, "X1: %d legs differ from the critic's" % nd)
    print("X1 legs equal to the critic's recursion: %d/%d (b = 0..4, k <= %d)" % (len(T) - nd, len(T), KMAX))
    n2 = 0
    for k in range(KMAX + 1):
        for b, j in ((0, 1), (4, 2), (3, 3)):
            den = 1
            for i in range(1, k + 1):
                den *= j + 4 * i
            pred = {(j + 4 * k, 0, 0): Fr((-1) ** k, den)}
            check(at_v0(NF(T[(b, k)])) == at_v0(NF(pred)), "X2 b=%d k=%d" % (b, k))
            n2 += 1
        den = 4
        for i in range(1, k + 1):
            den *= 4 + 4 * i
        pred = {(4 * k + 4, 0, 0): Fr((-1) ** (k + 1), den)}
        check(at_v0(NF(T[(2, k)])) == at_v0(NF(pred)), "X2 b=2 k=%d" % k)
        logpart = {key: c for key, c in at_v0(NF(T[(1, k)])).items() if key[0] == 'L'}
        check(logpart == {('L', 0, k): Fr((-1) ** k, factorial(k))}, "X2 b=1 k=%d log part" % k)
        n2 += 2
    print("X2 closed mod-v images: %d checks" % n2)
    n3 = 0
    for w in range(41):
        for c in range(w // 4 + 1):
            a = w - 4 * c
            G = {(a - c, c, 0): Fr(1)}
            sub = 0
            while any(j < 0 for (j, m, l) in G):
                jb = min(j for (j, m, l) in G)
                (j, m, l), cf = next((kk, cc) for kk, cc in G.items() if kk[0] == jb)
                check(m + jb >= 0, "X3 negative v-power needed for p^%d q^%d" % (a, c))
                G = add(G, mul({(0, m + jb, 0): cf}, xpow(-jb)), -1)
                sub += 1
            check(all(m >= 0 and j >= 0 for (j, m, l) in G), "X3 p^%d q^%d remainder" % (a, c))
            n3 += 1
    print("X3 monomials p^a q^c (a + 4c <= 40) reduced into Q[v][p] by subtracting Q[v][x]: %d" % n3)
    n4 = 0
    for k in range(KMAX + 1):          # X4: log part of Phi_k^1 = (-1)^k x^k/k! + (lower x-degree), at all v
        lp = {key: c for key, c in NF(T[(1, k)]).items() if key[0] == 'L'}
        top = max(n for (_, m, n) in lp)
        check(top == k and lp.get(('L', 0, k)) == Fr((-1) ** k, factorial(k))
              and all(m == 0 for (_, m, n) in lp if n == k), "X4 k=%d" % k)
        n4 += 1
    print("X4 log part of Phi_k^1 has x-degree k, leading coefficient (-1)^k/k!: %d checks" % n4)
    ok = L_test({(-1, 0): Fr(1)}, -1)
    print("C4 leg of tau_0(H^2)/v = 1/p in L:", ok, "(expected False)")
    check(not ok, "C4 control")
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS))


if __name__ == "__main__":
    main()
