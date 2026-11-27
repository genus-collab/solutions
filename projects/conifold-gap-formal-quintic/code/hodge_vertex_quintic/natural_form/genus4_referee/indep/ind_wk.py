# -*- coding: utf-8 -*-
"""ind_wk.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): Witten-Kontsevich intersection numbers.

wk(g, d) = <tau_d1 ... tau_dn>_g as a Fraction (0 outside the stable/dimension range).
Method: the Virasoro (DVV) recursion applied to the LARGEST index, for every k >= -1 (so the string and dilaton equations are
the cases k = -1, 0 of the same formula and are never coded separately):
   (2k+3)!! <tau_(k+1) tau_d>_g = Sum_j (2k+2d_j+1)!!/(2d_j-1)!! <tau_(d_j+k) tau_(d minus j)>_g
        + (1/2) Sum_(a+b=k-1) (2a+1)!! (2b+1)!! [ <tau_a tau_b tau_d>_(g-1) + Sum_(g1+g2=g, I u J = d) <tau_a tau_I>_g1 <tau_b tau_J>_g2 ],
with the two initial values <tau_0^3>_0 = 1, <tau_1>_1 = 1/24.  (fq_hodge.wk of the project removes the SMALLEST index and uses
the string equation separately; agreement of the two is a test of both.)
Self-test (py -I ind_wk.py): closed formulas (genus 0 multinomials, <tau_1^n>_1 = (n-1)!/24, <tau_(3g-2)>_g = 1/(24^g g!),
<tau_10>_4 = 1/7962624), classical genus-2 and genus-3 values, and a comparison with fq_hodge.wk on all (g, d) with g <= 4, n <= 4.
"""
import itertools, math, os, sys
from fractions import Fraction as Fr
from functools import lru_cache


def dfact(n):
    """n!! for odd n >= -1"""
    r = 1
    while n > 1: r *= n; n -= 2
    return r


def wk(g, d):
    return _wk(g, tuple(sorted(d)))


@lru_cache(maxsize=None)
def _wk(g, d):
    n = len(d)
    if g < 0 or n == 0 or 2 * g - 2 + n <= 0 or d[0] < 0 or sum(d) != 3 * g - 3 + n: return Fr(0)
    if g == 0 and n == 3: return Fr(1)
    if g == 1 and n == 1: return Fr(1, 24)
    k = d[-1] - 1; rest = d[:-1]; m = len(rest); tot = Fr(0)
    for j in range(m):
        if rest[j] + k < 0: continue
        r = rest[:j] + (rest[j] + k,) + rest[j + 1:]
        tot += Fr(dfact(2 * k + 2 * rest[j] + 1), dfact(2 * rest[j] - 1)) * wk(g, r)
    for a in range(0, k):
        b = k - 1 - a; c = Fr(dfact(2 * a + 1) * dfact(2 * b + 1), 2)
        acc = wk(g - 1, (a, b) + rest)
        for mask in range(1 << m):
            I = tuple(rest[i] for i in range(m) if (mask >> i) & 1); J = tuple(rest[i] for i in range(m) if not (mask >> i) & 1)
            # genus g1 of the first factor is fixed by the dimension: a + |I| = 3 g1 - 3 + (1 + #I)
            num = a + sum(I) + 2 - len(I)
            if num % 3: continue
            g1 = num // 3
            if g1 < 0 or g1 > g: continue
            x = wk(g1, (a,) + I)
            if x: acc += x * wk(g - g1, (b,) + J)
        tot += c * acc
    return tot / dfact(2 * k + 3)


def selftest(verbose=True):
    # closed formulas
    for n in range(3, 9):
        for d in itertools.combinations_with_replacement(range(n - 2), n):
            if sum(d) != n - 3: continue
            assert wk(0, d) == Fr(math.factorial(n - 3), math.prod(math.factorial(x) for x in d)), d
    for n in range(1, 8): assert wk(1, (1,) * n) == Fr(math.factorial(n - 1), 24)
    for g in range(1, 7): assert wk(g, (3 * g - 2,)) == Fr(1, 24 ** g * math.factorial(g)), g
    assert wk(4, (10,)) == Fr(1, 7962624)
    # classical values (Witten's and Faber's tables)
    known = {(2, (4,)): Fr(1, 1152), (2, (2, 3)): Fr(29, 5760), (2, (1, 4)): Fr(1, 384), (2, (2, 2, 2)): Fr(7, 240),
             (2, (1, 2, 3)): Fr(29, 1440), (2, (1, 1, 4)): Fr(1, 96), (2, (0, 5)): Fr(1, 1152),
             (3, (7,)): Fr(1, 82944), (3, (1, 7)): Fr(5, 82944), (3, (2, 6)): Fr(77, 414720), (3, (3, 5)): Fr(503, 1451520),
             (3, (4, 4)): Fr(607, 1451520), (3, (2, 2, 2, 2, 2, 2)): Fr(1225, 144)}
    for (g, d), v in known.items(): assert wk(g, d) == v, (g, d, wk(g, d), v)
    # comparison with the project's fq_hodge.wk (convention check only)
    here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(here, "..", "..", "..")))
    from fq_hodge import wk as wk_project
    cnt = 0
    for g in range(0, 5):
        for n in range(1, 5):
            if 2 * g - 2 + n <= 0: continue
            dim = 3 * g - 3 + n
            for d in itertools.combinations_with_replacement(range(dim + 1), n):
                if sum(d) != dim: continue
                assert wk(g, d) == wk_project(g, d), (g, d); cnt += 1
    for d in [(2, 2, 2, 2, 2, 2, 2, 2, 2), (2, 3, 4, 5), (0, 0, 2, 3, 3, 5), (3, 3, 3, 3), (2, 2, 3, 3, 3)]:
        g4 = (sum(d) - len(d) + 3)
        if g4 % 3 == 0: assert wk(g4 // 3, d) == wk_project(g4 // 3, d), d; cnt += 1
    if verbose: print("ind_wk ok: closed formulas, %d classical values, %d comparisons with fq_hodge.wk; <tau_10>_4 = %s" % (len(known), cnt, wk(4, (10,))))
    return True


if __name__ == "__main__":
    selftest()
