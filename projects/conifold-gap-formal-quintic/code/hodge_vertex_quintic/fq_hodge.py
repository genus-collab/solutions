# -*- coding: utf-8 -*-
"""fq_hodge.py -- Hodge integrals and QRR constants for the formal quintic.

Adapted (with light edits) from
  resolvedconifold/hodge_vertex_conifold/cy_torus_explore/exact_test/cylib.py
(functions wk, _wk, hi, _hi, partitions_odd, lam_from_ch, int_lam, class_s, VertexTable, h_coeffs,
series_exp), which is itself fresh code written from the Mumford/Faber-Pandharipande formulas.  Nothing is
imported from the old project; this file imports only fq_fp.

API (all exact; wk/hi return Fractions, VertexTable/hv work in an FpRing R):
  wk(g, d)            Witten-Kontsevich <tau_d1 ... tau_dn>_g (string/dilaton/DVV).
  hi(g, d, K)         Int_{Mbar_{g,n}} prod_{k in K} ch_k(E) prod_i psi_i^{d_i}, Mumford GRR recursion
                      (ch_even = 0 used through the cutoffs; K = tuple of odd k).
  int_lam(g, mono, d) Int prod_{i in mono} lambda_i prod psi^d via ch.
  class_s(R, w)       k -> s_k = -(k-1)! Sum_a w_a^(-k)  (mod p): the Hodge class of a weight multiset w is
                      Hodge_w := prod_a Lambda^vee(w_a)/w_a^g = exp(Sum_{k odd} s_k ch_k(E))   (FRAMEB A3.2).
  VertexTable(R, w)   .val(g, d): Int_{Mbar_{g,n}} Hodge_w prod psi^d  (mod p); w = None gives the point class.
  hv(R, w, zmax, sign=+1)  coefficients in z of h_w(z)^sign, h_w(z) = prod_a Theta(-z/a)
                      = exp(-Sum_k B_{2k}/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k))   (FRAMEB 0.1).
  wt(R, lams, i)      the 5-weight multiset wt_i = {lambda_i - lambda_j (j != i)} U {-5 lambda_i}
                      (numerator twist enters with a MINUS sign; loc_check.py conventions).
  W0(R, lams, i)      W_i^0 = prod_{j != i}(lambda_i - lambda_j)/(5 lambda_i)  (= zeta^(3i) at the zeta point).
"""
import itertools, math
from fractions import Fraction as Fr
from functools import lru_cache
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_fp import bernoulli


def df(n):
    r = 1
    while n > 0:
        r *= n; n -= 2
    return r


@lru_cache(maxsize=None)
def wk(g, d):
    d = tuple(sorted(d)); return _wk(g, d)


@lru_cache(maxsize=None)
def _wk(g, d):
    n = len(d)
    if g < 0 or n == 0 or d[0] < 0 or sum(d) != 3 * g - 3 + n or 2 * g - 2 + n <= 0: return Fr(0)
    if g == 0 and d == (0, 0, 0): return Fr(1)
    if g == 1 and d == (1,): return Fr(1, 24)
    if d[0] == 0:
        rest = d[1:]; tot = Fr(0)
        for j in range(len(rest)):
            if rest[j] > 0:
                r = list(rest); r[j] -= 1
                tot += wk(g, tuple(r))
        return tot
    k = d[0] - 1; rest = d[1:]; tot = Fr(0)
    for j in range(len(rest)):
        r = list(rest); r[j] += k
        tot += Fr(df(2 * k + 2 * rest[j] + 1), df(2 * rest[j] - 1)) * wk(g, tuple(r))
    for r_ in range(k):
        c = df(2 * r_ + 1) * df(2 * (k - 1 - r_) + 1)
        tot += Fr(c, 2) * wk(g - 1, (r_, k - 1 - r_) + rest)
        s = 0
        idx = range(len(rest))
        for m in range(len(rest) + 1):
            for S in itertools.combinations(idx, m):
                A = tuple(rest[i] for i in S); B = tuple(rest[i] for i in idx if i not in S)
                for g1 in range(g + 1):
                    x = wk(g1, (r_,) + A)
                    if x: s += x * wk(g - g1, (k - 1 - r_,) + B)
        tot += Fr(c, 2) * s
    return tot / df(2 * k + 3)


def dmax(g):
    if g <= 0: return 0
    if g == 1: return 1
    return min(g * (g + 1) // 2, 3 * g - 3)


def hi(g, d, K):
    return _hi(g, tuple(sorted(d)), tuple(sorted(K)))


@lru_cache(maxsize=None)
def _hi(g, d, K):
    """Int_{Mbar_{g,n}} prod ch_k (k in K) prod psi^d via Mumford's GRR formula (Faber-Pandharipande form)"""
    n = len(d)
    if g < 0 or 2 * g - 2 + n <= 0: return Fr(0)
    if sum(d) + sum(K) != 3 * g - 3 + n: return Fr(0)
    if not K: return wk(g, d) if n else Fr(0)
    if g == 0: return Fr(0)
    if sum(K) > dmax(g): return Fr(0)
    if K[-1] >= 2 * g: return Fr(0)
    k = K[-1]; Kr = K[:-1]; l = (k + 1) // 2
    tot = hi(g, d + (k + 1,), Kr)
    for i in range(n):
        dd = list(d); dd[i] += k
        tot -= hi(g, dd, Kr)
    b = Fr(0)
    for a in range(k):
        sgn = -1 if a % 2 else 1
        bb = hi(g - 1, d + (a, k - 1 - a), Kr)
        for m in range(n + 1):
            for S in itertools.combinations(range(n), m):
                A = tuple(d[i] for i in S); B = tuple(d[i] for i in range(n) if i not in S)
                for h in range(g + 1):
                    for mk in range(len(Kr) + 1):
                        for T in itertools.combinations(range(len(Kr)), mk):
                            K1 = tuple(Kr[i] for i in T); K2 = tuple(Kr[i] for i in range(len(Kr)) if i not in T)
                            x = hi(h, A + (a,), K1)
                            if x: bb += x * hi(g - h, B + (k - 1 - a,), K2)
        b += sgn * bb
    tot += b / 2
    return bernoulli(2 * l) / math.factorial(2 * l) * tot


def partitions_odd(D, maxpart):
    """multisets of odd integers <= maxpart with sum D"""
    def rec(D, mp):
        if D == 0:
            yield ()
            return
        for k in range(min(mp, D), 0, -1):
            if k % 2 == 1:
                for r in rec(D - k, k):
                    yield (k,) + r
    return list(rec(D, maxpart))


def lam_from_ch(g, nmax):
    """lambda_i as polynomials in ch_odd: c_t(E) = exp(sum (-1)^(k-1) (k-1)! ch_k t^k); dict i -> {K: Fraction}"""
    res = {0: {(): Fr(1)}}
    for i in range(1, nmax + 1):
        cur = {}
        for k in range(1, i + 1, 2):
            c = Fr(math.factorial(k), i)
            for K, v in res[i - k].items():
                K2 = tuple(sorted(K + (k,)))
                cur[K2] = cur.get(K2, Fr(0)) + c * v
        res[i] = cur
    return res


def int_lam(g, mono, d):
    """Int prod_{i in mono} lambda_i prod psi^d, via ch"""
    L = lam_from_ch(g, max(mono) if mono else 0)
    poly = {(): Fr(1)}
    for i in mono:
        new = {}
        for K, v in poly.items():
            for K2, w in L[i].items():
                KK = tuple(sorted(K + K2))
                new[KK] = new.get(KK, Fr(0)) + v * w
        poly = new
    return sum((v * hi(g, d, K) for K, v in poly.items()), Fr(0))


# ------------------------------------------------------------------ Hodge vertex class and QRR constants
def class_s(R, w):
    """k -> s_k = -(k-1)! Sum_a w_a^(-k) mod p, with Hodge_w = exp(Sum_{k odd} s_k ch_k).  w = None: point class."""
    if w is None: return lambda k: R.zero
    winv = [R.inv(R.of(x)) for x in w]
    def s(k):
        tot = R.zero
        for x in winv: tot = R.add(tot, R.pow(x, k))
        return R.mul(R.of(-math.factorial(k - 1)), tot)
    return s


class VertexTable(object):
    """Int_{Mbar_{g,n}} Hodge_w prod psi^d  (mod p), cached.  w = None gives the point class (WK numbers)."""
    def __init__(self, R, w):
        self.R = R
        self.w = None if w is None else [R.of(x) for x in w]
        self.s = class_s(R, self.w)
        self.cache = {}; self.sk = {}
    def sval(self, k):
        r = self.sk.get(k)
        if r is None:
            r = self.s(k); self.sk[k] = r
        return r
    def val(self, g, d):
        R = self.R
        d = tuple(sorted(d)); key = (g, d)
        r = self.cache.get(key)
        if r is not None: return r
        n = len(d); D = 3 * g - 3 + n - sum(d)
        if D < 0 or 2 * g - 2 + n <= 0: r = R.zero
        elif self.w is None: r = R.of(wk(g, d)) if D == 0 and n else R.zero
        elif D > dmax(g): r = R.zero
        else:
            r = R.zero
            for K in partitions_odd(D, max(2 * g - 1, 1)):
                x = hi(g, d, K)
                if x == 0: continue
                c = R.of(x)
                cnt = {}
                for k in K: cnt[k] = cnt.get(k, 0) + 1
                for k, m in cnt.items():
                    c = R.mul(c, R.mul(R.pow(self.sval(k), m), R.inv(R.of(math.factorial(m)))))
                r = R.add(r, c)
        self.cache[key] = r
        return r


def series_exp(R, ex, zmax):
    out = [R.zero] * (zmax + 1); out[0] = R.one
    for n in range(1, zmax + 1):
        t = R.zero
        for k in range(1, n + 1):
            if not R.iszero(ex[k]): t = R.add(t, R.mul(R.mul(R.of(k), ex[k]), out[n - k]))
        out[n] = R.mul(t, R.inv(R.of(n)))
    return out


def hv(R, w, zmax, sign=1):
    """coefficients [z^0..z^zmax] of h_w(z)^sign, h_w(z) = exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k)).
    w = None (or empty) gives 1."""
    out = [R.zero] * (zmax + 1); out[0] = R.one
    if not w: return out
    ex = [R.zero] * (zmax + 1)
    winv = [R.inv(R.of(x)) for x in w]
    k = 1
    while 2 * k - 1 <= zmax:
        c = R.of(bernoulli(2 * k) / (2 * k * (2 * k - 1)))
        tot = R.zero
        for x in winv: tot = R.add(tot, R.pow(x, 2 * k - 1))
        ex[2 * k - 1] = R.mul(R.of(-sign), R.mul(c, tot))
        k += 1
    return series_exp(R, ex, zmax)


def wt(R, lams, i):
    """wt_i = {lambda_i - lambda_j (j != i)} U {-5 lambda_i}"""
    return [R.sub(lams[i], lams[j]) for j in range(len(lams)) if j != i] + [R.mul(R.of(-5), lams[i])]


def W0(R, lams, i):
    """W_i^0 = prod_{j != i}(lambda_i - lambda_j)/(5 lambda_i)"""
    r = R.one
    for j in range(len(lams)):
        if j != i: r = R.mul(r, R.sub(lams[i], lams[j]))
    return R.div(r, R.mul(R.of(5), lams[i]))


if __name__ == "__main__":
    from fq_fp import FpRing
    assert wk(1, (1,)) == Fr(1, 24) and wk(2, (4,)) == Fr(1, 1152) and wk(2, (2, 3)) == Fr(29, 5760)
    assert hi(1, (0,), (1,)) == Fr(1, 24)                 # Int_{M_{1,1}} ch_1 = lambda_1 = 1/24
    assert int_lam(2, (1, 1, 1), ()) == Fr(1, 2880)       # Int_{M_2} lambda_1^3 = 1/2880
    assert int_lam(2, (2, 1), ()) == Fr(1, 5760)          # lambda_2 lambda_1 = 1/5760
    assert int_lam(3, (3, 2, 1), ()) == Fr(1, 1451520)    # lambda_3 lambda_2 lambda_1 = 1/1451520
    R = FpRing()
    lams = [R.z(i) for i in range(5)]
    for i in range(5):
        assert W0(R, lams, i) == R.z(3 * i)
    # constant-map genus-1 check: Int_{M_{1,1}} Hodge_w psi^0 = -(1/24) Sum 1/w  ... via exp(s_1 ch_1): s_1 = -Sum 1/w_a
    vt = VertexTable(R, wt(R, lams, 0))
    s1 = R.zero
    for a in wt(R, lams, 0): s1 = R.sub(s1, R.inv(a))
    assert vt.val(1, (0,)) == R.mul(s1, R.of(Fr(1, 24)))
    h0 = hv(R, wt(R, lams, 0), 7)
    print("h_wt0 coefficients:", [R.fmt(x) for x in h0])
    h1 = hv(R, wt(R, lams, 1), 7)
    assert all(R.eq(h1[n], R.mul(R.z(-n), h0[n])) for n in range(8)), "h_wt_i(z) = h_wt_0(zeta^-i z)"
    hm = hv(R, wt(R, lams, 0), 7, -1)
    prod = [sum(h0[a] * hm[n - a] for a in range(n + 1)) % R.p for n in range(8)]
    assert prod == [1] + [0] * 7
    print("fq_hodge ok")
