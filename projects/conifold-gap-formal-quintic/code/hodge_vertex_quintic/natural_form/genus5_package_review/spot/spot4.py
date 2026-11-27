# -*- coding: utf-8 -*-
"""spot4.py -- consistency of the shipped genus-5 point-frame results at the two primes (results/g5_89_id.json, g5_61_id.json, copies in pkg/).

(a) totals [u^0..u^7]: rational reconstruction at 2^89-1 alone (|num|, den <= sqrt(p/2)), reduced mod 2^61-1 and compared with the residue
    stored for 2^61-1 (and vice versa where 2^61-1 alone reconstructs); then CRT to p89*p61 and reconstruction there.
(b) the gap form checked with an INDEPENDENT t_c (referee module ref_rows.periods, not fq_hat.tc_series):
    [u^0..u^7] total == (1/660000) [u^0..u^7] (u/t_c)^8  at both primes.
(c) per outer graph (= inner graph in the point frame): for every one of the 4555 x 8 entries, does the 2^89-1 reconstruction reduce to the
    2^61-1 residue?  Entries are counted as consistent / not reconstructible at 2^89-1 alone / INCONSISTENT; the CRT reconstruction is then
    tried for the second class with a smoothness test of the denominator.
"""
import os, sys, json
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.join(HERE, "pkg"); OUT = os.path.join(HERE, "out")
REF = os.path.normpath(os.path.join(HERE, "..", "..", "referee"))
P89 = (1 << 89) - 1; P61 = (1 << 61) - 1


def ratrec(a, m, B=None):
    """x/y with x = a y mod m, |x|, y <= B (default sqrt(m/2)), or None"""
    a %= m; B = isqrt(m // 2) if B is None else B
    r0, r1 = m, a; s0, s1 = 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    f = Fr(r1, s1)
    if (f.numerator - a * f.denominator) % m: return None
    return f


def red(f, p): return f.numerator * pow(f.denominator, -1, p) % p


def crt(a, b): return (a + P89 * ((b - a) * pow(P89, -1, P61) % P61)) % (P89 * P61)


def smooth(n, bound=1000):
    for q in range(2, bound):
        while n % q == 0: n //= q
    return n == 1


def tc_series(p, K):
    sys.path.insert(0, REF)
    from ref_rows import periods
    from ref_fp import smul, sinv
    f0, f1, f2, g, h = periods(p, K + 3)
    return smul(g[:K + 2], sinv(f0[:K + 2], p), p)


def main():
    log_f = open(os.path.join(OUT, "spot4.log"), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); log_f.write(s + "\n")
    log("=" * 100)
    d89 = json.load(open(os.path.join(PKG, "results", "g5_89_id.json"))); d61 = json.load(open(os.path.join(PKG, "results", "g5_61_id.json")))
    assert d89['outer'] == d61['outer'], "outer lists differ between the primes"
    log("outer lists identical at both primes: %d graphs" % len(d89['outer']))
    # (a)
    t89 = d89['info']['total']; t61 = d61['info']['total']
    log("(a) totals")
    for k in range(8):
        r89 = ratrec(t89[k], P89); r61 = ratrec(t61[k], P61); c = ratrec(crt(t89[k], t61[k]), P89 * P61)
        s = "   u^%d: 2^89-1 alone -> %s%s; 2^61-1 alone -> %s%s; CRT -> %s (den smooth: %s)" % (
            k, r89, "" if r89 is None else (" (reduces to the 2^61-1 residue: %s)" % (red(r89, P61) == t61[k] % P61)),
            r61, "" if r61 is None else (" (reduces to the 2^89-1 residue: %s)" % (red(r61, P89) == t89[k] % P89)),
            c, None if c is None else smooth(c.denominator))
        log(s)
    # (b)
    log("(b) gap form with an independent t_c (referee ref_rows.periods)")
    for p, tot in ((P89, t89), (P61, t61)):
        tc = tc_series(p, 8)
        sys.path.insert(0, REF)
        from ref_fp import smul, sinv
        assert tc[0] % p == 0 and tc[1] % p == 1, "t_c = u + O(u^2) expected"
        r = sinv([tc[k + 1] for k in range(8)], p); w = [1] + [0] * 7
        for _ in range(8): w = smul(w, r, p, 8)
        c = red(Fr(1, 660000), p); gap = [c * x % p for x in w]
        log("   p = 2^%d-1: total == (1/660000) (u/t_c)^8 through u^7: %s" % (89 if p == P89 else 61, [x % p for x in tot] == gap))
    # (c) single-prime reconstruction at the balanced bound sqrt(p/2) succeeds for almost ANY residue, so it is no test on its own.
    # A reconstruction is accepted only when it is SMALL: |num| * den < 2^60 at one prime (a random residue does this with probability
    # about 2^-29), or |num| * den < 2^120 after CRT (probability about 2^-29 for the product modulus).
    def small(f, bits): return f is not None and abs(f.numerator) * f.denominator < (1 << bits) and smooth(f.denominator)
    T89 = d89['table']; T61 = d61['table']; cnt = dict(zero_both=0, zero_one=0, s89_ok=0, s89_bad=0, s61_ok=0, s61_bad=0, crt_small=0, undecided=0)
    bad = []; hist = {}
    for i in range(len(d89['outer'])):
        a = T89.get(str(i), [0] * 8); b = T61.get(str(i), [0] * 8)
        for k in range(8):
            x, y = a[k] % P89, b[k] % P61
            if x == 0 and y == 0: cnt['zero_both'] += 1; continue
            if (x == 0) != (y == 0): cnt['zero_one'] += 1; bad.append((i, k, 'zero at one prime')); continue
            r = ratrec(x, P89); r2 = ratrec(y, P61); c = ratrec(crt(x, y), P89 * P61)
            if small(c, 120) and not (small(r, 60) and r != c):         # CRT first: a small single-prime hit that differs from a small
                cnt['crt_small'] += 1; hist[k] = hist.get(k, 0) + 1    # smooth CRT rational is a coincidence at that prime
                if r == c: cnt['s89_ok'] += 1
                continue
            if small(r, 60):
                if red(r, P61) == y: cnt['s89_ok'] += 1
                else: cnt['s89_bad'] += 1; bad.append((i, k, str(r)))
                continue
            if small(r2, 45):
                if red(r2, P89) == x: cnt['s61_ok'] += 1
                else: cnt['s61_bad'] += 1; bad.append((i, k, str(r2)))
                continue
            cnt['undecided'] += 1
    log("(c) per outer graph (4555 x 8 entries): %s" % cnt)
    log("   meaning: crt_small = the CRT reconstruction is a small rational (|num| den < 2^120) with a 1000-smooth denominator (consistent;"
        " by power of u: %s; s89_ok of them are already reconstructed at 2^89-1 alone); s89_bad / s61_bad = small smooth rational at one prime"
        " that does NOT reduce to the other residue and no small CRT rational (a real inconsistency, or a coincidence); undecided = no small"
        " rational even after CRT (height > 2^120, or not rational)" % sorted(hist.items()))
    if bad: log("   inconsistent entries (first 20): %s" % bad[:20])
    log_f.close()


if __name__ == "__main__":
    main()
