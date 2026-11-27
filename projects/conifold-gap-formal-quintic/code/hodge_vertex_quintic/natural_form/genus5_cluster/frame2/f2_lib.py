# -*- coding: utf-8 -*-
"""f2_lib.py -- helpers for the SECOND-FRAME (alpha, beta) = (1/3, -2/7) audit of the genus-5 formal-quintic run.

Independent of g5_checks.checks: the part files are read and summed here, with this module's own code.  Nothing in
cluster_genus5_quintic/ is imported except fq_fp (the two primes, the field) and, for the t_c change of variables only,
fq_conifold / fq_hat / fq3_msum (no check logic lives there).

Part file format (see g5_checks.py):  {"tag", "pieces": {inner t: {outer q: [K residues]}}, "times": {t: seconds}}
A piece is stored only when at least one of its K residues is non-zero.  For the variant 'id' the only stored key is q = t
(the point-frame term of the inner graph t itself).
"""
import os, sys, json, re
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE)))), "cluster_genus5_quintic")
assert os.path.isdir(PKG), PKG
if PKG not in sys.path:
    sys.path.insert(0, PKG)

K = 8            # genus 5: u^0 .. u^7
NIN = 4555       # leg-free stable graphs of genus 5
AB = "1/3,-2/7"  # the second conifold frame
P89 = (1 << 89) - 1
P61 = (1 << 61) - 1


def prime(name):
    return P89 if str(name) == "89" else P61


def read_parts_dir(d, expect_tag=None):
    """merge every task<I>of<N>.json / local.json of the directory d; returns (names, pieces, times).

    An inner graph occurring in two files must carry identical content (asserted); it is counted once.
    """
    names = sorted(f for f in os.listdir(d) if re.match(r"^(task\d+of\d+|local)\.json$", f))
    assert names, "no part files in %s" % d
    pieces = {}
    times = {}
    for f in names:
        with open(os.path.join(d, f)) as fh:
            part = json.load(fh)
        if expect_tag is not None:
            assert part["tag"] == expect_tag, "%s carries the tag %s, expected %s" % (f, part["tag"], expect_tag)
        for t, res in part["pieces"].items():
            if t in pieces:
                assert pieces[t] == res, "inner graph %s differs between two part files in %s" % (t, d)
                continue
            pieces[t] = res
            times[t] = part["times"].get(t, 0.0)
    return names, pieces, times


def crt2(a1, p1, a2, p2):
    """the unique x mod p1 p2 with x = a1 mod p1, x = a2 mod p2"""
    a1 %= p1
    a2 %= p2
    x = (a1 + p1 * ((a2 - a1) * pow(p1, -1, p2) % p2)) % (p1 * p2)
    return x


def ratrec(a, m, bound=None):
    """rational reconstruction of a mod m by continued fractions: n/d with |n|, d <= bound (default sqrt(m/2)/8)."""
    if bound is None:
        bound = int((m // 2) ** 0.5) // 8
    a %= m
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
        if r0 <= bound and abs(s0) <= bound and s0 != 0:
            n, d = (r0, s0) if s0 > 0 else (-r0, -s0)
            if (n - d * a) % m == 0:
                return Fr(n, d)
    return None
