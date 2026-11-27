# -*- coding: utf-8 -*-
"""genus3_qzeta_rec.py -- reconstruction of elements of Q(zeta) from their residues at the two primes 2^61-1 and 2^89-1 (genus-3 O(u^3) data).

An element x = (n0 + n1 zeta + n2 zeta^2 + n3 zeta^3) / d of Q(zeta) with small integers is recovered from x mod p1, x mod p2 (with the
5th roots of unity zeta_1, zeta_2 used by FpRing at the two primes; CRT gives a ring map Z[zeta] -> Z/M, M = p1 p2 ~ 2^150) as a short
vector of the lattice {(d, n1, n2, n3, n0) : d x = n0 + n1 zeta + n2 zeta^2 + n3 zeta^3 mod M} (volume M, generic shortest vector ~ M^(1/5)
~ 2^30) by LLL (exact integer / rational arithmetic).  A vector with entries below ~2^22 is far below the generic length and is accepted
as the reconstruction ('strong'); anything longer is reported as not reconstructed.  The result is only as good as the assumption that the
two residues are images of the same element of Q(zeta) (true for canonically normalised solutions of the same Q(zeta)-rational system).
Applied to: the plane Pi (keys Pi_A, Pi_B of g3_cubic_<p>.json) and the isolated solution circ_b_solutions[0].
Run: py -I genus3_qzeta_rec.py
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2


def lll(B, delta=Fr(3, 4)):
    n = len(B); B = [list(b) for b in B]
    def gs():
        Bs = []; mu = [[Fr(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fr(x) for x in B[i]]
            for j in range(i):
                den = sum(x * x for x in Bs[j])
                mu[i][j] = sum(Fr(B[i][k]) * Bs[j][k] for k in range(n)) / den if den else Fr(0)
                v = [v[k] - mu[i][j] * Bs[j][k] for k in range(n)]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gs(); k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [B[k][t] - q * B[j][t] for t in range(n)]; Bs, mu = gs()
        nk = sum(x * x for x in Bs[k]); nk1 = sum(x * x for x in Bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1: k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]; Bs, mu = gs(); k = max(k - 1, 1)
    return B


def crt(a1, a2):
    M = P1 * P2
    return (a1 * P2 * pow(P2, P1 - 2, P1) + a2 * P1 * pow(P1, P2 - 2, P2)) % M


def reconstruct(x1, x2, z):
    """x1 mod P1, x2 mod P2, z = CRT of the two roots of unity; returns (d, [n0, n1, n2, n3], bits) of the shortest LLL vector"""
    M = P1 * P2; x = crt(x1, x2)
    B = [[1, 0, 0, 0, x], [0, 1, 0, 0, (-z) % M], [0, 0, 1, 0, (-z * z) % M], [0, 0, 0, 1, (-pow(z, 3, M)) % M], [0, 0, 0, 0, M]]
    R = lll(B); v = min(R, key=lambda r: sum(t * t for t in r))
    d, n1, n2, n3, n0 = v
    if d < 0: d, n1, n2, n3, n0 = -d, -n1, -n2, -n3, -n0
    return d, [n0, n1, n2, n3], max(abs(t) for t in v).bit_length()


def fmt(d, n):
    terms = []
    for k, c in enumerate(n):
        if c: terms.append(("%d" % c) + ("" if k == 0 else " zeta" if k == 1 else " zeta^%d" % k))
    s = " + ".join(terms).replace("+ -", "- ") or "0"
    return s if d == 1 else "(%s)/%d" % (s, d)


def main():
    R1 = FpRing(P1); R2 = FpRing(P2); z = crt(R1.zeta, R2.zeta)
    c1 = json.load(open(os.path.join(HERE, "g3_cubic_61.json"))); c2 = json.load(open(os.path.join(HERE, "g3_cubic_89.json")))
    names = c1['vis_names'] + ['c5_0', 'c5_1', 'c5_2', 'c6_1', 'c6_2']
    LOG = open(os.path.join(HERE, "logs", "genus3_qzeta_rec.log"), "w", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    # self-test: a known element
    t1 = (3 + 5 * R1.zeta - 7 * pow(R1.zeta, 3, P1)) * pow(1920, P1 - 2, P1) % P1; t2 = (3 + 5 * R2.zeta - 7 * pow(R2.zeta, 3, P2)) * pow(1920, P2 - 2, P2) % P2
    d, n, bits = reconstruct(t1, t2, z); log("self-test (3 + 5 zeta - 7 zeta^3)/1920 ->", fmt(d, n), "; bits", bits)
    assert (d, n) == (1920, [3, 5, 0, -7])
    out = {}
    for key, v1, v2 in (("Pi_A", c1['Pi_A'], c2['Pi_A']), ("Pi_B", c1['Pi_B'], c2['Pi_B']), ("isolated", c1['circ_b_solutions'][0], c2['circ_b_solutions'][0])):
        log("== %s" % key); out[key] = {}
        for j in range(15):
            if v1[j] == 0 and v2[j] == 0: log("   %-5s 0" % names[j]); out[key][names[j]] = "0"; continue
            d, n, bits = reconstruct(v1[j], v2[j], z)
            ok = bits <= 22
            log("   %-5s %s   [%d bits: %s]" % (names[j], fmt(d, n) if ok else "not reconstructed", bits, "strong" if ok else "generic length, rejected"))
            out[key][names[j]] = fmt(d, n) if ok else None
    json.dump(out, open(os.path.join(HERE, "g3_qzeta_rec.json"), "w"), indent=1); LOG.close()


if __name__ == "__main__":
    main()
