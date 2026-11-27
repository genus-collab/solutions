# -*- coding: utf-8 -*-
"""rc_t1_recon.py -- part (1) of the re-check: CRT + rational reconstruction of total and vertex,
from the residues stored in the two natural runs, with margins against the uniqueness bound.
"""
import os
import sys
import json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rc_rec

RES = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..",
                                   "cluster_genus5_quintic", "results"))


def load(tag):
    with open(os.path.join(RES, tag + ".json")) as f:
        return json.load(f)


ORCH = ["1/660000", "-7/825000", "199/9900000", "-717/27500000", "449357/22500000000",
        "-5593217/618750000000", "14686042363/6496875000000000", "-8370858511/32484375000000000"]


def fmt(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def reconstruct(r89, r61, label):
    out = []
    print("--- %s" % label)
    for k in range(len(r89)):
        v, n, d, cert, mg = rc_rec.recon_pair(r89[k], r61[k])
        out.append(v)
        print("   u^%d  %-34s  certified=%s  |n|=%d digits, d=%d digits, margin = mod/(2|n|d) = %.3e"
              % (k, fmt(v) if v is not None else "FAILED", cert,
                 len(str(abs(n))), len(str(d)), mg))
    return out


def main():
    n89 = load("g5_89_nat")
    n61 = load("g5_61_nat")
    print("modulus M = p89 * p61 has %d decimal digits; sqrt(M/2) has %d"
          % (len(str(rc_rec.M)), len(str(rc_rec.isqrt(rc_rec.M // 2)))))
    print()
    tot = reconstruct(n89['info']['total'], n61['info']['total'], "total (u^0..u^7), from info['total']")
    print()
    row89 = n89['table']['0']
    row61 = n61['table']['0']
    ver = reconstruct(row89, row61, "vertex graph row (table['0']), u^0..u^7")
    print()
    print("--- comparison with the orchestrator's series")
    ok = True
    for k, s in enumerate(ORCH):
        want = Fr(s)
        got = tot[k]
        same = (got == want)
        ok = ok and same
        print("   u^%d  orchestrator %-34s  mine %-34s  %s" % (k, s, fmt(got), "EQUAL" if same else "DIFFERENT"))
    print("   total series agrees with the orchestrator in all 8 coefficients:", ok)
    print()
    print("--- vertex row vs total")
    print("   identical:", ver == tot)
    print()
    print("--- single-prime reconstruction (p89 only), to show that two primes are needed")
    for k in range(8):
        n, d, cert = rc_rec.ratrec(n89['info']['total'][k], rc_rec.P89)
        print("   u^%d  %s  certified=%s  (two-prime value %s)"
              % (k, ("%d/%d" % (n, d)), cert, fmt(tot[k])))
    with open(os.path.join(HERE, "rc_total_rational.json"), "w") as f:
        json.dump(dict(total=[fmt(x) for x in tot], vertex=[fmt(x) for x in ver]), f, indent=1)
    print()
    print("written:", os.path.join(HERE, "rc_total_rational.json"))


if __name__ == "__main__":
    main()
