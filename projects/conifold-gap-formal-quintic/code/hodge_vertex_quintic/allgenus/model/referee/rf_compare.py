# -*- coding: utf-8 -*-
"""rf_compare.py -- (1) the one-point table h = 1..4 and the one-point law;  (2) EVERY line of the worker's ms_singlist_full.log
(value at v = 1 of each singular basis correlator) recomputed exactly with the referee engine;  (3) the three genus-dependent Hessians.
Reads the worker's log as DATA only.  Run: py -I -B rf_compare.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat as rf
import rf_tables as rt
import rf_gsum as gsm
WORKER = os.path.normpath(os.path.join(HERE, ".."))


def have(h, n):
    m1, m2 = rf.models()
    return all(os.path.exists(os.path.join(gsm.SCR, "rfamp_%d_1_%d_%d.json" % (m.p.bit_length(), h, n))) for m in (m1, m2))


def have61(h, n):
    m1, m2 = rf.models()
    return os.path.exists(os.path.join(gsm.SCR, "rfamp_61_1_%d_%d.json" % (h, n)))


def parse_ins(s):
    out = []
    for tok in s.split():
        k, a = tok[1:].split("(H")
        out.append((int(a.rstrip(")")), int(k)))
    return tuple(out)


def main():
    t0 = time.time()
    print("(1) one-point functions <tau_k(H^b)>_h (b + k = 1 mod 5), law PiC(X) F_h'(1), remainder")
    for h in (1, 2, 3, 4):
        if not have(h, 1): print("   genus %d: table missing" % h); continue
        for k in range(0, 3 * h - 1):
            for b in range(5):
                if (b + k) % 5 != 1: continue
                x = rt.exact_val(h, ((b, k),)); e = rf.vexp(h, ((b, k),)); pc = rf.pi(0, b, k)
                law = pc * rt.Fder(h, 1); rem = x - law
                pi5 = rf.pi(5, b, k)
                print("   h=%d <t%d(H%d)> e=%3s value %-22s law %-22s remainder %-22s Pi_5 %s inK %s"
                      % (h, k, b, e, x, law, rem, pi5, rf.in_K(b, k)))
    print("(2) the worker's list of singular correlators (ms_singlist_full.log), recomputed")
    nchk = nbad = nskip = n61 = 0
    with open(os.path.join(WORKER, "ms_singlist_full.log")) as f:
        lines = f.read().splitlines()[1:]
    for ln in lines:
        parts = [x.strip() for x in ln.split("|")]
        h = int(parts[0]); ins = parse_ins(parts[1]); val = Fr(parts[3])
        if not have(h, len(ins)):
            if have61(h, len(ins)):
                m1 = rf.models()[0]; r = rt.flat_table(m1, h, len(ins)).get(ins, 0); n61 += 1
                if r != m1.fr(val): nbad += 1; print("   MISMATCH (mod 2^61-1) h=%d %s" % (h, parts[1]))
            else: nskip += 1
            continue
        mine = rt.exact_val(h, ins); nchk += 1
        if mine != val:
            nbad += 1; print("   MISMATCH h=%d %s: worker %s, referee %s" % (h, parts[1], val, mine))
    print("   %d lines recomputed exactly (two primes), %d compared mod 2^61-1 only, %d mismatches, %d skipped" % (nchk, n61, nbad, nskip))
    print("(3) Hessians Hess_h(X,Y) = (<XY>_h - F_h'' PiC PiC)/F_h' for the three pairs reported genus dependent")
    for X, Y in (((1, 1), (3, 2)), ((1, 2), (3, 1)), ((1, 3), (3, 0))):
        row = []
        for h in (2, 3):
            if not have(h, 2): continue
            v = rt.exact_val(h, (X, Y)); e = rf.vexp(h, (X, Y))
            H = (v - rt.Fder(h, 2) * rf.pi(0, *X) * rf.pi(0, *Y)) / rt.Fder(h, 1)
            row.append("h=%d e=%s <XY>=%s Hess=%s" % (h, e, v, H))
        print("   t%d(H%d) t%d(H%d):  %s" % (X[1], X[0], Y[1], Y[0], "   ".join(row)))
    print("rf_compare done (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
