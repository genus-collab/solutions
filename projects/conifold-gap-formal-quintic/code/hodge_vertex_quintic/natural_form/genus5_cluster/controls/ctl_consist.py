# -*- coding: utf-8 -*-
"""ctl_consist.py -- an independent arithmetic test of the two genus-5 controls from the copied-back JSONs alone.

   The sum of the rows over ALL outer graphs is the point-frame total, which does not depend on M (it is the same
   series for the variants id, nat, M0 and nat_t7; the copied-back id and nat runs agree at both primes).  Hence
   the perturbation of the rows sums to zero:   Sum_Gamma delta row_Gamma = 0  mod u^8.
   The perturbation of one row has u-valuation at least k + (E - 1) for nat_tk and at least E for M0 (E = edges of
   the outer graph, vertex graph E = 0, k >= 1), so the sum rule truncates to the graphs the --maxout 2 controls
   computed exactly:
       nat_t7, order u^7 : only E <= 1    -> rows 0, 1, 2, 3
       M0,     order u^1 : only E <= 1    -> rows 0, 1, 2, 3
       M0,     order u^2 : only E <= 2    -> rows 0..14
   Each is a non-trivial identity between the control run and the natural run, with no new graph sum.

   py -I -B ctl_consist.py
"""
import os, sys, json
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
RES = os.path.join(PKG, "results")
sys.path.insert(0, PKG)
from fq_fp import P2


def rows(tag):
    d = json.load(open(os.path.join(RES, tag + ".json")))
    return {int(k): v for k, v in d['table'].items()}, d['outer']


def main():
    p = P2
    nat, outer = rows("g5_89_nat")
    print("=" * 100)
    print("sum rule  Sum_Gamma (row_Gamma(control) - row_Gamma(nat)) = 0  mod u^8   (2^89-1)")
    print("=" * 100)
    E = [o['E'] for o in outer]
    cases = [("g5_89_nat_t7_out2", 7, 1), ("g5_89_M0_out2", 1, 1), ("g5_89_M0_out2", 2, 2)]
    allok = True
    for tag, order, emax in cases:
        ctl, _ = rows(tag)
        idx = [i for i in range(len(E)) if E[i] <= emax]
        s = 0; terms = []
        for i in idx:
            d = (ctl.get(i, [0] * 8)[order] - nat.get(i, [0] * 8)[order]) % p
            s = (s + d) % p
            if d: terms.append((i, E[i], d))
        ok = (s == 0)
        allok = allok and ok
        print()
        print("%-20s order u^%d, graphs with E <= %d (%d rows: %s)" % (tag, order, emax, len(idx), idx))
        for i, e, d in terms:
            print("      row %2d (E=%d): delta = %d" % (i, e, d))
        print("   rows with a zero perturbation at this order: %d" % (len(idx) - len(terms)))
        print("   SUM = %d  ->  %s" % (s, "PASS (the perturbation cancels, as the M-independence of the total requires)" if ok else "FAIL"))
    print()
    print("ALL SUM RULES PASS" if allok else "SOME SUM RULE FAILED")


if __name__ == "__main__":
    main()
