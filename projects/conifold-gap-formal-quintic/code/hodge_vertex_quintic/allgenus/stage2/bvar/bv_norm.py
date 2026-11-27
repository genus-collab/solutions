# -*- coding: utf-8 -*-
"""bv_norm.py -- genus-4 normalisation response for a B-perturbed system: (stored - canonical)/m5 per outer graph (maxout 2 tables),
m5 = [z^5] P^stored_11(0,z) / P_11(0,0).  Predictions (univ/UNIVERSALITY.md s.2.3, valuation argument): the difference is confined
to graphs 0, 1, 2, 3, 4, 7, 8; its [u^4] column is the universal frozen response (vertex -1/8400000, bridge (1,3) -1/1008000,
bridge (2,2) 1/900000, others 0); the sum over the outer graphs vanishes at [u^4] and [u^5] (the total is normalisation independent).
Run: py -I -B bv_norm.py <dir> SYSTEM   (needs tab_<SYS>_g4_89_canon_none_out2.json and tab_<SYS>_g4_89_stored_none_out2.json)
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import bv_sys
from bv_run import gauge_of
from un_e4_cases import p11_series
from fq_fp import FpRing, P2
from fq3_rec import Rec


def main():
    d, name = sys.argv[1], sys.argv[2]
    R = FpRing(P2); p = R.p; rec = Rec(R)
    a = json.load(open(os.path.join(d, "tab_%s_g4_89_canon_none_out2.json" % name)))
    b = json.load(open(os.path.join(d, "tab_%s_g4_89_stored_none_out2.json" % name)))
    P11 = p11_series(gauge_of(name)['p']); m5 = P11[5] / P11[0]; im5 = R.inv(R.of(m5))
    print("system %s: m5 = %s" % (name, m5))
    ta = {int(i): c for i, c in a['table'].items()}; tb = {int(i): c for i, c in b['table'].items()}
    E = {i: (Ee, gv) for (i, Ee, gv, c) in a['rows']}; E[0] = (0, [4])
    tot = [0] * 6; diffs = {}
    for i in sorted(set(ta) | set(tb)):
        if i not in E or E[i][0] > 2: continue                    # graphs with > 2 edges are incomplete at maxout 2
        dd = [(tb.get(i, [0] * 6)[k] - ta.get(i, [0] * 6)[k]) * im5 % p for k in range(6)]
        if any(dd): diffs[i] = dd
        tot = [(x + y) % p for x, y in zip(tot, dd)]
    for i, dd in diffs.items():
        print("   outer graph %d (E=%d, gv=%s): (stored - canon)/m5 [u^0..u^5] = %s" % (i, E[i][0], E[i][1], [rec.fmt(x) for x in dd]))
    print("   graphs that move: %s (prediction: subset of [0, 1, 2, 3, 4, 7, 8])" % sorted(diffs))
    pred = {0: Fr(-1, 8400000), 1: Fr(-1, 1008000), 2: Fr(1, 900000)}
    ok4 = all(diffs.get(i, [0] * 6)[4] == R.of(pred.get(i, 0)) for i in set(diffs) | set(pred))
    print("   [u^4] column = universal frozen response (-1/8400000, -1/1008000, 1/900000 on graphs 0, 1, 2; 0 else): %s" % ("PASS" if ok4 else "FAIL"))
    print("   sum over outer graphs of (stored - canon)/m5: %s  -> %s" % ([rec.fmt(x) for x in tot], "total unchanged: PASS" if not any(tot) else "FAIL"))
    lo = all(not any(dd[:4]) for dd in diffs.values())
    print("   no difference below [u^4]: %s" % ("PASS" if lo else "FAIL"))


if __name__ == "__main__":
    main()
