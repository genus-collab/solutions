# -*- coding: utf-8 -*-
"""ctl_g4.py -- the same prediction tested on the FULL genus-4 tables of record (cluster_genus5_quintic/data/ref_g4.json:
   all 379 outer graphs through u^5 for nat, id, M0, nat_t5, 2^89-1).  Nothing is recomputed; nothing is written.

   (a) the row sum over all 379 outer graphs is the point-frame total and does not depend on M (id, nat, M0, nat_t5 agree);
   (b) val_u(row(M0) - row(nat)) = E exactly for every outer graph with 1 <= E <= 5 (and the vertex graph: valuation 1);
   (c) val_u(row(nat_t5) - row(nat)) >= 5 + E - 1, so only the vertex graph and the three one-edge graphs move at all.

   py -I -B ctl_g4.py
"""
import os, sys, json, collections
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
from fq_fp import P2

K = 6


def main():
    p = P2
    ref = json.load(open(os.path.join(PKG, "data", "ref_g4.json")))
    nat = ref["g4_89_nat"]['table']
    keys = sorted(nat, key=lambda k: (len(json.loads(k)[1]), k))
    E = {k: len(json.loads(k)[1]) for k in keys}
    print("=" * 100)
    print("genus 4, 2^89-1, full tables of record (379 outer graphs, u^0..u^5)")
    print("=" * 100)
    print("(a) row sum over all outer graphs (must be the M-independent point-frame total)")
    tots = {}
    for name in ("g4_89_id", "g4_89_nat", "g4_89_M0", "g4_89_nat_t5"):
        tab = ref[name]['table']
        tots[name] = [sum(tab[k][j] for k in tab) % p for j in range(K)]
        print("   %-14s %s" % (name, tots[name]))
    ok = len(set(tuple(v) for v in tots.values())) == 1
    print("   all four totals equal: %s" % ("PASS" if ok else "FAIL"))

    def val(d):
        for j in range(K):
            if d[j]: return j
        return K

    for name, k in (("g4_89_M0", 1), ("g4_89_nat_t5", 5)):
        tab = ref[name]['table']
        print()
        print("(%s) %s: valuation of row(control) - row(nat) against the prediction val >= %s"
              % ("b" if k == 1 else "c", name, "E (M0)" if k == 1 else "k + E - 1 = %d + E - 1" % k))
        bad = []
        rowsmoved = 0
        hist = collections.Counter()
        for key in keys:
            d = [(tab[key][j] - nat[key][j]) % p for j in range(K)]
            v = val(d)
            e = E[key]
            floor = max(e, 1) if k == 1 else (k + e - 1 if e >= 1 else k)
            hist[(e, v)] += 1
            if v < min(floor, K): bad.append((key, e, v, floor))
            if v < K: rowsmoved += 1
        print("   rows that move at all (within u^0..u^5): %d of %d" % (rowsmoved, len(keys)))
        print("   (E, valuation of the perturbation) histogram, valuation 6 = unchanged mod u^6:")
        for e in sorted(set(E.values())):
            row = {v: n for (ee, v), n in hist.items() if ee == e}
            print("      E=%d: %s" % (e, {v: row[v] for v in sorted(row)}))
        print("   predicted floor violated by: %s" % (bad if bad else "no graph (PASS)"))
        if k == 5:
            moved = [key for key in keys if any((tab[key][j] - nat[key][j]) % p for j in range(K))]
            print("   rows that move: %s" % [(json.loads(x)[0], json.loads(x)[1]) for x in moved])


if __name__ == "__main__":
    main()
