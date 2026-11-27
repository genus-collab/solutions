# -*- coding: utf-8 -*-
"""ctl_read.py -- read the copied-back genus-5 result JSONs and the stored genus-4 reference tables and tabulate
   exactly which outer graphs fail the regularity test in the two controls, at which orders, with E, gv, |Aut|.

   py -I -B ctl_read.py
   No cluster data beyond results/*.json is needed (logs/ and parts/ were not copied back).
"""
import os, sys, json, collections

PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
RES = os.path.join(PKG, "results")


def load(tag):
    return json.load(open(os.path.join(RES, tag + ".json")))


def gstr(o):
    return "gv=%s ed=%s |Aut|=%d" % (list(o['gv']), [list(e) for e in o['ed']], o['aut'])


def part1():
    print("=" * 110)
    print("(1) THE TWO GENUS-5 CONTROLS (both --maxout 2: exact for the vertex graph and for E <= 2)")
    print("=" * 110)
    nat = load("g5_89_nat")
    outer = nat['outer']
    hist = collections.Counter(o['E'] for o in outer)
    print("genus-5 leg-free stable graphs by number of edges E (from the 'outer' array, %d graphs):" % len(outer))
    print("   " + ", ".join("E=%d: %d" % (e, hist[e]) for e in sorted(hist)))
    print("   #(E=1) = %d, #(E=2) = %d, #(1 <= E <= 2) = %d" % (hist[1], hist[2], hist[1] + hist[2]))
    for tag in ("g5_89_M0_out2", "g5_89_nat_t7_out2"):
        d = load(tag)
        nr = {int(k): v for k, v in d['info']['nonregular'].items()}
        print()
        print("-" * 110)
        print("%s   failed checks: %s" % (tag, d['failed']))
        print("   rows stored in the table: %s (indices %s)" % (len(d['table']), sorted(int(k) for k in d['table'])))
        print("   non-regular outer graphs: %d" % len(nr))
        print("   %-5s %-3s %-34s %s" % ("idx", "E", "graph", "u-orders with a non-zero coefficient (regular = none of u^0..u^7)"))
        for i in sorted(nr):
            print("   %-5d %-3d %-34s %s" % (i, outer[i]['E'], gstr(outer[i]), nr[i]))
        tested = [i for i in range(len(outer)) if 1 <= outer[i]['E'] <= 2]
        ok = [i for i in tested if i not in nr]
        print("   tested (1 <= E <= 2): %d; failed: %d; still regular: %s" % (len(tested), len(nr), ok if ok else "none"))
        byE = collections.Counter(outer[i]['E'] for i in nr)
        print("   failures by E: %s" % sorted(byE.items()))
        first = {}
        for i in nr:
            first.setdefault(outer[i]['E'], set()).add(min(nr[i]))
        print("   first failing order by E: %s" % {e: sorted(s) for e, s in sorted(first.items())})
        V = d['info']['vertex']
        print("   vertex graph, polar part in t_c (coefficients of t_c^-8 .. t_c^-1), 0 = clean:")
        print("      %s" % ["0" if x == 0 else "NONZERO" for x in V])
        print("   [t_c^-8] vertex = %s ; c_vertex_equals_pattern = %s" % (V[0], d['info']['c_vertex_equals_pattern']))
    print()
    print("-" * 110)
    print("cross-check: [t_c^-8] of the vertex graph in all four 2^89-1 runs (must be identical: u^0 data is variant-independent)")
    for tag in ("g5_89_nat", "g5_89_nat_ab1o3_m2o7", "g5_89_M0_out2", "g5_89_nat_t7_out2"):
        d = load(tag)
        print("   %-24s %s   pattern=%s" % (tag, d['info']['vertex'][0], d['info']['c_vertex_equals_pattern']))
    print("   (frame (1/3,-2/7) is a different frame, so its t_c normalisation need not agree)")
    nz = [k for k in nat['table'] if any(nat['table'][k])]
    print("   g5_89_nat: %d stored rows, non-zero rows = %s (the vertex graph alone)" % (len(nat['table']), sorted(int(x) for x in nz)))
    e_stored = collections.Counter(outer[int(k)]['E'] for k in nat['table'])
    print("   g5_89_nat stored rows by E: %s  (no row with E >= 8: those vanish mod u^8 for free)" % sorted(e_stored.items()))


def part2():
    print()
    print("=" * 110)
    print("(2) THE GENUS-4 RECORD, recomputed from the stored tables data/ref_g4.json (orders u^0..u^5, 2^89-1)")
    print("=" * 110)
    ref = json.load(open(os.path.join(PKG, "data", "ref_g4.json")))
    print("runs stored: %s" % sorted(ref))
    base = ref["g4_89_nat"]['table']
    keys = list(base)
    E = {k: len(json.loads(k)[1]) for k in keys}
    hist = collections.Counter(E.values())
    print("genus-4 graphs by E: %s  (total %d)" % (sorted(hist.items()), len(keys)))
    print("   #(E=1) = %d, #(1 <= E <= 5) = %d" % (hist[1], sum(v for e, v in hist.items() if 1 <= e <= 5)))
    for name in ("g4_89_nat", "g4_89_M0", "g4_89_nat_t5"):
        tab = ref[name]['table']
        bad = {k: [j for j in range(6) if tab[k][j]] for k in keys if E[k] >= 1 and any(tab[k][:6])}
        byE = collections.Counter(E[k] for k in bad)
        print()
        print("   %-14s non-regular outer graphs with edges: %d; by E: %s" % (name, len(bad), sorted(byE.items())))
        if len(bad) <= 6:
            for k in sorted(bad, key=lambda x: (E[x], x)):
                gv, ed = json.loads(k)
                print("      gv=%s ed=%s  failing orders %s" % (gv, ed, bad[k]))
        else:
            fo = {}
            for k in bad:
                fo.setdefault(E[k], set()).add(min(bad[k]))
            print("      first failing order by E: %s" % {e: sorted(s) for e, s in sorted(fo.items())})
            allo = collections.Counter(tuple(v) for v in bad.values())
            print("      distinct failing-order patterns: %s" % dict(list(allo.items())[:8]))
            tested = [k for k in keys if 1 <= E[k] <= 5]
            print("      tested (1 <= E <= 5): %d; failed: %d; still regular: %d" % (len(tested), len(bad), len(tested) - len(bad)))


if __name__ == "__main__":
    part1(); part2()
