# -*- coding: utf-8 -*-
"""rc_t4_coherence.py -- parts (4) and (5) of the re-check: cross-run coherence of the seven copied JSONs,
the shipped point-frame part files, the controls, and what the stored 'failed' lists do and do not cover.
"""
import os
import sys
import json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rc_rec

PKG = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "cluster_genus5_quintic"))
RES = os.path.join(PKG, "results")
TAGS = ["g5_89_nat", "g5_61_nat", "g5_89_nat_ab1o3_m2o7", "g5_89_M0_out2", "g5_89_nat_t7_out2",
        "g5_89_id", "g5_61_id"]
P = {"89": rc_rec.P89, "61": rc_rec.P61}


def load(tag):
    with open(os.path.join(RES, tag + ".json")) as f:
        return json.load(f)


def main():
    D = {t: load(t) for t in TAGS}
    K = 8
    print("--- the 'outer' graph list is the same in all seven runs")
    ref = D["g5_89_nat"]['outer']
    for t in TAGS:
        o = D[t]['outer']
        print("    %-26s len=%d  identical to g5_89_nat: %s" % (t, len(o), o == ref))
    print("    (compared field by field: E, gv, ed, aut)")
    print()

    print("--- row sums vs info['total']")
    for t in TAGS:
        d = D[t]
        p = P[t.split("_")[1]]
        s = [sum(r[k] for r in d['table'].values()) % p for k in range(K)]
        tot = d['info'].get('total')
        print("    %-26s sum over the %4d stored rows == info['total']: %s"
              % (t, len(d['table']), "n/a (no total stored: --maxout run)" if tot is None else (s == tot)))
    print()

    print("--- natural total == point-frame total, prime by prime")
    for pr in ("89", "61"):
        a = D["g5_%s_nat" % pr]['info']['total']
        b = D["g5_%s_id" % pr]['info']['total']
        print("    p%s: g5_%s_nat total == g5_%s_id total: %s" % (pr, pr, pr, a == b))
    a = D["g5_89_nat_ab1o3_m2o7"]['info']['total']
    print("    second frame total differs from frame (0,0) at p89: %s" % (a != D["g5_89_nat"]['info']['total']))
    print()

    print("--- the shipped point-frame part files (parts/g5_<p>_id/local.json)")
    for pr in ("89", "61"):
        path = os.path.join(PKG, "parts", "g5_%s_id" % pr, "local.json")
        if not os.path.exists(path):
            print("    p%s: MISSING %s" % (pr, path))
            continue
        with open(path) as f:
            part = json.load(f)
        pieces = part['pieces']
        p = P[pr]
        miss = sorted(set(range(4555)) - set(int(k) for k in pieces))
        print("    p%s: tag=%s, %d inner graphs, missing %d" % (pr, part['tag'], len(pieces), len(miss)))
        # for variant id each inner graph contributes only its own index
        offdiag = [t for t, r in pieces.items() if set(r) - {t}]
        print("       entries off the diagonal (inner index != outer index): %d" % len(offdiag))
        tabfromparts = {}
        for t, r in pieces.items():
            for q, c in r.items():
                a = tabfromparts.setdefault(q, [0] * K)
                for k in range(K):
                    a[k] = (a[k] + c[k]) % p
        stored = D["g5_%s_id" % pr]['table']
        same = (set(tabfromparts) == set(stored)) and all(tabfromparts[q] == stored[q] for q in stored)
        print("       re-merged table equals results/g5_%s_id.json table: %s" % (pr, same))
        s = [sum(r[k] for r in tabfromparts.values()) % p for k in range(K)]
        print("       re-merged total equals the natural run's total: %s"
              % (s == D["g5_%s_nat" % pr]['info']['total']))
    print()

    print("--- controls")
    nat = D["g5_89_nat"]
    for t in ("g5_89_M0_out2", "g5_89_nat_t7_out2"):
        d = D[t]
        keys = sorted(int(k) for k in d['table'])
        Es = sorted(set(ref[i]['E'] for i in keys))
        allle2 = sorted(i for i, g in enumerate(ref) if g['E'] <= 2)
        print("    %-20s %d stored rows, edge counts present: %s" % (t, len(keys), Es))
        print("       stored rows are exactly the %d graphs with E <= 2: %s" % (len(allle2), keys == allle2))
        print("       failed = %s" % d['failed'])
        nr = d['info']['nonregular']
        print("       non-regular graphs: %d -> %s" % (len(nr), {k: (ref[int(k)]['E'], ref[int(k)]['gv'], v) for k, v in sorted(nr.items(), key=lambda x: int(x[0]))[:12]}))
        print("       info keys: %s (no 'total': --maxout runs skip the GAP check)" % sorted(d['info'].keys()))
        print("       c_vertex_equals_pattern = %s" % d['info']['c_vertex_equals_pattern'])
        agree = {}
        for i in keys:
            a = d['table'][str(i)]
            b = nat['table'].get(str(i), [0] * K)
            agree[i] = max([k for k in range(K) if a[k] == b[k] and all(a[j] == b[j] for j in range(k + 1))] or [-1])
        print("       agreement with g5_89_nat on these rows: first differing u-power per row")
        for i in keys:
            a = d['table'][str(i)]
            b = nat['table'].get(str(i), [0] * K)
            fd = next((k for k in range(K) if a[k] != b[k]), None)
            print("          graph %3d (E=%d, gv=%s): first difference at u^%s"
                  % (i, ref[i]['E'], ref[i]['gv'], "none" if fd is None else fd))
    print()

    print("--- info key sets (which checks the collect step could have stored)")
    for t in TAGS:
        print("    %-26s %s" % (t, sorted(D[t]['info'].keys())))


if __name__ == "__main__":
    main()
