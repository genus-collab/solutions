# -*- coding: utf-8 -*-
"""f2_seldetail.py -- where exactly do the two frames differ in which pieces are non-zero?

For the pieces found by f2_sel.py to be present in one frame only, or to have a different u-support, this prints the
u-support in each frame and the edge numbers E(t), E(q), so one can see whether the difference sits at the top of the
u-truncation (an artefact of the K = 8 window) or at a low order (a genuine structural difference).

usage:  py -I -B f2_seldetail.py <parts frame (0,0)> <parts frame2> <result json of either run> [--json OUT]
"""
import os, sys, json, argparse
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_lib import K, prime, read_parts_dir

P = prime("89")


def sup(c):
    return [k for k in range(K) if c[k] % P]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirA")
    ap.add_argument("dirB")
    ap.add_argument("resjson")
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    _, pA, _ = read_parts_dir(a.dirA, "g5_89_nat")
    _, pB, _ = read_parts_dir(a.dirB, "g5_89_nat_ab1o3_m2o7")
    with open(a.resjson) as f:
        E = [G["E"] for G in json.load(f)["outer"]]

    keys = set()
    for d in (pA, pB):
        for t, res in d.items():
            for q in res:
                keys.add((int(t), int(q)))
    rows = []
    for (t, q) in sorted(keys):
        cA = pA.get(str(t), {}).get(str(q))
        cB = pB.get(str(t), {}).get(str(q))
        sA = sup(cA) if cA else []
        sB = sup(cB) if cB else []
        if sA != sB:
            rows.append(dict(t=t, q=q, Et=E[t], Eq=E[q], supA=sA, supB=sB))
    print("pieces whose u-support differs between the frames: %d" % len(rows))
    print("lowest u-order at which any of them differs: %s"
          % min(min(set(r["supA"]) ^ set(r["supB"])) for r in rows))
    print("distribution of the lowest differing u-order: %s"
          % sorted(Counter(min(set(r["supA"]) ^ set(r["supB"])) for r in rows).items()))
    print("by (E(t), E(q)): %s" % sorted(Counter((r["Et"], r["Eq"]) for r in rows).items()))
    one = [r for r in rows if not r["supA"] or not r["supB"]]
    print("pieces stored in one frame only: %d; their supports: %s"
          % (len(one), sorted(Counter(tuple(r["supA"] or r["supB"]) for r in one).items())))
    for r in rows[:12]:
        print("   t=%4d (E=%2d) -> q=%4d (E=%d): support (0,0) %s, second frame %s" % (r["t"], r["Et"], r["q"], r["Eq"], r["supA"], r["supB"]))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(dict(n=len(rows), rows=rows), f, indent=1)
        print("written %s" % a.json)


if __name__ == "__main__":
    main()
