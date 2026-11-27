# -*- coding: utf-8 -*-
"""ig_more.py -- build the follow-up lists once the first p89 pass is stored.

  py -I -B ig_more.py --sup 12000 --p61 20000

--sup   budget (predicted CPU-s) for a SUPPLEMENT to the p89 natural run: the package-heaviest graphs that are still missing, then
        a cost-spread top-up of the strata E = 6..11 (the strata that are sampled, not complete).  Written to sample_sup89.txt.
--p61   budget for the SECOND PRIME: the same deliberate shape as the p89 sample (complete strata E <= 5 and E = 12, the
        package-heaviest, then a cost-spread sample of E = 6..11), restricted to graphs already done at p89 so that the two primes
        compare the same families.  Written to sample_nat61b.txt (and sample_id61b.txt, the same list for the point frame).
The graphs already stored in out/nat_89.jsonl are excluded from the supplement (ig_run would skip them anyway).
"""
import os, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RATE = 5e-6


def done(tag):
    p = os.path.join(OUT, "%s.jsonl" % tag)
    s = set()
    if os.path.exists(p):
        for line in open(p, encoding="utf-8"):
            if line.strip():
                s.add(json.loads(line)["t"])
    return s


def spread(pool, pred, budget):
    pool = sorted(pool, key=lambda t: pred[t])
    step = max(1, len(pool) // 200)
    order = []
    for off in range(step):
        order += pool[off::step]
    out = []
    spent = 0.0
    for t in order:
        if spent + pred[t] <= budget:
            out.append(t)
            spent += pred[t]
    return out, spent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sup", type=float, default=12000.0)
    ap.add_argument("--p61", type=float, default=20000.0)
    a = ap.parse_args()
    ix = json.load(open(os.path.join(OUT, "index.json")))
    G = {int(t): v for t, v in ix["graphs"].items()}
    ci = json.load(open(os.path.join(OUT, "cost_ind.json")))
    pred = {t: max(1.0, RATE * ci[str(t)][3]) for t in G}
    D = done("nat_89")
    byE = {}
    for t in G:
        byE.setdefault(G[t]["E"], []).append(t)
    heavy = sorted(G, key=lambda t: -G[t]["t89"])
    print("p89 already stored: %d inner graphs" % len(D))

    sup = []
    spent = 0.0
    for t in heavy[:20]:
        if t not in D and spent + pred[t] <= a.sup * 0.6:
            sup.append(t)
            spent += pred[t]
    print("   supplement: %d of the 20 package-heaviest graphs added (%.0f s)" % (len(sup), spent))
    for e in (6, 7, 8, 9, 10, 11):
        bud = (a.sup - spent) / (12 - e)
        got, s = spread([t for t in byE[e] if t not in D and t not in sup], pred, bud)
        sup += got
        spent += s
        print("   supplement E=%2d: %3d graphs (%.0f s)" % (e, len(got), s))
    print("SUPPLEMENT p89: %d graphs, predicted %.0f CPU-s" % (len(sup), spent))

    s61 = set()
    spent61 = 0.0
    for e in list(range(6)) + [12]:
        for t in byE[e]:
            if t in D:
                s61.add(t)
                spent61 += pred[t]
    print("p61: complete strata E <= 5 and E = 12 that are done at p89: %d graphs (%.0f s)" % (len(s61), spent61))
    for t in heavy[:20]:
        if t in D and t not in s61 and spent61 + pred[t] <= a.p61:
            s61.add(t)
            spent61 += pred[t]
    for e in (6, 7, 8, 9, 10, 11):
        bud = max(0.0, (a.p61 - spent61) / (12 - e))
        got, s = spread([t for t in byE[e] if t in D and t not in s61], pred, bud)
        s61 |= set(got)
        spent61 += s
    print("SAMPLE p61: %d graphs, predicted %.0f CPU-s; pieces %d"
          % (len(s61), spent61, sum(G[t]["npieces61"] for t in s61)))

    def write(name, items, note):
        with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
            f.write("# %s\n" % note)
            for t in sorted(items):
                f.write("%d   # E=%d V=%d pred=%.0fs\n" % (t, G[t]["E"], G[t]["V"], pred[t]))
        print("wrote %s (%d graphs)" % (name, len(items)))

    write("sample_sup89.txt", sup, "supplement to the p89 natural run")
    write("sample_nat61b.txt", sorted(s61), "second prime 2^61-1, natural M")
    write("sample_id61b.txt", sorted(s61), "second prime 2^61-1, point frame")


if __name__ == "__main__":
    main()
