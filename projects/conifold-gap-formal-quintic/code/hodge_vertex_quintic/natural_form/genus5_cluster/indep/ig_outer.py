# -*- coding: utf-8 -*-
"""ig_outer.py -- two kinds of extra bookkeeping on top of ig_cmp.

  py -I -B ig_outer.py --mode extra        write sample_extra89.txt: the inner graphs still missing from the p89 sample that carry
                                           a non-zero piece on one of the two cheapest-to-complete outer graphs (q = 1502, 839)
  py -I -B ig_outer.py --mode rows         for the outer graphs all of whose non-zero contributors are in the independent run,
                                           add up the INDEPENDENT pieces and compare with the cluster's own result table
                                           results/g5_89_nat.json (zero for every graph with an edge: the cancellation itself)

The "rows" mode is the only place where the independent code reproduces a whole outer-graph row rather than single pieces.  It is
complete only in this sense: every inner graph that the cluster records with a NON-ZERO piece on that outer graph has been run
independently, and for each of those inner graphs the independent run also confirms the absence of every piece the cluster does
not record.  That no further inner graph contributes to the row is the cluster's claim, checked independently only for the inner
graphs in the sample.
"""
import os, sys, json, argparse, glob

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
EV = os.environ.get("G5EV", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "g5ev"))
PKG = os.environ.get("G5PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic"))
TARGETS = [1502, 839]
PRIME = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def contributors():
    ix = json.load(open(os.path.join(OUT, "index.json")))
    G = {int(t): v for t, v in ix["graphs"].items()}
    c = {}
    for t, v in G.items():
        for q in v["qs"]:
            c.setdefault(q, set()).add(t)
    return G, c


def mode_extra():
    G, c = contributors()
    samp = set(json.load(open(os.path.join(OUT, "sample.json")))["id89"])
    ci = json.load(open(os.path.join(OUT, "cost_ind.json")))
    miss = set()
    for q in TARGETS:
        miss |= c[q] - samp
    pred = sum(max(1.0, 5e-6 * ci[str(t)][3]) for t in miss)
    with open(os.path.join(HERE, "sample_extra89.txt"), "w", encoding="utf-8") as f:
        f.write("# inner graphs added to complete the non-zero contributor sets of the outer graphs %s\n" % TARGETS)
        for t in sorted(miss):
            f.write("%d   # E=%d V=%d\n" % (t, G[t]["E"], G[t]["V"]))
    for q in TARGETS:
        print("outer graph q=%d (E=%d): %d inner graphs carry a non-zero piece there, %d of them were missing"
              % (q, G[q]["E"], len(c[q]), len(c[q] - samp)))
    print("wrote sample_extra89.txt: %d graphs, predicted %.0f CPU-s" % (len(miss), pred))


def mode_rows(prime, tag):
    p = PRIME[prime]
    G, c = contributors()
    rows = {}
    path = os.path.join(OUT, "%s.jsonl" % tag)
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        d = json.loads(line)
        for q, v in d["res"].items():
            a = rows.setdefault(int(q), [0] * len(v))
            for i, x in enumerate(v):
                a[i] = (a[i] + x) % p
    have = set()
    for line in open(path, encoding="utf-8"):
        if line.strip():
            have.add(json.loads(line)["t"])
    res = json.load(open(os.path.join(PKG, "results", "g5_%s_nat.json" % prime)))["table"]
    full = sorted(q for q in c if c[q] <= have)
    print("inner graphs run independently: %d; outer graphs whose non-zero contributors are all among them: %d"
          % (len(have), len(full)))
    bad = []
    for q in full:
        mine = [x % p for x in rows.get(q, [0] * 8)]
        theirs = [x % p for x in res[str(q)]] if str(q) in res else [0] * 8
        ok = mine == theirs
        print("   q=%4d E=%2d: %d non-zero contributors, independent sum %s cluster table row; row is %s"
              % (q, G[q]["E"], len(c[q]), "==" if ok else "!=", "zero (regular)" if not any(theirs) else "non-zero"))
        if not ok:
            bad.append(q)
    print("ROWS verdict: %s" % ("ALL EQUAL" if not bad else "FAILURES %s" % bad))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True)
    ap.add_argument("--prime", default="89")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    if a.mode == "extra":
        mode_extra()
    else:
        mode_rows(a.prime, a.tag or ("nat_%s" % a.prime))


if __name__ == "__main__":
    main()
