# -*- coding: utf-8 -*-
"""ig_index.py -- build an index of the genus-5 cluster evidence, for choosing the independent-verification sample.

Reads (read-only):
  <EV>/parts/g5_89_nat/task*of100.json      per-inner-graph pieces {t: {q: [8 residues]}} and wall times
  <EV>/parts/g5_61_nat/task*of20.json       the same at the second prime
  <EV>/parts/g5_89_id/local.json            the shipped point-frame pieces
  <PKG>/data/graphs_5_0_12.json             the package graph list [[gv, ed, legs], aut]
  <PKG>/data/cost_table_g5.json             est_split / est_dp per inner graph
Writes  out/index.json:
  graphs : per inner graph t -> dict(E, V, gv, aut, est_dp, est_split, t89, t61, npieces, qs)
  cancel : sorted list of outer graphs q with E(q) >= 1, E(q) <= 7 and at least one non-zero stored piece
           (these are exactly the outer graphs that are regular only by cancellation)
  zero_q : outer graphs with E <= 7 that receive no non-zero piece
  usage:  py -I -B ig_index.py
"""
import os, sys, json, glob

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
EV = os.environ.get("G5EV", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "g5ev"))
PKG = os.environ.get("G5PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic"))


def load_parts(tag):
    """{t: {q: [residues]}} and {t: seconds}, merged over the task files of one tag"""
    pieces = {}
    times = {}
    files = sorted(glob.glob(os.path.join(EV, "parts", tag, "task*of*.json")))
    if not files:
        files = sorted(glob.glob(os.path.join(EV, "parts", tag, "*.json")))
    for fn in files:
        d = json.load(open(fn))
        assert d["tag"] == tag, (fn, d["tag"])
        for t, fam in d["pieces"].items():
            assert int(t) not in pieces, "inner graph %s appears twice" % t
            pieces[int(t)] = {int(q): v for q, v in fam.items()}
        for t, s in d.get("times", {}).items():
            times[int(t)] = s
    return pieces, times, len(files)


def main():
    graphs = json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))
    cost = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))
    n = len(graphs)
    p89, t89, f89 = load_parts("g5_89_nat")
    p61, t61, f61 = load_parts("g5_61_nat")
    idp, _, fid = load_parts("g5_89_id")
    print("task files: g5_89_nat %d, g5_61_nat %d, g5_89_id %d" % (f89, f61, fid))
    print("inner graphs with a stored piece: 89 %d, 61 %d, id %d; graph list %d" % (len(p89), len(p61), len(idp), n))
    print("inner graphs with a time: 89 %d, 61 %d" % (len(t89), len(t61)))

    info = {}
    for t, (gk, aut) in enumerate(graphs):
        gv, ed, legs = gk
        assert not legs
        c = cost[str(t)]
        assert c["gv"] == gv and [list(e) for e in c["ed"]] == [list(e) for e in ed] and c["aut"] == aut, t
        fam = p89.get(t, {})
        info[t] = dict(E=len(ed), V=len(gv), gv=gv, aut=aut,
                       est_dp=c["est_dp"], est_split=c["est_split"],
                       t89=t89.get(t), t61=t61.get(t),
                       npieces=len(fam), qs=sorted(fam),
                       npieces61=len(p61.get(t, {})),
                       id_nonzero=bool(any(idp.get(t, {}).get(t, [0]))))
    # outer-graph totals from the pieces
    tot = {}
    for t, fam in p89.items():
        for q, v in fam.items():
            a = tot.setdefault(q, [0] * len(v))
            for i, x in enumerate(v):
                a[i] += x
    E = {t: info[t]["E"] for t in info}
    cancel = sorted(q for q in tot if E[q] >= 1 and E[q] <= 7)
    zero_q = sorted(q for q in range(n) if 1 <= E[q] <= 7 and q not in tot)
    print("outer graphs receiving a non-zero piece: %d (vertex graph included: %s)" % (len(tot), 0 in tot))
    print("  of these, with 1 <= E <= 7 (the cancellation set): %d" % len(cancel))
    print("  outer graphs with 1 <= E <= 7 and no non-zero piece: %d" % len(zero_q))
    print("  outer graphs with E >= 8: %d" % sum(1 for t in info if E[t] >= 8))
    byE = {}
    for t in info:
        byE.setdefault(info[t]["E"], []).append(t)
    print("inner graphs by E: " + ", ".join("E=%d:%d" % (e, len(byE[e])) for e in sorted(byE)))
    json.dump(dict(graphs={str(t): info[t] for t in info}, cancel=cancel, zero_q=zero_q),
              open(os.path.join(OUT, "index.json"), "w"))
    print("wrote %s" % os.path.join(OUT, "index.json"))


if __name__ == "__main__":
    main()
