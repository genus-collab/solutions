# -*- coding: utf-8 -*-
"""ig_lib.py -- worker side: drive the INDEPENDENT re-implementation (natural_form/genus4_referee/indep, imported read-only and
never modified) at genus 5.

What comes from the independent folder: the hatted rows and all edge/leaf/vertex tables (ind_tables.Tables, built on the referee
rows ref_fp/ref_rows/ref_hat), the Witten-Kontsevich numbers (ind_wk), the automorphism count and the canonical form and the
contraction (ind_graphs), and the momentum-space split sum (ind_sum.Summer.graph) with its colour-space brute force
(.graph_brute).  No fq_* / fq3_* / fq4_* / g5_* module is imported anywhere in this folder.

What comes from the package: only the genera and the edge lists of the 4555 inner graphs (data/graphs_5_0_12.json) and the
position in that list, which is the outer index used for the comparison.  |Aut|, the canonical forms and the assignment
piece -> outer graph are recomputed here (ig_canon.py, and ind_graphs.contract at run time).

Sharding: ShardSummer splits the decoration sum by striding the option list of vertex 0.  graph() is
    f * Sum_{decorations} (...)
with f a constant, so the shard results add up to the unsharded one; this is used for the few graphs that are too slow in one
process (and is itself checked on smaller graphs by --check-shards in ig_run.py).
"""
import os, sys, json, time, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
IND = os.path.normpath(os.path.join(HERE, "..", "..", "genus4_referee", "indep"))
PKG = os.environ.get("G5PKG", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic"))
if IND not in sys.path:
    sys.path.insert(0, IND)

_W = {}


def keystr(key):
    gv, ed = key
    return "%s|%s" % (",".join(str(x) for x in gv), ",".join("%d-%d" % e for e in ed))


class CanonIndex(dict):
    """{own canonical form (as returned by ind_graphs.canon / contract) -> package list position}"""

    def __init__(self, smap):
        super().__init__()
        self.smap = smap

    def __missing__(self, key):
        i = self.smap[keystr(key)]
        self[key] = i
        return i


def raw_graphs():
    g = json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))
    return [(tuple(k[0]), [tuple(e) for e in k[1]]) for k, a in g]


def setup(prime, variant):
    key = (prime, variant)
    if _W.get("key") != key:
        from ind_tables import Tables
        from ind_sum import Summer
        from ind_graphs import make
        raw = raw_graphs()
        smap = json.load(open(os.path.join(OUT, "canon_map.json")))["map"]
        T = Tables(prime, 5, 8, variant)

        class ShardSummer(Summer):
            shard = None

            def options(self, G):
                o = Summer.options(self, G)
                if self.shard is not None:
                    k, n = self.shard
                    o = [o[0][k::n]] + o[1:]
                return o

        S = ShardSummer(T, [dict(key=None) for _ in raw])
        S.index = CanonIndex(smap)
        S.raw = raw
        S.make = make
        _W.clear()
        _W.update(key=key, T=T, S=S)
    return _W["S"]


def job(args):
    """args = (prime, variant, t, shard or None, brute) -> (t, shard, dict)"""
    prime, variant, t, shard, brute = args
    t0 = time.time()
    S = setup(prime, variant)
    gv, ed = S.raw[t]
    S.graphs[t] = S.make(gv, ed)
    G = S.graphs[t]
    S.shard = tuple(shard) if shard else None
    r = S.graph(t)
    t1 = time.time()
    p = S.T.p
    out = dict(res={str(q): [x % p for x in v] for q, v in r.items() if any(x % p for x in v)},
               aut=G["aut"], sec=t1 - t0, V=G["V"], E=G["E"])
    if brute:
        b = S.graph_brute(t)
        out["brute_equal"] = ({str(q): [x % p for x in v] for q, v in b.items() if any(x % p for x in v)} == out["res"])
        out["sec_brute"] = time.time() - t1
    S.shard = None
    return t, shard, out


def merge_shards(rows, prime):
    """rows: list of out-dicts of the same t (shards); returns the merged {q: [8]} as strings -> list"""
    p = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}[prime]
    acc = {}
    for o in rows:
        for q, v in o["res"].items():
            a = acc.setdefault(q, [0] * len(v))
            for i, x in enumerate(v):
                a[i] = (a[i] + x) % p
    return {q: v for q, v in acc.items() if any(v)}
