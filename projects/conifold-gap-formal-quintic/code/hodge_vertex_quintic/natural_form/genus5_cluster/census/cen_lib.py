# -*- coding: utf-8 -*-
"""cen_lib.py -- independent reader/merger for the genus-5 cluster evidence (part files, result files, logs).

Nothing here imports the g5_* engine: everything is arithmetic on the stored residues plus my own re-implementation of the
longest-processing-time task assignment.  Written for the census of 2026-11-13; new file, nothing existing is touched.

Part file   : {"tag", "pieces": {inner t: {outer q: [K residues]}}, "times": {t: seconds}}
Result file : {"tag", "table": {outer q: [K residues]}, "info", "failed", "outer": [{E, gv, ed, aut} per graph]}
"""
import os, json, re

K = 8
NIN = 4555
PRIME = {"89": 2 ** 89 - 1, "61": 2 ** 61 - 1}
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
EV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "g5ev")
TAGS = ["g5_89_nat", "g5_61_nat", "g5_89_nat_ab1o3_m2o7", "g5_89_M0_out2", "g5_89_nat_t7_out2"]
IDTAGS = ["g5_89_id", "g5_61_id"]


def prime_of(tag):
    return PRIME[tag.split("_")[1]]


def partdir(tag):
    return os.path.join(EV, "parts", tag)


def logdir(tag):
    return os.path.join(EV, "logs", tag)


def part_files(tag):
    d = partdir(tag)
    names = [f for f in os.listdir(d) if f.endswith(".json")]
    def key(f):
        m = re.match(r"^task(\d+)of(\d+)\.json$", f)
        return (0, int(m.group(1))) if m else (1, f)
    return sorted(names, key=key)


def load_part(tag, name):
    with open(os.path.join(partdir(tag), name)) as f:
        d = json.load(f)
    assert d["tag"] == tag, "part file %s carries the tag %s" % (name, d["tag"])
    return d


def load_result(tag):
    with open(os.path.join(PKG, "results", "%s.json" % tag)) as f:
        return json.load(f)


def graph_meta():
    """per inner/outer graph index: E, gv, aut (identical in every result file -- checked by cen_merge)"""
    return load_result("g5_89_nat")["outer"]


def cost_table():
    with open(os.path.join(PKG, "data", "cost_table_g5.json")) as f:
        d = json.load(f)
    return {int(t): v for t, v in d.items()}


def my_lpt(costs, ntasks):
    """longest-processing-time-first, re-implemented from scratch: returns (slices, loads).

    Rule: take the inner graphs in order of decreasing cost (ties by increasing index) and give each to the task with the
    smallest current load (ties by smallest task index).  Costs are max(est, 0.01) of the column the variant selects.
    """
    order = sorted(costs.keys(), key=lambda t: (-costs[t], t))
    load = [0.0] * ntasks
    out = [[] for _ in range(ntasks)]
    for t in order:
        best = 0
        for j in range(1, ntasks):
            if load[j] < load[best]:
                best = j
        out[best].append(t)
        load[best] += costs[t]
    return out, load


def costs_for(variant):
    key = "est_split" if variant == "id" else "est_dp"
    c = cost_table()
    return {t: max(v[key], 0.01) for t, v in c.items()}


def variant_of(tag):
    """the --variant string the tag was produced with (tag = g5_<prime>_<variant>[_ab..][_out..])"""
    s = tag[len("g5_xx_"):]
    s = re.sub(r"_ab[0-9om_]+$", "", re.sub(r"_out\d+$", "", s))
    return s


def merge_pieces(pieces, p):
    """sum the stored pieces into the outer-graph table, modulo p (my own merge, independent of g5_checks.merge).

    Every outer graph that receives at least one stored piece gets a row, even when the sum is zero: those all-zero rows are
    exactly the cancellations, and the stored result tables keep them too.
    """
    tab = {}
    for t, res in pieces.items():
        for q, c in res.items():
            row = tab.get(int(q))
            if row is None:
                row = tab[int(q)] = [0] * K
            for k in range(K):
                row[k] = (row[k] + c[k]) % p
    return tab


def load_all_pieces(tag):
    """union of the part files: {inner t (int) -> {outer q (int) -> [K residues]}}, plus duplicate/omission bookkeeping"""
    pieces = {}
    owner = {}
    dup = []
    conflict = []
    times = {}
    per_file = {}
    for name in part_files(tag):
        d = load_part(tag, name)
        here = sorted(int(t) for t in d["pieces"])
        per_file[name] = dict(inner=here, cpu=sum(float(v) for v in d["times"].values()),
                              times={int(t): float(v) for t, v in d["times"].items()})
        for t, res in d["pieces"].items():
            ti = int(t)
            r = {int(q): list(v) for q, v in res.items()}
            if ti in pieces:
                dup.append((ti, owner[ti], name))
                if pieces[ti] != r:
                    conflict.append((ti, owner[ti], name))
                continue
            pieces[ti] = r
            owner[ti] = name
            times[ti] = float(d["times"].get(t, 0.0))
    return dict(pieces=pieces, owner=owner, dup=dup, conflict=conflict, times=times, per_file=per_file)
