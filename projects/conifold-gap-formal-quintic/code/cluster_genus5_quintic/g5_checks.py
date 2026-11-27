# -*- coding: utf-8 -*-
"""g5_checks.py -- tags, task assignment, part files, merging and the checks of the collect step (genus5_run.py).

Part file format (parts/<tag>/<name>.json, rewritten atomically after every inner graph):
  {"tag": tag, "pieces": {inner index t: {outer index: [K residues]}}, "times": {t: seconds}}
Every inner graph is stored separately, so part files from different task splittings can be merged (an inner graph present twice must
have identical pieces and is counted once) and the per-inner-graph additivity (sum over the outer graphs of its pieces = its point-frame
term) can be checked against an 'id' run of the same prime and frame.
"""
import os, sys, json, re, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)


def tag_of(genus, prime, variant, ab="0,0", maxout=None):
    s = "g%d_%s_%s" % (genus, prime, variant)
    if ab != "0,0": s += "_ab" + ab.replace("/", "o").replace("-", "m").replace(",", "_")
    if maxout is not None: s += "_out%d" % maxout
    return s


def lpt_assign(costs, ntasks):
    """longest processing time first: list over tasks of lists of inner indices (deterministic)"""
    order = sorted(costs, key=lambda t: (-costs[t], t)); load = [0.0] * ntasks; out = [[] for _ in range(ntasks)]
    for t in order:
        i = min(range(ntasks), key=lambda j: (load[j], j)); out[i].append(t); load[i] += costs[t]
    return out, load


def costs_for(genus, ninner, variant):
    if genus == 5:
        d = json.load(open(os.path.join(HERE, "data", "cost_table_g5.json")))
        key = 'est_split' if variant == 'id' else 'est_dp'
        c = {int(t): max(v[key], 0.01) for t, v in d.items()}
        assert len(c) == ninner, "cost table has %d entries, graph list %d" % (len(c), ninner)
        return c
    return {t: 1.0 for t in range(ninner)}


def write_part(path, tag, pieces, times):
    tmp = path + ".tmp%d" % os.getpid()
    json.dump(dict(tag=tag, pieces=pieces, times=times), open(tmp, "w")); os.replace(tmp, path)


def read_parts(tag):
    d = os.path.join(HERE, "parts", tag)
    if not os.path.isdir(d): return [], {}, {}
    names = sorted(f for f in os.listdir(d) if re.match(r"^(task\d+of\d+|local)\.json$", f))
    pieces = {}; times = {}
    for f in names:
        part = json.load(open(os.path.join(d, f)))
        assert part['tag'] == tag, "part file %s has the foreign tag %s" % (f, part['tag'])
        for t, res in part['pieces'].items():
            if t in pieces:
                assert pieces[t] == res, "inner graph %s occurs in two part files with different content" % t
                continue
            pieces[t] = res; times[t] = part['times'].get(t, 0.0)
    return names, pieces, times


def merge(pieces, K, p):
    tab = {}
    for t, res in pieces.items():
        for q, c in res.items():
            a = tab.setdefault(int(q), [0] * K)
            for k in range(K): a[k] = (a[k] + c[k]) % p
    return tab


def c_conj(g):
    """|B_2g| / (2g (2g-2) 5^(g-1)): 1/1200, 1/25200, 1/180000 for g = 2, 3, 4 (known); the value at g = 5 is the pattern's prediction"""
    from fq_fp import bernoulli
    return abs(Fr(bernoulli(2 * g))) / (2 * g * (2 * g - 2) * 5 ** (g - 1))


def checks(sp, tab, covered, maxout, log, rec, id_pieces=None, pieces=None):
    """returns (FAILED, info dict)"""
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3
    R = sp.R; p = sp.p; g = sp.g; K = sp.K; npol = 2 * g - 2; nin = len(sp.inner_graphs); FAILED = []; info = {}
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    missing = sorted(set(range(nin)) - covered)
    check("all %d inner graphs are covered by the part files" % nin, not missing, "%d missing, e.g. %s" % (len(missing), missing[:20]))
    Eof = [G['E'] for G in sp.inner_graphs]; tc = tc_series(sp.C, R, K + 4); cg = c_conj(g)
    def polar(row):
        T = tc_expansion3(row, tc, R, g)[:npol]; assert len(T) == npol, "t_c expansion too short"
        return T
    def fmt(T): return " + ".join("(%s) t_c^%d" % (rec.fmt(T[k]), k - npol) for k in range(npol))
    if sp.variant != 'id':
        check("every entry of E_out = E_R - E_M is O(u)", min(sp.valO.values()) >= 1)
        check("[u^0] E_R(a, b) = Eh_0(a, b) for all (a, b), all colours", bool(sp.frozen_ok))
        V = polar(tab.get(0, [0] * K)); info['vertex'] = V
        log("   vertex graph in t_c: %s" % fmt(V))
        check("VERTEX GRAPH: (c) t_c^-%d + 0 + ... + 0 t_c^-1 (exact for any maxout)" % npol, not any(V[1:]))
        log("   vertex graph leading coefficient c_%d = %s; pattern |B_2g|/(2g(2g-2)5^(g-1)) = %s: %s" % (g, rec.fmt(V[0]), cg, "EQUAL" if V[0] == R.of(cg) else "DIFFERENT"))
        info['c_vertex_equals_pattern'] = V[0] == R.of(cg)
        lim = maxout if maxout is not None else 10 ** 9
        test = [i for i in range(1, nin) if Eof[i] <= lim]
        bad = {i: [k for k in range(npol) if tab.get(i, [0] * K)[k]] for i in test if any(tab.get(i, [0] * K)[:npol])}
        nnz = {}
        for i in test:
            if i in tab: nnz[Eof[i]] = nnz.get(Eof[i], 0) + 1
        log("   outer graphs tested for regularity: %d (E <= %s); with a stored row (some non-zero piece), by E: %s" % (len(test), "all" if maxout is None else maxout, sorted(nnz.items())))
        check("REGULARITY: every tested outer graph with edges is O(u^%d)" % npol, not bad, {i: (Eof[i], sp.inner_graphs[i]['gv'], ks) for i, ks in list(bad.items())[:30]})
        info['nonregular'] = {str(i): ks for i, ks in bad.items()}
        for i in sorted(bad)[:40]: log("      non-regular outer graph %4d (E=%d, gv=%s): %s" % (i, Eof[i], sp.inner_graphs[i]['gv'], [rec.fmt(x) for x in tab[i]]))
    if sp.variant == 'id' or maxout is None:
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(nin)) % p for k in range(K)]; info['total'] = tot
        T = polar(tot); info['total_tc'] = T
        log("   total over all graphs: %s" % [rec.fmt(x) for x in tot]); log("   total in t_c: %s" % fmt(T))
        check("GAP: total = (c) t_c^-%d + 0 + ... + 0 t_c^-1" % npol, not any(T[1:]))
        log("   c_%d = %s; pattern |B_2g|/(2g(2g-2)5^(g-1)) = %s: %s" % (g, rec.fmt(T[0]), cg, "EQUAL" if T[0] == R.of(cg) else "DIFFERENT"))
        info['c_total_equals_pattern'] = T[0] == R.of(cg)
    if id_pieces is not None and pieces is not None and maxout is None and sp.variant != 'id':
        bad = []
        for t in range(nin):
            s = [0] * K
            for q, c in pieces.get(str(t), {}).items():
                for k in range(K): s[k] = (s[k] + c[k]) % p
            ref = id_pieces.get(str(t), {}).get(str(t), [0] * K)
            if s != ref: bad.append(t)
        check("ADDITIVITY: for every inner graph the pieces add up to its point-frame term (from the id run)", not bad, bad[:20])
    return FAILED, info
