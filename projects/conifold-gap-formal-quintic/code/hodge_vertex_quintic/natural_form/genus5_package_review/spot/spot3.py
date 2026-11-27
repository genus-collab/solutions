# -*- coding: utf-8 -*-
"""spot3.py -- genus-5 pieces from the INDEPENDENT re-implementation of the genus-4 referee (natural_form/genus4_referee/indep, imported
read-only, not modified) on a few small inner graphs, compared with the package's graph_dp pieces (out/spot1_89.json).

The independent code is genus-general except for its graph source (Mathematica StableGraph-g.bin, absent for g = 5) and its driver
(the gap constants stop at genus 4).  Here: the inner graphs are taken from pkg/data/graphs_5_0_12.json (genera and edge list only) and
rebuilt with ind_graphs.make (own |Aut| by depth-first search, own brute-force canonical form); contractions use ind_graphs.contract;
the outer index is found by matching the own canonical form against the package list (own canon of the candidates with the same
number of edges and the same multiset of genera).  Tables: ind_tables.Tables(prime, 5, K = 8, 'nat') (referee hatted rows, own E_R
division, own WK numbers, own Fourier transforms); sum: ind_sum.Summer.graph (momentum space).
  py -I -B spot3.py [--prime 89] [--list 1,3,..] [--nproc 6] [--brute 1,3]
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.join(HERE, "pkg"); OUT = os.path.join(HERE, "out")
IND = os.path.normpath(os.path.join(HERE, "..", "..", "genus4_referee", "indep"))
sys.path.insert(0, IND)
LIST = [1, 3, 5, 25, 28, 118, 383, 444, 1234]
_W = {}


class LazyIndex(dict):
    def __init__(self, raw):
        super().__init__(); self.raw = raw; self.cache = {}
    def __missing__(self, key):
        from ind_graphs import canon
        gv, ed = key; sig = (len(ed), tuple(sorted(gv)))
        for i, (g2, e2) in enumerate(self.raw):
            if (len(e2), tuple(sorted(g2))) != sig: continue
            if i not in self.cache: self.cache[i] = canon(tuple(g2), [tuple(e) for e in e2])
            if self.cache[i] == key: self[key] = i; return i
        raise KeyError(key)


def setup(prime):
    if _W.get('prime') != prime:
        from ind_tables import Tables
        from ind_sum import Summer
        from ind_graphs import make
        raw = [(tuple(k[0]), [tuple(e) for e in k[1]]) for k, a in json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))]
        T = Tables(prime, 5, 8, "nat")
        L = [dict(key=None) for _ in raw]
        S = Summer(T, L); S.index = LazyIndex(raw); S.raw = raw; S.make = make
        _W.update(prime=prime, T=T, S=S)
    return _W['S']


def job(args):
    prime, t, brute = args; t0 = time.time(); S = setup(prime)
    gv, ed = S.raw[t]; S.graphs[t] = S.make(gv, ed); G = S.graphs[t]
    r = S.graph(t); t1 = time.time()
    out = dict(res={str(q): [x % S.T.p for x in v] for q, v in r.items() if any(x % S.T.p for x in v)}, aut=G['aut'], t=t1 - t0)
    if brute:
        b = S.graph_brute(t); out['brute_equal'] = ({q: [x % S.T.p for x in v] for q, v in b.items() if any(x % S.T.p for x in v)} ==
                                                   {int(q): v for q, v in out['res'].items()}); out['tb'] = time.time() - t1
    return t, out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--list", default="")
    ap.add_argument("--nproc", type=int, default=6); ap.add_argument("--brute", default="1,3,5,25,28")
    a = ap.parse_args(); assert a.nproc <= 6
    lst = [int(x) for x in a.list.split(",")] if a.list else LIST; brute = set(int(x) for x in a.brute.split(",") if x)
    LOG = open(os.path.join(OUT, "spot3_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    T0 = time.time(); log("=" * 100)
    ref = json.load(open(os.path.join(OUT, "spot1_%s.json" % a.prime)))['full']
    pk = json.load(open(os.path.join(PKG, "data", "graphs_5_0_12.json")))
    import multiprocessing as mp
    res = {}
    with mp.Pool(processes=min(a.nproc, len(lst))) as pool:
        for t, o in pool.imap_unordered(job, [(a.prime, t, t in brute) for t in sorted(lst, reverse=True)]):
            res[t] = o; log("   done %d (%.0f s wall)" % (t, time.time() - T0))
    FAIL = []
    for t in lst:
        o = res[t]; mine = {int(q): v for q, v in o['res'].items()}; theirs = {int(q): v for q, v in ref[str(t)]['A'].items()}
        eq = (mine == theirs); autok = (o['aut'] == pk[t][1])
        log("   t=%4d gv=%s ed=%s: independent %.1f s, %d non-zero outer rows, package %d; EQUAL: %s; |Aut| own %d package %d: %s%s"
            % (t, pk[t][0][0], pk[t][0][1], o['t'], len(mine), len(theirs), eq, o['aut'], pk[t][1], autok,
               "" if 'brute_equal' not in o else "; momentum == colour-space brute force: %s (%.1f s)" % (o['brute_equal'], o['tb'])))
        if not eq or not autok or o.get('brute_equal') is False: FAIL.append(t)
    log("SPOT3 (%s) %s, wall %.0f s" % (a.prime, "ALL EQUAL" if not FAIL else "FAILURES %s" % FAIL, time.time() - T0)); LOG.close()


if __name__ == "__main__":
    main()
