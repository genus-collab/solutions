# -*- coding: utf-8 -*-
"""ms_chunk.py -- the amplitude AMP(g, n) of ms_engine computed in chunks of graphs (graphs from ms_graphs.legged_graphs2), so that no
single run exceeds the time budget; the merged result is written in ms_engine's cache format (scratch folder).
  py -I -B ms_chunk.py gen g n                 graph list -> scratch/graphs_g_n.json (with |Aut|)
  py -I -B ms_chunk.py run g n i0 i1           graphs [i0, i1), stopping after 1200 s -> scratch/part_g_n_i0_iend.json
  py -I -B ms_chunk.py drive g n [i0] [budget] the same in ONE process (not used for long jobs: a shell loop over `run` keeps each
                                                process below the time budget; `run` writes scratch/next_g_n.txt = 'end N')
  py -I -B ms_chunk.py merge g n [--compare]   check that the parts tile 0..N exactly once, sum them, write scratch/amp_g_n.json
                                                (--compare: compare with an existing amp_g_n.json instead of writing; regression)
"""
import json, os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr
import ms_engine as me
import ms_graphs as mg

SCR = mr.SCRATCH


def gpath(g, n): return os.path.join(SCR, "graphs_%d_%d.json" % (g, n))


def gen(g, n):
    t0 = time.time(); L = mg.legged_graphs2(g, n, 3 * g - 3 + n)
    with open(gpath(g, n), "w") as f:
        json.dump([dict(gv=list(x['gv']), ed=[list(e) for e in x['ed']], legs=list(x['legs']), aut=x['aut']) for x in L], f)
    print("gen (%d,%d): %d graphs, %.1f s" % (g, n, len(L), time.time() - t0), flush=True)


def load(g, n):
    with open(gpath(g, n)) as f: raw = json.load(f)
    out = []
    for x in raw:
        gv = tuple(x['gv']); ed = [tuple(e) for e in x['ed']]; legs = tuple(x['legs']); V = len(gv)
        H = [(e, s) for e in range(len(ed)) for s in (0, 1)]; hvtx = [ed[e][s] for e, s in H]
        out.append(dict(V=V, gv=gv, ed=ed, legs=legs, aut=x['aut'], H=H, E=len(ed),
                        hsAt=[[h for h, y in enumerate(hvtx) if y == v] for v in range(V)],
                        legsAt=[[l for l in range(n) if legs[l] == v] for v in range(V)]))
    return out


def enc(A):
    return {",".join(map(str, k[0])) + "|" + ",".join(map(str, k[1])): str(v) for k, v in A.items() if v}


def run(g, n, i0, i1, budget=1200.0):
    """graphs [i0, i1), stopping early (after a whole graph) once `budget` seconds are used; returns the end index reached"""
    L = load(g, n); t0 = time.time(); acc = {}; i = i0; i1 = min(i1, len(L))
    while i < i1 and time.time() - t0 < budget:
        me.graph_amp(L[i], n, acc); i += 1
    out = os.path.join(SCR, "part_%d_%d_%d_%d.json" % (g, n, i0, i))
    with open(out, "w") as f: json.dump(dict(i0=i0, i1=i, N=len(L), acc=enc(acc)), f)
    print("run (%d,%d) graphs [%d, %d): %d entries, %.1f s" % (g, n, i0, i, len(acc), time.time() - t0), flush=True)
    with open(os.path.join(SCR, "next_%d_%d.txt" % (g, n)), "w") as f: f.write("%d %d" % (i, len(L)) + chr(10))
    return i


def drive(g, n, i0, budget=1200.0):
    """one process: successive time-bounded parts from i0 to the end (each part < budget + one graph)"""
    N = len(load(g, n)); i = i0
    while i < N: i = run(g, n, i, N, budget)
    print("DRIVE (%d,%d) DONE at %d" % (g, n, i), flush=True)


def merge(g, n, compare=False):
    parts = []
    for fn in os.listdir(SCR):
        if fn.startswith("part_%d_%d_" % (g, n)) and fn.endswith(".json"):
            with open(os.path.join(SCR, fn)) as f: parts.append(json.load(f))
    parts.sort(key=lambda p: p['i0']); N = parts[0]['N'] if parts else -1
    pos = 0
    for p in parts:
        assert p['N'] == N and p['i0'] == pos, "parts do not tile: %s at %d" % (p['i0'], pos); pos = p['i1']
    assert pos == N, "parts cover %d of %d graphs" % (pos, N)
    tot = {}
    for p in parts:
        for k, v in p['acc'].items(): tot[k] = tot.get(k, Fr(0)) + Fr(v)
    tot = {k: str(v) for k, v in tot.items() if v}
    path = me._cache_path(g, n)
    if compare:
        with open(path) as f: ref = json.load(f)
        same = {k: Fr(v) for k, v in ref.items()} == {k: Fr(v) for k, v in tot.items()}
        print("merge (%d,%d): %d parts, %d graphs, %d entries; equal to the existing %s: %s" % (g, n, len(parts), N, len(tot), os.path.basename(path), same))
    else:
        assert not os.path.exists(path), "refusing to overwrite " + path
        with open(path, "w") as f: json.dump(tot, f)
        print("merge (%d,%d): %d parts, %d graphs, %d entries -> %s" % (g, n, len(parts), N, len(tot), path))


if __name__ == "__main__":
    a = sys.argv[1]; g, n = int(sys.argv[2]), int(sys.argv[3])
    if a == "gen": gen(g, n)
    elif a == "run": run(g, n, int(sys.argv[4]), int(sys.argv[5]))
    elif a == "drive": drive(g, n, int(sys.argv[4]) if len(sys.argv) > 4 else 0, float(sys.argv[5]) if len(sys.argv) > 5 else 1200.0)
    elif a == "merge": merge(g, n, "--compare" in sys.argv)
