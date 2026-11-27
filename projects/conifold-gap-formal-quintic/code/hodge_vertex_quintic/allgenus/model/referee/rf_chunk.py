# -*- coding: utf-8 -*-
"""rf_chunk.py -- time-bounded amplitude runs for large (g, n): the graph list is generated once (fq3_graphs, read only) and stored in
the scratch folder; graph index ranges are processed separately and merged into the same cache file rf_gsum.amp reads.
  py -I -B rf_chunk.py gen g n             -> scratch/rfgraphs_g_n.json
  py -I -B rf_chunk.py run g n p i0 i1     -> scratch/rfpart_p_g_n_i0_i1.json   (p = 61 or 89; stops cleanly after 25 min, reports where)
  py -I -B rf_chunk.py merge g n p         -> scratch/rfamp_p_1_g_n.json (checks that the parts cover all graphs exactly once)
"""
import glob, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as rc
import rf_gsum as gsm


def gpath(g, n): return os.path.join(gsm.SCR, "rfgraphs_%d_%d.json" % (g, n))


def main():
    cmd = sys.argv[1]; g = int(sys.argv[2]); n = int(sys.argv[3])
    if cmd == "gen":
        t0 = time.time(); L = gsm.legged_graphs(g, n, 3 * g - 3 + n)
        keep = [dict(V=G['V'], gv=list(G['gv']), ed=[list(e) for e in G['ed']], hsAt=G['hsAt'], legsAt=G['legsAt'], aut=G['aut']) for G in L]
        with open(gpath(g, n), "w") as f: json.dump(keep, f)
        print("gen (%d,%d): %d graphs, %.0f s" % (g, n, len(keep), time.time() - t0)); return
    pbits = int(sys.argv[4]); p = rc.P61 if pbits == 61 else rc.P89
    if cmd == "run":
        i0 = int(sys.argv[5]); i1 = int(sys.argv[6])
        with open(gpath(g, n)) as f: L = json.load(f)
        m = rc.Model(p, 16); gs = gsm.GSum(m); acc = {}; t0 = time.time(); i = i0
        for i in range(i0, min(i1, len(L))):
            Gm = L[i]; Gm['ed'] = [tuple(e) for e in Gm['ed']]
            gs.graph_amp(Gm, n, acc)
            if time.time() - t0 > 1500: i += 1; break
        else: i = min(i1, len(L))
        out = os.path.join(gsm.SCR, "rfpart_%d_%d_%d_%d_%d.json" % (pbits, g, n, i0, i))
        with open(out, "w") as f: json.dump([[list(k[0]), list(k[1]), v] for k, v in acc.items() if v], f)
        print("run (%d,%d) p=2^%d-1 graphs [%d, %d) done, %.0f s" % (g, n, pbits, i0, i, time.time() - t0)); return
    if cmd == "merge":
        with open(gpath(g, n)) as f: NG = len(json.load(f))
        parts = sorted(glob.glob(os.path.join(gsm.SCR, "rfpart_%d_%d_%d_*.json" % (pbits, g, n))))
        rng = sorted(tuple(map(int, os.path.basename(x)[:-5].split("_")[-2:])) for x in parts)
        cov = [0] * NG
        for a, b in rng:
            for i in range(a, b): cov[i] += 1
        assert all(c == 1 for c in cov), "parts do not cover the graphs exactly once"
        A = {}
        for x in parts:
            with open(x) as f:
                for c, e, v in json.load(f):
                    k = (tuple(c), tuple(e)); A[k] = (A.get(k, 0) + v) % p
        out = os.path.join(gsm.SCR, "rfamp_%d_1_%d_%d.json" % (pbits, g, n))
        with open(out, "w") as f: json.dump([[list(k[0]), list(k[1]), v] for k, v in A.items() if v], f)
        print("merge (%d,%d) p=2^%d-1: %d parts, %d graphs, %d entries" % (g, n, pbits, len(parts), NG, len(A)))


if __name__ == "__main__":
    main()
