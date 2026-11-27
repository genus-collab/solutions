# -*- coding: utf-8 -*-
"""g5_gen_check.py -- validate the fast canonical form of g5_lib against the original generator, then time genus 5.

1. genus 2, 3, 4 (leg-free) and a few legged types: the fast list, mapped through the ORIGINAL fq3_graphs.canon, equals the list of
   fq3_graphs.legged_graphs as a set of keys, with the same |Aut| graph by graph (genus 4: 379 graphs).
2. genus 5: generate with the fast canonical form (timed, per level), compare the total with the known 4555, and recompute |Aut| of
   every graph with the brute-force fq_graphs.aut_count (independent of the fast search).
Log: g5_gen_check.log (this folder).   Usage: py -I g5_gen_check.py
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g5_lib
from g5_lib import fast_levels, fast_graphs, canon_fast
from fq3_graphs import legged_graphs, canon
from fq_graphs import aut_count

LOG = open(os.path.join(HERE, "g5_gen_check.log"), "a", encoding="utf-8"); FAILED = []


def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


log("=" * 100); log("g5_gen_check")
for (g, n, Emax) in [(2, 0, 3), (3, 0, 6), (4, 0, 9), (2, 2, 5), (1, 4, 4), (3, 1, 5), (0, 6, 3)]:
    t0 = time.time(); old = legged_graphs(g, n, Emax); t_old = time.time() - t0
    t0 = time.time(); lev = fast_levels(g, n, Emax); t_new = time.time() - t0
    new = [x for l in lev for x in l]
    d_old = {G['key']: G['aut'] for G in old}; d_new = {canon(*k): a for k, a in new}
    check("type (%d, %d), E <= %d: %d graphs; fast list = original list (as canonical forms of fq3_graphs.canon) with equal |Aut|; original %.1f s, fast %.2f s"
          % (g, n, Emax, len(old), t_old, t_new), d_old == d_new and len(new) == len(old), (len(old), len(new)))
    # the fast canonical form is constant on isomorphism classes: relabel every graph by a few random permutations
    import random
    rng = random.Random(1); bad = 0
    for k, a in new:
        V = len(k[0])
        for _ in range(3):
            s = list(range(V)); rng.shuffle(s); gv2 = [0] * V
            for v in range(V): gv2[s[v]] = k[0][v]
            ed2 = [(s[u], s[w]) for u, w in k[1]]; rng.shuffle(ed2)
            if canon_fast(tuple(gv2), ed2, tuple(s[l] for l in k[2])) != k: bad += 1
    check("type (%d, %d): canon_fast is invariant under random relabellings (3 per graph)" % (g, n), bad == 0, bad)

t0 = time.time(); G5 = fast_graphs(5, 0, 12, log); t5 = time.time() - t0
byE = {}
for G in G5: byE[G['E']] = byE.get(G['E'], 0) + 1
log("genus 5, leg-free: %d graphs in %.1f s; per number of edges %s" % (len(G5), t5, [byE.get(E, 0) for E in range(13)]))
check("genus 5: total number of leg-free stable graphs = 4555 (known value)", len(G5) == 4555, len(G5))
byV = {}
for G in G5: byV[G['V']] = byV.get(G['V'], 0) + 1
log("   per number of vertices: %s" % sorted(byV.items()))
t0 = time.time(); bad = []
for i, G in enumerate(G5):
    V = G['V']; A = [[0] * V for _ in range(V)]
    for u, w in G['ed']:
        if u == w: A[u][u] += 1
        else: A[u][w] += 1; A[w][u] += 1
    if aut_count(G['gv'], A) != G['aut']: bad.append(i)
check("genus 5: |Aut| of all %d graphs equals the brute-force fq_graphs.aut_count (%.0f s)" % (len(G5), time.time() - t0), not bad, bad[:10])
# genus conservation and stability
ok = all(sum(G['gv']) + G['E'] - G['V'] + 1 == 5 and all(2 * G['gv'][v] - 2 + len(G['hsAt'][v]) > 0 for v in range(G['V'])) for G in G5)
check("genus 5: every graph is stable of genus 5", ok)
check("genus 5: keys pairwise distinct", len(set(G['key'] for G in G5)) == len(G5))
log("G5_GEN_CHECK %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
