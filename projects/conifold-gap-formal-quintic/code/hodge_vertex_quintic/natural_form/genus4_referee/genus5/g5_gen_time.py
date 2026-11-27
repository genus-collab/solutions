# -*- coding: utf-8 -*-
"""g5_gen_time.py -- time the ORIGINAL generator fq3_graphs (canon over all V! vertex permutations) on the leg-free genus-5 graphs.

Same level-by-level loop as fq3_graphs.legged_keys(5, 0, 12), with a wall-clock deadline (default 24 minutes) checked after every
canonical form; writes the counts per number of edges reached so far to g5_gen_time.log (this folder).  Nothing else is written.
Usage:  py -I g5_gen_time.py [genus] [deadline_seconds]
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
from fq3_graphs import canon, degenerations

g = int(sys.argv[1]) if len(sys.argv) > 1 else 5
DEADLINE = float(sys.argv[2]) if len(sys.argv) > 2 else 24 * 60.0
LOG = open(os.path.join(HERE, "g5_gen_time.log"), "a", encoding="utf-8")


def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


t0 = time.time(); Emax = 3 * g - 3
log("=" * 80); log("original generator, genus %d, Emax %d, deadline %.0f s" % (g, Emax, DEADLINE))
levels = [[canon((g,), [], ())]]; stopped = False
for E in range(1, Emax + 1):
    t1 = time.time(); new = set(); ncanon = 0; nsrc = 0
    for key in levels[-1]:
        for H in degenerations(key[0], list(key[1]), key[2]):
            new.add(canon(*H)); ncanon += 1
            if time.time() - t0 > DEADLINE: stopped = True; break
        if stopped: break
        nsrc += 1
    if stopped:
        log("   E = %2d: STOPPED at the deadline after %d of %d source graphs (%d canonical forms, %d distinct so far, %.0f s in this level)"
            % (E, nsrc, len(levels[-1]), ncanon, len(new), time.time() - t1))
        break
    levels.append(sorted(new))
    log("   E = %2d: %5d graphs (%d canonical forms, %.1f s; cumulative %.1f s)" % (E, len(new), ncanon, time.time() - t1, time.time() - t0))
log("counts per E (complete levels): %s; total %d; wall %.1f s; %s" % ([len(l) for l in levels], sum(len(l) for l in levels), time.time() - t0,
    "COMPLETE" if not stopped else "INCOMPLETE (deadline)"))
