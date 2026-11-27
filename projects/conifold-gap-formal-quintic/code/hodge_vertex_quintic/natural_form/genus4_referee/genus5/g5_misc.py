# -*- coding: utf-8 -*-
"""g5_misc.py -- small genus-5 checks for GENUS5_ESTIMATE.md (one process, a few minutes).

1. set-up of Split(p, 5, (0, 0), 8, 'nat') at the second prime 2^61-1 (time, memory, asserts, frozen_ok).
2. the collect step of genus4_split.py at genus 5: fq_hat.tc_series and fq3_msum.tc_expansion3 run with K = 8, g = 5 on a synthetic
   total  u^-8 G(u) = c t_c(u)^-8  (so the routine must return (c, 0, .., 0)); c = 1/660000 is the value that the genus 2, 3, 4
   constants (1/1200, 1/25200, 1/180000) = |B_2g| / (2g (2g-2) 5^(g-1)) suggest for g = 5 (a CONJECTURAL extrapolation, not computed).
3. time of the graph dictionaries of the ORIGINAL code (fq3_graphs._graph_dict -> fq_graphs.autos_of) for all 4555 genus-5 graphs.
4. time of the ORIGINAL canonical form on the contracted graphs of one E = 12, V = 8 inner graph (all 3302 subsets), against canon_fast.
Log: g5_misc.log.   Usage: py -I g5_misc.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g5_lib
from g5_lib import make_split, fast_graphs, canon_fast
from fq_fp import P1, P2

LOG = open(os.path.join(HERE, "g5_misc.log"), "a", encoding="utf-8")


def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


log("=" * 100); log("g5_misc, python %s" % sys.version.split()[0])
t0 = time.time(); sp = make_split(P1, 5, (0, 0), 8, 'nat')
log("1. PASS genus-5 set-up at 2^61-1: %.1f s, all asserts of Split.__init__ passed, frozen_ok = %s, valuations of E_out: %s"
    % (time.time() - t0, sp.frozen_ok, sorted(set(sp.valO.values()))))

from fq_hat import tc_series
from fq3_msum import tc_expansion3
from fq_series import Ser
R = sp.R; K = sp.K; g = 5; tc = tc_series(sp.C, R, K + 4); c5 = R.of(Fr(1, 660000))
# synthetic total: G(u) = c * (u / t_c(u))^8, K coefficients
ratio = tc.inv().pow(8)                      # t_c^-8 = u^-8 * (series)
coefs = [ratio.coef(k - 8) * c5 % sp.p for k in range(K)]
T = tc_expansion3(coefs, tc, R, g)[:8]
ok = (T[0] == c5 and not any(T[1:8]))
log("2. %s collect step at genus 5: tc_expansion3 with K = 8 returns (c, 0, .., 0) for the synthetic total c t_c^-8 (%d coefficients returned)" % ("PASS" if ok else "FAIL", len(T)))

import fq3_graphs
G5 = fast_graphs(5, 0, 12); t0 = time.time(); bad = 0; slow = (0.0, None)
for G in G5:
    t1 = time.time(); d = fq3_graphs._graph_dict(G['key'], 0); dt = time.time() - t1
    if d['aut'] != G['aut']: bad += 1
    if dt > slow[0]: slow = (dt, (G['E'], G['V']))
log("3. %s original fq3_graphs._graph_dict (autos_of) on all %d genus-5 graphs: %.0f s, slowest %.1f s at (E, V) = %s; |Aut| equal to the fast values (%d mismatches)"
    % ("PASS" if bad == 0 else "FAIL", len(G5), time.time() - t0, slow[0], slow[1], bad))

t = max(range(len(G5)), key=lambda i: (G5[i]['E'], G5[i]['V'])); G = G5[t]; E = G['E']
import fq4_split
subs = [S for S in range(1 << E) if bin(S).count("1") < 8]
fq4_split.canon = canon_fast; t0 = time.time(); a = [sp.contract(G, S) for S in subs]; tf = time.time() - t0
fq4_split.canon = fq3_graphs.canon; t0 = time.time(); b = [sp.contract(G, S) for S in subs]; to = time.time() - t0
fq4_split.canon = canon_fast
same = len(set(a)) == len(set(b)) and all(fq3_graphs.canon(*k) == kb for k, kb in zip(a, b))
log("4. %s contraction classes of the %d subsets of inner graph %d (E = %d, V = %d): original canon %.1f s, canon_fast %.2f s, %d outer graphs, same classes"
    % ("PASS" if same else "FAIL", len(subs), t, E, G['V'], to, tf, len(set(a))))
