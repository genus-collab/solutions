# -*- coding: utf-8 -*-
"""g5_calib4.py -- calibration of the sampling estimator on genus 4: sampled estimates (g5_sample_g4_89_nat_K6.jsonl) against the
measured per-inner-graph times of the real genus-4 run (parts/g4_89_nat_run.json, read only; 13 processes side by side).
Usage: py -I g5_calib4.py
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from g5_lib import fast_graphs, ROOT
from fq3_graphs import legged_graphs, canon

old = legged_graphs(4, 0, 9); old_pos = {G['key']: i for i, G in enumerate(old)}
new = fast_graphs(4, 0, 9); to_old = [old_pos[canon(*G['key'])] for G in new]
stored = json.load(open(os.path.join(ROOT, "parts", "g4_89_nat_run.json")))['times']
est = {}
for line in open(os.path.join(HERE, "g5_sample_g4_89_nat_K6.jsonl")):
    if line.strip(): r = json.loads(line); est[r['t']] = r
assert len(est) == 379
tot_e = sum(r['est'] for r in est.values()); tot_m = sum(stored.values())
print("genus 4: sum of sampled estimates %.0f s; measured CPU of the real run %.0f s; ratio measured/estimate %.2f" % (tot_e, tot_m, tot_m / tot_e))
byE = {}
for t, r in est.items():
    a = byE.setdefault(r['E'], [0.0, 0.0, 0]); a[0] += r['est']; a[1] += stored[str(to_old[t])]; a[2] += 1
for E in sorted(byE): a = byE[E]; print("   E = %d: %3d graphs, estimate %7.0f s, measured %7.0f s, ratio %.2f" % (E, a[2], a[0], a[1], a[1] / max(a[0], 1e-9)))
big = [t for t in est if stored[str(to_old[t])] > 5.0]
lr = [math.log(stored[str(to_old[t])] / est[t]['est']) for t in big]
m = sum(lr) / len(lr); sd = (sum((x - m) ** 2 for x in lr) / len(lr)) ** 0.5
print("   graphs with measured time > 5 s: %d; log(measured/estimate): mean %.2f (factor %.2f), s.d. %.2f (factor %.2f), min factor %.2f, max factor %.2f"
      % (len(big), m, math.exp(m), sd, math.exp(sd), math.exp(min(lr)), math.exp(max(lr))))
top = sorted(est, key=lambda t: -stored[str(to_old[t])])[:10]
for t in top:
    r = est[t]; print("   t=%3d E=%d V=%d val=%s: measured %.1f s, estimate %.1f s (ncombo %d, nsub %d, mean %.2f ms, %d samples)"
                      % (t, r['E'], r['V'], r['val'], stored[str(to_old[t])], r['est'], r['ncombo'], r['nsub'], 1e3 * r['mean'], r['samples']))
