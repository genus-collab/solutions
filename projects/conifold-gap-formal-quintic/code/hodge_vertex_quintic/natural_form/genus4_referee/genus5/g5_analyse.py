# -*- coding: utf-8 -*-
"""g5_analyse.py -- turn the genus-5 cost samples into the numbers of GENUS5_ESTIMATE.md (no heavy computation).

Inputs (this folder): g5_sample_g5_89_nat_K8.jsonl (all 4555 inner graphs, Split.graph cost sampled), g5_resample_* (the 150 heaviest,
independent larger samples: these replace the first-pass values), g5_full_* (graphs run completely), g5_dpsample_* / g5_dpresample_*
(the subset-DP prototype), if present.  Output: g5_analyse.log and g5_cost_table.json (per graph: gv, edges, estimates in seconds).
Task-balancing simulation: N array tasks filled (a) as genus4_split.py does now (round robin over its cost proxy), (b) round robin
over a structural proxy  ncombo * nsub * Sum_steps 5^(frontier+1) (1 + edges closed),  (c) longest-first into the least loaded task
using the sampled estimates; reported: the heaviest task in hours of the reference machine.
Usage: py -I g5_analyse.py
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from g5_lib import fast_graphs
import fq4_split

LOG = open(os.path.join(HERE, "g5_analyse.log"), "w", encoding="utf-8")


def log(s=""):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def load(name):
    path = os.path.join(HERE, name); out = {}
    if os.path.exists(path):
        for line in open(path):
            if line.strip(): r = json.loads(line); out[r['t']] = r
    return out


G5 = fast_graphs(5, 0, 12); n = len(G5); H = 3600.0
S = load("g5_sample_g5_89_nat_K8.jsonl"); RS = load("g5_resample_g5_89_nat_K8.jsonl"); F = load("g5_full_g5_89_nat_K8.jsonl")
D = load("g5_dpsample_g5_89_nat_K8.jsonl"); DR = load("g5_dpresample_g5_89_nat_K8.jsonl"); DF = load("g5_dpfull_g5_89_nat_K8.jsonl")
assert len(S) == n, (len(S), n)
est = {t: (RS[t]['est'] if t in RS else S[t]['est']) for t in S}
log("GENUS 5, Split.graph as it stands (2^89-1, natural M, K = 8), reference machine = this workstation, one core, CPython %s" % sys.version.split()[0])
log("   graphs: %d; sampled (decoration, subset) pairs: %d first pass + %d second pass on the %d heaviest" % (n, sum(r['samples'] for r in S.values()), sum(r['samples'] for r in RS.values()), len(RS)))
tot = sum(est.values()); log("   TOTAL estimated CPU: %.0f h   (first pass alone: %.0f h)" % (tot / H, sum(r['est'] for r in S.values()) / H))
if F:
    sf = sum(F[t]['full'] for t in F); se = sum(S[t]['est'] for t in F)
    lr = [math.log(F[t]['full'] / S[t]['est']) for t in F if S[t]['est'] > 5]
    log("   validation: %d graphs run completely: measured %.0f s, sampled estimate %.0f s, ratio %.2f; per graph (estimate > 5 s) the ratio ranges over %.2f .. %.2f"
        % (len(F), sf, se, sf / se, math.exp(min(lr)), math.exp(max(lr))))
log("   per number of edges:   E  graphs   sum [h]   mean [s]   max [h]")
for E in range(13):
    l = [est[t] for t in est if G5[t]['E'] == E]; log("                         %2d   %5d   %7.1f   %8.0f   %7.2f" % (E, len(l), sum(l) / H, sum(l) / len(l), max(l) / H))
for thr in (600, 1800, 3600, 7200):
    l = [x for x in est.values() if x > thr]; log("   graphs above %4.0f min: %4d  (they hold %.0f h = %.0f %% of the total)" % (thr / 60.0, len(l), sum(l) / H, 100 * sum(l) / tot))
log("   slowest graphs:")
for t in sorted(est, key=lambda t: -est[t])[:12]:
    G = G5[t]; r = RS.get(t, S[t])
    log("      t=%4d E=%2d V=%d gv=%s valences=%s edges=%s: %.2f h (%d decorations x %d subsets x %.0f ms, %d samples)"
        % (t, G['E'], G['V'], list(G['gv']), [len(h) for h in G['hsAt']], [list(e) for e in G['ed']], est[t] / H, r['ncombo'], r['nsub'], 1e3 * r['mean'], r['samples']))


# ---------------------------------------------------------------- proxies and task balancing
def proxy_now(G):
    E = G['E']; return (2 ** min(E, 6)) * (5 ** min(G['V'], 4)) * max(1, 12 - E) ** 2


def proxy_struct(t):
    G = G5[t]; steps = fq4_split.Split.plan(None, G); U = 0; fr = 0
    for v, conn, keep, loops in steps: U += 5 ** (fr + 1) * (1 + len(conn) + len(loops)); fr = len(keep)
    return S[t]['ncombo'] * S[t]['nsub'] * U


PS = {t: proxy_struct(t) for t in est}; PN = {t: proxy_now(G5[t]) for t in est}


def corr(P):
    xs = [math.log(P[t]) for t in est if est[t] > 1 and P[t] > 0]; ys = [math.log(est[t]) for t in est if est[t] > 1 and P[t] > 0]
    mx = sum(xs) / len(xs); my = sum(ys) / len(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys)); sxx = sum((x - mx) ** 2 for x in xs); syy = sum((y - my) ** 2 for y in ys)
    res = [y - my - (sxy / sxx) * (x - mx) for x, y in zip(xs, ys)]
    return sxy / math.sqrt(sxx * syy), sxy / sxx, math.exp(math.sqrt(sum(r * r for r in res) / len(res)))


for name, P in (("proxy of genus4_split.cost_order", PN), ("structural proxy ncombo*nsub*Sum 5^(f+1)(1+d)", PS)):
    c, slope, sc = corr(P); log("   %s: log-log correlation with the estimate %.3f, slope %.2f, residual scatter factor %.2f" % (name, c, slope, sc))


def round_robin(P, N, cost):
    order = sorted(cost, key=lambda t: -P[t]); load_ = [0.0] * N
    for r, t in enumerate(order): load_[r % N] += cost[t]
    return load_


def lpt(N, cost):
    load_ = [0.0] * N
    for t in sorted(cost, key=lambda t: -cost[t]):
        i = min(range(N), key=lambda j: load_[j]); load_[i] += cost[t]
    return load_


def balance(cost, title):
    log("   task balancing, %s (mean / heaviest task, hours on the reference machine):" % title)
    for N in (25, 50, 100, 200):
        a = round_robin(PN, N, cost); b = round_robin(PS, N, cost); c = lpt(N, cost)
        log("      N = %3d: mean %.2f | round robin, present proxy: %.2f | round robin, structural proxy: %.2f | longest-first with the estimates: %.2f"
            % (N, sum(cost.values()) / N / H, max(a) / H, max(b) / H, max(c) / H))


balance(est, "Split.graph as it stands")
table = {str(t): dict(gv=list(G5[t]['gv']), ed=[list(e) for e in G5[t]['ed']], aut=G5[t]['aut'], est_split=round(est[t], 2)) for t in est}
if len(D) == n:
    dest = {t: (DR[t]['est'] if t in DR else D[t]['est']) for t in D}; dtot = sum(dest.values())
    log(); log("SUBSET-DP PROTOTYPE (g5_dp.graph_dp), same machine: TOTAL estimated CPU %.1f h  (speed-up %.1f); peak memory of a worker %.0f MB" % (dtot / H, tot / dtot, max(r['peak_mb'] for r in list(D.values()) + list(DR.values()))))
    if DF:
        sf = sum(DF[t]['full'] for t in DF); se = sum(D[t]['est'] for t in DF)
        log("   validation: %d graphs run completely with graph_dp: measured %.0f s, sampled estimate %.0f s, ratio %.2f; peak memory %.0f MB" % (len(DF), sf, se, sf / se, max(r['peak_mb'] for r in DF.values())))
    log("   per number of edges:   E  graphs   sum [h]   mean [s]   max [h]")
    for E in range(13):
        l = [dest[t] for t in dest if G5[t]['E'] == E]; log("                         %2d   %5d   %7.2f   %8.1f   %7.2f" % (E, len(l), sum(l) / H, sum(l) / len(l), max(l) / H))
    for thr in (600, 1800, 3600):
        l = [x for x in dest.values() if x > thr]; log("   graphs above %4.0f min: %4d  (they hold %.1f h)" % (thr / 60.0, len(l), sum(l) / H))
    log("   slowest graphs:")
    for t in sorted(dest, key=lambda t: -dest[t])[:10]:
        G = G5[t]; r = DR.get(t, D[t])
        log("      t=%4d E=%2d V=%d gv=%s valences=%s: %.2f h (%d decorations x %.2f s, %d samples)" % (t, G['E'], G['V'], list(G['gv']), [len(h) for h in G['hsAt']], dest[t] / H, r['ncombo'], r['mean'], r['samples']))
    balance(dest, "subset-DP prototype")
    for t in dest: table[str(t)]['est_dp'] = round(dest[t], 2)
else:
    log("subset-DP samples incomplete: %d of %d" % (len(D), n))
json.dump(table, open(os.path.join(HERE, "g5_cost_table.json"), "w"))
log("wrote g5_cost_table.json (%d graphs)" % len(table))
