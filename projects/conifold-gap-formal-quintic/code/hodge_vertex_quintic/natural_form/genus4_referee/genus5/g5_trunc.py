# -*- coding: utf-8 -*-
"""g5_trunc.py -- prototype of ONE proposed change for genus 5: colour sums at the working precision K - |S|.

Every entry of E_out(a, b) is u * (series) (u-valuation exactly 1 at genus 5, see g5_setup.log), so a term with outer-edge set S is
u^|S| * (product of series), and only K - |S| coefficients of every factor are needed.  colour_sum_trunc does the same vertex
elimination as Split.colour_sum with E_out/u and all series cut to n = K - |S| coefficients; the result is shifted back by u^|S|.
This script (one process) draws graphs with probability proportional to their estimated cost (g5_sample_*.jsonl), draws random
(decoration, S) pairs, checks  colour_sum_trunc == Split.colour_sum  coefficient by coefficient (mod p) and compares the times.
Usage: py -I g5_trunc.py [number of graphs, default 60] [pairs per graph, default 6] [time limit s, default 900]
"""
import os, sys, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from g5_lib import make_split, SplitTools
from fq_fp import P2

LOG = open(os.path.join(HERE, "g5_trunc.log"), "a", encoding="utf-8")


def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def mulN(a, b, n, p):
    if a is None or b is None: return None
    out = [0] * n
    for i in range(n):
        x = a[i]
        if x:
            for j in range(n - i):
                y = b[j]
                if y: out[i + j] += x * y
    out = [c % p for c in out]
    return out if any(out) else None


def colour_sum_trunc(sp, steps, vecs, mats, n):
    """as Split.colour_sum, at precision n (all inputs are lists with at least n entries, or None)"""
    p = sp.p; cur = {(): [1] + [0] * (n - 1)}
    for v, conn, keep, loops in steps:
        new = {}; vec = vecs[v]
        for key, val in cur.items():
            for c in range(5):
                x = mulN(val, vec[c], n, p)
                for e in loops:
                    if x is None: break
                    x = mulN(x, mats[e][c][c], n, p)
                for e, pos, flip in conn:
                    if x is None: break
                    x = mulN(x, mats[e][c][key[pos]] if flip else mats[e][key[pos]][c], n, p)
                if x is None: continue
                full = key + (c,); nk = tuple(full[k] for k in keep); y = new.get(nk)
                new[nk] = x if y is None else [(a + b) % p for a, b in zip(x, y)]
        cur = new
        if not cur: return None
    r = cur.get(())
    return r if (r is None or any(r)) else None


ng = int(sys.argv[1]) if len(sys.argv) > 1 else 60; npairs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
limit = float(sys.argv[3]) if len(sys.argv) > 3 else 900.0
est = {}
for line in open(os.path.join(HERE, "g5_sample_g5_89_nat_K8.jsonl")):
    if line.strip(): r = json.loads(line); est[r['t']] = r['est']
sp = make_split(P2, 5, (0, 0), 8, 'nat'); tools = SplitTools(sp); K = sp.K; p = sp.p
assert all(v == 1 for v in sp.valO.values())
EOs = {key: [[(x[1:] if x is not None else None) for x in row] for row in M] for key, M in sp.EO.items()}     # E_out / u
rng = random.Random(11); ts = sorted(est); tot = sum(est.values())
chosen = rng.choices(ts, weights=[est[t] for t in ts], k=ng)
log("=" * 100); log("g5_trunc: %d graphs drawn with probability ~ estimated cost, %d (decoration, S) pairs each" % (ng, npairs))
t0 = time.time(); bad = 0; n_pairs = 0; ratios = []; T_std = 0.0; T_tr = 0.0
for t in chosen:
    if time.time() - t0 > limit: log("   time limit reached"); break
    P = tools.prepare(t); G = P['G']; V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; a_std = 0.0; a_tr = 0.0
    for _ in range(npairs):
        combo = [rng.choice(lst) for lst in P['opts']]; S, pos = rng.choice(P['subsets']); s = bin(S).count("1"); n = K - s; hp = [0] * (2 * E)
        for v in range(V):
            for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
        keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]; vecs = [combo[v][1] for v in range(V)]
        ta = time.time(); z = sp.colour_sum(P['steps'], vecs, [sp.EO[keys[e]] if (S >> e) & 1 else sp.EM[keys[e]] for e in range(E)]); tb = time.time()
        w = colour_sum_trunc(sp, P['steps'], vecs, [EOs[keys[e]] if (S >> e) & 1 else sp.EM[keys[e]] for e in range(E)], n); tc = time.time()
        zz = [x % p for x in z] if z is not None else [0] * K
        ww = [0] * s + ([x % p for x in w] if w is not None else [0] * n)
        if zz != ww: bad += 1
        a_std += tb - ta; a_tr += tc - tb; n_pairs += 1
    ratios.append(a_tr / a_std if a_std > 0 else 1.0); T_std += a_std; T_tr += a_tr
log("%s truncated colour sums equal the standard ones on %d random pairs (%d mismatches)" % ("PASS" if bad == 0 else "FAIL", n_pairs, bad))
if ratios:
    log("   cost-weighted mean of (truncated time / standard time) over %d graphs: %.3f  (speed-up factor %.2f); pooled time ratio %.3f"
        % (len(ratios), sum(ratios) / len(ratios), len(ratios) / sum(ratios), T_tr / T_std))
log("   wall %.0f s" % (time.time() - t0))
