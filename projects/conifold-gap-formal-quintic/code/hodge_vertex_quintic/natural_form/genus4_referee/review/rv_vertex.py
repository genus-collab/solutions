# -*- coding: utf-8 -*-
"""rv_vertex.py -- unit tests of the helpers of fq4_split.Split.vertex_table at the genus-4 parameters (2^89-1, K = 6, KZ = 9).

(1) partitions_min2(m, total_max, part_max) = all non-increasing m-tuples with parts in [2, part_max] and sum <= total_max (brute force);
(2) sorted_tuples(n, s) = all non-increasing n-tuples of non-negative integers with sum s (brute force);
(3) vertex_table(g, hs)[hp][j] against a brute-force sum over ORDERED leaf tuples with the factor 1/m!:
       Sum_m (1/m!) Sum_(e_1..e_m >= 2) WK(g; hp, e) (zeta^(4j) xhat)^(2g-2+hs+m) prod_i T_R(e_i)[j]     (own series arithmetic),
    for every ordered hp with the right total degree, for all vertex types (g, hs) that occur in the genus-4 graphs; a key missing
    from the table must have a zero brute-force value;
(4) the largest leaf exponent that occurs is 3g-2 = 10 = KZ + 1 (the cap of partitions_min2 is not active below it) and the largest
    power of xhat used is within the table xpow.
Usage: py -I rv_vertex.py   -> out/vertex.log
"""
import os, sys, itertools, math, time
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq_fp import P2
from fq_hodge import wk
import fq4_split as F
LOG = open(os.path.join(OUT, "vertex.log"), "w", encoding="utf-8"); FAILED = []
def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)
t0 = time.time()
ok = True
for m in range(0, 5):
    for tm in range(0, 13):
        for pm in range(2, 8):
            a = sorted(F.partitions_min2(m, tm, pm))
            b = sorted(t for t in itertools.product(range(2, pm + 1), repeat=m) if sum(t) <= tm and all(t[i] >= t[i + 1] for i in range(m - 1)))
            ok = ok and a == b and len(set(a)) == len(a)
check("(1) partitions_min2 = brute force for m <= 4, total <= 12, part_max <= 7 (no duplicates)", ok)
ok = True
for n in range(0, 6):
    for s in range(0, 9):
        a = sorted(F.sorted_tuples(n, s))
        b = sorted(t for t in itertools.product(range(s + 1), repeat=n) if sum(t) == s and all(t[i] >= t[i + 1] for i in range(n - 1)))
        ok = ok and a == b and len(set(a)) == len(a)
check("(2) sorted_tuples = brute force for n <= 5, s <= 8 (no duplicates)", ok)

sp = F.Split(P2, 4, (Fr(0), Fr(0)), 6, 'nat'); p = sp.p; K = sp.K; R = sp.R
def smul(a, b):
    out = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i):
                if b[j]: out[i + j] += a[i] * b[j]
    return [c % p for c in out]
def N(a): return [0] * K if a is None else list(a)
def comps(m, left):
    """ordered m-tuples of parts >= 2 with sum <= left"""
    if m == 0: yield (); return
    for x in range(2, left - 2 * (m - 1) + 1):
        for r in comps(m - 1, left - x): yield (x,) + r
types = sorted(set((G['gv'][v], len(G['hsAt'][v])) for G in sp.inner_graphs for v in range(G['V'])))
maxleaf = 0; maxn = 0; nvals = 0; bad = []
for (g, hs) in types:
    tab = sp.vertex_table(g, hs); dim0 = 3 * g - 3 + hs; brute = {}
    for m in range(0, dim0 + 1):
        dim = dim0 + m; n = 2 * g - 2 + hs + m; fm = R.inv(math.factorial(m))
        for leaves in comps(m, dim):
            rest = dim - sum(leaves)
            for hp in itertools.product(range(rest + 1), repeat=hs):
                if sum(hp) != rest: continue
                w = wk(g, tuple(hp) + tuple(leaves))
                if w == 0: continue
                maxleaf = max(maxleaf, max(leaves) if leaves else 0); maxn = max(maxn, n)
                c = R.of(w) * fm % p
                vec = brute.setdefault(hp, [[0] * K for _ in range(5)])
                for j in range(5):
                    b = [x * pow(sp.zph[j], n, p) % p * c % p for x in N(sp.xpow[n])]
                    for e in leaves: b = smul(b, N(sp.TR[e][j]))
                    vec[j] = [(x + y) % p for x, y in zip(vec[j], b)]
    for hp, vec in brute.items():
        key = tuple(sorted(hp, reverse=True)); t = tab.get(key); nvals += 1
        if [N(x) for x in (t or [None] * 5)] != vec: bad.append((g, hs, hp))
    extra = [k for k in tab if k not in set(tuple(sorted(h, reverse=True)) for h in brute)]
    if extra: bad.append((g, hs, "extra keys", extra[:3]))
check("(3) vertex_table = brute force over ordered leaves with 1/m! for the %d vertex types %s (%d ordered exponent tuples x 5 colours x %d u-orders)" % (len(types), types, nvals, K), not bad, bad[:5])
check("(4) largest leaf exponent with a non-zero KW number = %d = KZ + 1; largest xhat power %d <= %d" % (maxleaf, maxn, len(sp.xpow) - 1), maxleaf == sp.KZ + 1 and maxn <= len(sp.xpow) - 1)
log("RV_VERTEX %s  (%.0f s)" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
