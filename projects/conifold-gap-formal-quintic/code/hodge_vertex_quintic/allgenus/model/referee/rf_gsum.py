# -*- coding: utf-8 -*-
"""rf_gsum.py -- referee's graph sum for the frozen model theory in the NORMALISED CANONICAL frame (colours, zeta in F_p).

  <prod_l tau_(k_l)(H^(a_l))>_g = Sum_Gamma 1/|Aut| Sum_colourings Sum_exponents  prod_v W(g_v; slot exponents)[c_v]
                                   prod_edges E[(k, l)][c_u][c_v]  prod_legs (R^-1_j Psi^-1 H^a)_(c_v)  (leg exponent k_l + j),
  W(g; d)[i] = Sum_(leaf multisets) WK(g; d, leaves) prod (T_e[i])^m_e / m_e!  sqrt(Delta_i)^(2g - 2 + |d| + #leaves)
(Givental-Teleman sum with translation T(z) = z(1 - R^-1(z))1 = O(z^2), so p^* psi = psi on the vertex moduli spaces).
Stable graphs: fq3_graphs.legged_graphs (project module, read only).  Colour sums by variable elimination.  WK, R, E, T: rf_core.
amp(model, g, n) returns the amputated amplitude {(leg colours, leg exponents): value mod p}; cached as JSON in the scratch folder.
"""
import itertools, json, math, os, sys, time
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as rc
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(1, PROJ)
from fq3_graphs import legged_graphs
SCR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scratch")


def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for q in range(min(n, mx), 0, -1):
        for r in partitions(n - q, q): yield (q,) + r


class GSum(object):
    def __init__(self, model):
        self.m = model; self.p = model.p; self.vcache = {}

    def W(self, g, exps):
        key = (g, tuple(sorted(exps)))
        if key in self.vcache: return self.vcache[key]
        m = self.m; p = self.p; nv = len(exps); delta = 3 * g - 3 + nv - sum(exps)
        out = [0] * 5
        if delta >= 0 and 2 * g - 2 + nv > 0:
            for part in partitions(delta):
                leaves = tuple(q + 1 for q in part)
                w = rc.wk(g, tuple(sorted(tuple(exps) + leaves)))
                if w == 0: continue
                mult = {}
                for e in leaves: mult[e] = mult.get(e, 0) + 1
                den = 1
                for e, k in mult.items(): den *= math.factorial(k)
                c = m.fr(w / den)
                for i in range(5):
                    x = c * pow(m.D[i], 2 * g - 2 + nv + len(leaves), p) % p
                    for e, k in mult.items(): x = x * pow(m.T[e][i], k, p) % p
                    out[i] = (out[i] + x) % p
        self.vcache[key] = out
        return out

    def colour_sum(self, V, unary, pair, keep):
        """Sum over colours of all vertices not in keep of prod unary[v][c_v] prod pair[(a,b)][c_a][c_b]; table over keep"""
        p = self.p
        factors = [((v,), {(c,): unary[v][c] for c in range(5)}) for v in range(V)]
        for (a, b), M in pair.items():
            factors.append(((a, b), {(i, j): M[i][j] for i in range(5) for j in range(5)}))
        order = sorted([v for v in range(V) if v not in keep], key=lambda v: sum(1 for f in factors if v in f[0]))
        for x in order:
            inv_ = [f for f in factors if x in f[0]]; rest = [f for f in factors if x not in f[0]]
            nvars = tuple(sorted(set(y for f in inv_ for y in f[0]) - {x}))
            tab = {}
            for asg in itertools.product(range(5), repeat=len(nvars)):
                d = dict(zip(nvars, asg)); s = 0
                for cx in range(5):
                    d[x] = cx; pr = 1
                    for vs, t in inv_:
                        pr = pr * t[tuple(d[y] for y in vs)] % p
                        if not pr: break
                    s += pr
                tab[asg] = s % p
            factors = rest + [(nvars, tab)]
        ks = tuple(sorted(keep)); out = {}
        for asg in itertools.product(range(5), repeat=len(ks)):
            d = dict(zip(ks, asg)); pr = 1
            for vs, t in factors:
                pr = pr * t[tuple(d[y] for y in vs)] % p
                if not pr: break
            if pr: out[asg] = pr
        return ks, out

    def graph_amp(self, Gm, n, acc):
        p = self.p; m = self.m
        V = Gm['V']; gv = Gm['gv']; ed = Gm['ed']; hsAt = Gm['hsAt']; legsAt = Gm['legsAt']
        ia = pow(Gm['aut'], p - 2, p)
        opts = []
        for v in range(V):
            ns = len(hsAt[v]) + len(legsAt[v]); D = 3 * gv[v] - 3 + ns; lst = []
            for ex in itertools.product(range(D + 1), repeat=ns):
                if sum(ex) > D: continue
                w = self.W(gv[v], ex)
                if any(w): lst.append((ex, w))
            if not lst: return
            opts.append(lst)
        slot = {}
        for v in range(V):
            for r, h in enumerate(hsAt[v]): slot[h] = (v, r)
        legslot = {}
        for v in range(V):
            for r, l in enumerate(legsAt[v]): legslot[l] = (v, len(hsAt[v]) + r)
        keep = set(legslot[l][0] for l in range(n))
        for combo in itertools.product(*opts):
            unary = [list(combo[v][1]) for v in range(V)]; pair = {}; ok = True
            for e, (a, b) in enumerate(ed):
                va, ra = slot[2 * e]; vb, rb = slot[2 * e + 1]
                k = combo[va][0][ra]; l = combo[vb][0][rb]
                Em = m.E.get((k, l))
                if Em is None: ok = False; break
                if va == vb:
                    for c in range(5): unary[va][c] = unary[va][c] * Em[c][c] % p
                else:
                    key = (va, vb) if va < vb else (vb, va)
                    Mx = Em if va < vb else [[Em[j][i] for j in range(5)] for i in range(5)]
                    if key in pair:
                        old = pair[key]; pair[key] = [[old[i][j] * Mx[i][j] % p for j in range(5)] for i in range(5)]
                    else: pair[key] = Mx
            if not ok: continue
            ks, tab = self.colour_sum(V, unary, pair, keep)
            lexp = tuple(combo[legslot[l][0]][0][legslot[l][1]] for l in range(n))
            pos = [ks.index(legslot[l][0]) for l in range(n)]
            for asg, val in tab.items():
                key = (tuple(asg[q] for q in pos), lexp)
                acc[key] = (acc.get(key, 0) + ia * val) % p


_AMP = {}
def amp(model, g, n, verbose=False):
    key = (model.p, model.sign, g, n)
    if key in _AMP: return _AMP[key]
    path = os.path.join(SCR, "rfamp_%d_%d_%d_%d.json" % (model.p.bit_length(), 1 if model.sign > 0 else 0, g, n))
    if os.path.exists(path):
        with open(path) as f: raw = json.load(f)
        A = {(tuple(x[0]), tuple(x[1])): x[2] for x in raw}; _AMP[key] = A; return A
    assert model.K >= 3 * g - 1 + n, "model.K too small"
    gs = GSum(model); acc = {}; t0 = time.time()
    graphs = legged_graphs(g, n, 3 * g - 3 + n)
    for i, Gm in enumerate(graphs):
        gs.graph_amp(Gm, n, acc)
        if verbose and i % 100 == 0: print("   (%d,%d) graph %d/%d %.0f s" % (g, n, i, len(graphs), time.time() - t0), flush=True)
    A = {k: v for k, v in acc.items() if v}
    os.makedirs(SCR, exist_ok=True)
    with open(path, "w") as f: json.dump([[list(k[0]), list(k[1]), v] for k, v in A.items()], f)
    if verbose: print("   amp(%d,%d): %d graphs, %d entries, %.1f s" % (g, n, len(graphs), len(A), time.time() - t0), flush=True)
    _AMP[key] = A; return A


def corr(model, g, ins):
    """<prod tau_k(H^a)>_g at v = 1 mod p, ins = ((a, k), ...)"""
    p = model.p; n = len(ins); A = amp(model, g, n); tot = 0
    for (cols, exps), val in A.items():
        x = val
        for (a, k), c, d in zip(ins, cols, exps):
            j = d - k
            if j < 0: x = 0; break
            x = x * sum(model.Rinv[j][c][b] * model.flat[a % 5][b] for b in range(5)) % p
            if not x: break
        tot = (tot + x) % p
    return tot
