# -*- coding: utf-8 -*-
"""ms_engine.py -- correlators of the frozen (model) vertex theory V = M_0.omega^ in the FLAT basis, exact over Q.

Givental-Teleman graph sum in "momentum space" (ms_rows.py docstring): flat indices are traded for Z/5 charges, which are
conserved at the vertices, so a graph needs 5^(b_1) charge assignments instead of 5^V colourings and no root of unity.
  AMP(g, n)[(Q_1..Q_n), (e_1..e_n)]  = amputated amplitude: legs enter the vertices with charge Q_l and psi-exponent e_l
       = Sum_graphs 1/|Aut| Sum_charges Sum_exponents  prod_v U(g_v; exponents at v)  prod_edges ECOEF(k, l, q),
  U(g; e_1..e_m) = 5 Sum_{leaf multisets} WK(g; e, leaves)/prod(mult!) prod T(leaf),   Sum (leaf - 1) = 3g - 3 + m - Sum e;
  <tau_d1(h^a1) .. tau_dn(h^an)>_g = Sum_{k_l >= 0} prod_l (-1)^k_l g^(a_l)_(k_l)  AMP[(a_l + d_l - 1), (d_l + k_l)],
  <tau_d1(H^a1) ..>_g = prod N_(a_l) times that.
Graphs: fq3_graphs.legged_graphs (read only; legs labelled).  Witten-Kontsevich numbers: fq_hodge.wk (read only).
The theory is the LG model at v = 1 (H = d/dv); v-dependence: corr(v) = corr(1) v^vexp, see vexp().
Cache: one JSON per (g, n) in the scratch folder (never in the project).
"""
import itertools, json, math, os, sys, time
from fractions import Fraction as Fr
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr
sys.path.insert(1, mr.PROJ)
from fq_hodge import wk
from fq3_graphs import legged_graphs


def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for p in range(min(n, maxpart), 0, -1):
        for r in partitions(n - p, p): yield (p,) + r


@lru_cache(maxsize=None)
def U(g, exps):
    """vertex factor for sorted exponents exps (half-edges and legs); 0 if the dimension is exceeded"""
    m = len(exps); delta = 3 * g - 3 + m - sum(exps)
    if delta < 0 or 2 * g - 2 + m <= 0: return Fr(0)
    tot = Fr(0)
    for part in partitions(delta):
        leaves = tuple(p + 1 for p in part)
        w = wk(g, tuple(exps) + leaves)
        if w == 0: continue
        mult = {}
        for e in leaves: mult[e] = mult.get(e, 0) + 1
        c = w
        for e, k in mult.items(): c = c * mr.tleaf(e) ** k / math.factorial(k)
        tot += c
    return 5 * tot


def tuples_le(n, budget):
    if n == 0:
        yield (); return
    for f in range(budget + 1):
        for r in tuples_le(n - 1, budget - f): yield (f,) + r


def graph_amp(Gm, n, acc):
    """add the contribution of one legged graph to acc[(Q), (e)]"""
    V = Gm['V']; gv = Gm['gv']; ed = Gm['ed']; hsAt = Gm['hsAt']; legsAt = Gm['legsAt']; E = len(ed)
    inv_aut = Fr(1, Gm['aut'])
    # per-vertex options: exponents for the slots (half-edges of hsAt[v], then legs of legsAt[v])
    opts = []
    for v in range(V):
        ns = len(hsAt[v]) + len(legsAt[v]); lst = []
        for ex in tuples_le(ns, 3 * gv[v] - 3 + ns):
            u = U(gv[v], tuple(sorted(ex)))
            if u: lst.append((ex, u))
        if not lst: return
        opts.append(lst)
    # spanning tree (BFS from 0); loops and the other non-tree edges carry free charges
    parent = {0: None}; order = [0]; tree = set()
    for v in order:
        for e, (a, b) in enumerate(ed):
            if a == b or e in tree: continue
            w = b if a == v else (a if b == v else None)
            if w is None or w in parent: continue
            parent[w] = e; tree.add(e); order.append(w)
    assert len(order) == V
    loops = [e for e, (a, b) in enumerate(ed) if a == b]
    free = [e for e, (a, b) in enumerate(ed) if a != b and e not in tree]
    # charge assignments: for every (leg charges with sum 0, free charges) the tree charges
    assigns = {}
    for Q in itertools.product(range(5), repeat=n):
        if sum(Q) % 5: continue
        lst = []
        for fq in itertools.product(range(5), repeat=len(free)):
            qe = [0] * E
            for e, x in zip(free, fq): qe[e] = x
            net = [0] * V                                  # sum of the charges already known at v
            for v in range(V):
                for l in legsAt[v]: net[v] += Q[l]
            for e in free:
                a, b = ed[e]; net[a] += qe[e]; net[b] -= qe[e]
            for v in reversed(order[1:]):
                e = parent[v]; a, b = ed[e]
                s = (-net[v]) % 5                          # charge on v's end of its parent edge
                if a == v: qe[e] = s; net[b] -= s
                else: qe[e] = (-s) % 5; net[a] += qe[e]
            assert net[0] % 5 == 0
            lst.append(tuple(qe))
        assigns[Q] = lst
    nonloop = [e for e in range(E) if ed[e][0] != ed[e][1]]
    slot_of_h = {}
    for v in range(V):
        for r, h in enumerate(hsAt[v]): slot_of_h[h] = (v, r)
    leg_slot = {}
    for v in range(V):
        for r, l in enumerate(legsAt[v]): leg_slot[l] = (v, len(hsAt[v]) + r)
    for combo in itertools.product(*opts):
        base = inv_aut
        for v in range(V): base = base * combo[v][1]
        kl = []
        for e in range(E):
            v0, r0 = slot_of_h[2 * e]; v1, r1 = slot_of_h[2 * e + 1]
            kl.append((combo[v0][0][r0], combo[v1][0][r1]))
        for e in loops: base = base * mr.eloop(kl[e][0], kl[e][1])
        if not base: continue
        ev = tuple(combo[leg_slot[l][0]][0][leg_slot[l][1]] for l in range(n))
        for Q, lst in assigns.items():
            s = Fr(0)
            for qe in lst:
                x = Fr(1)
                for e in nonloop:
                    x = x * mr.ecoef(kl[e][0], kl[e][1], qe[e])
                    if not x: break
                s += x
            if s:
                key = (Q, ev); acc[key] = acc.get(key, 0) + base * s


def _cache_path(g, n):
    return os.path.join(mr.SCRATCH, "amp_%d_%d.json" % (g, n))


_AMP = {}
def amp(g, n, verbose=False):
    key = (g, n)
    if key in _AMP: return _AMP[key]
    path = _cache_path(g, n)
    if os.path.exists(path):
        with open(path) as f: raw = json.load(f)
        A = {(tuple(int(c) for c in k.split("|")[0].split(",") if c != ""), tuple(int(c) for c in k.split("|")[1].split(",") if c != "")): Fr(v) for k, v in raw.items()}
        _AMP[key] = A; return A
    t0 = time.time(); acc = {}
    graphs = legged_graphs(g, n, 3 * g - 3 + n)
    for i, Gm in enumerate(graphs):
        graph_amp(Gm, n, acc)
        if verbose and (i % 200 == 0): print("   amp(%d,%d): graph %d/%d, %.0f s" % (g, n, i, len(graphs), time.time() - t0), flush=True)
    A = {k: v for k, v in acc.items() if v}
    os.makedirs(mr.SCRATCH, exist_ok=True)
    with open(path, "w") as f:
        json.dump({",".join(map(str, k[0])) + "|" + ",".join(map(str, k[1])): str(v) for k, v in A.items()}, f)
    if verbose: print("   amp(%d,%d): %d graphs, %d entries, %.1f s" % (g, n, len(graphs), len(A), time.time() - t0), flush=True)
    _AMP[key] = A; return A


def corr_h(g, ins):
    """<prod tau_d(h^a)>_g for ins = ((a_1, d_1), ..) in the basis h^a = N_a H^a (model at v = 1)"""
    ins = tuple(ins); n = len(ins)
    if 2 * g - 2 + n <= 0: raise ValueError("unstable")
    D = 3 * g - 3 + n - sum(d for a, d in ins)
    if D < 0 or sum(a + d - 1 for a, d in ins) % 5: return Fr(0)
    A = amp(g, n); Q = tuple((a + d - 1) % 5 for a, d in ins); tot = Fr(0)
    for ks in tuples_le(n, D):
        x = A.get((Q, tuple(d + k for (a, d), k in zip(ins, ks))))
        if not x: continue
        for (a, d), k in zip(ins, ks): x = x * mr.legw(a, k)
        tot += x
    return tot


def corr(g, ins):
    """<prod tau_d(H^a)>_g, H = d/dv, model at v = 1"""
    s = 1
    for a, d in ins: s *= mr.N_SIGN[a % 5]
    return s * corr_h(g, tuple((a % 5, d) for a, d in ins))


def vexp(g, ins):
    """exponent of v: 2 - 2g - (1/5) Sum (4 - c_a - 4 d); an integer iff the charge selection rule holds"""
    return Fr(2 - 2 * g) - Fr(sum(4 - mr.C_EXP[a % 5] - 4 * d for a, d in ins), 5)


if __name__ == "__main__":
    for (g, n) in [(0, 3), (0, 4), (1, 1), (1, 2), (2, 0), (2, 1), (1, 3), (2, 2)]:
        amp(g, n, verbose=True)
    print("F_2 =", corr(2, ()), " <H>_2 =", corr(2, ((1, 0),)), " <tau_1(1)>_1 =", corr(1, ((0, 1),)), " <H>_1 =", corr(1, ((1, 0),)))
    print("<H,H,H>_0 =", corr(0, ((1, 0),) * 3), " <H,H>_2 =", corr(2, ((1, 0),) * 2))
