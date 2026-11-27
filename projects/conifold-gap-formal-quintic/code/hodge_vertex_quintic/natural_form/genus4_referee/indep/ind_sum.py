# -*- coding: utf-8 -*-
"""ind_sum.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): the split sum of one inner graph.

For a leg-free stable graph Gt (trivial-theory vertices, one colour each) and a labelled decoration (exponent on every half-edge;
the dilaton leaves are already summed inside the vertex tables) the point-frame term is
      (1/|Aut Gt|) Sum_colours prod_v V_v[c_v] prod_e E_e(a_e, b_e)[c_u][c_w],        E_e = E_M + E_out,
and the terms with the set S of edges taken from E_out belong to the outer graph Gt/(edges not in S).
Two implementations of the colour sum:
 (A) MOMENTUM SPACE (production; method different from the project's vertex elimination in colour space).  With the discrete
     Fourier transforms of ind_tables (vertex: one component kappa_v; edge: components (k, sigma_e - k)) the colour sum is
           5^V Sum_{k_e in Z/5} prod_e Eh_e[k_e] prod_v Vh_v    over the solutions of
           kappa_v + Sum_{e: u_e = v} k_e + Sum_{e: w_e = v} (sigma_e - k_e) = 0 mod 5 at every vertex v
     (5^(E-V+1) solutions: the momenta of the edges outside a spanning tree are free).  For the inner edges Eh is a number times
     a power of tau, so for each S the sum over the momenta of the inner edges is a sum of numbers, and only the momenta of the
     outer edges carry u-series (each divided by its u-valuation, which shortens the truncation).
 (B) COLOUR SPACE BRUTE FORCE (graph_brute; check of (A) on small graphs): plain sum over the 5^V colourings.
API: Summer(T, graphs).graph(t, max_outer=None) -> {outer index: [K coefficients]};  .graph_brute(t, max_outer=None).
max_outer = n keeps only the subsets S with at most n outer edges (the outer graphs with at most n edges are then complete).
"""
import itertools


class Summer(object):
    def __init__(self, T, graphs):
        self.T = T; self.graphs = graphs; self.index = {G['key']: i for i, G in enumerate(graphs)}
        self._struct = {}; self.stats = dict(decorations=0, charge_zero=0, subset_terms=0)

    # ------------------------------------------------------------ graph structure: spanning tree, momentum solution
    def structure(self, t):
        st = self._struct.get(t)
        if st is not None: return st
        from ind_graphs import contract
        G = self.graphs[t]; V = G['V']; E = G['E']; ed = G['ed']
        parent = {0: None}; order = [0]; i = 0
        while i < len(order):
            v = order[i]; i += 1
            for e, (u, w) in enumerate(ed):
                if u == w: continue
                o = w if u == v else (u if w == v else None)
                if o is not None and o not in parent: parent[o] = e; order.append(o)
        assert len(order) == V
        tree = set(parent[v] for v in order[1:]); free = [e for e in range(E) if e not in tree]
        elim = []                                                         # (v, parent edge, v is the u-end, [(edge, sign)] others)
        for v in reversed(order[1:]):
            te = parent[v]; others = []
            for e, (u, w) in enumerate(ed):
                if e == te: continue
                if u == v: others.append((e, 1, 0))                       # + k_e
                if w == v: others.append((e, -1, 1))                      # + sigma_e - k_e
            elim.append((v, te, ed[te][0] == v, others))
        root = [(e, s, c) for e, (u, w) in enumerate(ed) for (x, s, c) in ((u, 1, 0), (w, -1, 1)) if x == 0]
        # homogeneous solutions: momenta of all edges for every assignment of the free momenta (kappa = sigma = 0)
        cfgs = []
        for kf in itertools.product(range(5), repeat=len(free)):
            k = [0] * E
            for e, x in zip(free, kf): k[e] = x
            for v, te, uend, others in elim:
                c = sum(s * k[e] for e, s, cs in others) % 5
                k[te] = (-c) % 5 if uend else c
            cfgs.append(tuple(k))
        nmax = self.T.K - 1 if self.T.variant != "id" else E
        subsets = []
        for S in range(1 << E):
            ins = [e for e in range(E) if (S >> e) & 1]
            if self.T.variant == "id" and len(ins) != E: continue
            if len(ins) > nmax: continue
            subsets.append((S, self.index[contract(G, S)], ins, [e for e in range(E) if not (S >> e) & 1]))
        st = dict(order=order, elim=elim, root=root, cfgs=cfgs, subsets=subsets)
        self._struct[t] = st
        return st

    def options(self, G):
        """per vertex: list of (ordered exponents of its half-edges, kappa, Vh, [V[c]])"""
        out = []
        for v in range(G['V']):
            tab = self.T.vertex_table(G['gv'][v], len(G['hs'][v])); lst = []
            for hp, (kap, Vh, Vc) in tab.items():
                for perm in set(itertools.permutations(hp)): lst.append((perm, kap, Vh, Vc))
            out.append(lst)
        return out

    # ------------------------------------------------------------ (A) momentum space
    def graph(self, t, max_outer=None):
        T = self.T; p = T.p; K = T.K; G = self.graphs[t]; V = G['V']; E = G['E']; ed = G['ed']; hs = G['hs']
        st = self.structure(t); cfgs = st['cfgs']; mul = T.mul
        subsets = [s for s in st['subsets'] if max_outer is None or len(s[2]) <= max_outer or T.variant == "id"]
        acc = {}
        for combo in itertools.product(*self.options(G)):
            self.stats['decorations'] += 1
            hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hs[v]): hp[h] = combo[v][0][r]
            keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]; sig = [T.sig[k] for k in keys]; kap = [combo[v][1] for v in range(V)]
            # particular solution (free momenta 0), then the root constraint
            k0 = [0] * E
            for v, te, uend, others in st['elim']:
                c = (kap[v] + sum((k0[e] if s == 1 else sig[e] - k0[e]) for e, s, cs in others)) % 5
                k0[te] = (-c) % 5 if uend else (sig[te] + c) % 5
            if (kap[0] + sum((k0[e] if s == 1 else sig[e] - k0[e]) for e, s, cs in st['root'])) % 5:
                self.stats['charge_zero'] += 1; continue
            Vprod = combo[0][2]
            for v in range(1, V): Vprod = mul(Vprod, combo[v][2])
            if Vprod is None: continue
            # tables indexed by the homogeneous momentum b: actual momentum (b + k0) mod 5
            inner = [[T.FI[keys[e]][(b + k0[e]) % 5] for b in range(5)] for e in range(E)]
            outer = [[T.FO[keys[e]][(b + k0[e]) % 5] for b in range(5)] for e in range(E)]
            vals = [T.valO[k] for k in keys]; npow = [a + b + 1 for a, b in keys]; vt_cache = {}
            for S, pos, ins, outs in subsets:
                v0 = sum(vals[e] for e in ins)
                if v0 >= K: continue
                L = K - v0; part = {}
                for cfg in cfgs:
                    s = 1
                    for e in outs:
                        s = s * inner[e][cfg[e]] % p
                        if not s: break
                    if s:
                        kk = tuple([cfg[e] for e in ins]); part[kk] = part.get(kk, 0) + s
                C = [0] * L; hit = False
                for kk, s in part.items():
                    s %= p
                    if not s: continue
                    ser = [s]
                    for e, b in zip(ins, kk):
                        ser = mul(ser + [0] * (L - len(ser)), outer[e][b], L) if ser is not None else None
                        if ser is None: break
                    if ser is None: continue
                    hit = True
                    for i, x in enumerate(ser): C[i] += x
                if not hit: continue
                self.stats['subset_terms'] += 1
                N = sum(npow[e] for e in outs); pre = vt_cache.get(N)
                if pre is None: pre = mul(Vprod, T.tp[N]); vt_cache[N] = pre
                res = mul(pre, [x % p for x in C] + [0] * (K - L), L)
                if res is None: continue
                a = acc.setdefault(pos, [0] * K)
                for i in range(L): a[v0 + i] += res[i]
        f = pow(5, V, p) * pow(G['aut'], p - 2, p) % p
        return {pos: [x * f % p for x in a] for pos, a in acc.items()}

    # ------------------------------------------------------------ (B) colour space brute force
    def graph_brute(self, t, max_outer=None):
        T = self.T; p = T.p; K = T.K; G = self.graphs[t]; V = G['V']; E = G['E']; ed = G['ed']; hs = G['hs']
        st = self.structure(t); mul = T.mul; acc = {}
        subsets = [s for s in st['subsets'] if max_outer is None or len(s[2]) <= max_outer or T.variant == "id"]
        for combo in itertools.product(*self.options(G)):
            hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hs[v]): hp[h] = combo[v][0][r]
            keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
            for S, pos, ins, outs in subsets:
                mats = [(T.EO[keys[e]] if (S >> e) & 1 else T.EM[keys[e]]) for e in range(E)]
                for col in itertools.product(range(5), repeat=V):
                    x = [1] + [0] * (K - 1)
                    for v in range(V):
                        x = mul(x, combo[v][3][col[v]])
                        if x is None: break
                    if x is None: continue
                    for e, (u, w) in enumerate(ed):
                        x = mul(x, mats[e][col[u]][col[w]])
                        if x is None: break
                    if x is None: continue
                    a = acc.setdefault(pos, [0] * K)
                    for i in range(K): a[i] += x[i]
        f = pow(G['aut'], p - 2, p)
        return {pos: [x * f % p for x in a] for pos, a in acc.items()}
