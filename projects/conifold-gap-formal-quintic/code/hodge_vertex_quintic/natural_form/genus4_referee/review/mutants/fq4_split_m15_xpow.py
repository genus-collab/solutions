# -*- coding: utf-8 -*-
"""fq4_split.py -- the outer graph sum R'.(M.omega^) of the formal quintic at ANY genus by EDGE SPLITTING (no vertex tensors).

Composition law.  With R = R' M (hatted theory, colour = normalised canonical basis of the trivial rank-5 theory omega^) the edge
bivector of R splits as
        E_R(a, b) = E_M(a, b) + E_out(a, b),      E_out := M^-1 E_R' M^-T = E_R - E_M,
and the leaves of R are T_R = T_M + M^-1 T_R'.  Hence the point-frame sum over the leg-free stable graphs Gt of genus g, with every
edge marked INNER (factor E_M) or OUTER (factor E_out), is the outer sum: the terms with outer edge set S of Gt belong to the outer
graph Gt/(inner edges) (contract every inner edge), with weight 1/|Aut Gt| for every labelled S.  The term S = {} is the vertex graph.
All vertices are those of the trivial theory (Witten-Kontsevich numbers, one colour per vertex), so a term is a colour sum over the
vertices of Gt (done by vertex elimination) of products of u-series.
Ingredients (frame (alpha, beta), prime p, u-precision K = number of coefficients u^0..u^(K-1), zhat-order KZ = 3g - 3):
  E_R(a, b)[i][j], T_R(e)[i]   point-frame tables of the hatted rows (fq_sser.EdgeLeaf on fq3_msum.HatRowsU), series in u;
  E_M(a, b) = tau^(a+b+1) Eh_0(a, b)   for M = M_0(tau(u) zhat)  (Eh_0: fq3_mvertex formula for the symbol M_0 mod zhat^(KZ+1));
  variants: 'nat' tau = 1/2F1(1, 4/5; 9/5; u);  'M0' tau = 1 (constant M_0);  'id' M = 1 (E_M = 0: point frame, outer graph = Gt);
            'nat_t<k>' control: the natural tau with its u^k coefficient shifted by 1.
  vertex (g_v, half-edge exponents hp):  Sum_{m, leaves e_i >= 2} WK(g_v; hp, leaves)/m! (zeta^(4j) xhat)^(2g_v-2+n_v+m) prod T_R(e_i)[j].
A genus-g term is u^(-(2g-2)) G(u); an outer graph is regular at the conifold point iff [u^0..u^(2g-3)] G = 0.
Series are plain lists of K residues (None = zero).  API: Split(p, g, ab, K, variant); .inner_graphs (list of dicts, fq3_graphs format);
.outer_keys (canonical keys, position = outer index); .graph(t) -> {outer index: [K coefficients]} for the inner graph number t.
"""
import os, sys, itertools, math
POPMAX = 99
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing
from fq_hodge import wk
from fq3_graphs import legged_graphs, canon


def tau_pred(n):
    inv = [Fr(4, 5 * k + 4) for k in range(n)]; tau = [Fr(1)]
    for k in range(1, n): tau.append(-sum(inv[j] * tau[k - j] for j in range(1, k + 1)))
    return tau


def partitions_min2(m, total_max, part_max):
    """non-increasing m-tuples of parts in [2, part_max] with sum <= total_max"""
    def rec(k, cap, left):
        if k == 0: yield (); return
        for x in range(min(cap, left - 2 * (k - 1)), 1, -1):
            for r in rec(k - 1, x, left - x): yield (x,) + r
    return rec(m, part_max, total_max)


def sorted_tuples(n, s):
    """non-increasing n-tuples of non-negative integers with sum s"""
    def rec(k, cap, left):
        if k == 0:
            if left == 0: yield ()
            return
        for x in range(min(cap, left), -1, -1):
            if x * k < left: break
            for r in rec(k - 1, x, left - x): yield (x,) + r
    return rec(n, s, s)


class Split(object):
    def __init__(self, p, g, ab=(0, 0), K=None, variant='nat', log=None):
        from fq_rows import Rows
        from fq_conifold import Conifold
        from fq_sser import SRows, EdgeLeaf
        from fq_symbol import Symbol
        from fq3_msum import HatRowsU
        from fq3_mvertex import MVertex6
        self.g = g; self.KZ = KZ = 3 * g - 3; self.K = K = (2 * g - 2) if K is None else K; self.variant = variant
        R = self.R = FpRing(p); self.p = p; self.ab = ab
        rows = Rows(R, Nz=KZ); C = self.C = Conifold(R, Nu=K + 4, alpha=Fr(ab[0]), beta=Fr(ab[1])); S = SRows(rows, C, prec_abs=5 * K + 4)
        sym = Symbol(S); assert sym.val_ok, "symbol valuations"
        H = HatRowsU(S, K); assert H.val_ok and H.min_prec >= K, ("hatted rows", H.val_ok, H.min_prec)
        EL = EdgeLeaf(H, H.Q, None)
        def lst(s):
            assert s.prec >= K, "series precision %d < %d" % (s.prec, K)
            c = [s.coef(k) % p for k in range(K)]
            return c if any(c) else None
        self.ER = {(a, b): [[lst(x) for x in row] for row in EL.E(a, b)] for a in range(KZ) for b in range(KZ - a)}
        self.TR = {e: [lst(x) for x in EL.T(e)] for e in range(2, KZ + 2)}
        xs = lst(H.x); self.xpow = [[1] + [0] * (K - 1)]
        for k in range(1, 5 * g): self.xpow.append(self.mul(self.xpow[-1], xs))
        if variant == 'id':
            self.EM = {key: [[None] * 5 for _ in range(5)] for key in self.ER}; self.frozen_ok = None
        else:
            # M_0 = R~ mod zhat^(KZ+1): M_k = C_k (colour matrices of the symbol), M^-1_k = (-1)^k C_k^T (symplectic); checked here
            Cm = sym.colour_matrix(); assert len(Cm) >= KZ + 1
            Minv = [[[Cm[k][j][i] * (1 if k % 2 == 0 else p - 1) % p for j in range(5)] for i in range(5)] for k in range(KZ + 1)]
            for k in range(KZ + 1):
                for i in range(5):
                    for j in range(5):
                        assert sum(Cm[m][i][l] * Minv[k - m][l][j] for m in range(k + 1) for l in range(5)) % p == (1 if (k == 0 and i == j) else 0), "M_0 M_0^-1 != 1"
            class _M(object): pass
            M0 = _M(); M0.Minv = Minv
            Eh = MVertex6(R, M0, kmax=KZ).Eh
            self.frozen_ok = all((self.ER[key][i][j] or [0])[0] == Eh[key][i][j] for key in self.ER for i in range(5) for j in range(5))
            tau = [R.of(x) for x in tau_pred(K)] if variant.startswith('nat') else [1] + [0] * (K - 1)
            if variant.startswith('nat_t'):                       # CONTROL: tau_k -> tau_k + 1 (must destroy regularity)
                kk = int(variant[5:]); tau[kk] = (tau[kk] + 1) % p
            tp = [[1] + [0] * (K - 1)]
            for k in range(1, KZ + 1): tp.append(self.mul(tp[-1], tau))
            self.EM = {(a, b): [[(self.scale(tp[a + b + 1], Eh[(a, b)][i][j]) if Eh[(a, b)][i][j] else None) for j in range(5)] for i in range(5)] for (a, b) in self.ER}
        self.EO = {key: [[self.sub(self.ER[key][i][j], self.EM[key][i][j]) for j in range(5)] for i in range(5)] for key in self.ER}
        self.valO = {key: min(self.val(x) for row in self.EO[key] for x in row) for key in self.EO}
        # graph() keeps only the subsets with fewer than K outer edges: exact iff every outer edge is O(u)
        assert variant == 'id' or min(self.valO.values()) >= 1, "E_out is not O(u): the restriction to fewer than K outer edges would be wrong"
        self.zph = [R.z(4 * j) for j in range(5)]
        self.inner_graphs = legged_graphs(g, 0, KZ); self.outer_keys = [G['key'] for G in self.inner_graphs]
        self.outer_pos = {k: i for i, k in enumerate(self.outer_keys)}
        self._vt = {}

    # ---------------------------------------------------------------- series (lists of K residues, None = 0)
    def mul(self, a, b):
        if a is None or b is None: return None
        K = self.K; out = [0] * K
        for i in range(K):
            x = a[i]
            if x:
                for j in range(K - i):
                    y = b[j]
                    if y: out[i + j] += x * y
        p = self.p; out = [c % p for c in out]
        return out if any(out) else None

    def scale(self, a, c):
        if a is None or not c: return None
        p = self.p; return [x * c % p for x in a]

    def add(self, a, b):
        if a is None: return b
        if b is None: return a
        p = self.p; out = [(x + y) % p for x, y in zip(a, b)]
        return out if any(out) else None

    def sub(self, a, b):
        if b is None: return a
        p = self.p
        out = [(-y) % p for y in b] if a is None else [(x - y) % p for x, y in zip(a, b)]
        return out if any(out) else None

    def val(self, a):
        if a is None: return self.K
        for k, x in enumerate(a):
            if x: return k
        return self.K

    # ---------------------------------------------------------------- vertices of the trivial theory with the leaves of R
    def vertex_table(self, g, hs):
        """{sorted hp (non-increasing tuple of length hs): [5 series]}"""
        key = (g, hs); tab = self._vt.get(key)
        if tab is not None: return tab
        R = self.R; p = self.p; tab = {}; dim0 = 3 * g - 3 + hs
        for m in range(0, max(dim0, 0) + 1):
            dim = dim0 + m; n = 2 * g - 2 + hs + m
            base0 = [self.scale(self.xpow[n - 1], pow(self.zph[j], n, p)) for j in range(5)]
            for leaves in partitions_min2(m, dim, self.KZ + 1):
                mult = {}
                for e in leaves: mult[e] = mult.get(e, 0) + 1
                c0 = R.inv(math.prod(math.factorial(x) for x in mult.values()))
                base = []
                for j in range(5):
                    b = base0[j]
                    for e in leaves: b = self.mul(b, self.TR[e][j])
                    base.append(b)
                if all(b is None for b in base): continue
                for hp in sorted_tuples(hs, dim - sum(leaves)):
                    w = wk(g, tuple(hp) + tuple(leaves))
                    if w == 0: continue
                    c = R.of(w) * c0 % p; vec = tab.get(hp)
                    if vec is None: vec = [None] * 5; tab[hp] = vec
                    for j in range(5): vec[j] = self.add(vec[j], self.scale(base[j], c))
        self._vt[key] = tab
        return tab

    # ---------------------------------------------------------------- one inner graph
    def contract(self, G, S):
        """canonical key of G with the edges NOT in the bitmask S contracted"""
        V = G['V']; ed = G['ed']; par = list(range(V))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for e, (u, w) in enumerate(ed):
            if not (S >> e) & 1: par[find(u)] = find(w)
        roots = sorted(set(find(v) for v in range(V))); cid = {r: i for i, r in enumerate(roots)}
        gv = [0] * len(roots); nv = [0] * len(roots); ne = [0] * len(roots)
        for v in range(V): c = cid[find(v)]; gv[c] += G['gv'][v]; nv[c] += 1
        out = []
        for e, (u, w) in enumerate(ed):
            if (S >> e) & 1: out.append((cid[find(u)], cid[find(w)]))
            else: ne[cid[find(u)]] += 1
        gv = [gv[c] + ne[c] - nv[c] + 1 for c in range(len(roots))]
        return canon(tuple(gv), out, ())

    def plan(self, G):
        """vertex elimination order and, per step, (v, [(edge, frontier position, flip)], positions kept, loops at v)"""
        V = G['V']; ed = G['ed']; nb = [set() for _ in range(V)]
        for u, w in ed:
            if u != w: nb[u].add(w); nb[w].add(u)
        order = []; done = set(); frontier = []; steps = []
        for _ in range(V):
            cand = [v for v in range(V) if v not in done]
            v = min(cand, key=lambda x: (0 if (not frontier or nb[x] & set(frontier)) else 1, len(nb[x] - done - {x}), x))
            conn = []
            for e, (u, w) in enumerate(ed):
                if u == w: continue
                if u == v and w in frontier: conn.append((e, frontier.index(w), True))      # matrix index order [c_v][c_w]
                elif w == v and u in frontier: conn.append((e, frontier.index(u), False))   # [c_u][c_v]
            done.add(v); full = frontier + [v]
            keep = [k for k, f in enumerate(full) if nb[f] - done]
            loops = [e for e, (u, w) in enumerate(ed) if u == w == v]
            steps.append((v, conn, keep, loops)); frontier = [full[k] for k in keep]
        assert not frontier
        return steps

    def colour_sum(self, steps, vecs, mats):
        mul = self.mul; add = self.add
        cur = {(): [1] + [0] * (self.K - 1)}
        for v, conn, keep, loops in steps:
            new = {}; vec = vecs[v]
            for key, val in cur.items():
                for c in range(5):
                    x = mul(val, vec[c])
                    for e in loops:
                        if x is None: break
                        x = mul(x, mats[e][c][c])
                    for e, pos, flip in conn:
                        if x is None: break
                        x = mul(x, mats[e][c][key[pos]] if flip else mats[e][key[pos]][c])
                    if x is None: continue
                    full = key + (c,); nk = tuple(full[k] for k in keep)
                    new[nk] = add(new.get(nk), x)
            cur = new
            if not cur: return None
        return cur.get(())

    def graph(self, t):
        """{outer index: [K coefficients]} of the inner graph number t (all outer-edge subsets, all decorations), weight 1/|Aut|"""
        G = self.inner_graphs[t]; V = G['V']; ed = G['ed']; hsAt = G['hsAt']; E = len(ed); K = self.K; p = self.p
        opts = []
        for v in range(V):
            tab = self.vertex_table(G['gv'][v], len(hsAt[v])); lst = []
            for hp, vec in tab.items():
                for perm in set(itertools.permutations(hp)): lst.append((perm, vec))
            if not lst: return {}
            opts.append(lst)
        if self.variant == 'id': subsets = [((1 << E) - 1, t)]
        else: subsets = [(S, self.outer_pos[self.contract(G, S)]) for S in range(1 << E) if bin(S).count("1") <= POPMAX and bin(S).count("1") < K]
        steps = self.plan(G); acc = {}
        for combo in itertools.product(*opts):
            hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
            keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
            vecs = [combo[v][1] for v in range(V)]
            vals = [self.valO[k] for k in keys]; inner = [self.EM[k] for k in keys]; outer = [self.EO[k] for k in keys]
            for S, pos in subsets:
                if sum(vals[e] for e in range(E) if (S >> e) & 1) >= K: continue
                z = self.colour_sum(steps, vecs, [outer[e] if (S >> e) & 1 else inner[e] for e in range(E)])
                if z is None: continue
                a = acc.get(pos)
                if a is None: acc[pos] = list(z)
                else:
                    for k in range(K): a[k] += z[k]
        ia = self.R.inv(G['aut'])
        return {pos: [x * ia % p for x in a] for pos, a in acc.items()}
