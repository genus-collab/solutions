# -*- coding: utf-8 -*-
"""fq_msum.py -- the outer graph sum R'.(Omega^M) of the formal quintic at genus 2 (genus2_symplectic, step 3):
per stable graph Gamma (the 7 leg-free genus-2 graphs of fq_graphs.load_graphs), summed over ALL colourings of its
half-edges and dilaton leaves, as a Laurent series in s (conifold frame).

Term of a labelled decoration (psi-exponents a_h on the half-edges, ordered dilaton tuples (d_1..d_m) with d >= 2 at
each vertex, psi-DEFICIT allowed: Sum_h a_h + Sum d <= 3g_v-3+n_v+m_v) at a colouring of all half-edges and leaves:
  (1/|Aut Gamma|) prod_v [ (1/m_v!) Omega^M_{g_v, n_v+m_v}(e~_{i_h}..., e~_{i_l}...; a_h..., d_l...) ]  prod_e E'_{i_h i_h'}(a_h, a_h')  prod_l T'^{(i_l)}_{d_l}
with Omega^M_{g,n}(...) = x^{2g-2+n} OmegaHat (fq_mvertex), E', T' from fq_sser.EdgeLeaf.  The outer dilaton leaves are
treated as extra LEGS of the vertex tensor carrying the vector T'(psi) = Sum_d psi^d T'_d (Omega^M is multilinear, so the leaf
colour is summed against T'^{(i)}_d); this is the PPZ kappa-leaf term of R'.(Omega^M) integrated over Mbar_{g_v, n_v+m_v}
(no push-forward needed for an integral), exactly as gsum_q.term_value does for R.omega.  The sum over all labelled
decorations (not orbit representatives) with the 1/|Aut| weight is as in gsum_q.coloured_sum.
Contraction: tensors (labels = half-edge ids, flat list of Ser, index Sum_l i_l 5^l); vertex tensors are built from
OmegaHat (F_p) by contracting the leaf legs with T'_d and scaling by x^{2g-2+n}; then loops are contracted with E', bridges
applied as matrices and the next vertex merged over the common labels.  Precision is tracked by fq_series.Ser.
colmask (validation only, M = 1): restrict every vertex tensor to the entries with all legs of colour c_v, giving the
coloured graph term of gsum_q (the vertex colouring is then the colouring of its legs).
API: outer_graphs(); outer_decorations(G); GraphSummer(S, M).graph(G, colmask=None) -> Ser; .all_graphs() -> {idx: Ser};
     sum_series(list) -> Ser.
"""
import itertools, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_graphs import load_graphs, edges_of, aut_count, compositions, tuples_sum_le
from fq_sser import rows_prime, EdgeLeaf
from fq_mvertex import MVertex


def outer_graphs():
    out = []
    for gi, (gv, m, A) in enumerate(load_graphs(2)):
        V = len(gv); ed = edges_of(A)
        H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
        hvtx = [ed[e][s] for e, s in H]
        hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
        out.append(dict(idx=gi + 1, gv=gv, ed=ed, H=H, hsAt=hsAt, V=V, aut=aut_count(gv, A)))
    return out


def vertex_options_deficit(g, hs):
    """(hp, dils) with Sum hp + Sum dils <= 3g-3+hs+nd (psi-deficit allowed), dils ordered tuples with parts >= 2"""
    res = []
    for nd in range(0, max(0, 3 * g - 3 + hs) + 1):
        dim = 3 * g - 3 + hs + nd
        if dim < 2 * nd: continue
        dilss = [()] if nd == 0 else [c for s in range(2 * nd, dim + 1) for c in compositions(s, nd, 2)]
        for dils in dilss:
            for hp in tuples_sum_le(hs, dim - sum(dils)):
                res.append((tuple(hp), tuple(dils)))
    return res


def outer_decorations(G):
    """all labelled decorations: tuples over vertices of (hp, dils)"""
    opts = [vertex_options_deficit(G['gv'][v], len(G['hsAt'][v])) for v in range(G['V'])]
    return list(itertools.product(*opts))


# ------------------------------------------------------------------ tensors of Ser
def _pos(idx, n): return sum(i * 5 ** l for l, i in enumerate(idx))


def contract_pair(T, labels, h0, h1, E):
    """both labels in T: out[rest] = Sum_{i,j} T[.., i at h0, .., j at h1, ..] E[i][j]"""
    p0 = labels.index(h0); p1 = labels.index(h1); n = len(labels)
    rest = [l for l in labels if l not in (h0, h1)]
    out = {}
    for idx in itertools.product(range(5), repeat=n):
        t = T[_pos(idx, n)]
        if not t.c: continue
        key = tuple(idx[k] for k in range(n) if k not in (p0, p1))
        term = t.mul(E[idx[p0]][idx[p1]])
        out[key] = term if key not in out else out[key].add(term)
    res = []
    for idx in itertools.product(range(5), repeat=len(rest)):
        res.append(out.get(idx, _zero_like(T)))
    return res, rest


def apply_matrix(T, labels, h0, E, newlabel):
    """out[.., j at position of h0, ..] = Sum_i T[.., i, ..] E[i][j]; label h0 -> newlabel"""
    p0 = labels.index(h0); n = len(labels)
    out = [None] * len(T)
    for idx in itertools.product(range(5), repeat=n):
        acc = None
        for i in range(5):
            idx2 = idx[:p0] + (i,) + idx[p0 + 1:]
            t = T[_pos(idx2, n)]
            if not t.c: continue
            term = t.mul(E[i][idx[p0]])
            acc = term if acc is None else acc.add(term)
        out[_pos(idx, n)] = acc if acc is not None else _zero_like(T)
    return out, [newlabel if l == h0 else l for l in labels]


def merge(T1, l1, T2, l2):
    """product summed over the common labels; out labels = (l1 minus common) + (l2 minus common)"""
    common = [l for l in l1 if l in l2]
    r1 = [l for l in l1 if l not in common]; r2 = [l for l in l2 if l not in common]
    n1 = len(l1); n2 = len(l2)
    out = {}
    for idx1 in itertools.product(range(5), repeat=n1):
        t1 = T1[_pos(idx1, n1)]
        if not t1.c: continue
        cvals = {l: idx1[l1.index(l)] for l in common}
        k1 = tuple(idx1[l1.index(l)] for l in r1)
        for idx2r in itertools.product(range(5), repeat=len(r2)):
            idx2 = tuple(cvals[l] if l in common else idx2r[r2.index(l)] for l in l2)
            t2 = T2[_pos(idx2, n2)]
            if not t2.c: continue
            key = k1 + idx2r
            term = t1.mul(t2)
            out[key] = term if key not in out else out[key].add(term)
    rest = r1 + r2
    res = [out.get(idx, _zero_like(T1)) for idx in itertools.product(range(5), repeat=len(rest))]
    return res, rest


def _zero_like(T):
    prec = min(t.prec for t in T) if T else 10 ** 6
    return Ser(T[0].R, prec, [])


def sum_series(lst):
    acc = None
    for s in lst: acc = s if acc is None else acc.add(s)
    return acc


class GraphSummer(object):
    def __init__(self, S, M):
        self.S = S; self.M = M; self.R = S.R
        self.Qp = rows_prime(S, M)
        self.EL = EdgeLeaf(S, self.Qp, M)
        self.MV = MVertex(S.R, M)
        self.xpow = {0: Ser.const(S.R, 1, S.x.prec)}

    def _x(self, k):
        r = self.xpow.get(k)
        if r is None:
            r = self.S.x.pow(k); self.xpow[k] = r
        return r

    def vertex_tensor(self, g, hp, dils, cv=None):
        """Omega^M_{g, hs+m}(.. ; hp, dils) with the leaf legs contracted against T'_d: Ser tensor over the hs half-edge legs"""
        R = self.R; p = R.p
        hs = len(hp); m = len(dils); n = hs + m
        Om = self.MV.tensor(g, tuple(hp) + tuple(dils))
        if cv is not None:
            Om = [Om[_pos(idx, n)] if len(set(idx)) <= 1 and (n == 0 or idx[0] == cv) else 0 for idx in itertools.product(range(5), repeat=n)]
        xk = self._x(2 * g - 2 + n).scale(R.inv(math.factorial(m)))
        # contract leaves (last legs) one at a time
        cur = Om; is_const = True; nn = n
        for r in range(m - 1, -1, -1):
            vec = self.EL.T(dils[r])
            new = []
            for idx in itertools.product(range(5), repeat=nn - 1):
                acc = None
                for k in range(5):
                    t = cur[_pos(idx + (k,), nn)]
                    if is_const:
                        if not t: continue
                        term = vec[k].scale(t)
                    else:
                        if not t.c: continue
                        term = t.mul(vec[k])
                    acc = term if acc is None else acc.add(term)
                new.append(acc if acc is not None else Ser(R, vec[0].prec, []))
            cur = new; is_const = False; nn -= 1
        if is_const:
            cur = [xk.scale(t) if t else Ser(R, xk.prec, []) for t in cur]
        else:
            cur = [t.mul(xk) for t in cur]
        return cur

    def graph(self, G, colmask=None, verbose=False):
        R = self.R; V = G['V']; ed = G['ed']; hsAt = G['hsAt']; gv = G['gv']
        total = None
        for dec in outer_decorations(G):
            hp = [0] * (2 * len(ed))
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = dec[v][0][r]
            cur = None; labels = []; done_edges = set()
            for v in range(V):
                Tv = self.vertex_tensor(gv[v], dec[v][0], dec[v][1], None if colmask is None else colmask[v])
                lv = list(hsAt[v])
                if cur is None:
                    cur, labels = Tv, lv
                else:
                    # bridges from the current tensor to v: apply E' on the current side, relabel to the half-edge at v
                    for e in range(len(ed)):
                        u, w = ed[e]
                        if e in done_edges or u == w: continue
                        h0, h1 = 2 * e, 2 * e + 1
                        if h0 in labels and h1 in lv:
                            cur, labels = apply_matrix(cur, labels, h0, self.EL.E(hp[h0], hp[h1]), h1); done_edges.add(e)
                        elif h1 in labels and h0 in lv:
                            Et = self.EL.E(hp[h0], hp[h1]); EtT = [[Et[i][j] for i in range(5)] for j in range(5)]
                            cur, labels = apply_matrix(cur, labels, h1, EtT, h0); done_edges.add(e)
                    cur, labels = merge(cur, labels, Tv, lv)
                # loops (both half-edges now in cur)
                for e in range(len(ed)):
                    if e in done_edges: continue
                    h0, h1 = 2 * e, 2 * e + 1
                    if h0 in labels and h1 in labels:
                        cur, labels = contract_pair(cur, labels, h0, h1, self.EL.E(hp[h0], hp[h1])); done_edges.add(e)
            assert not labels and len(done_edges) == len(ed), (G['idx'], labels, done_edges)
            term = cur[0]
            total = term if total is None else total.add(term)
        return total.scale(R.inv(G['aut']))

    def all_graphs(self, graphs=None, verbose=None):
        if graphs is None: graphs = outer_graphs()
        out = {}
        for G in graphs:
            out[G['idx']] = self.graph(G)
            if verbose: verbose("   graph %d: val %s prec %d" % (G['idx'], out[G['idx']].v if out[G['idx']].c else None, out[G['idx']].prec))
        return out
