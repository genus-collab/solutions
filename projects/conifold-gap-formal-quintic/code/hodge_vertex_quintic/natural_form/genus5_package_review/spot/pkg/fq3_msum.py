# -*- coding: utf-8 -*-
"""fq3_msum.py -- the outer graph sum R'.(Omega^M) of the formal quintic at GENUS 3 in the HATTED theory, as u-series.

Hatted theory (fq_hat, FORMAL_QUINTIC s.3, OU_EXTENSION s.1), now in the variable u = s^5 directly: HatRowsU(S, K) turns the s-series
rows of fq_sser.SRows into u-series  Qhat[a][i][n] = Q[a][i][n] s^{4n+1-c_a} (asserted to be pure u-series), xhat = x/s, truncated at
O(u^K).  For a genus-g graph the un-hatted product sits at s^{-10(g-1)}, so with ANY symplectic M(zhat; u) the colour-summed term is
   T_Gamma = u^{-2(g-1)} G_Gamma(u),   G_Gamma = the sum computed here (genus 3: u^{-4} G_Gamma, regular iff [u^0..u^3] G_Gamma = 0).
GraphSummer6(rows, M, vt): fq_msum.GraphSummer with (i) rows_prime6 (M-orders <= 6, M^{-1} entries F_p or Ser), (ii) the vertex-tensor
provider vt (.tensor(g, d) -> flat list over 5^n of F_p ints (constant M: fq3_mvertex.MVertex6) or of Ser in u (u-dependent M: UVertex)),
(iii) tensor contractions contract_pair / merge and the leaf contraction in vertex_tensor rewritten with the index convention position = Sum_l i_l 5^l on OUTPUT as well as input (fq_msum's versions write in itertools.product order: harmless at genus 2, wrong at genus 3), (iv) the genus-3 outer graphs (fq3_graphs.outer_graphs3, 42 graphs) with ALL labelled decorations with psi-deficit (fq_msum.
vertex_options_deficit) and all colourings.  Edge and leaf tables: fq_sser.EdgeLeaf on the primed rows (unchanged).
UVertex(R, MUobj, K, nodes): vertex tensors of the u-dependent M(zhat; u) by interpolation in eps of MVertex6(M(eps)) at the nodes
eps = 0..NN-1 (the eps-degree of a weight-D entry is <= 3D <= 18 when M = M_0 exp(eps N_1 + eps^2 N_2 + eps^3 N_3); NN = 19 suffices,
NN = 21 also checks that the eps^19, eps^20 coefficients vanish); the coefficient of eps^e becomes [u^e], e < K.  The node tensors can be
supplied precomputed (node_mv: list of MVertex6 or of dicts {(g,d): tensor}) for parallel production (fq3_par).
frozen: fq_leadsolve.FrozenRows(sym) (constant rows R^, x = 1, precision O(u^1)) works unchanged as 'rows'.
API: HatRowsU; GraphSummer6(rows, M, vt).graph(G) -> Ser (u); .all_graphs(graphs); UVertex; hat_sums6(H, M, vt, graphs); frozen_sums6;
     ucoefs(ser, kmax); needed_keys(graphs) -> set of (g, d) vertex keys; tc_expansion3(acoefs, tc, R) -> t_c^{-4}.. coefficients.
"""
import itertools, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_sser import EdgeLeaf
from fq_msum import vertex_options_deficit, apply_matrix, _pos, _zero_like, sum_series
from fq_symbol import CEXP
from fq_hat import vander_inv
from fq3_graphs import outer_graphs3
from fq3_mvertex import MVertex6, KMAX
from fq3_m import rows_prime6
from gsum_q import is_u_series


# ------------------------------------------------------------------ tensor contractions with a CONSISTENT index convention
# NOTE: fq_msum.contract_pair / merge build their output lists in itertools.product order (last index fastest) while every reader
# (apply_matrix, contract_pair, merge, _pos) uses the convention position = Sum_l i_l 5^l (first index fastest).  For rank <= 1 or
# symmetric intermediate tensors this is invisible (all genus-2 graphs: validated against gsum_q), but at genus 3 asymmetric rank-2
# intermediates occur (e.g. graph 17, the chain g1-g1-g1 with the middle vertex last), so the outputs are written here with _pos.
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
    nr = len(rest); z = _zero_like(T); res = [z] * 5 ** nr
    for key, val in out.items(): res[_pos(key, nr)] = val
    return res, rest


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
    rest = r1 + r2; nr = len(rest); z = _zero_like(T1); res = [z] * 5 ** nr
    for key, val in out.items(): res[_pos(key, nr)] = val
    return res, rest


def s_to_u(ser, K=None):
    """pure u-series in s -> Ser in u (coefficients at s^{5k}); asserts purity"""
    R = ser.R
    assert is_u_series(ser), "not a u-series"
    if not ser.c: return Ser(R, (ser.prec + 4) // 5 if K is None else K, [])
    assert ser.v % 5 == 0
    c = [ser.c[5 * k] for k in range((len(ser.c) + 4) // 5)]
    out = Ser(R, ser.v // 5, c)
    return out.truncate(K) if K is not None else out


class HatRowsU(object):
    """duck-typed rows for GraphSummer6: Qhat[a][i][n] = Q s^{4n+1-c_a} and xhat = x/s as u-series, precision O(u^K)"""
    def __init__(self, S, K):
        R = S.R; self.R = R; self.Nz = S.Nz; self.eta_inv = S.eta_inv; self.S = S; self.K = K
        self.Q = [[[s_to_u(S.Q[a][i][n].shift(4 * n + 1 - CEXP[a]), K) for n in range(S.Nz + 1)] for i in range(5)] for a in range(5)]
        self.x = s_to_u(S.x.shift(-1), K)
        self.xinv = self.x.inv()
        allq = [self.Q[a][i][n] for a in range(5) for i in range(5) for n in range(S.Nz + 1)]
        self.min_prec = min(q.prec for q in allq + [self.x])
        self.val_ok = all(q.c and q.v == 0 for q in allq) and self.x.v == 0 and self.x.c[0] == R.one

    def symbol_check(self, sym):
        return all(self.Q[a][i][n].coef(0) == sym.Rh[a][i][n] for a in range(5) for i in range(5) for n in range(self.Nz + 1))


def outer_decorations(G):
    opts = [vertex_options_deficit(G['gv'][v], len(G['hsAt'][v])) for v in range(G['V'])]
    return list(itertools.product(*opts))


def needed_keys(graphs):
    keys = set()
    for G in graphs:
        for v in range(G['V']):
            for hp, dils in vertex_options_deficit(G['gv'][v], len(G['hsAt'][v])): keys.add((G['gv'][v], tuple(hp) + tuple(dils)))
    return keys


class GraphSummer6(object):
    def __init__(self, rows, M, vt):
        self.S = rows; self.M = M; self.R = rows.R; self.vt = vt
        self.Qp = rows_prime6(rows, M)
        self.EL = EdgeLeaf(rows, self.Qp, M)
        self.xpow = {0: Ser.const(rows.R, 1, rows.x.prec)}

    def _x(self, k):
        r = self.xpow.get(k)
        if r is None: r = self.S.x.pow(k); self.xpow[k] = r
        return r

    def vertex_tensor(self, g, hp, dils):
        R = self.R; p = R.p
        hs = len(hp); m = len(dils); n = hs + m
        Om = self.vt.tensor(g, tuple(hp) + tuple(dils))
        is_const = not (Om and isinstance(Om[0], Ser))
        xk = self._x(2 * g - 2 + n).scale(R.inv(math.factorial(m)))
        cur = Om; nn = n
        for r in range(m - 1, -1, -1):
            # contract the LAST leg (leaf dils[r]) with T'_{dils[r]}; the output is written with _pos (fq_msum.vertex_tensor appends in
            # itertools.product order, which is read back with _pos: harmless at genus 2 where the remaining legs are symmetric, wrong here)
            vec = self.EL.T(dils[r]); zero = Ser(R, vec[0].prec, []); new = [zero] * 5 ** (nn - 1)
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
                if acc is not None: new[_pos(idx, nn - 1)] = acc
            cur = new; is_const = False; nn -= 1
        if is_const: cur = [xk.scale(t) if t else Ser(R, xk.prec, []) for t in cur]
        else: cur = [t.mul(xk) for t in cur]
        return cur

    def graph(self, G):
        R = self.R; V = G['V']; ed = G['ed']; hsAt = G['hsAt']; gv = G['gv']
        total = None
        for dec in outer_decorations(G):
            hp = [0] * (2 * len(ed))
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = dec[v][0][r]
            cur = None; labels = []; done = set()
            for v in range(V):
                Tv = self.vertex_tensor(gv[v], dec[v][0], dec[v][1]); lv = list(hsAt[v])
                if cur is None: cur, labels = Tv, lv
                else:
                    for e in range(len(ed)):
                        u, w = ed[e]
                        if e in done or u == w: continue
                        h0, h1 = 2 * e, 2 * e + 1
                        if h0 in labels and h1 in lv:
                            cur, labels = apply_matrix(cur, labels, h0, self.EL.E(hp[h0], hp[h1]), h1); done.add(e)
                        elif h1 in labels and h0 in lv:
                            Et = self.EL.E(hp[h0], hp[h1]); EtT = [[Et[i][j] for i in range(5)] for j in range(5)]
                            cur, labels = apply_matrix(cur, labels, h1, EtT, h0); done.add(e)
                    cur, labels = merge(cur, labels, Tv, lv)
                for e in range(len(ed)):
                    if e in done: continue
                    h0, h1 = 2 * e, 2 * e + 1
                    if h0 in labels and h1 in labels:
                        cur, labels = contract_pair(cur, labels, h0, h1, self.EL.E(hp[h0], hp[h1])); done.add(e)
            assert not labels and len(done) == len(ed), (G['idx'], labels, done)
            term = cur[0]
            total = term if total is None else total.add(term)
        return total.scale(R.inv(G['aut']))

    def all_graphs(self, graphs=None, verbose=None):
        if graphs is None: graphs = outer_graphs3()
        out = {}
        for G in graphs:
            out[G['idx']] = self.graph(G)
            if verbose: verbose("   graph %d (E=%d): %s" % (G['idx'], G['E'], out[G['idx']].fmt(3)))
        return out


# ------------------------------------------------------------------ u-dependent vertex tensors by eps-interpolation
class UVertex(object):
    """vertex tensors of M(zhat; u) (fq3_m.MU) as Ser in u: interpolation of MVertex6(M(eps)) over eps = nodes"""
    def __init__(self, R, MUobj, K, NN=19, node_mv=None, degree_check=True):
        self.R = R; p = R.p; self.K = K; self.nodes = list(range(NN)); self.Vinv = vander_inv(self.nodes, p)
        self.mv = node_mv if node_mv is not None else [MVertex6(R, MUobj.at(e)) for e in self.nodes]
        self.degree_check = degree_check; self.maxdeg = 0; self.cache = {}

    def _node_tensor(self, t, g, d):
        mv = self.mv[t]
        return mv.tensor(g, d) if hasattr(mv, 'tensor') else mv[(g, d)]

    def tensor(self, g, d):
        d = tuple(d); key = (g, d); r = self.cache.get(key)
        if r is not None: return r
        R = self.R; p = R.p; NN = len(self.nodes); K = self.K
        Ts = [self._node_tensor(t, g, d) for t in range(NN)]
        size = len(Ts[0]); out = []
        for pos in range(size):
            vals = [Ts[t][pos] for t in range(NN)]
            if not any(vals): out.append(Ser(R, K, [])); continue
            coeffs = [sum(self.Vinv[e][t] * vals[t] for t in range(NN)) % p for e in range(NN)]
            deg = max([e for e, c in enumerate(coeffs) if c] or [0]); self.maxdeg = max(self.maxdeg, deg)
            if self.degree_check and NN > 3 * KMAX + 1: assert deg <= 3 * KMAX, "eps-degree %d > %d at %s" % (deg, 3 * KMAX, key)
            out.append(Ser(R, 0, coeffs[:K]))
        self.cache[key] = out
        return out


def hat_sums6(H, M, vt, graphs=None, verbose=None):
    return GraphSummer6(H, M, vt).all_graphs(graphs, verbose)


def frozen_sums6(frozen, M, graphs):
    """{idx: F_p} frozen (u = 0) sums for a constant M (MConst)"""
    GS = GraphSummer6(frozen, M, MVertex6(frozen.R, M)); out = {}
    for G in graphs:
        s = GS.graph(G); out[G['idx']] = s.coef(0) if s.prec > 0 else 0
    return out


def ucoefs(ser, kmax):
    return [ser.coef(k) if k < ser.prec else None for k in range(kmax + 1)]


def tc_expansion3(acoefs, tc, R, g=3):
    """u^{-2(g-1)} Sum_k a_k u^k rewritten in t_c: coefficients of t_c^{-2(g-1)}, .., as far as the data allow"""
    G = Ser(R, 0, [R.of(a) for a in acoefs]); u_of_t = tc.reverse()
    T = G.compose(u_of_t).mul(u_of_t.inv().pow(2 * (g - 1)))
    return [T.coef(k) for k in range(-2 * (g - 1), T.prec)]
