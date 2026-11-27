# -*- coding: utf-8 -*-
"""fq_mvertex.py -- the colour-mixing vertex CohFT Omega^M := M.omega of the formal quintic (genus2_symplectic, step 2).

omega = the TFT in the normalised canonical basis e~_i (eta = identity): omega_{g,n}(e~_i^{(x)n}) = Delta_i^{(2g-2+n)/2},
Delta_i^{1/2} = zeta^{4i} x with x = s f_0(s^5) (fq_sser), unit 1 = Sum_k Delta_k^{-1/2} e~_k.  M(z) = exp(m_1 z + m_2 z^2 + m_3 z^3)
constant and symplectic (fq_sser.SympM).  PPZ action (same conventions as the Teleman R-action of gsum_q / FRAMEB):
  Omega^M_{g,n}(e~_{i_1},..,e~_{i_n}) prod psi_l^{d_l}  =  Sum_{stable graphs Gamma' of genus g with n legs} (1/|Aut Gamma'|) Sum_{decorations}
     prod_v [ omega_{g_v, n_v+m_v}(e~_{j_v}...) WK(g_v; exponents at v) / m_v! ]
     prod_legs (M^{-1}_{k_l})_{j_{v(l)} i_l}  [leg l carries psi^{d_l + k_l}]
     prod_edges Ehat_{j j'}(a, b)  [= [psi^a psi'^b] (delta - M^{-1}(psi) M^{-1}(psi')^T)_{jj'}/(psi+psi'), half-edges psi^a, psi'^b]
     prod_kappa-leaves That_{j, e}  [= [z^e] z(1 - M^{-1}(z) 1)_j = -(M^{-1}_{e-1} 1)_j, e >= 2, the leaf carries psi^e],
  with Sum of exponents at v = 3g_v - 3 + n_v + m_v (TFT vertices = point class, fq_hodge.wk).
  Weight bookkeeping: Sum of all M-orders = D := 3g-3+n - Sum d_l (the psi-deficit), so inner graphs have <= D edges.
HOMOGENEITY: Sum_v (2g_v-2+n_v+m_v) = 2g-2+n+m, each kappa-leaf carries Delta^{-1/2} = x^{-1} zeta^{-4k}: every term is x^{2g-2+n}
  times an F_p constant.  So  Omega^M_{g,n}(e~_{i_1}..e~_{i_n}; d) = x^{2g-2+n} OmegaHat(i_1..i_n; d)  with OmegaHat computed over
  F_p from the constant TFT Delta^hat_j^{1/2} = zeta^{4j}, unit 1^hat = Sum_k zeta^{-4k} e~_k.  MVertex(R, M).tensor(g, d) returns
  OmegaHat as a flat list of length 5^n (index Sum_l i_l 5^l), cached per (g, d).
Legged stable graphs: legged_graphs(g, n, Emax) enumerates isomorphism classes (vertex genera, edge multiset with loops, leg ->
  vertex map; legs labelled) with at most Emax edges, by canonical forms under S_V; |Aut| = number of (vertex perm, edge perm,
  loop flips) preserving everything and fixing the legs (fq_graphs.autos_of filtered).
Checks (driver): M = 1 gives omega (diagonal, WK); M = diag(exp(a1_i z + a2_i z^3)) gives zeta^{4i(2g-2+n)} PhiVertexTable(fq_phi);
  symmetry; string and dilaton equations (flat unit of M.omega).
"""
import itertools, math, os, sys
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_hodge import wk
from fq_graphs import autos_of, compositions, tuples_sum_le


def _connected(V, ed):
    adj = {v: set() for v in range(V)}
    for u, v in ed: adj[u].add(v); adj[v].add(u)
    seen = {0}; stack = [0]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen: seen.add(w); stack.append(w)
    return len(seen) == V


def _canon(V, gv, ed, legs):
    best = None
    for sig in itertools.permutations(range(V)):
        gv2 = [0] * V
        for v in range(V): gv2[sig[v]] = gv[v]
        ed2 = tuple(sorted(tuple(sorted((sig[u], sig[v]))) for u, v in ed))
        legs2 = tuple(sig[l] for l in legs)
        key = (tuple(gv2), ed2, legs2)
        if best is None or key < best: best = key
    return best


@lru_cache(maxsize=None)
def legged_graphs(g, n, Emax):
    """list of dicts: V, gv, ed (list of (u,v), u <= v, loops allowed), legs (vertex of each leg), aut, hsAt, legsAt, H"""
    reps = {}
    for E in range(Emax + 1):
        for V in range(1, E + 2):
            G = g - E + V - 1
            if G < 0: continue
            pairs = [(u, v) for u in range(V) for v in range(u, V)]
            for ed in itertools.combinations_with_replacement(pairs, E):
                if not _connected(V, ed): continue
                for gv in compositions(G, V, 0):
                    for legs in itertools.product(range(V), repeat=n):
                        nv = [0] * V
                        for u, v in ed: nv[u] += 1; nv[v] += 1
                        for l in legs: nv[l] += 1
                        if any(2 * gv[v] - 2 + nv[v] <= 0 for v in range(V)): continue
                        key = _canon(V, gv, ed, legs)
                        if key not in reps: reps[key] = True
    out = []
    for (gv, ed, legs) in sorted(reps):
        V = len(gv); ed = list(ed)
        A = [[0] * V for _ in range(V)]
        for u, v in ed:
            if u == v: A[u][u] += 1
            else: A[u][v] += 1; A[v][u] += 1
        auts = [a for a in autos_of(gv, A) if all(a[0][l] == l for l in legs)]
        H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
        hvtx = [ed[e][s] for e, s in H]
        hsAt = [[h for h, x in enumerate(hvtx) if x == v] for v in range(V)]
        legsAt = [[l for l in range(n) if legs[l] == v] for v in range(V)]
        out.append(dict(V=V, gv=gv, ed=ed, legs=legs, aut=len(auts), H=H, hsAt=hsAt, legsAt=legsAt, E=len(ed)))
    return out


class MVertex(object):
    """OmegaHat^M_{g,n}(i_1..i_n; d) over F_p (see module docstring)"""
    def __init__(self, R, M):
        self.R = R; self.M = M; p = R.p; self.p = p
        self.Minv = M.Minv                                        # [z^k] M^-1, k = 0..3
        self.zh = [R.z(4 * j) for j in range(5)]                  # Delta^hat_j^{1/2}
        self.unit = [R.z(-4 * k) for k in range(5)]               # 1^hat components
        self.cache = {}
        # inner edge coefficients Ehat[(a,b)][j][j'] for a + b + 1 <= 3, inner leaves That[e][j] for e - 1 <= 3
        Gh = {}
        for Mo in range(4):
            for No in range(4 - Mo):
                if Mo == 0 and No == 0: Gh[(0, 0)] = [[0] * 5 for _ in range(5)]; continue
                Gh[(Mo, No)] = [[(-sum(self.Minv[Mo][j][k] * self.Minv[No][jj][k] for k in range(5))) % p for jj in range(5)] for j in range(5)]
        self.Eh = {}
        for a in range(3):
            for b in range(3 - a):
                tab = [[0] * 5 for _ in range(5)]
                for t in range(b + 1):
                    Gt = Gh[(a + 1 + t, b - t)]; sg = 1 if t % 2 == 0 else p - 1
                    for j in range(5):
                        for jj in range(5): tab[j][jj] = (tab[j][jj] + sg * Gt[j][jj]) % p
                self.Eh[(a, b)] = tab
        self.Th = {e: [(-sum(self.Minv[e - 1][j][k] * self.unit[k] for k in range(5))) % p for j in range(5)] for e in range(2, 5)}

    def _vertex_options(self, g, hs, legd, D):
        """list of (m, leaves, ks, hp, wkval): leaves ordered tuple (e >= 2), ks M-orders of the legs, hp half-edge exponents,
        with Sum(legd + ks) + Sum hp + Sum leaves = 3g-3+hs+nlegs+m and total M-weight Sum ks + Sum(e-1) <= D"""
        R = self.R; p = self.p
        nl = len(legd); out = []
        dim0 = 3 * g - 3 + hs + nl
        for m in range(0, D + 1):
            dim = dim0 + m; rem = dim - sum(legd)
            if rem < 2 * m or rem < 0: continue
            leaf_opts = [()] if m == 0 else [c for s in range(2 * m, min(rem, 2 * m + D) + 1) for c in compositions(s, m, 2)]
            for leaves in leaf_opts:
                if sum(e - 1 for e in leaves) > D: continue
                left = rem - sum(leaves)
                for ks in tuples_sum_le(nl, min(left, D - sum(e - 1 for e in leaves))):
                    left2 = left - sum(ks)
                    for hp in tuples_sum_le(hs, left2):
                        if sum(hp) != left2: continue
                        exps = tuple(d + k for d, k in zip(legd, ks)) + tuple(hp) + tuple(leaves)
                        w = wk(g, exps)
                        if w == 0: continue
                        c = R.of(w) * R.inv(math.factorial(m)) % p
                        out.append((m, leaves, ks, hp, c))
        return out

    def tensor(self, g, d):
        """OmegaHat_{g,n}(i_1..i_n; d), flat list of length 5^n, index Sum_l i_l 5^l"""
        d = tuple(d); key = (g, d)
        r = self.cache.get(key)
        if r is not None: return r
        R = self.R; p = self.p; n = len(d)
        D = 3 * g - 3 + n - sum(d)
        assert D <= 3, "psi-deficit %d > 3: M-orders beyond z^3 would be needed (genus-2 programme only)" % D
        size = 5 ** n
        out = [0] * size
        if D >= 0 and 2 * g - 2 + n > 0:
            for Gm in legged_graphs(g, n, D):
                self._add_graph(Gm, d, D, out)
        self.cache[key] = out
        return out

    def _add_graph(self, Gm, d, D, out):
        R = self.R; p = self.p; n = len(d)
        V = Gm['V']; gv = Gm['gv']; ed = Gm['ed']; hsAt = Gm['hsAt']; legsAt = Gm['legsAt']; H = Gm['H']
        opts = [self._vertex_options(gv[v], len(hsAt[v]), tuple(d[l] for l in legsAt[v]), D) for v in range(V)]
        if any(not o for o in opts): return
        inv_aut = R.inv(Gm['aut'])
        strides = [5 ** l for l in range(n)]
        for dec in itertools.product(*opts):
            # total weight must be D (automatic, but check cheaply)
            hp = [0] * len(H)
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = dec[v][3][r]
            wtot = sum(sum(dec[v][2]) + sum(e - 1 for e in dec[v][1]) for v in range(V)) + sum(hp) + len(ed)
            if wtot != D: continue
            for col in itertools.product(range(5), repeat=V):
                val = inv_aut
                for v in range(V):
                    m, leaves, ks, hpv, c = dec[v]
                    j = col[v]
                    val = val * c % p * R.pow(self.zh[j], 2 * gv[v] - 2 + len(hsAt[v]) + len(legsAt[v]) + m) % p
                    for e in leaves: val = val * self.Th[e][j] % p
                    if not val: break
                if not val: continue
                for e in range(len(ed)):
                    u, v = ed[e]
                    val = val * self.Eh[(hp[2 * e], hp[2 * e + 1])][col[u]][col[v]] % p
                    if not val: break
                if not val: continue
                # legs: outer product of the rows Minv[k_l][j_v][.]
                vecs = []
                for v in range(V):
                    for r, l in enumerate(legsAt[v]):
                        vecs.append((l, self.Minv[dec[v][2][r]][col[v]]))
                if n == 0:
                    out[0] = (out[0] + val) % p; continue
                vecs.sort()
                self._outer_add(out, val, [vc for (l, vc) in vecs], strides)

    def _outer_add(self, out, val, vecs, strides):
        p = self.p; n = len(vecs)
        for idx in itertools.product(range(5), repeat=n):
            c = val
            for l in range(n):
                c = c * vecs[l][idx[l]] % p
                if not c: break
            if c:
                pos = sum(idx[l] * strides[l] for l in range(n))
                out[pos] = (out[pos] + c) % p


def tensor_get(T, idx):
    return T[sum(i * 5 ** l for l, i in enumerate(idx))]


def tensor_permute(T, n, perm):
    """T'(i_{perm[0]}, ..., i_{perm[n-1]}) = T(i_0..i_{n-1}): returns the tensor with legs reordered by perm"""
    out = [0] * len(T)
    for idx in itertools.product(range(5), repeat=n):
        new = tuple(idx[perm[l]] for l in range(n))
        out[sum(i * 5 ** l for l, i in enumerate(new))] = T[sum(i * 5 ** l for l, i in enumerate(idx))]
    return out
