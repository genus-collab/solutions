# -*- coding: utf-8 -*-
"""rf_s15_sum.py -- referee engine (own code): the flat P-action graph sum of P.Omega^mod_{v(u)}, graph by graph, computed as ONE
Givental sum over composite graphs Gamma' whose edges are 2-coloured (inner / outer), in the idempotent frame of the FROZEN model.

Why this is the flat sum (standard composition of R-actions, arXiv:1303.1043 s.2): the vertex of the flat sum is a model correlator,
which is itself a Givental sum (Teleman, arXiv:0712.0160) with R = R_frozen; expanding every flat vertex gives composite graphs
(Gamma', S), S = the outer edges, with weight 1/|Aut(Gamma', S)|; contracting the inner edges gives the outer graph Gamma.
  inner edge  E_B(psi, psi') = (eta_c^-1 - B^-1(psi) eta_c^-1 B^-1(psi')^T)/(psi + psi'),   B = frozen R, eta_c^-1 = diag(Delta);
  outer edge  Eo = v^2 (L(psi) x L(psi')) Ehat,  L(psi) = B^-1(psi) Psi^-1,  Ehat = v-weighted flat edge (rf_s15_outer);
  vertex      Delta_i^(h-1) Sum_{leaf multisets} WK(h; exps, leaves) prod X_i(k)/mult!,   X = z(1 - B^-1(z) 1) + L(z) (d vhat T)  (rf_s15_outer);
  global      v^(2-2g) d^(2g-2)   (dilaton factor d^(2h-2+n) multiplied over the vertices; v^(2-2h) of the frozen scaling).
Only series with non-negative u-offset occur before the global factor, so plain truncated power series suffice.
"""
import itertools, math, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_s15_base import wk, smul, inv, fr, stable_graphs, canon2, aut2, contract


def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield (); return
    for q in range(min(n, maxpart), 0, -1):
        for r in partitions(n - q, q): yield (q,) + r


def tuples_le(n, s):
    if n == 0:
        yield (); return
    for k in range(s + 1):
        for rest in tuples_le(n - 1, s - k): yield (k,) + rest


class Engine(object):
    def __init__(self, Fz, N, X=None, Eo=None, nodil_leaf=False):
        """Fz: Frozen; N: series length; X[(i, k)]: total translation series (k >= 2) or None (frozen only);
        Eo[((m,i),(n,j))]: outer edge series"""
        self.F = Fz; self.p = Fz.p; self.N = N; p = self.p; K = Fz.K
        Bi = Fz.Ri; D = Fz.Delta
        # inner edges (constants)
        f = {}
        for r in range(K + 1):
            for s in range(K + 1 - r):
                if r == 0 and s == 0: continue
                f[(r, s)] = [[(-sum(Bi[r][i][l] * D[l] * Bi[s][j][l] for l in range(5))) % p for j in range(5)] for i in range(5)]
        self.EB = {}
        for m in range(K):
            for n in range(K - m):
                for i in range(5):
                    for j in range(5):
                        val = sum((-1) ** t * f[(m + 1 + t, n - t)][i][j] for t in range(n + 1)) % p
                        if val: self.EB[((m, i), (n, j))] = val
        # symplecticity check of the division: Sum_t (-1)^t f_(M-t, t) = 0
        for Mt in range(1, K + 1):
            for i in range(5):
                for j in range(5):
                    assert sum((-1) ** t * f[(Mt - t, t)][i][j] for t in range(Mt + 1)) % p == 0, "E_B not divisible"
        TB = {}
        for k in range(2, K + 1):
            for i in range(5):
                val = (-sum(Bi[k - 1][i][l] for l in range(5))) % p
                if val: TB[(i, k)] = val
        self.X = {}
        keys = set(TB) | set(X or {})
        for key in keys:
            s = [0] * N; s[0] = TB.get(key, 0)
            if X and key in X: s = [(a + b) % p for a, b in zip(s, X[key])]
            if any(s): self.X[key] = s
        self.Eo = Eo or {}
        self.vcache = {}

    def vertex(self, h, exps, i):
        key = (h, tuple(sorted(exps)), i)
        if key in self.vcache: return self.vcache[key]
        p = self.p; N = self.N; n = len(exps); budget = 3 * h - 3 + n - sum(exps); out = None
        if budget >= 0:
            for part in partitions(budget):
                leaves = [q + 1 for q in part]
                if 2 * h - 2 + n + len(leaves) <= 0: continue
                if any((i, k) not in self.X for k in leaves): continue
                w = wk(h, tuple(sorted(tuple(exps) + tuple(leaves))))
                if not w: continue
                c = fr(w, p)
                for k in set(leaves): c = c * inv(math.factorial(leaves.count(k)), p) % p
                s = [c] + [0] * (N - 1)
                for k in leaves: s = smul(s, self.X[(i, k)], p, N)
                out = s if out is None else [(a + b) % p for a, b in zip(out, s)]
        if out is not None:
            dh = pow(self.F.Delta[i], h - 1, p) if h >= 1 else inv(self.F.Delta[i], p)
            out = [x * dh % p for x in out]
            if not any(out): out = None
        self.vcache[key] = out
        return out

    def composite_value(self, h, ed, outer):
        """sum over decorations and colourings of the labelled composite graph (h, ed) with outer edge set 'outer'"""
        p = self.p; N = self.N; V = len(h)
        hs = [[] for _ in range(V)]
        for e, (a, b) in enumerate(ed): hs[a].append(2 * e); hs[b].append(2 * e + 1)
        decs = [list(tuples_le(len(hs[v]), 3 * h[v] - 3 + len(hs[v]))) for v in range(V)]
        total = [0] * N
        for dec in itertools.product(*decs):
            ex = {}
            for v in range(V):
                for r, hh in enumerate(hs[v]): ex[hh] = dec[v][r]
            for col in itertools.product(range(5), repeat=V):
                const = 1; ser = []; ok = True
                for v in range(V):
                    x = self.vertex(h[v], dec[v], col[v])
                    if x is None: ok = False; break
                    ser.append(x)
                if not ok: continue
                for e, (a, b) in enumerate(ed):
                    k = ((ex[2 * e], col[a]), (ex[2 * e + 1], col[b]))
                    if e in outer:
                        x = self.Eo.get(k)
                        if x is None: ok = False; break
                        ser.append(x)
                    else:
                        x = self.EB.get(k)
                        if x is None: ok = False; break
                        const = const * x % p
                if not ok: continue
                ser.sort(key=lambda s: sum(1 for y in s if y))
                acc = [const] + [0] * (N - 1)
                for s in ser:
                    acc = smul(acc, s, p, N)
                total = [(a + b) % p for a, b in zip(total, acc)]
        return total

    def outer_sums(self, g, only=None, log=None):
        """{outer canonical key: series} (before the global factor), plus bookkeeping"""
        p = self.p; N = self.N; res = {}; seen = set(); info = {}
        self.full = {}      # regrouping by the UNCONTRACTED graph Gamma' = the full R_tot-action term of Gamma' (PPZ Prop. 2.4)
        for (h, ed) in stable_graphs(g):
            E = len(ed); fkey = canon2(h, ed, ())
            for mask in range(1 << E):
                outer = frozenset(e for e in range(E) if mask >> e & 1)
                ein = [ed[e] for e in range(E) if e not in outer]; eout = [ed[e] for e in outer]
                key2 = canon2(h, ein, eout)
                if key2 in seen: continue
                seen.add(key2)
                ho, eo = contract(h, ed, [e for e in range(E) if e not in outer])
                okey = canon2(ho, eo, ())
                if only is not None and okey not in only: continue
                if outer and not self.Eo: continue
                val = self.composite_value(h, ed, outer)
                a = aut2(h, ein, eout); ia = inv(a, p)
                val = [x * ia % p for x in val]
                res[okey] = [(x + y) % p for x, y in zip(res.get(okey, [0] * N), val)]
                self.full[fkey] = [(x + y) % p for x, y in zip(self.full.get(fkey, [0] * N), val)]
                info[okey] = info.get(okey, 0) + 1
        return res, info
