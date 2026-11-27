# -*- coding: utf-8 -*-
"""un_engine.py -- the outer graph sum R'.(M.omega^) by edge splitting (fq4_split.Split, imported read-only) for an ARBITRARY table
of hatted rows of the full R = R' M, with M = M_0(tau(u) zhat) the natural vertex factor of the formal quintic, and with an optional
formal parameter eps kept to order eps^(L-1) (dual numbers for L = 2: exact first-order responses).

Var(base, Qh, L): base = fq4_split.Split(p, g, ab, K, 'nat') (supplies E_M, the inner graphs, the Witten-Kontsevich vertices and the
vertex-elimination colour sum); Qh = rows as (u, eps)-series (un_core); the tables E_R, T_R, xhat are rebuilt from Qh,
E_out = E_R - E_M.  Var.graph(t) -> {outer index: series}; nothing of base is modified (its dict is copied).
Differences from Split.graph: series have KU*L entries; the subsets S of outer edges are restricted by a valuation bound that
allows eps^(>=1) parts of E_out of valuation 0 (Split assumes every outer edge is O(u)); option maxout.

run(spec) / run_pool(spec, nproc): spec = (prime, g, (alpha, beta) strings, K, L, Qh or None, maxout); Qh None = natural rows.
"""
import os, sys, itertools, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq4_split import Split
from un_core import Arith, PRIMES, nat_rows, lift_rows, tables, xhat_of


class Var(Split):
    def __init__(self, base, Qh, L=1, maxout=None):
        self.__dict__.update(base.__dict__)
        self.KU = base.K; self.L = L; self.K = self.KU * L; self.maxout = maxout
        ar = self.ar = Arith(self.p, self.KU, L); KZ = self.KZ; R = self.R
        self._vt = {}
        self.ER, self.TR, self.Gtab = tables(ar, Qh, KZ, self.eta_inv_tab)
        x = xhat_of(ar, Qh, R)
        self.xpow = [ar.one()]
        for k in range(1, 5 * self.g): self.xpow.append(ar.mul(self.xpow[-1], x))
        self.EM = {key: [[(ar.from_u(c) if c is not None else None) for c in row] for row in tab] for key, tab in base.EM.items()}
        self.EO = {key: [[ar.sub(self.ER[key][i][j], self.EM[key][i][j]) for j in range(5)] for i in range(5)] for key in self.ER}
        # valuations per eps-order
        self.v0 = {key: min(ar.val(x_, 0) for row in self.EO[key] for x_ in row) for key in self.EO}
        self.vhi = {key: (min(ar.val(x_, e) for row in self.EO[key] for x_ in row for e in range(1, L)) if L > 1 else self.KU) for key in self.EO}

    # series arithmetic (overrides of Split)
    def mul(self, a, b): return self.ar.mul(a, b)
    def add(self, a, b): return self.ar.add(a, b)
    def sub(self, a, b): return self.ar.sub(a, b)
    def scale(self, a, c): return self.ar.scale(a, c)

    def bound(self, keys, S, E):
        """lower bound of the u-valuation of the product of the outer edges in S (any eps-order < L)"""
        s = 0; gains = []
        for e in range(E):
            if (S >> e) & 1:
                k = keys[e]; s += self.v0[k]; gains.append(max(0, self.v0[k] - self.vhi[k]))
        if self.L > 1 and gains:
            gains.sort(reverse=True); s -= sum(gains[:self.L - 1])
        return s

    def graph(self, t):
        G = self.inner_graphs[t]; V = G['V']; ed = G['ed']; hsAt = G['hsAt']; E = len(ed); N = self.K; p = self.p; KU = self.KU
        opts = []
        for v in range(V):
            tab = self.vertex_table(G['gv'][v], len(hsAt[v])); lst = []
            for hp, vec in tab.items():
                for perm in set(itertools.permutations(hp)): lst.append((perm, vec))
            if not lst: return {}
            opts.append(lst)
        subsets = []
        for S in range(1 << E):
            if self.maxout is not None and bin(S).count("1") > self.maxout: continue
            subsets.append((S, self.outer_pos[self.contract(G, S)]))
        steps = self.plan(G); acc = {}
        for combo in itertools.product(*opts):
            hp = [0] * (2 * E)
            for v in range(V):
                for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
            keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
            vecs = [combo[v][1] for v in range(V)]
            inner = [self.EM[k] for k in keys]; outer = [self.EO[k] for k in keys]
            for S, pos in subsets:
                if self.bound(keys, S, E) >= KU: continue
                z = self.colour_sum(steps, vecs, [outer[e] if (S >> e) & 1 else inner[e] for e in range(E)])
                if z is None: continue
                a = acc.get(pos)
                if a is None: acc[pos] = list(z)
                else:
                    for k in range(N): a[k] += z[k]
        ia = self.R.inv(G['aut'])
        return {pos: [x * ia % p for x in a] for pos, a in acc.items()}


_BASE = {}; _VAR = {}; _NAT = {}


def get_base(prime, g, ab, K):
    key = (prime, g, tuple(ab), K); b = _BASE.get(key)
    if b is None:
        b = Split(PRIMES[prime], g, (Fr(ab[0]), Fr(ab[1])), K, 'nat')
        nr = get_nat(prime, g, ab, K); b.eta_inv_tab = nr['eta_inv']
        _BASE[key] = b
    return b


def get_nat(prime, g, ab, K):
    key = (prime, g, tuple(ab), K); n = _NAT.get(key)
    if n is None: n = nat_rows(prime, g, ab, K); _NAT[key] = n
    return n


def build(spec, cache=True):
    prime, g, ab, K, L, Qh, maxout = spec
    base = get_base(prime, g, ab, K)
    if Qh is None:
        Qh = lift_rows(Arith(base.p, K, L), get_nat(prime, g, ab, K)['Qh'])
    return Var(base, Qh, L, maxout)


def _job(args):
    sid, spec, t = args
    v = _VAR.get(sid)
    if v is None: v = build(spec); _VAR.clear(); _VAR[sid] = v
    t0 = time.time(); res = v.graph(t)
    return t, res, time.time() - t0


def collect(var, results):
    tab = {}; p = var.p; N = var.K
    for res in results:
        for pos, c in res.items():
            a = tab.setdefault(pos, [0] * N)
            for k in range(N): a[k] = (a[k] + c[k]) % p
    return tab


def run(spec, var=None, log=None):
    """single process: {outer index: series} summed over all inner graphs"""
    v = build(spec) if var is None else var
    return v, collect(v, [v.graph(t) for t in range(len(v.inner_graphs))])


def run_pool(spec, nproc=4, log=None, sid="job"):
    """multiprocessing over the inner graphs (heaviest first); returns (var built in the parent, table)"""
    import multiprocessing as mp
    v = build(spec); n = len(v.inner_graphs)
    order = sorted(range(n), key=lambda t: -(2 ** min(v.inner_graphs[t]['E'], 6)) * 5 ** min(v.inner_graphs[t]['V'], 4))
    results = []; t0 = time.time(); cpu = 0.0
    with mp.Pool(processes=min(nproc, 5)) as pool:
        for k, (t, res, dt) in enumerate(pool.imap_unordered(_job, [(sid, spec, t) for t in order])):
            results.append(res); cpu += dt
            if log and (k + 1) % max(1, n // 10) == 0: log("      %d/%d inner graphs (%.0f s wall, %.0f s CPU)" % (k + 1, n, time.time() - t0, cpu))
    return v, collect(v, results)


def outer_info(var):
    """[(outer index, E, gv)]"""
    return [(i, G['E'], tuple(G['gv'])) for i, G in enumerate(var.inner_graphs)]


def polar_report(var, tab, e=0, npol=None):
    """{outer index: [u^0..u^(npol-1)] of the eps^e part} for the outer graphs with a non-zero polar coefficient"""
    KU = var.KU; npol = min(KU, 2 * var.g - 2) if npol is None else npol; out = {}
    for i, c in tab.items():
        part = c[e * KU:e * KU + npol]
        if any(part): out[i] = part
    return out
