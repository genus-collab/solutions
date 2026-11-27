# -*- coding: utf-8 -*-
"""rf_indep.py -- REFEREE (univ): graph sums of ARBITRARY hatted rows with the module-disjoint engine of the genus-4 referee
(natural_form/genus4_referee/indep: ind_tables, ind_graphs, ind_sum; imported read-only).  That engine reads its rows from
natural_form/referee/ref_hat.HatRows; here the class is replaced IN MEMORY (in every worker process) by RowsH, which serves rows
computed by the referee's own code (rf_rows).  Nothing on disk is modified.  It shares no module with fq4_split / un_engine.
Its graph lists: genus 4 in the project's outer order (index 0 = vertex graph), genus 2/3 in the Mathematica order.
run(rows_file, g, K, max_outer, nproc, deadline, parts_file) -> ({outer index: [K residues]}, graph list, complete?)
Rows file (JSON): p, zeta, K, NZ, Q[a][i][n] = list of K residues.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
IND = os.path.join(PROJ, "natural_form", "genus4_referee", "indep"); REFD = os.path.join(PROJ, "natural_form", "referee")
_W = {}


class RowsH(object):
    def __init__(self, d):
        p = self.p = d['p']; self.zeta = d['zeta']; self.K = d['K']; self.NZ = d['NZ']; self.Q = d['Q']; K = self.K
        q0 = self.Q[0][0][0]; assert q0[0] % p
        x = [0] * K; x[0] = pow(q0[0], p - 2, p)                       # xhat = 1/Qhat_0^(0)[zhat^0]
        for n in range(1, K): x[n] = (-sum(q0[j] * x[n - j] for j in range(1, n + 1)) * x[0]) % p
        self.C0 = x
        i5 = pow(5, p - 2, p)
        self.Rt = {(n, 0): [[sum(self.Q[(3 - a) % 5][j][0][0] * self.Q[a][i][n][0] for a in range(5)) * i5 % p for i in range(5)]
                            for j in range(5)] for n in range(self.NZ + 1)}
        self.checks = {"rows supplied by the referee (rf_rows)": True}
        self.tc = [0, 1] + [0] * K


def setup(rows_file, g):
    key = (rows_file, g)
    if _W.get('key') == key: return _W['S']
    for d in (IND, REFD):
        if d not in sys.path: sys.path.insert(0, d)
    import ind_tables
    from ind_graphs import graphs
    from ind_sum import Summer
    d = json.load(open(rows_file))
    assert d['NZ'] == 3 * g - 3
    ind_tables.HatRows = lambda p, NZ, K, alpha=0, beta=0, zeta=None: RowsH(d)
    T = ind_tables.Tables(d['prime'], g, d['K'], "nat")
    assert min(T.valO.values()) >= 1, "E_out not O(u): the engine's subset bound would be wrong"
    L = graphs(g); S = Summer(T, L)
    _W.update(key=key, S=S)
    return S


def _job(args):
    rows_file, g, t, mo = args; t0 = time.time(); S = setup(rows_file, g)
    return t, {str(k): v for k, v in S.graph(t, mo).items()}, time.time() - t0


def run(rows_file, g, max_outer=None, nproc=2, deadline=1500.0, parts_file=None, log=print):
    t00 = time.time(); S = setup(rows_file, g); L = S.graphs; n = len(L); K = S.T.K; p = S.T.p
    res = {}
    if parts_file and os.path.exists(parts_file):
        dd = json.load(open(parts_file))
        if dd.get('rows') == rows_file and dd.get('max_outer') == max_outer: res = dd['res']
    todo = sorted([t for t in range(n) if str(t) not in res], key=lambda t: -(2 ** min(L[t]['E'], 7)) * 5 ** min(L[t]['V'], 4))
    def save():
        if parts_file:
            json.dump(dict(rows=rows_file, max_outer=max_outer, res=res), open(parts_file + ".new", "w")); os.replace(parts_file + ".new", parts_file)
    if todo and nproc <= 1:
        for t in todo:
            res[str(t)] = {str(k): v for k, v in S.graph(t, max_outer).items()}
            if time.time() - t00 > deadline: save(); break
    elif todo:
        import multiprocessing as mp
        pool = mp.Pool(processes=min(nproc, 2)); it = pool.imap_unordered(_job, [(rows_file, g, t, max_outer) for t in todo]); done = 0; last = time.time()
        try:
            while done < len(todo):
                left = deadline - (time.time() - t00)
                if left <= 0: raise mp.TimeoutError()
                t, r, dt = it.next(timeout=left); res[str(t)] = r; done += 1
                if time.time() - last > 60: save(); last = time.time(); log("      %d/%d inner graphs (%.0f s)" % (done, len(todo), time.time() - t00))
            pool.close(); pool.join()
        except mp.TimeoutError:
            pool.terminate(); pool.join()
    save()
    complete = len(res) == n
    tab = {}
    for t, r in res.items():
        for pos, c in r.items():
            acc = tab.setdefault(int(pos), [0] * K)
            for k in range(K): acc[k] = (acc[k] + c[k]) % p
    return tab, L, complete


def series_inv(a, p, K):
    x = [0] * K; x[0] = pow(a[0], p - 2, p)
    for n in range(1, K): x[n] = (-sum(a[j] * x[n - j] for j in range(1, n + 1)) * x[0]) % p
    return x


def smul(a, b, p, K):
    out = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i): out[i + j] += a[i] * b[j]
    return [x % p for x in out]


def polar_in_T(G, T1, g, p, K):
    """G = [u^0..u^(K-1)] of u^(2g-2) x (genus-g term), T = u T1(u) (T1(0) = 1): coefficients c_j of T^(j+2-2g), j < K,
    i.e. G(u) u^(2-2g) = Sum_j c_j T^(j+2-2g) + O(u^(K+2-2g)).  Own reversion: u = T s(T)."""
    # reversion of T = u T1(u): s with s = 1/T1(T s), fixed point iteration
    s = [1] + [0] * (K - 1); iT1 = series_inv(T1, p, K)
    for _ in range(K + 1):
        us = [0] + s[:K - 1]                                            # u = T s(T) as a series in T
        pw = [1] + [0] * (K - 1); comp = [0] * K
        for c in iT1:
            comp = [(x + c * y) % p for x, y in zip(comp, pw)]; pw = smul(pw, us, p, K)
        s = comp
    us = [0] + s[:K - 1]; pw = [1] + [0] * (K - 1); Gc = [0] * K
    for c in G:
        Gc = [(x + c * y) % p for x, y in zip(Gc, pw)]; pw = smul(pw, us, p, K)
    e = 2 * g - 2; spow = [1] + [0] * (K - 1); isv = series_inv(s, p, K)
    for _ in range(e): spow = smul(spow, isv, p, K)
    return smul(Gc, spow, p, K)                                          # G(u) u^-e = G(Ts) T^-e s^-e
