# -*- coding: utf-8 -*-
"""rf_indep2.py -- the referee's B-perturbed rows (rf_bsys/rf_run.make_rows) summed on the MODULE-DISJOINT engine
natural_form/genus4_referee/indep (ind_tables, ind_graphs, ind_sum; shares no module with fq4_split / un_engine), imported read-only.
Two in-memory patches (nothing on disk is changed): ind_tables.HatRows serves our rows (as in univ/referee/rf_indep.py), and
ind_tables.Fr returns, for the call Fr(4, 5k+4) of the 'nat' branch, the k-th coefficient of 1/(tau/tau(0)) of OUR system, so the
engine builds E_M = (tau/tau0)^(a+b+1) E_M0 itself (E_M0 is read off our rows at u = 0, which already contain M_0(tau0 zhat)) and
runs its own asserts ([u^0] E_R = E_M0, E_M = division of M_0^-1(tau z), Z/5 charges, symplectic division, E_out = O(u)).
Compares with the referee's un_engine tables rf_<tag>.json: genus 4 graph by graph (same outer order); genus 2/3 (different graph
order) as multisets of per-graph series, plus regularity counts, the vertex law and the total.
Usage: py -I -B rf_indep2.py SYS g prime [--maxout m] [--nproc n] [--out DIR]
"""
import os, sys, json, time, argparse
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_run
from rf_run import PR, SCR, CG
from rf_ser import fr_mod
PROJ = rf_run.PROJ
IND = os.path.join(PROJ, "natural_form", "genus4_referee", "indep"); REFD = os.path.join(PROJ, "natural_form", "referee")
_W = {}


class RowsH(object):
    def __init__(self, d):
        p = self.p = d['p']; self.zeta = d['zeta']; self.K = d['K']; self.NZ = d['NZ']; self.Q = d['Q']; K = self.K
        q0 = self.Q[0][0][0]; assert q0[0] % p
        x = [0] * K; x[0] = pow(q0[0], p - 2, p)
        for n in range(1, K): x[n] = (-sum(q0[j] * x[n - j] for j in range(1, n + 1)) * x[0]) % p
        self.C0 = x; i5 = pow(5, p - 2, p)
        self.Rt = {(n, 0): [[sum(self.Q[(3 - a) % 5][j][0][0] * self.Q[a][i][n][0] for a in range(5)) * i5 % p for i in range(5)]
                            for j in range(5)] for n in range(self.NZ + 1)}
        self.checks = {"rows supplied by the referee (rf_run.make_rows)": True}; self.tc = [0, 1] + [0] * K


def setup(fn, g):
    if _W.get('fn') == fn: return _W['S']
    for d in (IND, REFD):
        if d not in sys.path: sys.path.insert(0, d)
    import ind_tables
    from ind_graphs import graphs
    from ind_sum import Summer
    d = json.load(open(fn)); itau = [Fraction(x) for x in d['invtau']]
    def Fr(a=0, b=1):
        if a == 4 and b >= 4 and (b - 4) % 5 == 0: return itau[(b - 4) // 5]
        return Fraction(a, b)
    ind_tables.Fr = Fr
    ind_tables.HatRows = lambda p, NZ, K, alpha=0, beta=0, zeta=None: RowsH(d)
    T = ind_tables.Tables(d['prime'], g, d['K'], "nat")
    assert min(T.valO.values()) >= 1, "E_out not O(u)"
    S = Summer(T, graphs(g)); _W.update(fn=fn, S=S)
    return S


def _job(a):
    fn, g, t, mo = a; S = setup(fn, g)
    return {str(k): v for k, v in S.graph(t, mo).items()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("sys"); ap.add_argument("g", type=int); ap.add_argument("prime")
    ap.add_argument("--maxout", type=int, default=None); ap.add_argument("--nproc", type=int, default=1); ap.add_argument("--out", default=SCR)
    a = ap.parse_args(); g = a.g; p = PR[a.prime]; K = 2 * g - 2 + (0 if g == 4 else 2); t0 = time.time()
    S_, pn, inf = rf_run.system_gauge(a.sys, "canon"); Qh, nr = rf_run.make_rows(S_, pn, a.prime, g, K)
    from fq_fp import FpRing
    zeta = FpRing(p).z(1) % p; tau = S_['tau']; t0c = tau.co(0)
    it = tau.sc(1 / t0c); inv = [Fraction(1)]
    tl = [it.co(k) for k in range(K)]
    for k in range(1, K): inv.append(-sum(tl[j] * inv[k - j] for j in range(1, k + 1)))
    tag = "%s_g%d_%s%s" % (a.sys, g, a.prime, "" if a.maxout is None else "_out%d" % a.maxout)
    fn = os.path.join(a.out, "rows_%s.json" % tag)
    json.dump(dict(prime=a.prime, p=p, zeta=zeta, K=K, NZ=3 * g - 3, Q=Qh, invtau=[str(x) for x in inv]), open(fn, "w"))
    S = setup(fn, g); L = S.graphs; n = len(L)
    if a.nproc <= 1: res = [_job((fn, g, t, a.maxout)) for t in range(n)]
    else:
        import multiprocessing as mp
        order = sorted(range(n), key=lambda t: -(2 ** min(L[t]['E'], 7)) * 5 ** min(L[t]['V'], 4))
        with mp.Pool(processes=min(a.nproc, 3)) as pool: res = list(pool.imap_unordered(_job, [(fn, g, t, a.maxout) for t in order]))
    tab = {}
    for r in res:
        for pos, c in r.items():
            acc = tab.setdefault(int(pos), [0] * K)
            for k in range(K): acc[k] = (acc[k] + c[k]) % p
    npol = 2 * g - 2; rng = [i for i in range(n) if L[i]['E'] > 0 and (a.maxout is None or L[i]['E'] <= a.maxout)]
    sing = [i for i in rng if any(tab.get(i, [0] * K)[:npol])]
    iv = [i for i in range(n) if L[i]['E'] == 0]; assert len(iv) == 1; iv = iv[0]
    V = rf_run.in_T(tab.get(iv, [0] * K), S_['T'].sc(inf['c0']), g, K, p)
    okv = V[0] == fr_mod(CG[g], p) and not any(V[1:])
    msg = "indep engine %s: %d edge graphs in range, singular %s; vertex in T'' %s (%s)" % (tag, len(rng), sing or "none", [rf_run.rat(x, p) for x in V], "law holds" if okv else "FAILS")
    mine = os.path.join(a.out, "rf_%s_g%d_%s_canon_none%s_k%d.json" % (a.sys, g, a.prime, "" if a.maxout is None else "_out%d" % a.maxout, K))
    if os.path.exists(mine):
        J = json.load(open(mine)); TJ = {int(i): c for i, c in J['table'].items()}; outer = J['outer']
        sel_ind = [i for i in range(n) if a.maxout is None or L[i]['E'] <= a.maxout]
        sel_un = [i for i in range(len(outer)) if a.maxout is None or outer[i][0] <= a.maxout]
        A = sorted(tuple(tab.get(i, [0] * K)) for i in sel_ind); B = sorted(tuple(x % p for x in TJ.get(i, [0] * K)) for i in sel_un)
        same_ms = (A == B)
        same_idx = (g == 4 and all(tab.get(i, [0] * K) == [x % p for x in TJ.get(i, [0] * K)] for i in sel_ind))
        msg += "; vs un_engine run: multiset of per-graph series %s%s" % ("IDENTICAL" if same_ms else "DIFFERENT",
                                                                       ("; graph by graph %s" % ("IDENTICAL" if same_idx else "DIFFERENT")) if g == 4 else "")
    if g == 4:     # the worker's genus-4 table (same outer order), if present
        wt = os.path.join(rf_run.WORKER_TAB, "tab_%s_g4_%s_canon_none%s.json" % (a.sys, a.prime, "" if a.maxout is None else "_out%d" % a.maxout))
        if os.path.exists(wt):
            W = {int(i): c for i, c in json.load(open(wt))['table'].items()}
            sel = [i for i in range(n) if a.maxout is None or L[i]['E'] <= a.maxout]
            dif = [i for i in sel if tab.get(i, [0] * K) != [x % p for x in W.get(i, [0] * K)[:K]]]
            msg += "; vs WORKER table %s: %s" % (os.path.basename(wt), "IDENTICAL on all %d graphs" % len(sel) if not dif else "DIFFERS on %s" % dif[:10])
    print(msg + "; %.0f s" % (time.time() - t0), flush=True)
    with open(os.path.join(a.out, "rf_indep2.log"), "a", encoding="utf-8") as f: f.write(msg + "\n")


if __name__ == "__main__":
    main()
