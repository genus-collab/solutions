# -*- coding: utf-8 -*-
"""rf_run.py -- referee's independent graph-sum runs for B-perturbed systems (rows from rf_bsys, E_M = tau^(a+b+1) Eh_0 with Eh_0 taken
from a fq4_split.Split base of variant 'M0' (tau = 1), NOT from the u^0 part of the 'nat' base as the worker does).
Stage-1 / project modules used read-only: un_engine.Var + collect (the outer graph sum), un_core (Arith, lift_rows, symplectic_check,
left_dress), un_engine.get_nat (frozen model rows Rh, eta^-1, project rows), fq4_split.Split, fq_fp.  No bv_* module is imported.
Usage: py -I -B rf_run.py SYS g prime [--norm canon|stored] [--ctrl none|staleM|staleHat|dress0] [--maxout m] [--kext k] [--nproc n]
       [--out DIR]   -> DIR/rf_<tag>.log, DIR/rf_<tag>.json
"""
import os, sys, json, time, argparse, copy
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); AG = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
PROJ = os.path.normpath(os.path.join(AG, ".."))
for d in (PROJ, os.path.join(AG, "univ"), HERE):
    if d not in sys.path: sys.path.insert(0, d)
import rf_bsys as RB
from rf_ser import L, fr_mod, pmul, pinv, ppow
from fq_fp import P1, P2
PR = {"89": P2, "61": P1}
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}
SCR = os.path.join(HERE, "scratch")
WORKER_TAB = os.path.normpath(os.path.join(HERE, "..", "tables"))
_GA = {}


def system_gauge(name, norm):
    key = (name, norm)
    if key not in _GA:
        S = RB.build(name); p, info = RB.gauge(S, 12); pn, C, sig, c0, m5 = RB.normalise(p, norm)
        _GA[key] = (S, pn, dict(C=C, sig=sig, c0=c0, m5=m5, info=info))
    return _GA[key]


def resid(s, K, p):
    assert s.zero() or s.v >= 0, "pole in a hatted entry"
    assert s.P >= K, "precision O(u^%d) < %d" % (s.P, K)
    return [fr_mod(s.co(n), p) for n in range(K)]


def make_rows(S, pn, prime, g, K, stale=None):
    """Qh[a][i][n] (plain residue lists) = Sum_k Ghat_(k,b) tau^(n-k) Rh[b][i][n-k], a = b + k"""
    from un_engine import get_nat
    p = PR[prime]; nr = get_nat(prime, g, ("0", "0"), K); Rh = nr['Rh']; KZ = 3 * g - 3
    hat = RB.hatted(S, pn, stale_tau=stale)
    taur = resid(S['tau'], K, p)
    tp = [[1] + [0] * (K - 1)]
    for _ in range(KZ): tp.append(pmul(tp[-1], taur, K, p))
    Qh = [[[[0] * K for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (k, b), s in hat.items():
        if k > KZ: continue
        h = resid(s, K, p); a = (b + k) % 5
        if not any(h): continue
        for n in range(k, KZ + 1):
            ht = pmul(h, tp[n - k], K, p)
            for i in range(5):
                c = Rh[b][i][n - k] % p
                if c:
                    q = Qh[a][i][n]
                    for j in range(K): q[j] = (q[j] + ht[j] * c) % p
    return Qh, nr


_BASE = {}


def em_base(prime, g, K, taur):
    """Split('M0') base with E_M = tau^(a+b+1) Eh_0 (Eh_0 = the constant M0 table)"""
    from fq4_split import Split
    from un_engine import get_nat
    key = (prime, g, K)
    if key not in _BASE:
        b = Split(PR[prime], g, (Fr(0), Fr(0)), K, 'M0'); b.eta_inv_tab = get_nat(prime, g, ("0", "0"), K)['eta_inv']; _BASE[key] = b
    b0 = _BASE[key]; p = PR[prime]; b2 = copy.copy(b0); KZ = 3 * g - 3
    tp = [[1] + [0] * (K - 1)]
    for _ in range(KZ + 1): tp.append(pmul(tp[-1], taur, K, p))
    EM = {}
    for (a, bb), tab in b0.EM.items():
        out = []
        for row in tab:
            r2 = []
            for x in row:
                if x is None: r2.append(None); continue
                assert not any(x[1:]), "M0 base E_M not constant"
                r2.append([x[0] * t % p for t in tp[a + bb + 1]])
            out.append(r2)
        EM[(a, bb)] = out
    b2.EM = EM; b2._vt = {}
    return b2


_V = {}


def make_var(spec):
    from un_engine import Var
    from un_core import Arith, lift_rows
    prime, g, K, taur, Qh, maxout = spec
    base = em_base(prime, g, K, taur); QL = lift_rows(Arith(PR[prime], K, 1), Qh)
    return Var(base, QL, 1, maxout)


def _job(args):
    sid, spec, t = args
    if sid not in _V: _V.clear(); _V[sid] = make_var(spec)
    return _V[sid].graph(t)


def graph_sum(spec, nproc, sid):
    from un_engine import collect
    v = make_var(spec); n = len(v.inner_graphs)
    if nproc <= 1: return v, collect(v, [v.graph(t) for t in range(n)])
    import multiprocessing as mp
    order = sorted(range(n), key=lambda t: -(2 ** min(v.inner_graphs[t]['E'], 6)) * 5 ** min(v.inner_graphs[t]['V'], 4))
    with mp.Pool(processes=min(nproc, 3)) as pool:
        res = list(pool.imap_unordered(_job, [(sid, spec, t) for t in order]))
    return v, collect(v, res)


def in_T(G, Tq, g, K, p):
    """u^(2-2g) G(u) = Sum_{j=2-2g}^{-1} d_j T^j + O(1); returns [d_{2-2g}, ..., d_{-1}] (residues)"""
    X = [fr_mod(Tq.co(n + 1), p) for n in range(K)]                     # T = u X
    rem = list(G[:K]); out = []
    for j in range(2 - 2 * g, 0):
        e = j + 2 * g - 2; Xj = ppow(X, j, K, p); d = rem[e] * pow(Xj[0], p - 2, p) % p; out.append(d)
        for t in range(e, K): rem[t] = (rem[t] - d * Xj[t - e]) % p
    return out


def rat(x, p):
    """small rational reconstruction for printing"""
    from fq_fp import FpRing
    from fq3_rec import Rec
    return Rec(FpRing(p)).fmt(x)


def run(a):
    name, g, prime = a.sys, a.g, a.prime; p = PR[prime]; K = 2 * g - 2 + a.kext; npol = 2 * g - 2
    tag = "%s_g%d_%s_%s_%s%s_k%d" % (name, g, prime, a.norm, a.ctrl, "" if a.maxout is None else "_out%d" % a.maxout, K)
    os.makedirs(a.out, exist_ok=True); lf = open(os.path.join(a.out, "rf_%s.log" % tag), "w", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); lf.write(s + "\n"); lf.flush()
    t0 = time.time(); S, pn, inf = system_gauge(name, a.norm)
    log("rf_run %s: C = %s, sigma = %s, c0 = %s, m5 = %s" % (tag, inf['C'], inf['sig'], inf['c0'], inf['m5']))
    stale = None
    if a.ctrl == "staleHat": stale = RB.build("R0")['w']
    Qh, nr = make_rows(S, pn, prime, g, K, stale=stale)
    taur = [fr_mod(S['tau'].co(n), p) for n in range(K)]
    if a.ctrl == "staleM": taur = [fr_mod(RB.build("R0")['tau'].co(n), p) for n in range(K)]
    from un_core import Arith, lift_rows, symplectic_check, left_dress
    ar = Arith(p, K, 1)
    if a.ctrl in ("dress0", "dress2"):
        QL = left_dress(ar, lift_rows(ar, Qh), {(2, 1, 1): ar.mono(fr_mod(Fr(1, 7), p), int(a.ctrl[-1]))}, 3 * g - 3)
        Qh = [[[(x if x is not None else [0] * K) for x in col] for col in row] for row in QL]
    sb = symplectic_check(ar, lift_rows(ar, Qh), 3 * g - 3, nr['eta_inv'])
    log("   rows symplectic: %s;  R0 regression (rows = project rows): %s" % (not sb, all(Qh[x][i][n] == nr['Qh'][x][i][n] for x in range(5) for i in range(5) for n in range(3 * g - 2)) if name == "R0" else "n/a"))
    v, tab = graph_sum((prime, g, K, taur, Qh, a.maxout), a.nproc, tag)
    outer = [(G['E'], list(G['gv'])) for G in v.inner_graphs]
    rng = [i for i in range(1, len(outer)) if a.maxout is None or outer[i][0] <= a.maxout]
    bad = {i: tab.get(i, [0] * K)[:npol] for i in rng if any(tab.get(i, [0] * K)[:npol])}
    log("   edge graphs in range: %d; singular: %s" % (len(rng), sorted(bad) if bad else "none"))
    for i in sorted(bad)[:20]: log("      graph %d %s: polar %s" % (i, outer[i], [rat(x, p) for x in bad[i]]))
    c0 = inf['c0']; Tq = S['T'].sc(c0); V = in_T(tab.get(0, [0] * K), Tq, g, K, p)
    okv = V[0] == fr_mod(CG[g], p) and not any(V[1:])
    log("   vertex in T'' = c0 T: %s -> %s" % ([rat(x, p) for x in V], "c_g T''^(2-2g) + 0" if okv else "FAILS"))
    if a.maxout is None:
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(len(outer))) % p for k in range(K)]
        log("   total in T'': %s" % [rat(x, p) for x in in_T(tot, Tq, g, K, p)])
    # comparison with the worker's table (same system, g, prime, norm, ctrl, maxout, K) if present
    wt = "tab_%s_g%d_%s_%s_%s%s%s.json" % (name, g, prime, a.norm, a.ctrl, "" if a.maxout is None else "_out%d" % a.maxout,
                                           "" if g == 4 else "_k%d" % K)
    wp = os.path.join(WORKER_TAB, wt)
    if os.path.exists(wp):
        W = {int(i): c for i, c in json.load(open(wp))['table'].items()}
        dif = [i for i in range(len(outer)) if (a.maxout is None or outer[i][0] <= a.maxout or i == 0)
               and [x % p for x in tab.get(i, [0] * K)] != [x % p for x in W.get(i, [0] * K)[:K]]]
        log("   worker table %s: %s" % (wt, "IDENTICAL on all %d graphs in range (all K = %d coefficients)" % (len(rng) + 1, K) if not dif else "DIFFERS on %s" % dif[:10]))
    else: log("   worker table %s: not present" % wt)
    json.dump(dict(tag=tag, table={str(i): c for i, c in tab.items()}, singular=sorted(bad), vertex_ok=okv, outer=outer),
              open(os.path.join(a.out, "rf_%s.json" % tag), "w"))
    log("RF_RUN %s: singular %s; vertex law %s; %.0f s" % (tag, sorted(bad) if bad else "none", "holds" if okv else "FAILS", time.time() - t0))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("sys"); ap.add_argument("g", type=int); ap.add_argument("prime")
    ap.add_argument("--norm", default="canon"); ap.add_argument("--ctrl", default="none"); ap.add_argument("--maxout", type=int, default=None)
    ap.add_argument("--kext", type=int, default=2); ap.add_argument("--nproc", type=int, default=1); ap.add_argument("--out", default=SCR)
    run(ap.parse_args())
