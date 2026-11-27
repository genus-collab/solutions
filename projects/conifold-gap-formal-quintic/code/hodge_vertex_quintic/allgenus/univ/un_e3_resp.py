# -*- coding: utf-8 -*-
"""un_e3_resp.py -- E3, step 1: the COMPLETE first-order response map of the outer graph sum to a left dressing of the natural R'.

For every Z/5-graded, infinitesimally symplectic constant hatted monomial
        h(k, b) = zhat^k [ E_{b+k,b} - (-1)^k E_{b*+k,b*} ],   b* = 3 - b - k mod 5      (one term, k odd, if b* = b)
(matrices acting on the flat row index of the hatted rows; symplectic: p(z) eta^-1 + eta^-1 p(-z)^T = 0  <=>  p_{k,b} = -(-1)^k p_{k,b*})
the engine un_engine.Var is run with dual numbers (L = 2) on rows -> (1 + eps h) rows, M = M_0(tau zhat) fixed, and the eps^1 part of
every outer graph is stored as a u-series mod u^(2g-2) (genus 2 and 3).  A flat-regular left dressing
        p(u, z) = Sum_{k,b} p_{k,b}(u) z^k E_{b+k,b}     is, in hatted variables,     Sum p_{k,b}(u) u^mu(b,k) zhat^k E_{b+k,b},
so its first-order response is Sum_{orbits} p_{k,b}(u) u^mu(b,k) rho_{k,b}(Gamma; u): the map is determined by these runs.
Output: un_resp_<prime>.json = {"g": {"k,b": {outer index: [coefficients]}}, "outer": {...}, "orbits": [[k, b, b*]]}.
Run: py -I -B un_e3_resp.py 89|61 [nproc<=4]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)


def orbits(kmax):
    """[(k, b, b*)] orbit representatives (b <= b* in the cyclic sense of first occurrence); self-dual orbits only for k odd"""
    out = []
    for k in range(kmax + 1):
        seen = set()
        for b in range(5):
            if b in seen: continue
            bs = (3 - b - k) % 5; seen.add(b); seen.add(bs)
            if bs == b and k % 2 == 0: continue
            out.append((k, b, bs))
    return out


def monomial(ar, k, b, bs):
    p = ar.p
    X = {((b + k) % 5, b, k): ar.mono(1, 0, 1)}
    if bs != b: X[((bs + k) % 5, bs, k)] = ar.mono(p - 1 if k % 2 == 0 else 1, 0, 1)
    return X


def job(args):
    prime, g, K, k, b, bs = args
    from un_core import Arith, lift_rows, left_dress
    from un_engine import Var, get_base, get_nat, run
    base = get_base(prime, g, ("0", "0"), K); nr = get_nat(prime, g, ("0", "0"), K)
    ar = Arith(nr['p'], K, 2); Q0 = lift_rows(ar, nr['Qh'])
    t0 = time.time(); v = Var(base, left_dress(ar, Q0, monomial(ar, k, b, bs), base.KZ), 2); _, tab = run(None, v)
    res = {str(i): ar.part(c, 1) for i, c in tab.items() if any(ar.part(c, 1))}
    nat0 = {str(i): ar.part(c, 0) for i, c in tab.items() if any(ar.part(c, 0))}
    return (g, k, b), res, nat0, time.time() - t0


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"; nproc = min(4, int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    LOG = open(os.path.join(HERE, "un_e3_resp_%s.log" % prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    import multiprocessing as mp
    from un_engine import get_base
    log("=" * 100); log("un_e3_resp: p = 2^%s-1" % prime)
    jobs = []; out = {}; outer = {}
    for g, K in ((2, 2), (3, 4)):
        base = get_base(prime, g, ("0", "0"), K)
        outer[str(g)] = [dict(E=G['E'], gv=list(G['gv']), ed=[list(e) for e in G['ed']]) for G in base.inner_graphs]
        out[str(g)] = {}
        for k, b, bs in orbits(base.KZ): jobs.append((prime, g, K, k, b, bs))
    t0 = time.time(); nat = {}
    pool = mp.Pool(processes=nproc) if nproc > 1 else None          # nproc 1: in-process (no worker process)
    try:
        for (g, k, b), res, nat0, dt in (pool.imap_unordered(job, jobs) if pool else map(job, jobs)):
            out[str(g)]["%d,%d" % (k, b)] = res
            if str(g) in nat: assert nat[str(g)] == nat0, "eps^0 parts differ between runs"
            nat[str(g)] = nat0
            log("   genus %d, h(%d, %d): %d outer graphs respond (%.0f s)" % (g, k, b, len(res), dt))
    finally:
        if pool: pool.close(); pool.join()
    for g in ("2", "3"):
        ok = set(nat[g]) == {"0"}
        log(("PASS " if ok else "FAIL ") + "genus %s: the eps^0 part of every run is the natural run (only the vertex graph is non-zero in the window)" % g)
    json.dump(dict(resp=out, outer=outer, nat=nat, orbits={g: [list(o) for o in orbits(3 * int(g) - 3)] for g in ("2", "3")}),
              open(os.path.join(HERE, "un_resp_%s.json" % prime), "w"))
    log("UN_E3_RESP done: %d runs, %.0f s wall" % (len(jobs), time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
