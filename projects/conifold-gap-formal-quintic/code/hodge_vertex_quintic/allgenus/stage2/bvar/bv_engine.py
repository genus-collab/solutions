# -*- coding: utf-8 -*-
"""bv_engine.py -- the outer graph sum for rows of a system with VARIED B: univ/un_engine.Var (imported read-only) on a copy of the
fq4_split.Split base whose E_M = tau^(a+b+1) Eh_0 uses the system's own tau (bv_sys.base_with_tau).  un_engine.run_pool cannot be used
as is (its workers rebuild the base with the quintic tau), hence the small pool below (same job order and caching as un_engine).
spec = (prime, g, K, tau list (K residues) or None for the quintic tau, rows (lifted, L = 1), maxout).
API: make_var(spec);  run_single(spec);  run_pool(spec, nproc, log, sid) -> (var, {outer index: [K residues]}).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import bv_sys  # noqa: F401  (sets sys.path for the project and univ modules)
from un_engine import Var, get_base, collect

_VAR = {}


def make_var(spec):
    prime, g, K, tl, QL, maxout = spec
    base = get_base(prime, g, ("0", "0"), K)
    b2 = bv_sys.base_with_tau(base, tl) if tl is not None else base
    return Var(b2, QL, 1, maxout)


def _job(args):
    sid, spec, t = args
    v = _VAR.get(sid)
    if v is None:
        v = make_var(spec); _VAR.clear(); _VAR[sid] = v
    t0 = time.time(); res = v.graph(t)
    return t, res, time.time() - t0


def run_single(spec, var=None):
    v = make_var(spec) if var is None else var
    return v, collect(v, [v.graph(t) for t in range(len(v.inner_graphs))])


def run_pool(spec, nproc=4, log=None, sid="bv"):
    import multiprocessing as mp
    v = make_var(spec); n = len(v.inner_graphs)
    order = sorted(range(n), key=lambda t: -(2 ** min(v.inner_graphs[t]['E'], 6)) * 5 ** min(v.inner_graphs[t]['V'], 4))
    results = []; t0 = time.time(); cpu = 0.0
    with mp.Pool(processes=min(nproc, 4)) as pool:
        for k, (t, res, dt) in enumerate(pool.imap_unordered(_job, [(sid, spec, t) for t in order])):
            results.append(res); cpu += dt
            if log and (k + 1) % max(1, n // 4) == 0:
                log("      %d/%d inner graphs (%.0f s wall, %.0f s CPU)" % (k + 1, n, time.time() - t0, cpu))
    return v, collect(v, results), cpu
