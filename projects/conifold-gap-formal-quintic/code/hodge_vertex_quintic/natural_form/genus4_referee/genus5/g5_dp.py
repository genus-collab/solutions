# -*- coding: utf-8 -*-
"""g5_dp.py -- PROTOTYPE of the main change proposed for genus 5: one vertex elimination per decoration for ALL outer-edge subsets.

Split.graph repeats the whole colour sum for each of the (up to 3302 at genus 5) outer-edge subsets S.  Here the subset is carried
through the elimination as part of the state: after eliminating some vertices the state is (colours of the frontier vertices, the
inner/outer marks of the edges closed so far) -> series; an edge is closed when its second end is eliminated and then branches into
"inner" (factor E_M) and "outer" (factor E_out, only while |S| < K - 1 before the branch, i.e. |S| <= K - 1 after it: the same subsets
as Split.graph).  The work is  Sum_steps 2^(closed edges) 5^(frontier + 1)  instead of  2^E Sum_steps 5^(frontier + 1).
Same tables, same vertex-elimination plan (Split.plan), same output format as Split.graph.
Modes:
  check4            genus 4, 2^89-1, natural M: graph_dp on ALL 379 inner graphs = stored tables parts/g4_89_nat_run.json (read only)
  check5 --list ..  genus 5: graph_dp == Split.graph on the listed (cheap) inner graphs
  sample [--slice i/n]  genus 5 cost sampling of graph_dp (random decorations, one full subset-DP each) -> g5_dpsample_<tag>.jsonl
  full --list ..    graph_dp on listed graphs, measured seconds and peak memory -> g5_dpfull_<tag>.jsonl
This file holds the functions; the command-line driver is g5_dp_run.py:   py -I g5_dp_run.py check4 --nproc 5
"""
import os, sys, json, time, random, argparse, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

_SP = {}
POPC = [bin(i).count("1") for i in range(1 << 12)]


def get(spec):
    sp = _SP.get(spec)
    if sp is None:
        from g5_lib import make_split
        from fq_fp import P1, P2
        prime, g, K, variant = spec
        sp = make_split(P2 if prime == "89" else P1, g, (0, 0), K, variant); _SP[spec] = sp
    return sp


def peak_mb():
    try:
        import ctypes
        from ctypes import wintypes
        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL("kernel32"); k32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi = ctypes.WinDLL("psapi"); psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return pmc.PeakWorkingSetSize / 2.0 ** 20
    except Exception:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        except Exception:
            return float("nan")


def combo_dp(sp, steps, vecs, inner, outer):
    """{S (bitmask of outer edges, |S| < K): series} for one decoration; inner/outer = per-edge 5x5 tables E_M / E_out"""
    mul = sp.mul; add = sp.add; cap = sp.K - 1
    cur = {(): {0: [1] + [0] * (sp.K - 1)}}
    for v, conn, keep, loops in steps:
        new = {}; vec = vecs[v]
        for key, sub in cur.items():
            for c in range(5):
                vc = vec[c]
                if vc is None: continue
                facs = [(1 << e, inner[e][c][c], outer[e][c][c]) for e in loops]
                for e, pos, flip in conn:
                    k = key[pos]
                    facs.append((1 << e, inner[e][c][k], outer[e][c][k]) if flip else (1 << e, inner[e][k][c], outer[e][k][c]))
                items = []
                for S, val in sub.items():
                    x = mul(val, vc)
                    if x is not None: items.append((S, x))
                for bit, fi, fo in facs:
                    nxt = []
                    for S, x in items:
                        a = mul(x, fi)
                        if a is not None: nxt.append((S, a))
                        if POPC[S] < cap:
                            b = mul(x, fo)
                            if b is not None: nxt.append((S | bit, b))
                    items = nxt
                if not items: continue
                full = key + (c,); nk = tuple(full[k] for k in keep); dst = new.get(nk)
                if dst is None: dst = {}; new[nk] = dst
                for S, x in items:
                    y = dst.get(S); dst[S] = x if y is None else add(y, x)
        cur = new
        if not cur: return {}
    return cur.get((), {})


def prepare(sp, t):
    G = sp.inner_graphs[t]; V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; opts = []
    for v in range(V):
        tab = sp.vertex_table(G['gv'][v], len(hsAt[v])); lst = []
        for hp, vec in tab.items():
            for perm in set(itertools.permutations(hp)): lst.append((perm, vec))
        opts.append(lst)
    assert sp.variant != 'id'
    pos = {S: sp.outer_pos[sp.contract(G, S)] for S in range(1 << E) if POPC[S] < sp.K}
    return G, opts, pos, sp.plan(G)


def one(sp, G, combo, steps):
    V = G['V']; E = len(G['ed']); hsAt = G['hsAt']; hp = [0] * (2 * E)
    for v in range(V):
        for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
    keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
    return combo_dp(sp, steps, [combo[v][1] for v in range(V)], [sp.EM[k] for k in keys], [sp.EO[k] for k in keys])


def graph_dp(sp, t):
    """same output as Split.graph(t) (variant != 'id')"""
    G, opts, pos, steps = prepare(sp, t); K = sp.K; p = sp.p; acc = {}
    if any(not lst for lst in opts): return {}
    for combo in itertools.product(*opts):
        for S, z in one(sp, G, combo, steps).items():
            if z is None: continue
            a = acc.get(pos[S])
            if a is None: acc[pos[S]] = list(z)
            else:
                for k in range(K): a[k] += z[k]
    ia = sp.R.inv(G['aut'])
    return {q: [x * ia % p for x in a] for q, a in acc.items()}


def job_check4(args):
    spec, t, ref = args; sp = get(spec); p = sp.p; t0 = time.time(); res = graph_dp(sp, t)
    return t, {str(k): v for k, v in res.items() if any(v)}, time.time() - t0


def job_check5(args):
    spec, t = args; sp = get(spec); t0 = time.time(); a = graph_dp(sp, t); t1 = time.time(); b = sp.graph(t); t2 = time.time()
    a = {k: v for k, v in a.items() if any(v)}; b = {k: v for k, v in b.items() if any(v)}
    G = sp.inner_graphs[t]
    return dict(t=t, E=G['E'], V=G['V'], val=[len(h) for h in G['hsAt']], equal=(a == b), dp=t1 - t0, std=t2 - t1, nouter=len(b))


def job_sample(args):
    spec, t, budget, mins, maxs, seed = args; sp = get(spec); t0 = time.time(); G, opts, pos, steps = prepare(sp, t); tp = time.time() - t0
    ncombo = 1
    for lst in opts: ncombo *= len(lst)
    out = dict(t=t, E=G['E'], V=G['V'], gv=list(G['gv']), val=[len(h) for h in G['hsAt']], ncombo=ncombo, nsub=len(pos), t_prep=tp)
    if ncombo == 0: out.update(samples=0, mean=0.0, tmax=0.0, est=tp, peak_mb=peak_mb()); return out
    rng = random.Random(seed * 100003 + t); n = 0; tsum = 0.0; tmax = 0.0; t1 = time.time()
    while n < maxs and (n < mins or time.time() - t1 < budget):
        combo = [rng.choice(lst) for lst in opts]; ta = time.time(); one(sp, G, combo, steps); dt = time.time() - ta
        n += 1; tsum += dt; tmax = max(tmax, dt)
    out.update(samples=n, mean=tsum / n, tmax=tmax, est=tp + ncombo * tsum / n, peak_mb=peak_mb())
    return out


def job_full(args):
    spec, t = args; sp = get(spec); G = sp.inner_graphs[t]; t0 = time.time(); res = graph_dp(sp, t); dt = time.time() - t0
    return dict(t=t, E=G['E'], V=G['V'], gv=list(G['gv']), val=[len(h) for h in G['hsAt']], full=dt, nouter=len(res), peak_mb=peak_mb())
