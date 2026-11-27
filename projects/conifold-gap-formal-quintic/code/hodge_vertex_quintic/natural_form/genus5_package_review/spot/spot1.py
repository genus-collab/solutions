# -*- coding: utf-8 -*-
"""spot1.py -- genus-5 spot checks of g5_engine.graph_dp against fq4_split.Split.graph (slow path), on a COPY of the package (spot/pkg).

  py -I -B spot1.py [--nproc 6] [--prime 89]
Per inner graph t of LIST (natural M, K = 8):
  A = graph_dp(sp, t)             (package fast path, all outer-edge subsets with <= K-1 outer edges)
  B = sp.graph(t)                 (original slow path of fq4_split, Split built by g5_engine.make_split)
  compare A and B (rows that are zero mod p dropped);  sum_q A[q], sum_q B[q] against parts/g5_<prime>_id/local.json piece of t.
Per inner graph t of LIST2:
  C = graph_dp(sp, t, maxout=2)   D = split_graph_restricted(sp, t, 2) (copy of Split.graph, outer-edge sets of size <= 2)
  compare C and D; and C against A on the outer graphs with <= 2 edges (exactness claim of maxout).
Output: out/spot1_<prime>.json and out/spot1_<prime>.log.
"""
import os, sys, json, time, itertools, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.join(HERE, "pkg"); OUT = os.path.join(HERE, "out")
sys.path.insert(0, PKG)

LIST = [1, 3, 5, 25, 28, 118, 119, 149, 245, 383, 444, 1234, 2318, 2539, 2668, 3638, 3847, 4356, 4484, 4552]
LIST2 = [383, 2318, 3638, 3847]
_SP = {}


def split(prime):
    if prime not in _SP:
        from fq_fp import P1, P2
        from g5_engine import make_split
        _SP[prime] = make_split(P2 if prime == "89" else P1, 5, (0, 0), 8, 'nat')
    return _SP[prime]


def split_graph_restricted(sp, t, m):
    """copy of fq4_split.Split.graph with the subsets restricted to at most m outer edges (and < K)"""
    G = sp.inner_graphs[t]; V = G['V']; ed = G['ed']; hsAt = G['hsAt']; E = len(ed); K = sp.K; p = sp.p
    opts = []
    for v in range(V):
        tab = sp.vertex_table(G['gv'][v], len(hsAt[v])); lst = []
        for hp, vec in tab.items():
            for perm in set(itertools.permutations(hp)): lst.append((perm, vec))
        if not lst: return {}
        opts.append(lst)
    subsets = [(S, sp.outer_pos[sp.contract(G, S)]) for S in range(1 << E) if bin(S).count("1") < K and bin(S).count("1") <= m]
    steps = sp.plan(G); acc = {}
    for combo in itertools.product(*opts):
        hp = [0] * (2 * E)
        for v in range(V):
            for r, h in enumerate(hsAt[v]): hp[h] = combo[v][0][r]
        keys = [(hp[2 * e], hp[2 * e + 1]) for e in range(E)]
        vecs = [combo[v][1] for v in range(V)]
        vals = [sp.valO[k] for k in keys]; inner = [sp.EM[k] for k in keys]; outer = [sp.EO[k] for k in keys]
        for S, pos in subsets:
            if sum(vals[e] for e in range(E) if (S >> e) & 1) >= K: continue
            z = sp.colour_sum(steps, vecs, [outer[e] if (S >> e) & 1 else inner[e] for e in range(E)])
            if z is None: continue
            a = acc.get(pos)
            if a is None: acc[pos] = list(z)
            else:
                for k in range(K): a[k] += z[k]
    ia = sp.R.inv(G['aut'])
    return {pos: [x * ia % p for x in a] for pos, a in acc.items()}


def norm(d, p):
    out = {}
    for q, v in d.items():
        w = [x % p for x in v]
        if any(w): out[int(q)] = w
    return out


def job(args):
    kind, prime, t = args
    from g5_engine import graph_dp
    sp = split(prime); p = sp.p
    if kind == "full":
        t0 = time.time(); A = norm(graph_dp(sp, t), p); t1 = time.time(); B = norm(sp.graph(t), p); t2 = time.time()
        return kind, t, dict(A=A, B=B, tA=t1 - t0, tB=t2 - t1)
    t0 = time.time(); C = norm(graph_dp(sp, t, 2), p); t1 = time.time(); D = norm(split_graph_restricted(sp, t, 2), p); t2 = time.time()
    return kind, t, dict(C=C, D=D, tC=t1 - t0, tD=t2 - t1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nproc", type=int, default=6); ap.add_argument("--prime", default="89")
    ap.add_argument("--list", default=""); ap.add_argument("--list2", default=None)
    a = ap.parse_args(); assert a.nproc <= 6
    lst = [int(x) for x in a.list.split(",")] if a.list else LIST
    lst2 = LIST2 if a.list2 is None else [int(x) for x in a.list2.split(",") if x]
    os.makedirs(OUT, exist_ok=True); tag = "spot1_%s" % a.prime; LOG = open(os.path.join(OUT, tag + ".log"), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    T0 = time.time(); sp = split(a.prime); p = sp.p; L = sp.inner_graphs
    log("=" * 100); log("spot1 %s: python %s, %d processes, set-up %.1f s, %d inner graphs" % (tag, sys.version.split()[0], a.nproc, time.time() - T0, len(L)))
    idp = json.load(open(os.path.join(PKG, "parts", "g5_%s_id" % a.prime, "local.json")))
    assert idp['tag'] == "g5_%s_id" % a.prime
    jobs = [("full", a.prime, t) for t in lst] + [("out2", a.prime, t) for t in lst2]
    import multiprocessing as mp
    res = []
    with mp.Pool(processes=min(a.nproc, len(jobs))) as pool:
        for r in pool.imap_unordered(job, sorted(jobs, key=lambda j: -j[2])):
            res.append(r); log("   done %s %4d (%.0f s wall)" % (r[0], r[1], time.time() - T0))
    full = {t: d for k, t, d in res if k == "full"}; out2 = {t: d for k, t, d in res if k == "out2"}
    FAIL = []
    zero = [0] * sp.K
    log("-- (1) graph_dp vs Split.graph, all subsets, and (2) additivity against the shipped id piece")
    for t in lst:
        d = full[t]; G = L[t]; A = d['A']; B = d['B']
        eq = (A == B)
        sA = [sum(v[k] for v in A.values()) % p for k in range(sp.K)]; sB = [sum(v[k] for v in B.values()) % p for k in range(sp.K)]
        ref = [x % p for x in idp['pieces'][str(t)].get(str(t), zero)]
        addA = (sA == ref); addB = (sB == ref)
        byE = {}
        for q in A: byE[L[q]['E']] = byE.get(L[q]['E'], 0) + 1
        log("   t=%4d E=%2d V=%d gv=%s aut=%d: dp %.1f s, split %.1f s; %d non-zero outer rows (by E_out %s); dp == split: %s; sum = id piece: dp %s, split %s; id piece non-zero: %s"
            % (t, G['E'], G['V'], list(G['gv']), G['aut'], d['tA'], d['tB'], len(A), sorted(byE.items()), eq, addA, addB, any(ref)))
        if not (eq and addA and addB): FAIL.append(("full", t, eq, addA, addB))
    log("-- maxout = 2: graph_dp(maxout=2) vs restricted Split.graph, and vs the full graph_dp on the outer graphs with <= 2 edges")
    for t in lst2:
        d = out2[t]; C = d['C']; D = d['D']; eq = (C == D)
        small = lambda X: {q: v for q, v in X.items() if L[q]['E'] <= 2}
        ex = None
        if t in full: ex = (small(C) == small(full[t]['A']))
        extra = sorted(set(q for q in C if L[q]['E'] > 2))
        log("   t=%4d E=%2d: dp(maxout 2) %.1f s, restricted split %.1f s; %d rows; equal: %s; rows with <= 2 edges equal the full dp: %s; rows with > 2 edges (must be none): %d"
            % (t, L[t]['E'], d['tC'], d['tD'], len(C), eq, ex, len(extra)))
        if not eq or ex is False or extra: FAIL.append(("out2", t, eq, ex, len(extra)))
    json.dump(dict(prime=a.prime, full={str(t): dict(A={str(q): v for q, v in full[t]['A'].items()}, tA=full[t]['tA'], tB=full[t]['tB'],
                                                    equal=full[t]['A'] == full[t]['B']) for t in full},
                   out2={str(t): dict(C={str(q): v for q, v in out2[t]['C'].items()}, equal=out2[t]['C'] == out2[t]['D']) for t in out2}, failed=FAIL),
              open(os.path.join(OUT, tag + ".json"), "w"))
    log("SPOT1 (%s) %s, wall %.0f s" % (tag, "ALL EQUAL" if not FAIL else "FAILURES %s" % FAIL, time.time() - T0)); LOG.close()


if __name__ == "__main__":
    main()
