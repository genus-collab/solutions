# -*- coding: utf-8 -*-
"""ind_run.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): driver.

  py -I ind_run.py --genus 4 --variant nat --prime 89 [--K 6] [--max-outer n] [--nproc 6] [--deadline 1500] [--tasks i,j,..]
Runs ind_sum.Summer.graph on the inner graphs (process pool of at most 6), stores every finished inner graph in
out/<tag>_parts.json (a rerun resumes), stops cleanly at the deadline (seconds of wall time, default 1500 < 30 minutes), and when
all inner graphs are present collects the table {outer index: [u^0..u^(K-1)]}, compares it with the project's stored tables and
writes out/<tag>.json and out/<tag>.log.  tag = g<genus>_<prime>_<variant>[_E<max-outer>][_K<K>].
Reference tables: genus 4: <project>/g4_split_g4_<prime>_<variant>.json (missing index = 0);
  genus 3: g3_state_<prime>.json (A_id, A_M0, six coefficients), g3_nat_<prime>_nat.json (A3: natural M of the OLD vertex-tensor
  code, valid through u^3); genus 2: the values of the notes.  --max-outer n: only subsets with <= n outer edges, so only the
outer graphs with <= n edges are complete (and compared).  Other modes: --brute i,j,.. compares graph and graph_brute on the
listed inner graphs;  --cost prints the cost proxy.
"""
import os, sys, json, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out")
_W = {}


def setup(spec):
    if _W.get('spec') != spec:
        from ind_tables import Tables
        from ind_graphs import graphs
        from ind_sum import Summer
        prime, g, K, variant = spec
        T = Tables(prime, g, K, variant); L = graphs(g)
        _W.update(spec=spec, T=T, L=L, S=Summer(T, L))
    return _W['S']


def job(args):
    spec, t, max_outer = args; t0 = time.time(); S = setup(spec)
    res = S.graph(t, max_outer)
    return t, {str(k): v for k, v in res.items()}, time.time() - t0


def proxy(S, t):
    G = S.graphs[t]; st = S.structure(t); n = 1
    for o in S.options(G): n *= len(o)
    return n * len(st['cfgs']) * len(st['subsets'])


def tc_coeffs(T, g):
    """first K coefficients of (u/t_c)^(2g-2)"""
    from ref_fp import sinv, smul
    p = T.p; K = T.K; r = sinv([T.H.tc[k + 1] for k in range(K)], p); out = [1] + [0] * (K - 1)
    for _ in range(2 * g - 2): out = smul(out, r, p, K)
    return out


def reference(g, prime, variant, proj):
    """{outer position: coefficients}, number of coefficients to compare, description"""
    if g == 4:
        path = os.path.join(proj, "g4_split_g4_%s_%s.json" % (prime, variant))
        if not os.path.exists(path): return None, 0, "no reference file"
        d = json.load(open(path)); return {int(k): v for k, v in d['table'].items()}, 6, os.path.basename(path)
    if g == 3:
        if variant == "nat":
            d = json.load(open(os.path.join(proj, "g3_nat_%s_nat.json" % prime)))['A3']; n = 4; name = "g3_nat_%s_nat.json:A3 (through u^3)" % prime
        elif variant in ("id", "M0"):
            d = json.load(open(os.path.join(proj, "g3_state_%s.json" % prime)))['A_' + variant]; n = 6; name = "g3_state_%s.json:A_%s" % (prime, variant)
        else: return None, 0, "no reference"
        return {int(k) - 1: v for k, v in d.items()}, n, name
    if g == 2:
        nat = {0: [Fr(1, 1200), Fr(-7, 6000), Fr(287, 864000)], 1: [0, 0, Fr(479, 756000)], 2: [0, 0, Fr(1, 2250)], 3: [0, 0, Fr(-127, 2016000)],
               4: [0, 0, 0], 5: [0, 0, 0], 6: [0, 0, 0]}
        if variant == "nat": return nat, 3, "values of the notes (OU_EXTENSION, genus 2)"
    return None, 0, "no reference"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", type=int, default=4); ap.add_argument("--prime", default="89"); ap.add_argument("--variant", default="nat")
    ap.add_argument("--K", type=int, default=None); ap.add_argument("--max-outer", type=int, default=None); ap.add_argument("--nproc", type=int, default=6)
    ap.add_argument("--deadline", type=float, default=1500.0); ap.add_argument("--tasks", default=""); ap.add_argument("--brute", default="")
    ap.add_argument("--cost", action="store_true")
    a = ap.parse_args(); g = a.genus; K = a.K if a.K is not None else max(2 * g - 2, 4); t00 = time.time()
    assert a.nproc <= 6, "at most 6 processes"
    tag = "g%d_%s_%s%s%s" % (g, a.prime, a.variant, "" if a.max_outer is None else "_E%d" % a.max_outer, "" if a.K is None else "_K%d" % K)
    os.makedirs(OUT, exist_ok=True); LOG = open(os.path.join(OUT, tag + ".log"), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    spec = (a.prime, g, K, a.variant); S = setup(spec); T = S.T; p = T.p; L = S.graphs; nin = len(L)
    from ref_fp import of, fmt
    F = lambda x: fmt(x, p)
    log("=" * 100); log("ind_run %s: %d inner graphs, K = %d, NZ = %d, zeta rows checks %s, setup %.1f s" % (tag, nin, K, T.NZ, T.H.checks, time.time() - t00))
    if a.cost:
        c = sorted(((proxy(S, t), t) for t in range(nin)), reverse=True)
        log("   cost proxy total %d; top: %s" % (sum(x for x, t in c), c[:12])); return
    if a.brute:
        bad = []
        for t in [int(x) for x in a.brute.split(",")]:
            t1 = time.time(); A = S.graph(t, a.max_outer); t2 = time.time(); B = S.graph_brute(t, a.max_outer); t3 = time.time()
            ok = all(A.get(k, [0] * K) == B.get(k, [0] * K) for k in set(A) | set(B))
            log("   inner graph %3d (E=%d, V=%d, gv=%s): momentum %.2f s, brute %.2f s, %d outer graphs, equal: %s" % (t, L[t]['E'], L[t]['V'], L[t]['gv'], t2 - t1, t3 - t2, len(B), ok))
            if not ok: bad.append(t)
        log("BRUTE-FORCE CHECK (%s) %s" % (tag, "ALL EQUAL" if not bad else "MISMATCH %s" % bad)); return
    # ---------------------------------------------------------------- run / resume
    ppath = os.path.join(OUT, tag + "_parts.json"); res = {}; times = {}
    if os.path.exists(ppath): d = json.load(open(ppath)); res = d['res']; times = d['times']
    want = [int(x) for x in a.tasks.split(",")] if a.tasks else list(range(nin))
    todo = sorted([t for t in want if str(t) not in res], key=lambda t: -proxy(S, t))
    def save():
        json.dump(dict(tag=tag, res=res, times=times), open(ppath + ".new", "w")); os.replace(ppath + ".new", ppath)
    if todo:
        import multiprocessing as mp
        log("   %d inner graphs to do with %d processes, deadline %.0f s" % (len(todo), min(a.nproc, len(todo)), a.deadline))
        pool = mp.Pool(processes=min(a.nproc, len(todo))); it = pool.imap_unordered(job, [(spec, t, a.max_outer) for t in todo]); done = 0; last = time.time()
        try:
            while done < len(todo):
                left = a.deadline - (time.time() - t00)
                if left <= 0: raise mp.TimeoutError()
                t, r, dt = it.next(timeout=left); res[str(t)] = r; times[str(t)] = dt; done += 1
                if time.time() - last > 60 or done == len(todo):
                    save(); last = time.time(); log("      %d/%d done (%.0f s wall)" % (done, len(todo), time.time() - t00))
            pool.close(); pool.join()
        except mp.TimeoutError:
            pool.terminate(); pool.join(); save()
            log("DEADLINE reached after %.0f s: %d of %d inner graphs stored; rerun the same command to resume" % (time.time() - t00, len(res), nin)); return
        save()
        log("   run: wall %.0f s, CPU (this and earlier runs) %.0f s, slowest %s" % (time.time() - t00, sum(times.values()), sorted(((round(v, 1), int(k)) for k, v in times.items()), reverse=True)[:6]))
    if len(res) < nin and not a.tasks: log("INCOMPLETE: %d of %d" % (len(res), nin)); return
    # ---------------------------------------------------------------- collect and compare
    tab = {}
    for t, r in res.items():
        for pos, c in r.items():
            acc = tab.setdefault(int(pos), [0] * K)
            for k in range(K): acc[k] = (acc[k] + c[k]) % p
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    full = (len(res) == nin)
    complete = [i for i in range(nin) if full and (a.max_outer is None or a.variant == "id" or L[i]['E'] <= a.max_outer)]
    proj = os.path.normpath(os.path.join(HERE, "..", "..", "..")); ref, ncmp, name = reference(g, a.prime, a.variant, proj); ncmp = min(ncmp, K)
    zero = [0] * K; npol = 2 * g - 2
    log("   complete outer graphs: %d of %d (%s)" % (len(complete), nin, "all" if len(complete) == nin else "those with at most %s edges" % a.max_outer))
    if ref is not None and complete:
        mism = [i for i in complete if [of(x, p) if isinstance(x, Fr) else x % p for x in ref.get(i, zero)[:ncmp]] != tab.get(i, zero)[:ncmp]]
        check("%d outer graphs equal the reference %s, coefficients u^0..u^%d" % (len(complete), name, ncmp - 1), not mism, mism[:30])
        byE = {}
        for i in complete: byE[L[i]['E']] = byE.get(L[i]['E'], 0) + 1
        log("      compared, by number of outer edges: %s; non-zero in my table: %d" % (sorted(byE.items()), sum(1 for i in complete if any(tab.get(i, zero)))))
    if 0 in tab: log("   vertex graph: %s" % [F(x) for x in tab[0]])
    for i in sorted(tab):
        if a.variant != "id" and i and i in complete and any(tab[i][:min(K, npol)]): log("      outer graph %3d (E=%d, gv=%s, ed=%s): %s" % (i, L[i]['E'], L[i]['gv'], L[i]['ed'], [F(x) for x in tab[i]]))
    cg = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}[g]; gap = [of(cg, p) * x % p for x in tc_coeffs(T, g)][:min(K, npol)]
    if a.variant != "id" and complete:
        bad = [i for i in complete if i and any(tab.get(i, zero)[:min(K, npol)])]
        if a.variant != "nat":
            log("INFO (control variant %s, regularity NOT expected): %d of the %d complete outer graphs with edges are not O(u^%d); vertex graph has the pure gap form: %s"
                % (a.variant, len(bad), len(complete) - 1, min(K, npol), tab.get(0, zero)[:len(gap)] == gap))
        else:
            check("REGULARITY (%s): the %d complete outer graphs with edges are O(u^%d)" % (a.variant, len(complete) - 1, min(K, npol)), not bad, bad[:30])
            check("VERTEX GRAPH = (%s) t_c^-%d + 0 + .. + 0, i.e. [u^0..u^%d] G = (%s) (u/t_c)^%d" % (cg, npol, len(gap) - 1, cg, npol), tab.get(0, zero)[:len(gap)] == gap)
    tot = None
    if len(complete) == nin:
        tot = [sum(tab.get(i, zero)[k] for i in range(nin)) % p for k in range(K)]
        log("   total: %s" % [F(x) for x in tot])
        check("GAP: total = (%s) t_c^-%d + 0 + .. + 0" % (cg, npol), tot[:len(gap)] == gap)
        if a.variant != "id":
            ipath = os.path.join(OUT, "g%d_%s_id%s.json" % (g, a.prime, "" if a.K is None else "_K%d" % K))
            if os.path.exists(ipath):
                idres = json.load(open(ipath)); check("total of this variant = total of my point-frame run (variant id)", idres['total'] == tot)
                # per inner graph: the split pieces add up to the point-frame term
                ipp = json.load(open(ipath.replace(".json", "_parts.json")))['res']; bad = []
                for t in range(nin):
                    mine = [sum(c[k] for c in res[str(t)].values()) % p for k in range(K)]
                    if mine != ipp[str(t)].get(str(t), zero): bad.append(t)
                check("for each of the %d inner graphs the split pieces add up to my point-frame term" % nin, not bad, bad[:20])
    json.dump(dict(tag=tag, table={str(i): tab[i] for i in sorted(tab)}, total=tot, complete=complete, failed=FAILED,
                   outer=[dict(E=G['E'], gv=list(G['gv']), ed=[list(e) for e in G['ed']], aut=G['aut']) for G in L]), open(os.path.join(OUT, tag + ".json"), "w"))
    log("IND_RUN (%s) %s" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))


if __name__ == "__main__":
    main()
