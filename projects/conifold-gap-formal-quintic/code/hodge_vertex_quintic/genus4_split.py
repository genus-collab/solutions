# -*- coding: utf-8 -*-
"""genus4_split.py -- driver for fq4_split.Split: the outer graph sum of the formal quintic by edge splitting, genus 2, 3, 4.

Modes (one of):
  --regress            genus 2 and 3 against the stored tables of the old vertex-tensor code (g3_state_<p>.json: A_id, A_M0;
                       g3_nat_<p>_nat.json: A3; genus-2 values of OU_EXTENSION_genus2/3.md): the validation of this code.
  --run                all inner graphs (or --tasks i,j,..) with a local process pool (--nproc); writes the part file and collects.
  --task I --ntasks N  cluster mode: the inner graphs of slice I of N (cost-balanced round robin) -> parts/<tag>_task<I>of<N>.json
  --collect            merge all part files of <tag>, run the checks, write g4_split_<tag>.json and the log.
  --timing             time the inner graphs listed in --tasks (one process), print seconds per graph.
Common options: --genus g (default 4)  --prime 89|61  --variant nat|M0|id  --ab alpha,beta  --K (default 2g-2)
tag = g<genus>_<prime>_<variant>[_ab..][_K..].  Results: {outer index: [u^0..u^(K-1)]}; a genus-g term is u^(-(2g-2)) G(u).
"""
import os, sys, json, argparse, time, re
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import P1, P2

_CACHE = {}


def get_split(spec):
    sp = _CACHE.get(spec)
    if sp is None:
        from fq4_split import Split
        prime, g, a, b, K, variant = spec
        sp = Split(P2 if prime == "89" else P1, g, (Fr(a), Fr(b)), K, variant); _CACHE[spec] = sp
    return sp


def job(args):
    spec, t = args; t0 = time.time(); sp = get_split(spec); res = sp.graph(t)
    return t, {str(k): v for k, v in res.items()}, time.time() - t0


def tag_of(a):
    ab = "" if a.ab == "0,0" else "_ab" + a.ab.replace("/", "o").replace("-", "m").replace(",", "_")
    return "g%d_%s_%s%s%s" % (a.genus, a.prime, a.variant, ab, "" if a.K is None else "_K%d" % a.K)


def cost_order(sp):
    """inner graph numbers sorted by a crude cost proxy (heaviest first)"""
    def proxy(G):
        E = G['E']; return (2 ** min(E, 6)) * (5 ** min(G['V'], 4)) * max(1, 12 - E) ** 2
    return sorted(range(len(sp.inner_graphs)), key=lambda t: -proxy(sp.inner_graphs[t]))


def merge(parts, K, p):
    """sum of the per-inner-graph pieces; an inner graph present in several part files must have identical pieces and is counted once"""
    tab = {}; seen = {}
    for part in parts:
        for t, res in part['res'].items():
            if t in seen:
                assert seen[t] == res, "inner graph %s occurs in two part files with different content" % t
                continue
            seen[t] = res
            for pos, c in res.items():
                a = tab.setdefault(int(pos), [0] * K)
                for k in range(K): a[k] = (a[k] + c[k]) % p
    return tab


def checks(sp, tab, log, rec, covered):
    R = sp.R; p = sp.p; g = sp.g; K = sp.K; FAILED = []; nin = len(sp.inner_graphs)
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    check("all %d inner graphs are covered by the part files" % nin, covered == set(range(nin)), sorted(set(range(nin)) - covered)[:20])
    if sp.frozen_ok is not None:
        check("[u^0] E_R(a, b) = Eh_0(a, b) for all (a, b), all colours", sp.frozen_ok)
        check("every entry of E_out = E_R - E_M is O(u) (needed for the restriction to fewer than K outer edges)", min(sp.valO.values()) >= 1)
    npol = 2 * g - 2; tot = [sum(tab.get(i, [0] * K)[k] for i in range(nin)) % p for k in range(K)]
    log("   total over all outer graphs: %s" % [rec.fmt(x) for x in tot])
    Eof = [G['E'] for G in sp.inner_graphs]
    for i in sorted(tab):
        c = tab[i]
        if i == 0 or any(c[k] for k in range(min(K, npol))): log("      outer graph %3d (E=%d, gv=%s): %s" % (i, Eof[i], sp.inner_graphs[i]['gv'], [rec.fmt(x) for x in c]))
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3
    tc = tc_series(sp.C, R, K + 4)
    cg = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}.get(g)
    T = tc_expansion3(tot, tc, R, g)[:npol]
    log("   total in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(T[k]), k - npol) for k in range(len(T))))
    if K >= npol and cg is not None: check("GAP: total = (%s) t_c^-%d + 0 + ... + 0 t_c^-1" % (cg, npol), len(T) == npol and T[0] == R.of(cg) and not any(T[1:npol]))
    if sp.variant != 'id':
        bad = {i: [k for k in range(min(K, npol)) if tab[i][k]] for i in tab if i != 0 and any(tab[i][k] for k in range(min(K, npol)))}
        check("REGULARITY (variant %s): every outer graph with edges is O(u^%d)" % (sp.variant, npol), not bad, {i: (Eof[i], sp.inner_graphs[i]['gv'], ks) for i, ks in bad.items()})
        V = tc_expansion3(tab.get(0, [0] * K), tc, R, g)[:npol]
        log("   vertex graph in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(V[k]), k - npol) for k in range(len(V))))
        if K >= npol and cg is not None: check("VERTEX GRAPH carries the polar part (%s) t_c^-%d + 0 + ... + 0 (equivalent to REGULARITY given GAP)" % (cg, npol), len(V) == npol and V[0] == R.of(cg) and not any(V[1:npol]))
    return FAILED, tot


def regress(a, log):
    from fq3_graphs import outer_graphs
    from fq3_rec import Rec
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    for prime in ("89", "61"):
        g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % prime))); nat = json.load(open(os.path.join(HERE, "g3_nat_%s_nat.json" % prime)))
        for variant, ref, kcmp in (("id", g3['A_id'], 6), ("M0", g3['A_M0'], 6), ("nat", nat['A3'], 4)):
            t0 = time.time(); sp = get_split((prime, 3, "0", "0", 6, variant)); tab = {}
            for t in range(len(sp.inner_graphs)):
                for pos, c in sp.graph(t).items():
                    acc = tab.setdefault(pos, [0] * 6)
                    for k in range(6): acc[k] = (acc[k] + c[k]) % sp.p
            mism = []
            for G in outer_graphs(3):
                mine = tab.get(sp.outer_pos[G['key']], [0] * 6)
                if [x % sp.p for x in ref[str(G['idx'])][:kcmp]] != mine[:kcmp]: mism.append(G['idx'])
            check("genus 3, 2^%s-1, variant %s: all 42 outer graphs equal the stored table through u^%d (%.0f s)" % (prime, variant, kcmp - 1, time.time() - t0), not mism, mism)
        R = sp.R; rec = Rec(R)
        for variant, exp in (("nat", {1: [Fr(1, 1200), Fr(-7, 6000), Fr(287, 864000)], 2: [0, 0, Fr(479, 756000)], 3: [0, 0, Fr(1, 2250)], 4: [0, 0, Fr(-127, 2016000)], 5: [0, 0, 0], 6: [0, 0, 0], 7: [0, 0, 0]}),
                             ("M0", {2: [0, Fr(1, 750)], 4: [0, Fr(-13, 8000)], 3: [0, 0], 5: [0, 0], 6: [0, 0], 7: [0, 0]})):
            sp = get_split((prime, 2, "0", "0", 4, variant)); tab = {}
            for t in range(len(sp.inner_graphs)):
                for pos, c in sp.graph(t).items():
                    acc = tab.setdefault(pos, [0] * 4)
                    for k in range(4): acc[k] = (acc[k] + c[k]) % sp.p
            mism = []
            for G in outer_graphs(2):
                if G['idx'] not in exp: continue
                mine = tab.get(sp.outer_pos[G['key']], [0] * 4); e = exp[G['idx']]
                if [R.of(x) for x in e] != mine[:len(e)]: mism.append((G['idx'], [rec.fmt(x) for x in mine]))
            check("genus 2, 2^%s-1, variant %s: outer graphs equal the values of the notes" % (prime, variant), not mism, mism)
    log("GENUS4_SPLIT REGRESSION %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", type=int, default=4); ap.add_argument("--prime", default="89"); ap.add_argument("--variant", default="nat")
    ap.add_argument("--ab", default="0,0"); ap.add_argument("--K", type=int, default=None); ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--regress", action="store_true"); ap.add_argument("--run", action="store_true"); ap.add_argument("--collect", action="store_true")
    ap.add_argument("--timing", action="store_true"); ap.add_argument("--task", type=int, default=None); ap.add_argument("--ntasks", type=int, default=None)
    ap.add_argument("--tasks", default="")
    a = ap.parse_args(); tag = tag_of(a); os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); os.makedirs(os.path.join(HERE, "parts"), exist_ok=True)
    LOG = open(os.path.join(HERE, "logs", "genus4_split_%s.log" % ("regress" if a.regress else tag)), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    if a.regress: log("=" * 100); log("genus4_split regression"); regress(a, log); return
    al, be = a.ab.split(","); spec = (a.prime, a.genus, al, be, a.K if a.K is not None else 2 * a.genus - 2, a.variant)
    log("=" * 100); log("genus4_split: %s" % tag)
    t0 = time.time(); sp = get_split(spec); order = cost_order(sp); nin = len(sp.inner_graphs)
    log("   setup %.0f s; %d inner graphs, K = %d, KZ = %d" % (time.time() - t0, nin, sp.K, sp.KZ))
    if a.timing:
        for t in [int(x) for x in a.tasks.split(",")]:
            t1 = time.time(); res = sp.graph(t); G = sp.inner_graphs[t]
            log("   inner graph %3d (E=%d, V=%d, gv=%s): %.1f s, %d outer graphs" % (t, G['E'], G['V'], G['gv'], time.time() - t1, len(res)))
        return
    if a.task is not None:
        mine = [t for r, t in enumerate(order) if r % a.ntasks == a.task]; res = {}; times = {}
        path = os.path.join(HERE, "parts", "%s_task%dof%d.json" % (tag, a.task, a.ntasks))
        if os.path.exists(path): prev = json.load(open(path)); res = prev['res']; times = prev['times']
        for t in mine:
            if str(t) in res: continue
            _, r, dt = job((spec, t)); res[str(t)] = r; times[str(t)] = dt
            json.dump(dict(tag=tag, res=res, times=times), open(path + ".tmp", "w")); os.replace(path + ".tmp", path)
            log("   task %d/%d: inner graph %d done in %.0f s" % (a.task, a.ntasks, t, dt))
        log("TASK %d of %d (%s) COMPLETE: %d inner graphs" % (a.task, a.ntasks, tag, len(mine))); return
    if a.run:
        import multiprocessing as mp
        todo = [int(x) for x in a.tasks.split(",")] if a.tasks else order; res = {}; times = {}
        with mp.Pool(processes=min(a.nproc, len(todo))) as pool:
            for k, (t, r, dt) in enumerate(pool.imap_unordered(job, [(spec, t) for t in todo])):
                res[str(t)] = r; times[str(t)] = dt
                if (k + 1) % max(1, len(todo) // 10) == 0: log("      %d/%d inner graphs done (%.0f s)" % (k + 1, len(todo), time.time() - t0))
        json.dump(dict(tag=tag, res=res, times=times), open(os.path.join(HERE, "parts", "%s_run.json" % tag), "w"))
        log("   run: %d inner graphs, wall %.0f s, CPU %.0f s, slowest %s" % (len(todo), time.time() - t0, sum(times.values()), sorted(((round(v), int(k)) for k, v in times.items()), reverse=True)[:5]))
    if a.collect or a.run:
        from fq3_rec import Rec
        pat = re.compile("^" + re.escape(tag) + r"_(run|task\d+of\d+)\.json$")          # exact tag: 'g4_89_nat' must not pick up 'g4_89_nat_t5_run.json'
        names = [f for f in sorted(os.listdir(os.path.join(HERE, "parts"))) if pat.match(f)]
        parts = [json.load(open(os.path.join(HERE, "parts", f))) for f in names]
        assert all(part['tag'] == tag for part in parts), "part file with a foreign tag"
        log("   part files: %s" % names)
        covered = set(int(t) for part in parts for t in part['res']); tab = merge(parts, sp.K, sp.p)
        cpu = sum(v for part in parts for v in part['times'].values())
        log("   collect: %d part files, %d inner graphs, total CPU %.0f s" % (len(parts), len(covered), cpu))
        FAILED, tot = checks(sp, tab, log, Rec(sp.R), covered)
        json.dump(dict(tag=tag, table={str(i): tab[i] for i in sorted(tab)}, total=tot, outer=[dict(E=G['E'], gv=list(G['gv']), ed=[list(e) for e in G['ed']], aut=G['aut']) for G in sp.inner_graphs], failed=FAILED),
                  open(os.path.join(HERE, "g4_split_%s.json" % tag), "w"))
        log("GENUS4_SPLIT (%s) %s" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))


if __name__ == "__main__":
    main()
