# -*- coding: utf-8 -*-
"""rv_g4mut.py -- genus-4 runs of a mutated copy of fq4_split.Split (rv_common.MUTANTS), judged by the checks of genus4_split.py.

Usage: py -I rv_g4mut.py <mutant> <variant> <prime 89|61> <popmax | full> [nproc (default 5)] [max minutes (default 25)]
  full      : all outer-edge subsets (the original algorithm); the verdict is the ORIGINAL genus4_split.checks (GAP, REGULARITY, VERTEX GRAPH).
  popmax = m: only the subsets with at most m outer edges (module global POPMAX of the copy).  This is EXACT for the outer graphs with at
              most m edges (an outer graph with E edges only receives subsets with |S| = E) and costs a few per cent of a full run.
              Verdict: (i) REGULARITY of the outer graphs with 1..m edges through u^(K-1); (ii) the vertex graph in t_c = gap;
              (iii) equality of these rows with the stored table of record g4_split_g4_<prime>_nat.json (project folder, read only).
A mutant is DETECTED when a check fails.  The run is abandoned (pool terminated, verdict INCOMPLETE) after the time limit.
Output: out/g4mut_<mutant>_<variant>_<prime>_<mode>.json and .log (nothing is written outside this folder).
"""
import os, sys, json, time
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT, MUTANTS, write_mutants, mutant_module
from fractions import Fraction as Fr

_SP = {}


def get_sp(name, prime, variant, popmax):
    key = (name, prime, variant, popmax); sp = _SP.get(key)
    if sp is None:
        from fq_fp import P1, P2
        mod = mutant_module(name); mod.POPMAX = popmax
        sp = mod.Split(P2 if prime == "89" else P1, 4, (Fr(0), Fr(0)), 6, variant); _SP[key] = sp
    return sp


def job(args):
    name, prime, variant, popmax, t = args; t0 = time.time()
    try:
        res = get_sp(name, prime, variant, popmax).graph(t)
        return t, {str(k): v for k, v in res.items()}, time.time() - t0, None
    except Exception as e:
        return t, {}, time.time() - t0, "%s: %s" % (type(e).__name__, str(e)[:200])


def main():
    name, variant, prime, mode = sys.argv[1:5]
    nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 5; maxmin = float(sys.argv[6]) if len(sys.argv) > 6 else 25.0
    popmax = 99 if mode == "full" else int(mode)
    write_mutants([name])
    tag = "g4mut_%s_%s_%s_%s" % (name, variant, prime, "full" if popmax == 99 else "pop%d" % popmax)
    LOG = open(os.path.join(OUT, tag + ".log"), "w", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    import multiprocessing as mp
    import genus4_split as drv
    from fq3_rec import Rec
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3
    t0 = time.time(); log("=" * 100); log("rv_g4mut %s: %s" % (tag, MUTANTS[name][1]))
    try:
        sp = get_sp(name, prime, variant, popmax)
    except Exception as e:
        log("CRASH in the constructor: %s: %s" % (type(e).__name__, str(e)[:300])); log("VERDICT %s: DETECTED (crash)" % tag)
        json.dump(dict(tag=tag, verdict="DETECTED (crash in the constructor)", failed=["constructor"]), open(os.path.join(OUT, tag + ".json"), "w")); return
    R = sp.R; p = sp.p; K = sp.K; nin = len(sp.inner_graphs); order = drv.cost_order(sp); Eof = [G['E'] for G in sp.inner_graphs]
    log("   setup %.0f s; %d inner graphs; min val E_out = %d; frozen_ok = %s; %d processes" % (time.time() - t0, nin, min(sp.valO.values()), sp.frozen_ok, nproc))
    tab = {}; covered = set(); cpu = 0.0; crashes = []; incomplete = False
    pool = mp.Pool(processes=nproc)
    try:
        it = pool.imap_unordered(job, [(name, prime, variant, popmax, t) for t in order])
        for k in range(nin):
            left = maxmin * 60 - (time.time() - t0)
            if left <= 0: incomplete = True; break
            try: t, res, dt, err = it.next(timeout=left)
            except mp.TimeoutError: incomplete = True; break
            cpu += dt
            if err: crashes.append((t, err)); continue
            covered.add(t)
            for pos, c in res.items():
                a = tab.setdefault(int(pos), [0] * K)
                for j in range(K): a[j] = (a[j] + c[j]) % p
            if (k + 1) % 76 == 0: log("      %d/%d inner graphs (%.0f s)" % (k + 1, nin, time.time() - t0))
    finally:
        pool.terminate(); pool.join()
    log("   run: wall %.0f s, CPU %.0f s, covered %d/%d, crashes %d %s" % (time.time() - t0, cpu, len(covered), nin, len(crashes), crashes[:3]))
    rec = Rec(R); FAILED = []
    if incomplete: FAILED.append("INCOMPLETE (time limit %.0f min)" % maxmin)
    if crashes: FAILED.append("crash in %d inner graphs" % len(crashes))
    if not incomplete and not crashes:
        if popmax == 99:
            F, tot = drv.checks(sp, tab, log, rec, covered); FAILED += F
        else:
            def check(nm, ok, detail=""):
                log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail else ""))
                if not ok: FAILED.append(nm)
            I = [i for i in range(nin) if Eof[i] <= popmax]; zero = [0] * K
            bad = {i: [k for k in range(K) if tab.get(i, zero)[k]] for i in I if i != 0 and any(tab.get(i, zero))}
            check("REGULARITY restricted: the %d outer graphs with 1..%d edges are O(u^%d)" % (len(I) - 1, popmax, K), not bad,
                  {i: (Eof[i], sp.inner_graphs[i]['gv'], ks) for i, ks in sorted(bad.items())[:12]})
            log("      number of non-regular outer graphs: %d of %d; by E: %s" % (len(bad), len(I) - 1, {E: sum(1 for i in bad if Eof[i] == E) for E in range(1, popmax + 1)}))
            if variant != 'id':
                tc = tc_series(sp.C, R, K + 4); V = tc_expansion3(tab.get(0, zero), tc, R, 4)[:6]
                log("   vertex graph: %s" % [rec.fmt(x) for x in tab.get(0, zero)])
                log("   vertex graph in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(V[k]), k - 6) for k in range(6)))
                check("VERTEX GRAPH carries the polar part (1/180000) t_c^-6 + 0 + ... + 0", V[0] == R.of(Fr(1, 180000)) and not any(V[1:6]))
            ref = os.path.join(PROJ, "g4_split_g4_%s_nat.json" % prime)
            if variant == 'nat' and os.path.exists(ref):
                st = json.load(open(ref))['table']; diff = [i for i in I if [x % p for x in st.get(str(i), zero)] != tab.get(i, zero)]
                check("rows of the %d outer graphs with <= %d edges equal the stored table of record" % (len(I), popmax), not diff, diff[:12])
    verdict = "NOT DETECTED (all checks pass)" if not FAILED else "DETECTED: %s" % FAILED
    if name == 'baseline': verdict = "control ok" if not FAILED else "CONTROL BROKEN: %s" % FAILED
    log("VERDICT %s: %s  (wall %.0f s)" % (tag, verdict, time.time() - t0))
    json.dump(dict(tag=tag, desc=MUTANTS[name][1], verdict=verdict, failed=FAILED, table={str(i): tab[i] for i in sorted(tab)}, cpu=cpu, wall=time.time() - t0),
              open(os.path.join(OUT, tag + ".json"), "w"))


if __name__ == "__main__":
    main()
