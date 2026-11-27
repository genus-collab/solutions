# -*- coding: utf-8 -*-
"""g5_validate.py -- validation and calibration of the genus-5 package (genus5_run.py --validate; V_validate.sh on the cluster).

 (1) interpreter and graph lists: Python >= 3.8; 7, 42, 379, 4555 leg-free stable graphs of genus 2..5; |Aut| of the genus 2-4 lists equal
     to the brute-force count fq_graphs.aut_count; the cost table of genus 5 refers to the same graphs (genera and edges, index by index).
 (2) genus 2 and 3, both primes: variants nat, M0 (graph_dp) and id (point frame) against the stored tables of the OLD vertex-tensor code
     (data/ref_g23.json: genus 3 through u^5 for id, M0, through u^3 for nat; genus 2 values of the notes).
 (3) genus 4: eight runs against the tables of record (data/ref_g4.json, all 379 outer graphs through u^5): nat and id at both primes, the
     controls M0 and nat_t5 at 2^89-1, the second conifold frame (1/3, -2/7) at both primes; the collect checks on nat at 2^89-1.
 (4) genus 5 set-up at both primes: all assertions of fq4_split.Split, E_out = O(u) for all keys, frozen symbol, tau(u) Sum 4u^n/(5n+4) = 1.
 (5) calibration: three genus-5 inner graphs (808, 1109, 4511) timed against the reference workstation; prints the machine factor f and the
     expected wall time per array task for several array sizes.
Log: logs/validate.log; last line "G5_VALIDATE ALL CHECKS PASSED" or the list of failures.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REF_DP = {808: 10.81, 1109: 20.92, 4511: 3.49}          # graph_dp seconds on the reference workstation (feasibility study, 2^89-1, natural M)


def key_str(G): return json.dumps([list(G['gv']), [list(e) for e in G['ed']]])


def main(nproc):
    import genus5_run as D, g5_checks as C, g5_graphs
    from fq_graphs import aut_count
    from fq3_rec import Rec
    log, f = D.opener(os.path.join(HERE, "logs", "validate.log")); FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    T0 = time.time(); log("=" * 90); log("g5_validate: python %s, %d processes, host %s" % (sys.version.split()[0], nproc, __import__("socket").gethostname()))
    # (1)
    check("Python >= 3.8", sys.version_info >= (3, 8), sys.version.split()[0])
    cnt = {g: len(g5_graphs.graphs(g, 3 * g - 3)) for g in (2, 3, 4, 5)}
    check("graph lists: 7, 42, 379, 4555 graphs of genus 2..5", cnt == {2: 7, 3: 42, 4: 379, 5: 4555}, cnt)
    bad = 0
    for g in (2, 3, 4):
        for G in g5_graphs.graphs(g, 3 * g - 3):
            V = G['V']; A = [[0] * V for _ in range(V)]
            for u, w in G['ed']:
                if u == w: A[u][u] += 1
                else: A[u][w] += 1; A[w][u] += 1
            if aut_count(list(G['gv']), A) != G['aut']: bad += 1
    check("|Aut| of the 428 graphs of genus 2-4 = brute-force count", bad == 0, bad)
    ct = json.load(open(os.path.join(HERE, "data", "cost_table_g5.json"))); G5 = g5_graphs.graphs(5, 12)
    check("cost table of genus 5 refers to the same graphs (genera, edges) index by index",
          all(tuple(ct[str(t)]['gv']) == G5[t]['gv'] and [tuple(e) for e in ct[str(t)]['ed']] == list(G5[t]['ed']) for t in range(4555)))
    # (2)
    ref23 = json.load(open(os.path.join(HERE, "data", "ref_g23.json")))
    for g, K in ((2, 4), (3, 6)):
        for prime in ("89", "61"):
            for variant in ("nat", "M0", "id"):
                name = "g%d_%s_%s" % (g, prime, variant)
                if name not in ref23: continue
                t0 = time.time(); spec = (prime, g, "0,0", K, variant); sp = D.get_split(spec); tab = {}
                for t, res, dt in D.pool_map([(spec, t, None) for t in range(len(sp.inner_graphs))], 1):
                    for q, c in res.items():
                        a = tab.setdefault(int(q), [0] * K)
                        for k in range(K): a[k] = (a[k] + c[k]) % sp.p
                ref = ref23[name]; kc = ref['K']; mine = {key_str(sp.inner_graphs[i]): tab.get(i, [0] * K)[:kc] for i in range(len(sp.inner_graphs))}
                mism = [k for k, v in ref['table'].items() if mine.get(k) != v]
                check("genus %d, 2^%s-1, %s: %d outer graphs equal the stored table through u^%d (%.0f s)" % (g, prime, variant, len(ref['table']), kc - 1, time.time() - t0), not mism, mism[:5])
    # (3)
    ref4 = json.load(open(os.path.join(HERE, "data", "ref_g4.json")))
    runs = [("89", "nat", "0,0", "g4_89_nat"), ("61", "nat", "0,0", "g4_61_nat"), ("89", "id", "0,0", "g4_89_id"), ("61", "id", "0,0", "g4_61_id"),
            ("89", "M0", "0,0", "g4_89_M0"), ("89", "nat_t5", "0,0", "g4_89_nat_t5"), ("89", "nat", "1/3,-2/7", "g4_89_nat_ab1o3_m2o7"),
            ("61", "nat", "1/3,-2/7", "g4_61_nat_ab1o3_m2o7")]
    for prime, variant, ab, name in runs:
        t0 = time.time(); spec = (prime, 4, ab, 6, variant); sp = D.get_split(spec); pieces = {}
        for t, res, dt in D.pool_map([(spec, t, None) for t in range(379)], nproc): pieces[str(t)] = res
        tab = C.merge(pieces, 6, sp.p)
        if variant == 'id': tab = {i: v for i, v in tab.items()}
        mine = {key_str(sp.inner_graphs[i]): tab.get(i, [0] * 6) for i in range(379)}
        mism = [k for k, v in ref4[name]['table'].items() if mine.get(k) != v]
        check("genus 4, %s: all 379 outer graphs equal the table of record through u^5 (%.0f s wall)" % (name, time.time() - t0), not mism, mism[:3])
        if name == "g4_89_nat":
            Fc, info = C.checks(sp, tab, set(range(379)), None, log, Rec(sp.R))
            check("genus 4, g4_89_nat: the collect checks pass and c_4 = 1/180000 matches the Bernoulli pattern", not Fc and info.get('c_total_equals_pattern') and info.get('c_vertex_equals_pattern'))
    # (4)
    from g5_engine import tau_identity_ok
    for prime in ("89", "61"):
        t0 = time.time(); sp = D.get_split((prime, 5, "0,0", 8, "nat"))
        check("genus 5 set-up, 2^%s-1 (%.0f s): Split assertions, frozen symbol, E_out = O(u) for all %d keys, tau identity" % (prime, time.time() - t0, len(sp.valO)),
              bool(sp.frozen_ok) and min(sp.valO.values()) >= 1 and tau_identity_ok(sp))
    # (5)
    spec = ("89", 5, "0,0", 8, "nat"); t0 = time.time()
    res = D.pool_map([(spec, t, None) for t in sorted(REF_DP)], min(nproc, 3))
    meas = {t: dt for t, r, dt in res}; fac = sum(meas.values()) / sum(REF_DP.values())
    log("   calibration (wall %.0f s): %s; machine factor f = %.2f (cluster time / reference time)" % (time.time() - t0, {t: "%.1f s (ref %.1f s)" % (meas[t], REF_DP[t]) for t in sorted(meas)}, fac))
    costs = C.costs_for(5, 4555, "nat")
    for n in (25, 50, 100, 200):
        _, load = C.lpt_assign(costs, n)
        log("   NTASKS = %3d: heaviest task %.2f h, mean %.2f h on this machine (estimate x f)" % (n, max(load) * fac / 3600, sum(load) * fac / n / 3600))
    log("   total CPU per full natural-M run on this machine: %.0f h (estimate x f)" % (sum(costs.values()) * fac / 3600))
    log("   validation wall time %.0f s" % (time.time() - T0))
    log("G5_VALIDATE ALL CHECKS PASSED" if not FAILED else "G5_VALIDATE FAILURES: %s" % FAILED); f.close()
    return not FAILED
