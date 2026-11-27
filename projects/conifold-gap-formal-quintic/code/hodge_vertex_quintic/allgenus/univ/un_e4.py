# -*- coding: utf-8 -*-
"""un_e4.py -- E4/U4: genus-4 universality of the regular-gauge factorisation, WITH and WITHOUT the canonical normalisation.

For a case of un_e4_cases (cyclic system + normalisation of its Theorem-U gauge P to the model) the hatted rows Qhat = Ghat Rhat(tau zhat)
are formed (zhat^0..9, u^0..5) and the outer graph sum R'.(M.omega^), M = M_0(tau zhat), is computed over all 379 genus-4 inner graphs
by un_engine.Var (edge splitting, fq4_split.Split imported read-only), prime 2^89-1.  A genus-4 term is u^-6 G(u); an outer graph is
regular iff [u^0..u^5] G = 0.  --maxout m keeps the outer-edge subsets of size <= m: exact for the vertex graph and the outer graphs
with <= m edges (the others are then incomplete and are NOT reported).
Modes:  prep CASES            gauge checks ((G), regularity), rows, symplecticity, row regressions (quintic cases)
        run CASE [--maxout m] [--nproc n]   graph sum; table -> e4_tab_<case>_<prime>[_out<m>].json; log e4.log
        report                compare the stored tables (a vs the natural run of record, b - a, d - c1, ...)
Run: py -I -B un_e4.py prep a,b,c1,c2,d  |  py -I -B un_e4.py run a --nproc 5  |  py -I -B un_e4.py report
"""
import os, sys, json, time, argparse, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)

LOGF = os.path.join(HERE, "e4.log"); FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True)
    with open(LOGF, "a", encoding="utf-8") as f: f.write(s + "\n")


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def tab_path(case, prime, maxout):
    return os.path.join(HERE, "e4_tab_%s_%s%s.json" % (case, prime, "" if maxout is None else "_out%d" % maxout))


def prep(case, prime, g=4):
    import un_e4_cases as C
    res, info = C.build(case)
    log("-" * 110); log("case %s: %s" % (case, info['desc']))
    log("   P^stored_11(0, z) = %s;  multiplier m(z) = %s" % (info['P11_stored'], info['mult']))
    C.gauge_checks(res, info, check, log)
    Qh, nat = C.rows(res, prime, g)
    QL = C.row_checks(Qh, nat, info, prime, check, log, g)
    return res, info, Qh, QL, nat


def report_table(tab, outer, g, maxout, rec, R, res, info, nat, prime):
    from fq_series import Ser
    from fq3_msum import tc_expansion3
    from fq_hat import tc_series
    from un_engine import get_nat
    K = 2 * g - 2; npol = K
    shown = [i for i in range(len(outer)) if maxout is None or outer[i]['E'] <= maxout]
    bad = {i: tab[i][:npol] for i in shown if i != 0 and i in tab and any(tab[i][:npol])}
    nE = sum(1 for i in shown if i != 0)
    first = {}
    for i, c in bad.items():
        n0 = next(n for n in range(npol) if c[n]); first.setdefault((outer[i]['E'], n0), []).append(i)
    if not bad: log("   RESULT: all %d outer graphs with edges%s are O(u^%d): REGULAR" % (nE, "" if maxout is None else " and <= %d edges" % maxout, npol))
    else:
        log("   RESULT: %d of %d outer graphs with edges%s are SINGULAR; (E, first non-zero u-order): graphs = %s" % (len(bad), nE, "" if maxout is None else " (<= %d edges)" % maxout, {k: v for k, v in sorted(first.items())}))
        for i in sorted(bad): log("      outer graph %3d (E=%d, gv=%s): [u^0..u^5] = %s" % (i, outer[i]['E'], outer[i]['gv'], [rec.fmt(x) for x in bad[i]]))
    v0 = tab.get(0, [0] * K)
    log("   vertex graph [u^0..u^5] = %s" % [rec.fmt(x) for x in v0])
    tcq = tc_series(get_nat(prime, g, info['frame'], K)['C'], R, K + 4)
    Vt = tc_expansion3(v0, tcq, R, g)[:npol]
    log("   vertex graph in the quintic t_c (frame %s): %s" % (info['frame'], " + ".join("(%s) t_c^%d" % (rec.fmt(Vt[k]), k - npol) for k in range(npol))))
    Ts = res['T'].scale(1 / info['lam']); T = Ser(R, Ts.v, [R.of(x) for x in Ts.c]).truncate(K + 4)
    V = tc_expansion3(v0, T, R, g)[:npol]
    okv = V[0] == R.of(Fr(1, 180000)) and not any(V[1:])
    log("   vertex graph in T/lambda (T = Int f, lambda = %s): %s  -> %s" % (info['lam'], " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)), "(1/180000) (T/lambda)^-6 + 0 + .. + 0" if okv else "NOT c_4 (T/lambda)^-6 + 0"))
    if maxout is None:
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(len(outer))) % R.p for k in range(K)]
        Tt = tc_expansion3(tot, T, R, g)[:npol]
        log("   total over all outer graphs, in T/lambda: %s" % " + ".join("(%s) T^%d" % (rec.fmt(Tt[k]), k - npol) for k in range(npol)))
    return bad, okv


def run(case, prime, maxout, nproc):
    from un_engine import run_pool
    from fq3_rec import Rec
    g = 4; K = 6
    log("=" * 110); log("un_e4 run: case %s, p = 2^%s-1, genus %d, maxout %s, %d worker processes   [%s]" % (case, prime, g, maxout, nproc, datetime.datetime.now().isoformat(timespec='seconds')))
    res, info, Qh, QL, nat = prep(case, prime, g)
    t0 = time.time(); spec = (prime, g, ("0", "0"), K, 1, QL, maxout)
    v, tab = run_pool(spec, nproc, log, sid="e4_%s_%s_%s" % (case, prime, maxout))
    wall = time.time() - t0; R = nat['R']; rec = Rec(R)
    outer = [dict(E=G['E'], gv=list(G['gv'])) for G in v.inner_graphs]
    log("   graph sum: %d inner graphs, wall %.0f s" % (len(v.inner_graphs), wall))
    json.dump(dict(case=case, prime=prime, maxout=maxout, info={k: (str(x) if isinstance(x, Fr) else x) for k, x in info.items()},
                   table={str(i): c for i, c in sorted(tab.items())}, outer=outer, wall=wall), open(tab_path(case, prime, maxout), "w"))
    report_table(tab, outer, g, maxout, rec, R, res, info, nat, prime)
    log("UN_E4 RUN %s %s" % (case, "checks passed" if not FAILED else "FAILURES: %s" % FAILED))


def load(case, prime, maxout):
    path = tab_path(case, prime, maxout)
    if not os.path.exists(path): return None
    d = json.load(open(path)); d['tab'] = {int(i): c for i, c in d['table'].items()}
    return d


def report(prime):
    from fq_fp import FpRing
    from fq3_rec import Rec
    from un_core import PRIMES
    R = FpRing(PRIMES[prime]); rec = Rec(R); p = R.p; K = 6
    log("=" * 110); log("un_e4 report, p = 2^%s-1   [%s]" % (prime, datetime.datetime.now().isoformat(timespec='seconds')))
    ref = json.load(open(os.path.join(PROJ, "g4_split_g4_%s_nat.json" % prime))); reft = {int(i): c for i, c in ref['table'].items()}
    for maxout in (None, 2):
        a = load("a", prime, maxout)
        if a is None: continue
        sel = [i for i in range(len(a['outer'])) if maxout is None or a['outer'][i]['E'] <= maxout]
        diff = [i for i in sel if [x % p for x in a['tab'].get(i, [0] * K)] != [x % p for x in reft.get(i, [0] * K)]]
        check("REGRESSION (a, maxout %s): all %d outer graphs in range equal g4_split_g4_%s_nat.json (u^0..u^5)" % (maxout, len(sel), prime), not diff, diff[:10])
        for x, y in (("b", "a"), ("d", "c1")):
            X = load(x, prime, maxout); Y = load(y, prime, maxout)
            if X is None or Y is None: continue
            P11 = X['info']['P11_stored']; m5 = Fr(P11['5']) / Fr(P11['0']); im5 = R.inv(R.of(m5))
            dd = {i: [(X['tab'].get(i, [0] * K)[k] - Y['tab'].get(i, [0] * K)[k]) % p for k in range(K)] for i in sel}
            dd = {i: c for i, c in dd.items() if any(c)}
            log("   %s - %s (maxout %s): %d outer graphs differ; difference / m5, m5 = %s (the rows differ by the scalar factor 1 + m5 u^4 zhat^5):" % (x, y, maxout, len(dd), m5))
            for i in sorted(dd): log("      outer graph %3d (E=%d, gv=%s): [u^0..u^5] = %s" % (i, a['outer'][i]['E'], a['outer'][i]['gv'], [rec.fmt(c * im5 % p) for c in dd[i]]))
            hi = [i for i in dd if a['outer'][i]['E'] > 2 or any(dd[i][:4]) or (a['outer'][i]['E'] == 2 and dd[i][4])]
            check("   %s - %s is confined to E <= 2, to [u^4], [u^5] for E <= 1 and to [u^5] for E = 2 (valuation argument)" % (x, y), not hi, hi[:10])
            if 0 in dd:
                from fq_hat import tc_series
                from fq3_msum import tc_expansion3
                from un_engine import get_nat
                tc = tc_series(get_nat(prime, 4, ("0", "0"), K)['C'], R, K + 4)
                Vd = tc_expansion3([c * im5 % p for c in dd[0]], tc, R, 4)[:K]
                log("      vertex graph difference / m5 in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(Vd[k]), k - K) for k in range(K)))
        for c in ("c1", "c2", "c3", "c4", "e"):
            X = load(c, prime, maxout)
            if X is None: continue
            sel = [i for i in range(len(X['outer'])) if maxout is None or X['outer'][i]['E'] <= maxout]
            bad = [i for i in sel if i != 0 and any(X['tab'].get(i, [0] * K))]
            check("   case %s (maxout %s): all %d outer graphs with edges in range are regular" % (c, maxout, len(sel) - 1), not bad, bad[:10])
    log("UN_E4 REPORT %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("mode"); ap.add_argument("cases", nargs="?", default="a")
    ap.add_argument("--prime", default="89"); ap.add_argument("--maxout", type=int, default=None); ap.add_argument("--nproc", type=int, default=5)
    a = ap.parse_args()
    if a.mode == "prep":
        log("=" * 110); log("un_e4 prep, p = 2^%s-1   [%s]" % (a.prime, datetime.datetime.now().isoformat(timespec='seconds')))
        for c in a.cases.split(","): prep(c, a.prime)
        log("UN_E4 PREP %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
    elif a.mode == "run": run(a.cases, a.prime, a.maxout, min(a.nproc, 5))
    elif a.mode == "report": report(a.prime)


if __name__ == "__main__":
    main()
