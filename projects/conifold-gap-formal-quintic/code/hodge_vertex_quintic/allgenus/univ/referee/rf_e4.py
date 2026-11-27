# -*- coding: utf-8 -*-
"""rf_e4.py -- REFEREE (univ): E2/E4 re-derived by an independent route.  Rows from the referee's own PF solution, own regular gauge,
own symbol (rf_pf, rf_gauge, rf_rows); graph sums with the module-disjoint engine natural_form/genus4_referee/indep (rf_indep).
Cases: q = quintic (0,0); qab = quintic (1/3,-2/7); c1, c2 = worker's perturbed systems; r1, r2 = NEW systems (rf_g_run.CASES).
Normalisations: canon; shift = canon x (1 + z^5) (m = 1; at genus 4 only the z^5 term is visible).  Reports: regularity of every
outer graph with edges, vertex graph and total re-expanded in the system's own vanishing coordinate T'' = c0 Int f, and, when
the worker's table exists, a residue-by-residue comparison (e4_tab_<case>_89[_out2].json; shift: (d - c1)/m5, (b - a)/m5).
Run: py -I -B rf_e4.py CASE NORM GENUS [max_outer|full] [nproc]     -> rf_e4.log (appended)
"""
import os, sys, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SCR = os.environ.get("RF_SCRATCH", os.path.join(os.path.dirname(os.path.abspath(__file__)), "scratch"))
WORK = os.path.normpath(os.path.join(HERE, ".."))
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}
ALIAS = {"q": "quintic", "qab": "quintic_ab"}


def make_rows(case, norm, g, K, prime=None):
    prime = prime or os.environ.get("RF_PRIME", "89")
    from rf_g_run import CASES
    from rf_pf import system
    from rf_gauge import gauge, times_c, P11
    from rf_rows import hatted_gauge, rows, prime_zeta
    name = ALIAS.get(case, case); Fp, Gp, ab, _ = CASES[name]
    KZ = 3 * g - 3; S = system(30, Fp, Gp, *ab); pg, info = gauge(S['b'], S['bt'], KZ + 3, "canon")
    c0 = P11(pg)[0]
    if norm == "shift": pg = times_c(pg, {0: Fr(1), 5: Fr(1)})
    p, z = prime_zeta(prime); Gh = hatted_gauge(pg, S['tau'], K, KZ, p); Q = rows(Gh, S['tau'], K, KZ, p, z)
    T = S['f'].integ().scale(c0)                                        # T'' = c0 Int_0^u f
    T1 = [(T.coef(n + 1).numerator * pow(T.coef(n + 1).denominator, p - 2, p)) % p for n in range(K)]
    path = os.path.join(SCR, "rows_%s_%s_g%d_%s.json" % (case, norm, g, prime))
    json.dump(dict(prime=prime, p=p, zeta=z, K=K, NZ=KZ, Q=Q), open(path, "w"))
    return path, T1, p


def fmtr(x, p):
    sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "natural_form", "referee")))
    from ref_fp import fmt
    return fmt(x, p)


def main():
    case, norm, g = sys.argv[1], sys.argv[2], int(sys.argv[3])
    mo = None if (len(sys.argv) <= 4 or sys.argv[4] == "full") else int(sys.argv[4]); nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    K = 2 * g - 2; LOG = open(os.path.join(HERE, "rf_e4.log"), "a", encoding="ascii")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    from rf_indep import run, polar_in_T
    t0 = time.time(); path, T1, p = make_rows(case, norm, g, K)
    PR = os.environ.get("RF_PRIME", "89")
    log("=" * 100); log("rf_e4: case %s, normalisation %s, genus %d, max_outer %s, p = 2^%s-1   [%s]" % (case, norm, g, mo, PR, time.strftime("%Y-%m-%dT%H:%M:%S")))
    tag = "%s_%s_g%d_%s%s" % (case, norm, g, "full" if mo is None else "E%d" % mo, "" if PR == "89" else "_p" + PR)
    tab, L, complete = run(path, g, mo, nproc, 1600.0, os.path.join(SCR, "parts_%s.json" % tag), log)
    log("   rows built and graph sum %s in %.0f s (%d inner graphs)" % ("COMPLETE" if complete else "INCOMPLETE (rerun to resume)", time.time() - t0, len(L)))
    if not complete: return
    shown = [i for i in range(len(L)) if mo is None or L[i]['E'] <= mo]; v0 = [i for i in shown if L[i]['E'] == 0][0]
    bad = {i: tab[i] for i in shown if L[i]['E'] > 0 and i in tab and any(tab[i])}
    if not bad: log("   RESULT: all %d outer graphs with edges%s are O(u^%d): REGULAR" % (len(shown) - 1, "" if mo is None else " (<= %d edges)" % mo, K))
    else:
        log("   RESULT: %d of %d outer graphs with edges are SINGULAR" % (len(bad), len(shown) - 1))
        for i in sorted(bad): log("      outer graph %d (E=%d, gv=%s): %s" % (i, L[i]['E'], list(L[i]['gv']), [fmtr(x, p) for x in bad[i]]))
    V = tab.get(v0, [0] * K); log("   vertex graph [u^0..u^%d] = %s" % (K - 1, [fmtr(x, p) for x in V]))
    c = polar_in_T(V, T1, g, p, K); okv = c[0] == CG[g].numerator * pow(CG[g].denominator, p - 2, p) % p and not any(c[1:])
    log("   vertex graph in T'' = c0 Int f: %s  -> %s" % (" + ".join("(%s) T^%d" % (fmtr(c[j], p), j + 2 - 2 * g) for j in range(K)), "c_g T''^(2-2g) + 0" if okv else "NOT c_g T''^(2-2g) + 0"))
    if mo is None:
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(len(L))) % p for k in range(K)]
        ct = polar_in_T(tot, T1, g, p, K)
        log("   total in T'': %s" % " + ".join("(%s) T^%d" % (fmtr(ct[j], p), j + 2 - 2 * g) for j in range(K)))
    if g == 4:
        wcase = {"q": "a", "c1": "c1", "c2": "c2"}.get(case)
        wpath = os.path.join(WORK, "e4_tab_%s_89%s.json" % (wcase, "" if mo is None else "_out%d" % mo)) if wcase else None
        if norm == "canon" and wpath and os.path.exists(wpath) and PR == "89":
            W = {int(k): v for k, v in json.load(open(wpath))['table'].items()}
            diff = [i for i in shown if [x % p for x in tab.get(i, [0] * K)] != [x % p for x in W.get(i, [0] * K)]]
            log("   %s worker's %s: %d of %d outer rows differ %s" % ("PASS" if not diff else "FAIL", os.path.basename(wpath), len(diff), len(shown), diff[:8]))
        if norm == "shift":
            bp = os.path.join(SCR, "tab_%s.json" % tag.replace("_shift_", "_canon_"))
            base = json.load(open(bp)) if os.path.exists(bp) else None
            if base:
                B = {int(k): v for k, v in base.items()}
                for i in shown:
                    dlt = [(tab.get(i, [0] * K)[k] - B.get(i, [0] * K)[k]) % p for k in range(K)]
                    if any(dlt): log("      (shift - canon) per unit m, outer graph %d (E=%d, gv=%s): %s" % (i, L[i]['E'], list(L[i]['gv']), [fmtr(x, p) for x in dlt]))
    json.dump({str(i): c for i, c in tab.items()}, open(os.path.join(SCR, "tab_%s.json" % tag), "w"))
    log("   (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
