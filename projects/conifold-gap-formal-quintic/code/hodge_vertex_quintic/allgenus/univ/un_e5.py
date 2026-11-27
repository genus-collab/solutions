# -*- coding: utf-8 -*-
"""un_e5.py -- E5: are the two first-order regular directions of genus 3 that VIOLATE annihilation (un_e3_lin: ker Resp is 2
dimensions larger than ker Ann at genus 3, both extra directions at the top hatted order u^3) a window artefact?

For each extra direction w (flat coefficients p_{k,b,j}, representative from the reduced echelon form; hatted entry
x_{k,b} = Sum_j p_{k,b,j} u^(mu(b,k) + j), all at hatted order u^3):
  genus 3: D = exp(X), X = w  (X has valuation 3, so D = 1 + X mod u^4: the first-order statement, now as a finite dressing);
  genus 4: X completed by column-1 corrections at hatted orders u^4, u^5 only, so that every annihilation condition of the genus-4
           window holds EXCEPT the genus-3 top-order ones (hatted order u^3), which w violates by construction; then the genus-4
           graph sum (maxout 2: vertex graph and the outer graphs with <= 2 edges, exact for those).
Run: py -I -B un_e5.py [--nproc n]   (log un_e5.log)
"""
import os, sys, json, time, argparse, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress, symplectic_check, mu
from un_lin import variables, response_rows, annihilation_rows, yc_series, nullspace, rref
from un_ann import yhat_series, expX, delta, set_entry, lie_symplectic_defect, flat_irregular, teff


def extra_directions(p, R, ys):
    d = json.load(open(os.path.join(HERE, "un_resp_89.json")))
    g = 3; VARS = variables(g); nv = len(VARS)
    Rr, _ = response_rows(d['resp']['3'], d['outer']['3'], g, VARS, p); Ar, Al = annihilation_rows(g, VARS, ys, R)
    KR = nullspace(Rr, p, nv)
    vals = [[sum(r[c] * v[c] for c in range(nv)) % p for r in Ar] for v in KR]
    piv, M2 = rref([vals[i] + KR[i] for i in range(len(KR))], p, len(Ar) + nv)
    out = []
    for row in M2:
        a = row[:len(Ar)]; v = row[len(Ar):]
        if any(a): out.append(({Al[j]: a[j] for j in range(len(Ar)) if a[j]}, {VARS[c]: v[c] for c in range(nv) if v[c]}))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nproc", type=int, default=1); a = ap.parse_args()
    from fq_fp import FpRing
    from fq_series import Ser
    from fq3_rec import Rec
    from fq3_msum import tc_expansion3
    from un_engine import Var, get_base, get_nat, run, run_pool
    import ur_gauge
    LOGF = os.path.join(HERE, "un_e5.log")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True)
        with open(LOGF, "a", encoding="utf-8") as f: f.write(s + "\n")
    p = PRIMES["89"]; R = FpRing(p); rec = Rec(R)
    log("=" * 110); log("un_e5: the genus-3 top-order regular directions outside ker Ann, tested at genus 3 and 4   [%s]" % datetime.datetime.now().isoformat(timespec='seconds'))
    bq = ur_gauge.coefficients(30)[0]; ys = yc_series(bq, 14, 24)
    dirs = extra_directions(p, R, ys)
    for nd, (annv, w) in enumerate(dirs):
        log("-" * 110); log("direction %d: flat coefficients %s;  annihilation values (m, u^i): %s" % (nd + 1, {k: rec.fmt(v) for k, v in w.items()}, {k: rec.fmt(v) for k, v in annv.items()}))
        for g in (3, 4):
            K = 2 * g - 2; KZ = 3 * g - 3; ar = Arith(p, K, 1); yh = yhat_series(ys, R, K)
            X = {}
            for (k, b, j), c in w.items():
                cur = X.get(((b + k) % 5, b, k)); set_entry(ar, X, k, b, ar.add(cur, ar.mono(c, mu(b, k) + j)))
            if g == 4:                                   # column-1 corrections at hatted orders >= 4 only
                for it in range(12):
                    D = expX(ar, X, KZ); dl = {m: delta(ar, D, yh, m, KZ) for m in range(1, KZ + 1)}
                    hi = {m: [0] * 4 + ar.part(v, 0)[4:] for m, v in dl.items() if v is not None and any(ar.part(v, 0)[4:])}
                    if not hi: break
                    for m, c in hi.items():
                        if (2 - m) % 5 == 1 and m % 2 == 0: continue           # self-dual even column-1 entry: cannot be used
                        cur = X.get(((m + 1) % 5, 1, m)); set_entry(ar, X, m, 1, ar.add(cur, ar.from_u(c)))
            D = expX(ar, X, KZ); dl = {m: delta(ar, D, yh, m, KZ) for m in range(1, KZ + 1)}
            viol = {m: [k for k, c in enumerate(ar.part(v, 0)) if c] for m, v in dl.items() if v is not None}
            log("   genus %d: X symplectic %s, flat-irregular %s; annihilation violated at (m: hatted orders) %s" % (g, not lie_symplectic_defect(ar, X), flat_irregular(ar, X), viol))
            nr = get_nat("89", g, ("0", "0"), K); Dm1 = dict(D)
            for i in range(5): Dm1[(i, i, 0)] = ar.sub(Dm1.get((i, i, 0)), ar.one())
            Dm1 = {key: v for key, v in Dm1.items() if v is not None}
            Qh = left_dress(ar, lift_rows(ar, nr['Qh']), Dm1, KZ); maxout = None if g == 3 else 2; t0 = time.time()
            if g == 3 or a.nproc <= 1:
                base = get_base("89", g, ("0", "0"), K); v = Var(base, Qh, 1, maxout); _, tab = run(None, v)
            else:
                v, tab = run_pool(("89", g, ("0", "0"), K, 1, Qh, maxout), a.nproc, None, sid="e5_%d" % nd)
            npol = K; shown = [i for i, G in enumerate(v.inner_graphs) if maxout is None or G['E'] <= maxout]
            bad = {i: tab[i][:npol] for i in shown if i != 0 and i in tab and any(tab[i][:npol])}
            log("   genus %d%s: %s   (%.0f s)" % (g, "" if maxout is None else " (maxout %d)" % maxout, "ALL %d edge graphs REGULAR" % (len(shown) - 1) if not bad else "%d of %d edge graphs SINGULAR" % (len(bad), len(shown) - 1), time.time() - t0))
            for i in sorted(bad)[:6]: log("      outer graph %d (E=%d, gv=%s): %s" % (i, v.inner_graphs[i]['E'], v.inner_graphs[i]['gv'], [rec.fmt(x) for x in bad[i]]))
            Tu, d0 = teff(ar, D, yh, KZ, R); T = Ser(R, 1, list(Tu[:K])); V = tc_expansion3(tab.get(0, [0] * K)[:K], T, R, g)[:npol]
            log("      vertex graph in T_eff: %s" % " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)))
    log("UN_E5 done")


if __name__ == "__main__":
    main()
