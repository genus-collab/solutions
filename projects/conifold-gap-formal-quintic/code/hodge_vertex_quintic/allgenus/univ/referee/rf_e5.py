# -*- coding: utf-8 -*-
"""rf_e5.py -- REFEREE (univ): the two genus-3 'surplus' directions w1, w2 of the worker (un_e5.log, flat coefficients p_{k,b,j}
copied from that log as DATA) and the controlled flat-irregular case C1c, rebuilt with the referee's own code.
D = 1 + X for w1, w2 (X = O(u^3) hatted, so D = 1 + X is exact and symplectic mod u^6); C1c: own random annihilating A plus
(3/7) u^2 zhat^3 h(3,2) + (2/5) u^1 zhat^2 h(2,2) (each irregular by one order), annihilation re-imposed (Cayley, rf_e3core).
Genus 3 full; genus 4 with max_outer 1 (exact for the vertex graph and the one-edge graphs): the [u^3] coefficients of the one-edge
graphs depend only on the u^3 part of X (valuation argument), so the worker's completion at u^4, u^5 is not needed.
Run: py -I -B rf_e5.py w1|w2|w12|C1c GENUS [max_outer|full]  -> rf_e5.log
"""
import os, sys, json, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_e4 import make_rows, fmtr, SCR, CG
from rf_e3 import yhat
from rf_e3core import set_orbit, checks, random_X, solve_ann, cayley, madd, ident
from rf_rows import dress, red, mu

W = {"w1": {(3, 0, 1): Fr(130, 27), (3, 1, 1): Fr(70, 27), (4, 0, 0): Fr(-260, 27), (4, 1, 0): Fr(-16, 9)},
     "w2": {(3, 0, 1): Fr(-460, 27), (3, 1, 1): Fr(-40, 27), (4, 0, 0): Fr(866, 27), (4, 1, 0): Fr(40, 9)}}
W["w12"] = {k: 5 * W["w1"].get(k, 0) + 2 * W["w2"].get(k, 0) for k in W["w1"]}


def X_of(coefs, K, p):
    X = {}
    for (k, b, j), c in coefs.items():
        e = mu(b, k) + j
        if e >= K: continue
        ser = [0] * K; ser[e] = red(c, p); set_orbit(X, k, b, ser, p)
    return X


def main():
    case, g = sys.argv[1], int(sys.argv[2]); mo = None if (len(sys.argv) <= 3 or sys.argv[3] == "full") else int(sys.argv[3])
    K = 2 * g - 2; KZ = 3 * g - 3; LOG = open(os.path.join(HERE, "rf_e5.log"), "a", encoding="ascii")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    from rf_indep import run, polar_in_T
    t0 = time.time(); base, T1q, p = make_rows("q", "canon", g, K); d = json.load(open(base)); yh = yhat(g, K, p)
    log("=" * 100); log("rf_e5: case %s, genus %d, max_outer %s, p = 2^89-1   [%s]" % (case, g, mo, time.strftime("%Y-%m-%dT%H:%M:%S")))
    if case.startswith("w"):
        X = X_of(W[case], K, p); D = madd(ident(K), X, p)
    else:
        X = random_X(random.Random(2026), KZ, K, p); X, D, it, skip = solve_ann(X, yh, KZ, p, K)
        for (k, b, j, c) in ((3, 2, 2, Fr(3, 7)), (2, 2, 1, Fr(2, 5))):
            cur = X.get(((b + k) % 5, b, k), [0] * K); add = [0] * K
            if j < K: add[j] = red(c, p)
            set_orbit(X, k, b, [(x + y) % p for x, y in zip(cur, add)], p)
        X, D, it, skip = solve_ann(X, yh, KZ, p, K)
    sd, irr, ann, d0, Tu = checks(D, yh, KZ, p, K)
    log("   D: symplectic %s; flat-irregular entries %s; annihilation violated at m = %s" % (not sd, sorted(irr), sorted(ann)))
    Q = dress(d['Q'], D, K, KZ, p); path = os.path.join(SCR, "rows_e5_%s_g%d.json" % (case, g))
    json.dump(dict(prime="89", p=p, zeta=d['zeta'], K=K, NZ=KZ, Q=Q), open(path, "w"))
    tag = "e5_%s_g%d_%s" % (case, g, "full" if mo is None else "E%d" % mo)
    tab, L, complete = run(path, g, mo, 1, 1600.0, os.path.join(SCR, "parts_%s.json" % tag), log)
    if not complete: log("   INCOMPLETE"); return
    shown = [i for i in range(len(L)) if mo is None or L[i]['E'] <= mo]; v0 = [i for i in shown if L[i]['E'] == 0][0]
    bad = {i: tab[i] for i in shown if L[i]['E'] > 0 and i in tab and any(tab[i])}
    log("   %d of %d outer graphs with edges%s SINGULAR" % (len(bad), len(shown) - 1, "" if mo is None else " (<= %d edges)" % mo))
    for i in sorted(bad)[:8]: log("      outer graph %d (E=%d, gv=%s): %s" % (i, L[i]['E'], list(L[i]['gv']), [fmtr(x, p) for x in bad[i]]))
    c = polar_in_T(tab.get(v0, [0] * K), Tu, g, p, K)
    log("   vertex graph in T_eff: %s" % " + ".join("(%s) T^%d" % (fmtr(c[j], p), j + 2 - 2 * g) for j in range(K)))
    log("   (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
