# -*- coding: utf-8 -*-
"""un_e3_g4.py -- E3 at GENUS 4: a random annihilating, flat-regular, symplectic, NON-gauge left dressing D = exp(X) of R' (case A at
genus 4), M = M_0(tau zhat) fixed, frame (0, 0), p = 2^89-1.

Construction: X = un_ann.random_X (flat-regular, graded, infinitesimally symplectic, every hatted entry O(u), column-1 orbits left
out), then the column-1 orbits are solved for delta_m(exp X) = 0 mod u^6, m = 1..9, EXCEPT m = 6: the column-1 entry of delta_6 is the
self-dual even entry (6, 1), which must vanish, and every other entry of delta_6 inside the window belongs to a column-1 orbit already
fixed by delta_7..delta_9.  delta_6 is then CHECKED (it vanishes for every seed tried: symplecticity makes the system consistent).
Witness of 'not a gauge': the z^0 coefficient delta_0(u) e_1 of D Y^C_Q is not constant.
Reported: polar coefficients of every outer graph with edges (all of them, or those with <= maxout edges), vertex graph in
T_eff = -[z^-1 e_0] D Y^C_Q.  Table -> e3g4_tab_<seed>_<prime>[_out<m>].json, log -> un_e3_g4.log.
Run: py -I -B un_e3_g4.py [--seed 20261113] [--maxout m] [--nproc n]   (nproc 1: single process, no pool)
"""
import os, sys, json, time, random, argparse, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress, symplectic_check, mu
from un_lin import yc_series
from un_ann import yhat_series, random_X, expX, delta, set_entry, flat_irregular, lie_symplectic_defect, teff


def solve_ann4(ar, X, yh, KZ, skip):
    K = ar.KU; X = dict(X)
    for it in range(4 * K + 4):
        D = expX(ar, X, KZ); dl = {m: delta(ar, D, yh, m, KZ) for m in range(1, KZ + 1)}
        if all(v is None for m, v in dl.items() if m not in skip): break
        dX = {}
        for m in range(KZ, 0, -1):
            if m in skip: continue
            acc = dl[m]
            for n in range(1, len(yh)):
                k = m + n
                if k > KZ: break
                s = dX.get(((m + 1) % 5, (1 - n) % 5, k))
                if s is not None: acc = ar.add(acc, ar.mul(s, ar.from_u(yh[n])))
            if acc is not None: set_entry(ar, dX, m, 1, acc)
        for key, v in dX.items(): X[key] = ar.add(X.get(key), v)
        X = {key: v for key, v in X.items() if v is not None}
    else:
        raise AssertionError("no convergence")
    return X, it


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=20261113); ap.add_argument("--maxout", type=int, default=None)
    ap.add_argument("--nproc", type=int, default=1); ap.add_argument("--prime", default="89"); a = ap.parse_args()
    from fq_fp import FpRing
    from fq_series import Ser
    from fq3_rec import Rec
    from fq3_msum import tc_expansion3
    from un_engine import Var, get_base, get_nat, run, run_pool
    import ur_gauge
    LOGF = os.path.join(HERE, "un_e3_g4.log")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True)
        with open(LOGF, "a", encoding="utf-8") as f: f.write(s + "\n")
    p = PRIMES[a.prime]; R = FpRing(p); rec = Rec(R); g = 4; K = 2 * g - 2; KZ = 3 * g - 3; ar = Arith(p, K, 1)
    log("=" * 110); log("un_e3_g4: case A at genus 4, seed %d, maxout %s, p = 2^%s-1   [%s]" % (a.seed, a.maxout, a.prime, datetime.datetime.now().isoformat(timespec='seconds')))
    bq = ur_gauge.coefficients(K + 16)[0]; ys = yc_series(bq, KZ + 2, K + 14); yh = yhat_series(ys, R, K)
    skip = [m for m in range(1, KZ + 1) if (2 - m) % 5 == 1 and m % 2 == 0]
    X = random_X(ar, R, KZ, random.Random(a.seed), flat_regular=True, with_k0=(a.seed % 2 == 0))
    X, it = solve_ann4(ar, X, yh, KZ, skip); D = expX(ar, X, KZ)
    dl = {m: delta(ar, D, yh, m, KZ) for m in range(1, KZ + 1)}
    irr = flat_irregular(ar, X); lsd = lie_symplectic_defect(ar, X); Tu, d0 = teff(ar, D, yh, KZ, R)
    log("   X: %d entries, zhat-orders %s; infinitesimally symplectic: %s; flat-regular: %s; solved in %d iterations (skipped m = %s)" % (len(X), sorted(set(k for (_, _, k) in X)), not lsd, not irr, it, skip))
    ok_ann = all(v is None for v in dl.values())
    log("   annihilation (D Y^C_Q)_+ = 0 mod u^%d for m = 1..%d, INCLUDING the unsolved m = %s: %s" % (K, KZ, skip, ok_ann))
    log("   z^0 coefficient of D Y^C_Q: (%s) e_1   [not constant: D is not a gauge between cyclic systems]" % " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(d0) if c))
    log("   T_eff = u (%s)" % " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(Tu) if c))
    assert ok_ann and not irr and not lsd
    nr = get_nat(a.prime, g, ("0", "0"), K); Dm1 = dict(D)
    for i in range(5): Dm1[(i, i, 0)] = ar.sub(Dm1.get((i, i, 0)), ar.one())
    Dm1 = {key: v for key, v in Dm1.items() if v is not None}
    Qh = left_dress(ar, lift_rows(ar, nr['Qh']), Dm1, KZ); sb = symplectic_check(ar, Qh, KZ, nr['eta_inv'])
    log("   dressed rows symplectic mod (zhat^%d, u^%d): %s" % (KZ + 1, K, not sb))
    t0 = time.time(); spec = (a.prime, g, ("0", "0"), K, 1, Qh, a.maxout)
    if a.nproc <= 1:
        base = get_base(a.prime, g, ("0", "0"), K); v = Var(base, Qh, 1, a.maxout); _, tab = run(None, v)
    else:
        v, tab = run_pool(spec, a.nproc, log, sid="e3g4_%d" % a.seed)
    outer = [dict(E=G['E'], gv=list(G['gv'])) for G in v.inner_graphs]; npol = K
    json.dump(dict(seed=a.seed, maxout=a.maxout, table={str(i): c for i, c in sorted(tab.items())}, outer=outer, wall=time.time() - t0),
              open(os.path.join(HERE, "e3g4_tab_%d_%s%s.json" % (a.seed, a.prime, "" if a.maxout is None else "_out%d" % a.maxout)), "w"))
    shown = [i for i in range(len(outer)) if a.maxout is None or outer[i]['E'] <= a.maxout]
    bad = {i: tab[i][:npol] for i in shown if i != 0 and i in tab and any(tab[i][:npol])}
    if not bad: log("   RESULT: all %d outer graphs with edges%s are O(u^6): REGULAR   (%.0f s)" % (len(shown) - 1, "" if a.maxout is None else " (<= %d edges)" % a.maxout, time.time() - t0))
    else:
        first = {}
        for i, c in bad.items(): first.setdefault((outer[i]['E'], next(n for n in range(npol) if c[n])), []).append(i)
        log("   RESULT: %d of %d outer graphs with edges are SINGULAR: (E, first order): %s   (%.0f s)" % (len(bad), len(shown) - 1, {k: sorted(l) for k, l in sorted(first.items())}, time.time() - t0))
        for i in sorted(bad)[:8]: log("      outer graph %d (E=%d, gv=%s): %s" % (i, outer[i]['E'], outer[i]['gv'], [rec.fmt(x) for x in bad[i]]))
    T = Ser(R, 1, list(Tu[:K])); V = tc_expansion3(tab.get(0, [0] * K)[:K], T, R, g)[:npol]
    okv = V[0] == R.of(Fr(1, 180000)) and not any(V[1:])
    log("   vertex graph in T_eff: %s  -> %s" % (" + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)), "c_4 T_eff^-6 + 0 + .. + 0" if okv else "NOT c_4 T_eff^-6 + 0"))
    log("UN_E3_G4 done")


if __name__ == "__main__":
    main()
