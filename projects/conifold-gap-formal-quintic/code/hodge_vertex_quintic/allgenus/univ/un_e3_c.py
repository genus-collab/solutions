# -*- coding: utf-8 -*-
"""un_e3_c.py -- E3, non-linear, second set: robustness of case A (other seeds) and controlled violations of FLAT REGULARITY that keep
the annihilation (D Y^C_Q)_+ = 0 (replaces case C1 of un_e3_nl.py, whose random flat-irregular X cannot be made annihilating inside
the zhat-window: delta_6 sits on the self-dual column-1 entry (6, 1), which must vanish, and the entries that would compensate it
have k >= 7, outside the genus-3 window; for flat-regular X they are O(u^4) and harmless).

Annihilation-neutral orbits (first order): an entry (k, b) enters delta_m with m = k - n0(b), n0(b) = (1 - b) mod 5, so the orbits
(1,0,2), (1,3,4), (2,2,4), (3,2,3) (and k = 0) never enter the positive part linearly.  Their flat regularity thresholds are
mu = 1, 1, 2, 3.  Cases (all built on A1 = the random annihilating flat-regular X of un_e3_nl.py, seed 20261113):
  A3, A4   A1-type with seeds 1, 2, X multiplied by 23/7 resp. -17/5          -> robustness of case A
  C1a      A1 + (1/3) u^0 zhat h(1,3)     (Pi-type H^4(x)H^3 + H^3(x)H^4; irregular by one order; X = O(u^0))
  C1b      A1 + (3/7) u^1 zhat^3 h(3,2)   (irregular by two orders)
  C1c      A1 + (3/7) u^2 zhat^3 h(3,2) + (2/5) u^1 zhat^2 h(2,2)   (irregular by one order each)
After adding the terms, solve_ann re-imposes (D Y^C_Q)_+ = 0 with column-1 corrections (non-linear terms only).
Run: py -I -B un_e3_c.py 89 [cases]
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress, symplectic_check, mu
from un_engine import Var, get_base, get_nat, run
from un_lin import yc_series
from un_ann import yhat_series, random_X, solve_ann, expX, delta, flat_irregular, lie_symplectic_defect, teff, set_entry
from un_e3_nl import restrict, CG
from fq_fp import FpRing
from fq_series import Ser
from fq3_rec import Rec
from fq3_msum import tc_expansion3


def random_big(ar, R, KZ, seed):
    rng = random.Random(seed)
    X = random_X(ar, R, KZ, rng, flat_regular=True, with_k0=(seed % 2 == 0))
    s = R.of({1: Fr(23, 7), 2: Fr(-17, 5)}[seed])                     # ONE common factor: keeps X infinitesimally symplectic
    return {key: ar.scale(v, s) for key, v in X.items()}


def build(name, ar, R, KZ, yh, log):
    if name in ("A3", "A4"):
        X = random_big(ar, R, KZ, {"A3": 1, "A4": 2}[name])
    else:
        X = random_X(ar, R, KZ, random.Random(20261113), flat_regular=True, with_k0=False)
        X, bad = solve_ann(ar, X, yh, KZ, log); assert not bad
        adds = {"C1a": [(1, 3, 0, Fr(1, 3))], "C1b": [(3, 2, 1, Fr(3, 7))], "C1c": [(3, 2, 2, Fr(3, 7)), (2, 2, 1, Fr(2, 5))]}[name]
        for k, b, j, c in adds:
            cur = X.get(((b + k) % 5, b, k))
            set_entry(ar, X, k, b, ar.add(cur, ar.mono(R.of(c), j)))
        X = {key: v for key, v in X.items() if v is not None}
    X, bad = solve_ann(ar, X, yh, KZ, log)
    return X, bad


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"; cases = sys.argv[2].split(",") if len(sys.argv) > 2 else ["A3", "A4", "C1a", "C1b", "C1c"]
    LOG = open(os.path.join(HERE, "un_e3_c_%s.log" % prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    p = PRIMES[prime]; R = FpRing(p); rec = Rec(R); gmax = 3; KB = 2 * gmax - 2; KZB = 3 * gmax - 3; arB = Arith(p, KB, 1)
    import ur_gauge, datetime
    bq = ur_gauge.coefficients(KB + 16)[0]; ys = yc_series(bq, KZB + 2, KB + 14); yhB = yhat_series(ys, R, KB)
    log("=" * 110); log("un_e3_c: p = 2^%s-1, D built mod (zhat^%d, u^%d)   [%s]" % (prime, KZB + 1, KB, datetime.datetime.now().isoformat(timespec='seconds')))
    for name in cases:
        log("-" * 110); log("case %s" % name)
        X, bad = build(name, arB, R, KZB, yhB, log); D = expX(arB, X, KZB)
        irr = flat_irregular(arB, X); lsd = lie_symplectic_defect(arB, X)
        dl = {m: delta(arB, D, yhB, m, KZB) for m in range(1, KZB + 1)}
        viol = {m: arB.val(v) for m, v in dl.items() if v is not None}
        Tu, d0 = teff(arB, D, yhB, KZB, R)
        log("   X: %d entries; symplectic: %s; flat-irregular entries (a, b, k): hatted valuation vs mu: %s" % (len(X), not lsd, {key: (arB.val(X[key]), mu(key[1], key[2])) for key in sorted(irr)}))
        log("   annihilation in the window: %s%s" % (not viol, "" if not viol else "  violated: %s" % viol))
        log("   z^0 coefficient of D Y^C_Q: (%s) e_1" % " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(d0) if c))
        for g in (2, 3):
            t0 = time.time(); K = 2 * g - 2; base = get_base(prime, g, ("0", "0"), K); nr = get_nat(prime, g, ("0", "0"), K); KZ = base.KZ
            ar = Arith(p, K, 1); Dg = restrict(D, ar, KZ, arB)
            for a in range(5): Dg[(a, a, 0)] = ar.sub(Dg.get((a, a, 0)), ar.one())
            Dg = {key: v for key, v in Dg.items() if v is not None}
            Qh = left_dress(ar, lift_rows(ar, nr['Qh']), Dg, KZ); sb = symplectic_check(ar, Qh, KZ, nr['eta_inv'])
            v = Var(base, Qh, 1); _, tab = run(None, v); npol = 2 * g - 2
            bad = {i: c[:npol] for i, c in tab.items() if i != 0 and any(c[:npol])}
            first = {}
            for i, c in bad.items(): first.setdefault((v.inner_graphs[i]['E'], next(n for n in range(npol) if c[n])), []).append(i)
            lam = Tu[0]; il = R.inv(lam); T = Ser(R, 1, [x * il % p for x in Tu[:K]])            # T_eff/lam = u + O(u^2)
            V = tc_expansion3(tab[0][:K], T, R, g)[:npol]
            V = [V[k] * pow(il, (k - npol) % (p - 1), p) % p for k in range(npol)]                  # coefficients of T_eff^(k - npol)
            okv = V[0] == R.of(CG[g]) and not any(V[1:])
            log("   genus %d: rows symplectic %s; %s; vertex graph in T_eff: %s %s   (%.0f s)" % (g, not sb, "ALL %d edge graphs REGULAR" % (len(v.inner_graphs) - 1) if not bad else "%d edge graphs SINGULAR, (E, first order): graphs %s" % (len(bad), {k: sorted(l) for k, l in sorted(first.items())}),
                                                                            " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)), "= c_g T_eff^(2-2g)" if okv else "NOT c_g T_eff^(2-2g)", time.time() - t0))
            for i in sorted(bad)[:3]: log("        e.g. outer graph %d (E=%d, gv=%s): %s" % (i, v.inner_graphs[i]['E'], v.inner_graphs[i]['gv'], [rec.fmt(x) for x in bad[i]]))
    LOG.close()


if __name__ == "__main__":
    main()
