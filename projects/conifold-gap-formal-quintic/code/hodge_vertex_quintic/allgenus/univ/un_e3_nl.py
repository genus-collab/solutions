# -*- coding: utf-8 -*-
"""un_e3_nl.py -- E3, NON-LINEAR: the outer graph sum with R' replaced by D R', D = exp(X) a finite left dressing (un_ann), M fixed.

Cases (X random, exact rational coefficients reduced mod p; one D per case, built at the precision of the highest genus requested and
restricted to the lower genera):
  A1  flat-regular, symplectic, (D Y^C_Q)_+ = 0 exactly in the window, k >= 1          -> hypothesis (A), not a gauge
  A2  as A1, another seed, with u-dependent diagonal k = 0 terms                       -> hypothesis (A), not a gauge
  B1  flat-regular, symplectic, column-1 entries random instead of solved              -> violates annihilation
  B2  A1 plus the single flat-regular term c u^2 zhat H(x)H                            -> minimal violation
  C1  symplectic, annihilating exactly, hatted-regular but NOT flat-regular            -> annihilation without flat regularity
Witness that A1, A2 are not Theorem-U gauges to any cyclic system: the z^0 coefficient delta_0(u) of the transported covector
D Y^C_Q would have to be constant (first equation of th Y = A' Y for a pure-1/z solution); it is printed.
Reported: polar coefficients of all outer graphs with edges; vertex graph re-expanded in T_eff = -[z^-1 e_0](D Y^C_Q).
Run: py -I -B un_e3_nl.py 89|61 [maxgenus=3] [cases=A1,A2,B1,B2,C1] [nproc]     (genus 4: pool of <= 4 processes, one case at a time)
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress, symplectic_check, mu
from un_engine import Var, get_base, get_nat, run, run_pool
from un_lin import yc_series
from un_ann import yhat_series, random_X, solve_ann, expX, delta, flat_irregular, lie_symplectic_defect, teff, set_entry
from fq_fp import FpRing
from fq_series import Ser
from fq3_rec import Rec
from fq3_msum import tc_expansion3
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}


def build_case(name, ar, R, KZ, yh, log):
    if name in ("A1", "B1", "B2"):
        X = random_X(ar, R, KZ, random.Random(20261113), flat_regular=True, with_k0=False)
        if name == "B1":
            X.update(random_X(ar, R, KZ, random.Random(77), flat_regular=True, with_k0=False, exclude_col1=False))
            X2 = random_X(ar, R, KZ, random.Random(20261113), flat_regular=True, with_k0=False)
            for key, v in X2.items(): X[key] = v                      # same non-column-1 entries as A1, random column-1 entries
            return X
        X, bad = solve_ann(ar, X, yh, KZ, log); assert not bad, bad
        if name == "B2":
            c = ar.mono(R.of(Fr(1, 7)), 2); X[(2, 1, 1)] = ar.add(X.get((2, 1, 1)), c)
        return X
    if name == "A2":
        X = random_X(ar, R, KZ, random.Random(4242), flat_regular=True, with_k0=True)
        X, bad = solve_ann(ar, X, yh, KZ, log); assert not bad, bad
        return X
    if name == "C1":
        X = random_X(ar, R, KZ, random.Random(999), flat_regular=False, with_k0=False)
        X, bad = solve_ann(ar, X, yh, KZ, log)
        return X
    raise ValueError(name)


def restrict(D, arS, KZS, arB):
    out = {}
    for (a, b, k), v in D.items():
        if k > KZS or v is None: continue
        s = arS.from_u(arB.part(v, 0)[:arS.KU])
        if s is not None: out[(a, b, k)] = s
    return out


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"; gmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    cases = sys.argv[3].split(",") if len(sys.argv) > 3 else ["A1", "A2", "B1", "B2", "C1"]
    nproc = min(4, int(sys.argv[4])) if len(sys.argv) > 4 else 4
    genera = [g for g in (2, 3, 4) if g <= gmax] if gmax < 4 else [4]
    LOG = open(os.path.join(HERE, "un_e3_nl_%s%s.log" % (prime, "_g4" if gmax == 4 else "")), "a", encoding="utf-8"); FAILED = []
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    p = PRIMES[prime]; R = FpRing(p); rec = Rec(R)
    KB = 2 * gmax - 2; KZB = 3 * gmax - 3; arB = Arith(p, KB, 1)
    import ur_gauge
    bq = ur_gauge.coefficients(KB + 16)[0]; ys = yc_series(bq, KZB + 2, KB + 14); yhB = yhat_series(ys, R, KB)
    log("=" * 110); log("un_e3_nl: p = 2^%s-1, D built mod (zhat^%d, u^%d), genera %s" % (prime, KZB + 1, KB, genera))
    for name in cases:
        log("-" * 110); log("case %s" % name)
        X = build_case(name, arB, R, KZB, yhB, log); D = expX(arB, X, KZB)
        irr = flat_irregular(arB, X); lsd = lie_symplectic_defect(arB, X)
        dl = {m: delta(arB, D, yhB, m, KZB) for m in range(1, KZB + 1)}
        viol = {m: arB.val(v) + mu(1, m) for m, v in dl.items() if v is not None and arB.val(v) + mu(1, m) < KB + mu(1, m)}
        Tu, d0 = teff(arB, D, yhB, KZB, R)
        log("   X: %d non-zero entries (zhat-orders %s); infinitesimally symplectic: %s; flat-regular: %s%s" % (len(X), sorted(set(k for (_, _, k) in X)), not lsd, not irr, "" if not irr else " (irregular entries (a,b,k): %s)" % sorted(irr)[:8]))
        log("   annihilation (D Y^C_Q)_+ = 0 in the window: %s%s" % (not viol, "" if not viol else "  violated at zhat^m, first hatted u-order: %s" % {m: arB.val(dl[m]) for m in viol}))
        log("   z^0 coefficient of D Y^C_Q: (%s) e_1   [constant <=> could be a gauge]" % " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(d0) if c))
        log("   T_eff = -[z^-1 e_0] D Y^C_Q = u (%s)" % " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(Tu) if c))
        for g in genera:
            t0 = time.time(); K = 2 * g - 2; base = get_base(prime, g, ("0", "0"), K); nr = get_nat(prime, g, ("0", "0"), K); KZ = base.KZ
            ar = Arith(p, K, 1); Dg = restrict(D, ar, KZ, arB)
            Dm1 = dict(Dg)
            for a in range(5): Dm1[(a, a, 0)] = ar.sub(Dm1.get((a, a, 0)), ar.one())
            Dm1 = {key: v for key, v in Dm1.items() if v is not None}
            Qh = left_dress(ar, lift_rows(ar, nr['Qh']), Dm1, KZ); sb = symplectic_check(ar, Qh, KZ, nr['eta_inv'])
            log("   genus %d: dressed rows symplectic mod (zhat^%d, u^%d): %s" % (g, KZ + 1, K, not sb))
            if sb: FAILED.append("%s genus %d rows not symplectic" % (name, g))
            if g == 4:
                spec = (prime, 4, ("0", "0"), K, 1, [[[q for q in col] for col in row] for row in Qh], None)
                v, tab = run_pool(spec, nproc, log, sid="%s_%s" % (name, prime))
            else:
                v = Var(base, Qh, 1); _, tab = run(None, v)
            npol = 2 * g - 2
            bad = {i: c[:npol] for i, c in tab.items() if i != 0 and any(c[:npol])}
            nE = len(v.inner_graphs) - 1
            if not bad: log("   genus %d: ALL %d outer graphs with edges are O(u^%d): REGULAR" % (g, nE, npol))
            else:
                first = {}
                for i, c in bad.items():
                    n0 = next(n for n in range(npol) if c[n]); first.setdefault((v.inner_graphs[i]['E'], n0), []).append(i)
                log("   genus %d: %d of %d outer graphs with edges are SINGULAR; (edges, first non-zero u-order): number of graphs = %s" % (g, len(bad), nE, {key: len(l) for key, l in sorted(first.items())}))
                for i in sorted(bad)[:4]: log("        e.g. outer graph %d (E=%d, gv=%s): %s" % (i, v.inner_graphs[i]['E'], v.inner_graphs[i]['gv'], [rec.fmt(x) for x in bad[i]]))
            T = Ser(R, 1, list(Tu[:K])); V = tc_expansion3(tab[0][:K], T, R, g)[:npol]
            okv = V[0] == R.of(CG[g]) and not any(V[1:])
            log("   genus %d: vertex graph in T_eff: %s   -> %s" % (g, " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)), "c_g T_eff^%d + 0" % (-npol) if okv else "NOT a pure power of T_eff"))
            log("   RESULT %s genus %d: edge graphs %s, vertex graph %s   (%.0f s)" % (name, g, "REGULAR" if not bad else "SINGULAR", "= c_g T_eff^(2-2g)" if okv else "not c_g T_eff^(2-2g)", time.time() - t0))
    log("UN_E3_NL done%s" % ("" if not FAILED else "; FAILURES: %s" % FAILED)); LOG.close()


if __name__ == "__main__":
    main()
