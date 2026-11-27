# -*- coding: utf-8 -*-
"""un_e3_model.py -- Conjecture A far from the quintic: annihilating dressings of the PURE MODEL theory (no quintic gauge at all).

Rows: P = identity in ur_gauge's hatted form, Qhat = diag(tau^((1-c_b)/4)) Rhat(tau zhat) = the model x^4 + v/x at v(u) = u tau^(-5/4) in
the quintic's hatted frame (then E_R = E_M exactly and every edge graph vanishes; the vertex graph is the model potential, which by
homogeneity is c_g v^(2-2g)).  Vanishing covector: the pure-1/z solution of the PULLED-BACK model system bt = v'(u) bmod(v(u))
(ur_gauge.coefficients), y_0 = -e_1, y_1 = -v e_0/z, ...  Cases (genus 2 and 3):
  M0   D = 1 (control: edge graphs 0, vertex graph c_g T^(2-2g) with T = v)
  MA1, MA2   D = exp(X), X random flat-regular symplectic (seeds 7, 8; MA2 with k = 0 terms), column-1 orbits solved for (D Y^C)_+ = 0
  MB   MA1 + (1/7) u^2 zhat H(x)H   (minimal violation of annihilation)
Run: py -I -B un_e3_model.py [89|61]   (log un_e3_model_<p>.log)
"""
import os, sys, time, random, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress, symplectic_check, mu
from un_engine import Var, get_base, get_nat, run
from un_lin import yc_series
from un_ann import yhat_series, random_X, solve_ann, expX, delta, flat_irregular, lie_symplectic_defect, teff
from un_gauge import rows_from_gauge
from fq_fp import FpRing, QRing
from fq_series import Ser
from fq3_rec import Rec
from fq3_msum import tc_expansion3
CG = {2: Fr(1, 1200), 3: Fr(1, 25200)}


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"; p = PRIMES[prime]; R = FpRing(p); rec = Rec(R); Q = QRing()
    LOG = open(os.path.join(HERE, "un_e3_model_%s.log" % prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    import ur_gauge
    gmax = 3; KB = 2 * gmax - 2; KZB = 3 * gmax - 3; arB = Arith(p, KB, 1)
    bq, bt, tau, v, C = ur_gauge.coefficients(KB + 16)
    ys = yc_series(bt, KZB + 2, KB + 14); yhB = yhat_series(ys, R, KB)
    log("=" * 110); log("un_e3_model: dressings of the pure model, p = 2^%s-1   [%s]" % (prime, datetime.datetime.now().isoformat(timespec='seconds')))
    log("   model vanishing covector: y_1 = %s (= -v), y_2 = %s (= v^2/2)" % ([str(c) for c in ys[1][:4]], [str(c) for c in ys[2][:4]]))
    NU = KB + 12; ident = [[Ser.const(Q, 1 if (k == 0) else 0, NU) for b in range(5)] for k in range(KZB + 1)]
    res = dict(p=ident, tau=tau)
    for name in ("M0", "MA1", "MA2", "MB"):
        log("-" * 110); log("case %s" % name)
        if name == "M0": X = {}
        else:
            seed = 8 if name == "MA2" else 7
            X = random_X(arB, R, KZB, random.Random(seed), flat_regular=True, with_k0=(name == "MA2"))
            X, bad = solve_ann(arB, X, yhB, KZB, None); assert not bad
            if name == "MB": X[(2, 1, 1)] = arB.add(X.get((2, 1, 1)), arB.mono(R.of(Fr(1, 7)), 2))
        D = expX(arB, X, KZB); dl = {m: delta(arB, D, yhB, m, KZB) for m in range(1, KZB + 1)}
        Tu, d0 = teff(arB, D, yhB, KZB, R)
        log("   X: %d entries; symplectic %s; flat-regular %s; annihilation in the window %s; z^0 coeff of D Y^C: (%s) e_1"
            % (len(X), not lie_symplectic_defect(arB, X), not flat_irregular(arB, X), all(x is None for x in dl.values()),
               " + ".join("%s u^%d" % (rec.fmt(c), j) for j, c in enumerate(d0) if c)))
        for g in (2, 3):
            t0 = time.time(); K = 2 * g - 2; KZ = 3 * g - 3; ar = Arith(p, K, 1); nr = get_nat(prime, g, ("0", "0"), K)
            rows = lift_rows(ar, rows_from_gauge(res, nr['Rh'], R, KZ, K))
            Dg = {}
            for (a_, b_, k), s in D.items():
                if k > KZ: continue
                x = ar.from_u(arB.part(s, 0)[:K]) if s is not None else None
                if a_ == b_ and k == 0: x = ar.sub(x, ar.one())
                if x is not None: Dg[(a_, b_, k)] = x
            Qh = left_dress(ar, rows, Dg, KZ); sb = symplectic_check(ar, Qh, KZ, nr['eta_inv'])
            base = get_base(prime, g, ("0", "0"), K); vv = Var(base, Qh, 1); _, tab = run(None, vv); npol = K
            bad = {i: c[:npol] for i, c in tab.items() if i != 0 and any(c[:npol])}
            T = Ser(R, 1, list(Tu[:K])); V = tc_expansion3(tab.get(0, [0] * K)[:K], T, R, g)[:npol]
            okv = V[0] == R.of(CG[g]) and not any(V[1:])
            first = {}
            for i, c in bad.items(): first.setdefault((vv.inner_graphs[i]['E'], next(n for n in range(npol) if c[n])), []).append(i)
            log("   genus %d: rows symplectic %s; %s; vertex graph in T_eff: %s -> %s   (%.0f s)" % (g, not sb,
                "ALL %d edge graphs REGULAR" % (len(vv.inner_graphs) - 1) if not bad else "%d SINGULAR, (E, first order): %s" % (len(bad), {k: sorted(l) for k, l in sorted(first.items())}),
                " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)), "c_g T_eff^(2-2g)" if okv else "NOT c_g T_eff^(2-2g)", time.time() - t0))
    log("UN_E3_MODEL done"); LOG.close()


if __name__ == "__main__":
    main()
