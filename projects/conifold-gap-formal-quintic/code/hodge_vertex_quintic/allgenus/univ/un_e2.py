# -*- coding: utf-8 -*-
"""un_e2.py -- E2 (universality): the full NON-LINEAR outer graph sum at genus 2 and 3 for fictitious self-dual cyclic systems.

For b = (f, B/(f^2 g^2), f, g, g), f = t_c' F(u), g = -(varpi_0/(1-u)) Gm(u), B = -1/(u(1-u)^5) (un_gauge.system) the Theorem-U gauge P'
to the model is computed over Q, its rows Qhat' = Ghat' . Rhat(tau zhat) are formed (Rhat = the frozen symbol rows of the project) and
the engine un_engine.Var is run with these rows and the natural M = M_0(tau zhat).  Everything else (Psi, xhat, leaves) is read off
the rows:  xhat = 1/Qhat'_0^(0)[zhat^0],  leaves = unit row.  Reported per case:
   regularity of the gauge, hypothesis (G) (P' Y^C = P'_11(0,z) Y^C_system, z^-6..z^6, mod u^6), symplecticity of the rows,
   the polar coefficients [u^0..u^(2g-3)] of every outer graph with edges, and the vertex graph re-expanded in T/lambda,
   T = Int_0^u f, lambda = F(0)  (the system's own vanishing coordinate; (G) gives P' Y^C = -(1/lambda)(e_1 + T e_0/z + ..)).
Regressions: F = Gm = 1 reproduces the project rows in frame (0, 0) and in frame (1/3, -2/7) (the latter also re-run).
The normalisation P'_11(0, z) = 1 is irrelevant here (z^5 = u^4 zhat^5 is outside the genus-3 window); see un_e4.py.
Run: py -I -B un_e2.py [89|61|both]
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from un_core import Arith, lift_rows, symplectic_check
from un_engine import Var, get_base, get_nat, run
from un_gauge import pert_gauge, rows_from_gauge, check_G
from fq_series import Ser
from fq3_rec import Rec
from fq3_msum import tc_expansion3

F = Fr
CASES = [   # name, F, Gm, frame, regression against project rows?
    ("quintic, frame (0,0)  [regression]", (1,), (1,), ("0", "0"), True),
    ("quintic, frame (1/3,-2/7)  [regression]", (1,), (1,), ("1/3", "-2/7"), True),
    ("f -> f (1 + u/3)", (1, F(1, 3)), (1,), ("0", "0"), False),
    ("g -> g (1 - 3u/5)", (1,), (1, F(-3, 5)), ("0", "0"), False),
    ("f -> f (1 + u/3 - 2u^2/7 + 5u^3/11), g -> g (1 - 3u/5 + u^2/2 + 2u^3/3)", (1, F(1, 3), F(-2, 7), F(5, 11)), (1, F(-3, 5), F(1, 2), F(2, 3)), ("0", "0"), False),
    ("f -> f (1 - 7u/4 + 3u^2 - u^3/9), g -> g (1 + 11u/6 - 5u^2/3 + 4u^3/13)", (1, F(-7, 4), 3, F(-1, 9)), (1, F(11, 6), F(-5, 3), F(4, 13)), ("0", "0"), False),
    ("lambda = f(0) = 3/2: f -> f (3/2 + u/3), g -> g (1 - 3u/5 + u^2/2)", (F(3, 2), F(1, 3)), (1, F(-3, 5), F(1, 2)), ("0", "0"), False),
    ("perturbed, in frame (1/3,-2/7): f -> f (1 + 2u/5 - u^2/3), g -> g (1 + u/7 + 3u^2/4)", (1, F(2, 5), F(-1, 3)), (1, F(1, 7), F(3, 4)), ("1/3", "-2/7"), False),
]
CG = {2: Fr(1, 1200), 3: Fr(1, 25200)}


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    LOG = open(os.path.join(HERE, "un_e2.log"), "a", encoding="utf-8"); FAILED = []
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    gauges = {}
    log("=" * 110); log("un_e2: gauges over Q (k <= 12, working precision O(u^24))")
    for name, Fp, Gp, ab, regr in CASES:
        t0 = time.time(); res = pert_gauge(12, 24, Fp, Gp, ab); gauges[name] = res
        hol = all(x is None or x >= 0 for inf in res['info'] for x in inf['vals'])
        ok, c, bad, nchk, rng = check_G(res, 12, 6)
        log("-" * 110); log("%s   (%.1f s)" % (name, time.time() - t0))
        log("   T = Int f = %s;  P'_11(0, z) = %s" % (res['T'].fmt(5), " + ".join("(%s) z^%d" % (c[k], k) for k in sorted(c))))
        check("   Theorem U: every p_{k,b}, k <= 12, is a power series in u", hol)
        check("   (G): P' Y^C_model(v(u)) = P'_11(0, z) Y^C_system on z^%d..z^%d mod u^6 (%d coefficients of z, exact over Q)" % (rng[0], rng[1], nchk), ok, bad[:3])
    for prime in (("89", "61") if which == "both" else (which,)):
        log("=" * 110); log("un_e2: graph sums, p = 2^%s-1" % prime)
        for name, Fp, Gp, ab, regr in CASES:
            res = gauges[name]; lam = Fr(Fp[0]); log("-" * 110); log(name)
            for g, K in ((2, 2), (3, 4)):
                t0 = time.time(); base = get_base(prime, g, ("0", "0"), K); nr = get_nat(prime, g, ("0", "0"), K); R = nr['R']; p = nr['p']; rec = Rec(R)
                KZ = base.KZ; ar = Arith(p, K, 1)
                Qh = rows_from_gauge(res, nr['Rh'], R, KZ, K); QL = lift_rows(ar, Qh)
                sb = symplectic_check(ar, QL, KZ, nr['eta_inv'])
                check("   genus %d: the rows Ghat' Rhat(tau zhat) are symplectic mod (zhat^%d, u^%d)" % (g, KZ + 1, K), not sb, sb[:4])
                if regr:
                    ref = get_nat(prime, g, ab, K)['Qh']
                    same = all(Qh[a][i][n] == ref[a][i][n] for a in range(5) for i in range(5) for n in range(KZ + 1))
                    check("   genus %d REGRESSION: rows from the gauge = the project rows in frame %s (all 25 x %d x %d residues)" % (g, ab, KZ + 1, K), same)
                v = Var(base, QL, 1); _, tab = run(None, v); npol = 2 * g - 2
                bad = {i: [rec.fmt(x) for x in c[:npol]] for i, c in tab.items() if i != 0 and any(c[:npol])}
                nE = len(v.inner_graphs) - 1
                check("   genus %d: all %d outer graphs with edges are O(u^%d) (regular)" % (g, nE, npol), not bad, {i: (v.inner_graphs[i]['E'], bad[i]) for i in sorted(bad)[:6]})
                log("      vertex graph [u^0..]: %s" % [rec.fmt(x) for x in tab[0][:npol]])
                Ts = res['T'].scale(1 / lam); T = Ser(R, Ts.v, [R.of(x) for x in Ts.c]).truncate(K + 4)
                V = tc_expansion3(tab[0][:K], T, R, g)[:npol]
                log("      vertex graph in T/lambda: %s" % " + ".join("(%s) T^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)))
                check("   genus %d: vertex graph = (%s) (T/lambda)^-%d + 0 + .. + 0 (T = Int f, lambda = f(0)/t_c'(0) = %s)" % (g, CG[g], npol, lam), V[0] == R.of(CG[g]) and not any(V[1:]), [rec.fmt(x) for x in V])
                if not regr:
                    tcq = res['C'].tc.truncate(K + 4); tc = Ser(R, tcq.v, [R.of(x) for x in tcq.c]); V2 = tc_expansion3(tab[0][:K], tc, R, g)[:npol]
                    log("      control: the same vertex graph in the QUINTIC t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(V2[k]), k - npol) for k in range(npol)))
                log("      (%.0f s)" % (time.time() - t0))
    log("UN_E2 %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); LOG.close()


if __name__ == "__main__":
    main()
