# -*- coding: utf-8 -*-
"""rf_ml.py -- independent test of the Model Lemma (ML) on K with the referee engine and the referee's own K-generators.

(ML1) <H^n>_h = d_v^n F_h at v = 1.   (ML2) every correlator of K-generators (tau_0(H), kappa(X,k) = [S(z) X z^k]_+, X != H) with at least
one kappa and v-exponent e < 0 vanishes.  A generator tuple is tested if its charge is 0 and e < 0 (all k with e < 0 are enumerated:
e increases with k).  Zero tests at BOTH primes; ML1 exact.  CONTROL (must fail): the same test with the calibration of the wrong sign,
kappa^-(X,k) = Sum_j (-1)^j S_j X z^(k-j), and with bare basis insertions tau_k(X) in place of kappa(X,k).
Run: py -I -B rf_ml.py h,n [h,n ...]
"""
import itertools, os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat as rf
import rf_tables as rt
ONE = "--one" in sys.argv


def gens(kmax, variant):
    out = [("H", ((Fr(1), (1, 0)),), (1, 0))]
    for k in range(kmax + 1):
        for b in (0, 2, 3, 4):
            if variant == "kappa": ex = tuple(rf.kappa(b, k))
            elif variant == "wrongsign": ex = tuple(((-1) ** j * c, ins) for j, (c, ins) in enumerate(rf.kappa(b, k)))
            else: ex = ((Fr(1), (b, k)),)
            out.append((("k", b, k), ex, (b, k)))
    return out


def value(tabs, h, gl):
    """multilinear expansion with pruning by the total psi-degree (<= 3h - 3 + n)"""
    D = 3 * h - 3 + len(gl); vals = []
    for m, tab in tabs:
        comps = [[(m.fr(cf), ins) for cf, ins in g[1]] for g in gl]
        tot = 0; stack = [(0, 1, (), 0)]
        while stack:
            i, c, ins, deg = stack.pop()
            if i == len(comps):
                x = tab.get(ins)
                if x: tot = (tot + c * x) % m.p
                continue
            for cf, x in comps[i]:
                if deg + x[1] <= D: stack.append((i + 1, c * cf % m.p, ins + (x,), deg + x[1]))
        vals.append(tot)
    return vals


def singular_tuples(G, h, n):
    """multisets of n generators (indices into G sorted by weight, descending) with v-exponent e < 0, by pruned recursion"""
    order = sorted(range(len(G)), key=lambda i: -rf.weight(*G[i][2]))
    w = [rf.weight(*G[i][2]) for i in order]; need = 2 - 2 * h          # e < 0  <=>  Sum w > 2 - 2h
    out = []
    def rec(start, chosen, sw):
        r = n - len(chosen)
        if r == 0:
            if sw > need: out.append(tuple(order[i] for i in chosen))
            return
        for i in range(start, len(order)):
            if sw + r * w[i] <= need: break
            rec(i, chosen + [i], sw + w[i])
    rec(0, [], Fr(0))
    return out


def run(h, n, variant):
    m1, m2 = rf.models()
    tabs = [(m1, rt.flat_table(m1, h, n))] + ([] if ONE else [(m2, rt.flat_table(m2, h, n))])
    G = gens(3 * h + n + 8, variant)
    ntest = nbad = 0; ml1 = []; ex = []
    for combo in singular_tuples(G, h, n):
        gl = [G[i] for i in combo]; lead = [g[2] for g in gl]
        if not rf.charge_ok(lead): continue
        e = rf.vexp(h, lead); assert e < 0
        nk = sum(1 for g in gl if g[0] != "H")
        if nk == 0:
            if h == 0 and n < 3: continue
            ml1.append(((rt.exact_val(h, tuple(lead)) if not ONE else (m1.fr(rt.Fder(h, n)) == tabs[0][1].get(tuple(lead), 0) and rt.Fder(h, n))), rt.Fder(h, n))); continue
        vals = value(tabs, h, gl); ntest += 1
        if any(vals):
            nbad += 1
            if len(ex) < 3: ex.append(" ".join("H" if g[0] == "H" else "k(H%d,%d)" % (g[0][1], g[0][2]) for g in gl) + " e=%s" % e)
    if not ml1 and not (h == 0 and n < 3) and not ONE: ml1.append((rt.exact_val(h, ((1, 0),) * n), rt.Fder(h, n)))
    return ntest, nbad, ml1, ex


def main():
    t0 = time.time(); T = {"kappa": [0, 0], "wrongsign": [0, 0], "bare": [0, 0]}
    for arg in [a for a in sys.argv[1:] if a != "--one"]:
        h, n = map(int, arg.split(","))
        for variant in ("kappa", "wrongsign", "bare"):
            nt, nb, ml1, ex = run(h, n, variant)
            T[variant][0] += nt; T[variant][1] += nb
            line = "(%d,%d) %-9s: %5d singular generator correlators with >= 1 kappa, %5d NON-ZERO" % (h, n, variant, nt, nb)
            if variant == "kappa":
                okml1 = all(a == b for a, b in ml1)
                line += ";  ML1 <H^%d>_%d: %s" % (n, h, ("ok " + str(ml1[0][0])) if ml1 and okml1 else ("FAIL %s" % ml1 if ml1 else "n/a"))
            if ex and variant != "kappa": line += "   e.g. " + "; ".join(ex[:2])
            if ex and variant == "kappa": line += "   VIOLATIONS e.g. " + "; ".join(ex)
            print(line, flush=True)
    for k, (nt, nb) in T.items(): print("TOTAL %-9s: %d tested, %d non-zero" % (k, nt, nb))
    print("rf_ml done (%.0f s)%s" % (time.time() - t0, "  [zero tests at 2^61-1 only]" if ONE else ""))


if __name__ == "__main__":
    main()
