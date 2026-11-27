# -*- coding: utf-8 -*-
"""ms_ktest.py -- the Model Lemma on the subspace K (MODEL_SINGULAR.md s.3).

K = {q in H[[z]] : Pi_i(q) = 0 for i >= 1},  Pi_i(q) = (1/5)[z^i] eta(S(-z)H^2, q(z)) = the descendant time t^H_i of the ancestor vector q.
K is spanned by tau_0(H) (PiC = 1) and kappa(X, k) = [S(z) X z^k]_+ (X in {1, H^2, H^3, H^4}, k >= 0; PiC = 0), the ancestor images of
the descendant insertions tau_k(X).  Claim (ML): for q in K the singular part of Fbar_h(q) is F_h(v + PiC(q)), i.e.
   (ML1)  <H^n>_h = F_h^(n)(v)                                       (primary H-sector),
   (ML2)  every correlator of K-generators with at least one kappa has v-exponent >= 0 or vanishes.
This script evaluates every admissible correlator of K-generators (kappa(X, k) expanded with ms_model.kappa, all at v = 1): all of
them with Sum k <= 3h - 3 + n, and in addition every SINGULAR one with k up to 3h + 3 + n (kappa(X, k) has components of lower psi-degree,
so k beyond the dimension can contribute; by weight, e < 0 needs k < (10h - 11 + 5n)/4, inside this range), and reports the singular
ones (e < 0) that do not vanish.  Run: py -I -B ms_ktest.py
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc
import ms_model as mm
import ms_engine as me


def gens(kmax):
    """K-generators: ('H',) and ('k', b, k); each with its expansion [(coef, (a, k))] and its (charge, weight-class insertion)"""
    out = [(("H", 1, 0), [(Fr(1), (1, 0))])]
    for k in range(kmax + 1):
        for b in (0, 2, 3, 4):
            out.append((("k", b, k), [(c, ins) for c, ins in mm.kappa(b, k)]))
    return out


def corr_gen(h, glist):
    """<g_1 .. g_n>_h for generators given by their expansions"""
    tot = Fr(0)
    for combo in itertools.product(*[g[1] for g in glist]):
        c = Fr(1); ins = []
        for cf, x in combo: c *= cf; ins.append(x)
        if sum(k for a, k in ins) > 3 * h - 3 + len(ins): continue
        if mc.charge(ins) != 0: continue
        val = mc.corr(h, tuple(sorted(ins)))
        if val: tot += c * val
    return tot


def F0der(n):
    """n-th derivative (n >= 3) of F_0 = -(5/2) v^2 log v at v = 1 (F_0''' = -5/v)"""
    c = Fr(-5)
    for j in range(n - 3): c *= -(j + 1)
    return c


def run():
    plan = [(0, 3), (0, 4), (0, 5), (0, 6), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3),
            (4, 1), (4, 2)]
    if len(sys.argv) > 1: plan = [tuple(map(int, x.split(","))) for x in sys.argv[1:]]
    ntest = nsing = nbad = 0; ml1 = []; stats = {}
    for h, n in plan:
        if not os.path.exists(me._cache_path(h, n)): print("skip (%d,%d)" % (h, n)); continue
        D = 3 * h - 3 + n; G = gens(D + 6); c0 = (ntest, nsing, nbad, len(ml1))   # kappa(X, k) with k > D can be non-zero
        for combo in itertools.combinations_with_replacement(range(len(G)), n):
            gl = [G[i] for i in combo]
            lead = [(g[0][1], g[0][2]) for g in gl]                 # the tau_k(H^b) of the same weight and charge
            if mc.charge(lead) != 0: continue
            e = mc.vexp(h, tuple(lead))
            if e >= 0 and sum(k for a, k in lead) > D: continue           # beyond the old range only the singular ones are evaluated
            val = corr_gen(h, gl); ntest += 1
            nk = sum(1 for g in gl if g[0][0] == "k")
            if nk and val != 0:
                towers = "H2only" if all(g[0][0] == "H" or g[0][1] == 2 for g in gl) else "mixed"
                st = stats.setdefault((h, n, towers), {}); st[e] = st.get(e, 0) + 1
            if nk == 0:
                pred = mc.Fder(h, n) if h >= 1 else F0der(n)
                ml1.append((h, n, val, pred))
                continue
            if e < 0:
                nsing += 1
                if val != 0:
                    nbad += 1
                    print("   NONZERO singular: h=%d  %s  e=%d  value %s" % (h, " ".join(
                        ("H" if g[0][0] == "H" else "kappa(H%d,%d)" % (g[0][1], g[0][2])) for g in gl), e, val))
        print("(%d,%d): %d correlators of K-generators with >= 1 kappa, %d singular, %d singular non-zero"
              % (h, n, ntest - c0[0] - (len(ml1) - c0[3]), nsing - c0[1], nbad - c0[2]), flush=True)
    for key in sorted(stats):
        pass
    for key in sorted(stats):
        print("   non-zero correlators with >= 1 kappa, (h, n, towers) = %s: v-exponent histogram %s" % (key, dict(sorted(stats[key].items()))))
    for h, n, val, pred in ml1:
        print("(ML1) <H^%d>_%d = %s   F_h^(n)(1) = %s   %s" % (n, h, val, pred, "ok" if val == pred else "DIFFERENT"))
    print("(ML2) %d correlators of K-generators with >= 1 kappa; %d of them singular (e < 0); %d singular ones NON-ZERO"
          % (ntest - len(ml1), nsing, nbad))
    print("ML on K: %s" % ("CONSISTENT with all data" if nbad == 0 and all(v == p for h, n, v, p in ml1) else "VIOLATED"))


if __name__ == "__main__":
    run()
