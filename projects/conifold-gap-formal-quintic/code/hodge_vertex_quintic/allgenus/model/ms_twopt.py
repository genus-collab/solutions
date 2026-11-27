# -*- coding: utf-8 -*-
"""ms_twopt.py -- task 2(b), two points.  For F_h(V(q)) with dV = PiC the two-point singular part is
      <X Y>_h = F_h'' PiC(X) PiC(Y) + F_h' Hess(X, Y)      (+ regular),
so  Hess_h(X, Y) := (<X Y>_h - F_h'' PiC(X) PiC(Y)) / F_h'  must be independent of h (h >= 2; h = 1 with F_1 = -(1/12) log v) on
every pair with negative v-exponent.  Prints, for every admissible pair (all psi-degrees), e and Hess_h for h = 1..hmax, and the
string-only prediction  Hess_str(X, Y) = [X = t0(1)] Pi_1(Y) + [Y = t0(1)] Pi_1(X),  Pi_1(t_k(H^b)) = PiC(t_(k-1)(H^b)).
Run: py -I -B ms_twopt.py [hmax] [--all]
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc
import ms_engine as me


def pi1(x):
    a, k = x
    return mc.pic((a, k - 1)) if k >= 1 else Fr(0)


def hess_str(x, y):
    return (pi1(y) if x == (0, 0) else 0) + (pi1(x) if y == (0, 0) else 0)


def hess_h(h, x, y):
    val = mc.corr(h, (x, y))
    return (val - mc.Fder(h, 2) * mc.pic(x) * mc.pic(y)) / mc.Fder(h, 1)


def run(hmax, show_all):
    pairs = set()
    for h in range(1, hmax + 1):
        if not os.path.exists(me._cache_path(h, 2)): continue
        pairs.update(mc.all_insertions(h, 2))
    pairs = sorted(pairs, key=lambda p: (p[0][1] + p[1][1], p))
    nsing = nbad_h = nbad_str = 0; rows = []
    for x, y in pairs:
        line = []; hs = {}
        for h in range(1, hmax + 1):
            if not os.path.exists(me._cache_path(h, 2)) or x[1] + y[1] > 3 * h - 1: continue
            e = mc.vexp(h, (x, y)); H = hess_h(h, x, y); hs[h] = (e, H)
        sing = [h for h, (e, H) in hs.items() if e < 0 and h >= 1]
        if not sing:
            if not show_all: continue
        st = hess_str(x, y)
        vals_sing = set(hs[h][1] for h in sing)
        ok_h = len(vals_sing) <= 1
        ok_s = all(hs[h][1] == st for h in sing)
        nsing += len(sing); nbad_h += (not ok_h); nbad_str += sum(1 for h in sing if hs[h][1] != st)
        txt = "  ".join("h%d:e=%d,H=%s" % (h, e, H) for h, (e, H) in sorted(hs.items()))
        print("%-22s str=%-10s %s %s%s" % (mc.fmt_ins((x, y)), st, txt, "" if ok_h else " [h-DEPENDENT]", "" if ok_s else " [!=str]"))
    print("singular (pair, genus) instances: %d; pairs whose Hess_h differs between singular genera: %d; singular instances != Hess_str: %d"
          % (nsing, nbad_h, nbad_str))


if __name__ == "__main__":
    hmax = int(sys.argv[1]) if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else 3
    run(hmax, "--all" in sys.argv)
