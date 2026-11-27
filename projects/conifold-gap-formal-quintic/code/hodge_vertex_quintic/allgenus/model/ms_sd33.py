# -*- coding: utf-8 -*-
"""ms_sd33.py -- integrity check of an amplitude table that was NOT produced in this session's runs (amp(3,3), finished by a process
of the first attempt): every (h, n) correlator with a tau_0(1) or tau_1(1) insertion is compared with the string / dilaton equation
applied to the (h, n-1) table (recomputed here), and every correlator is compared with its permutations.
  string:   <tau_0(1) X_1 .. X_m>_h = Sum_j <.. tau_(k_j - 1)(X_j) ..>_h
  dilaton:  <tau_1(1) X_1 .. X_m>_h = (2h - 2 + m) <X_1 .. X_m>_h
Run: py -I -B ms_sd33.py [h n]   (default 3 3)
"""
import itertools, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc
import ms_engine as me


def main(h, n):
    ns = nd = bad_s = bad_d = bad_p = 0
    for ins in mc.all_insertions(h, n):
        val = mc.corr(h, ins)
        for perm in set(itertools.permutations(ins)):
            if me.corr(h, tuple(perm)) != val: bad_p += 1
        if (0, 0) in ins:
            rest = list(ins); rest.remove((0, 0)); rhs = Fr(0)
            for j, (a, k) in enumerate(rest):
                if k > 0:
                    r2 = rest[:j] + [(a, k - 1)] + rest[j + 1:]
                    rhs += mc.corr(h, tuple(sorted(r2)))
            ns += 1; bad_s += (val != rhs)
        if (0, 1) in ins:
            rest = list(ins); rest.remove((0, 1))
            rhs = (2 * h - 2 + len(rest)) * mc.corr(h, tuple(sorted(rest)))
            nd += 1; bad_d += (val != rhs)
    print("(%d,%d): string equation on %d correlators, %d failures; dilaton on %d, %d failures; permutation failures %d"
          % (h, n, ns, bad_s, nd, bad_d, bad_p))
    print("PASS" if bad_s == bad_d == bad_p == 0 else "FAIL")


if __name__ == "__main__":
    h = int(sys.argv[1]) if len(sys.argv) > 1 else 3; n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    main(h, n)
