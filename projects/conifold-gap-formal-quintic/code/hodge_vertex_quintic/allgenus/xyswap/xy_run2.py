# -*- coding: utf-8 -*-
"""xy_run2.py -- closed formula of arXiv:2212.00320 Thm 7.1 (xy_jets, exact over Q at rational points) against own direct TR
(xy_tr) for n >= 2.  Usage: py -I -B xy_run2.py [qz|fp89|fp61] gmax
"""
import os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from xy_field import QZ, make_fp, P89, P61
from xy_tr import TR
from xy_jets import closed_value


def main():
    fld = sys.argv[1] if len(sys.argv) > 1 else "qz"
    gmax = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    F = QZ if fld == "qz" else make_fp(P89 if fld == "fp89" else P61)
    order = [(0, 3), (1, 1)]
    if gmax >= 1: order += [(0, 4), (1, 2), (2, 1)]
    if gmax >= 2: order += [(0, 5), (1, 3), (2, 2)]
    tr = TR(F, 2 if gmax >= 1 else 1, max(n for (_, n) in order)); t0 = time.time()
    for (g, n) in order: tr.step(g, n)
    print("field", F.name, "TR done %.1fs" % (time.time() - t0))
    rng = random.Random(20261113); v = 4
    tests = [(0, 3), (1, 2), (0, 4)] + ([(2, 2), (1, 3)] if gmax >= 2 else [])
    for (g, n) in tests:
        ok = True; npts = 3 if n <= 3 else 2
        for _ in range(npts):
            b = tuple(Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(n))
            while len(set(b)) < n or any(x == 0 or x ** 5 == 1 for x in b):
                b = tuple(Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(n))
            cl = closed_value(g, b, v) * (-1) ** n
            t = tr.value(g, n, b)
            same = (t - F.frac(cl)).is_zero()
            ok = ok and same
            extra = ""
            if fld == "qz": extra = " TR value rational: %s" % (t.rational() is not None)
            print("  (%d,%d) at %s: closed(EO) = %s ; equal %s%s" % (g, n, tuple(str(x) for x in b), cl, same, extra))
        print("CHECK omega_{%d,%d}: closed formula (Thm 7.1, EO sign) = direct TR: %s  (%.1fs)" % (g, n, ok, time.time() - t0))


if __name__ == "__main__":
    main()
