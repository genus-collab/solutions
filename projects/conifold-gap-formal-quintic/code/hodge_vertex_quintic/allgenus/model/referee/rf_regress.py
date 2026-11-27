# -*- coding: utf-8 -*-
"""rf_regress.py -- fix the sign convention of the referee engine by the numbers ON RECORD (not by the worker's tables):
F_2 = 1/1200, F_3 = 1/25200 (NATURAL_FORM.md s.0), <tau_1(1)>_1 = 5/24, <H>_1 = -1/12 (UNIFORM_REDUCTION s.6.3), <H,H,H>_0 = -5
(F_0 = -(5/2) v^2 log v), and the primary H-sector d_v^n F_h.  Both R(z) (sign +1) and R(-z) (sign -1) are run; exactly one must pass.
Run: py -I -B rf_regress.py 61|89 [sign]
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as rc
import rf_gsum as gsm


def record(h, n):
    """d_v^n F_h at v = 1 (record values: c_2 = 1/1200, c_3 = 1/25200; F_1 = -(1/12) log v; F_0''' = -5/v)"""
    if h == 0:
        c = Fr(-5)
        for j in range(n - 3): c *= -(j + 1)
        return c
    if h == 1:
        c = Fr(-1, 12)
        for j in range(1, n): c *= -j
        return c
    c = {2: Fr(1, 1200), 3: Fr(1, 25200)}[h]; e = 2 - 2 * h
    for j in range(n): c *= (e - j)
    return c


def main():
    p = rc.P61 if sys.argv[1] == "61" else rc.P89
    signs = [int(sys.argv[2])] if len(sys.argv) > 2 else [1, -1]
    for sg in signs:
        m = rc.Model(p, 14, sign=sg); t0 = time.time(); fails = 0; checks = []
        checks.append(("F_2", gsm.corr(m, 2, ()), Fr(1, 1200)))
        checks.append(("F_3", gsm.corr(m, 3, ()), Fr(1, 25200)))
        checks.append(("<tau_1(1)>_1", gsm.corr(m, 1, ((0, 1),)), Fr(5, 24)))
        for (h, n) in [(0, 3), (0, 4), (1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (3, 1)]:
            checks.append(("<H^%d>_%d" % (n, h), gsm.corr(m, h, ((1, 0),) * n), record(h, n)))
        for name, val, ref in checks:
            ok = val == m.fr(ref); fails += (not ok)
            print("   sign %+d  %-14s %s (record %s)" % (sg, name, "PASS" if ok else "FAIL", ref))
        print("SIGN %+d: %s  (%.0f s)" % (sg, "ALL PASS" if fails == 0 else "%d FAIL" % fails, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
