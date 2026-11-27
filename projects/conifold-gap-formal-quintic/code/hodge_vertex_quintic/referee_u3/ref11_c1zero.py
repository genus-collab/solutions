# -*- coding: utf-8 -*-
"""ref11_c1zero.py -- claim C6, symbolically: the stratum c1 = 0 of {m_1 = 0, genus-2 rows} (variables t = c2_1, w = c3_2, high).
Run: py -I -B ref11_c1zero.py --prime 61|89"""
import sys, os, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_loc import *
VN = ['z', 'x', 'y', 't', 'w', 'c4_1', 'c4_2', 'c5_0', 'c5_1', 'c5_2', 'c6_1', 'c6_2']


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); a = ap.parse_args()
    S = load(a.prime); p = S.p; Y = Sys(S); P = load_polys(a.prime); FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""), flush=True)
        if not ok: FAIL.append(name)
    print("=" * 100); print("ref11_c1zero, p = 2^%s-1" % a.prime)
    Z = {k: msubs(f, {0: 0, 1: 0, 2: 0}, p) for k, f in P.items()}
    nzc = {k: sorted((tuple(VN[i] + ("^%d" % e[i] if e[i] > 1 else "") for i in range(12) if e[i])) for e in f) for k, f in Z.items() if f and k.startswith('c_')}
    print("   c1 = 0: non-zero cubics and their monomials: %s" % nzc)
    check("c1 = 0: the only non-zero cubic is c_9 = kappa t^3 (so c2_1 = 0, hence c2 = 0 by m_1 = 0)", list(nzc) == ['c_9'] and nzc['c_9'] == [('t^3',)])
    Z2 = {k: msubs(f, {3: 0}, p) for k, f in Z.items()}
    nzq = {k: sorted((tuple(VN[i] + ("^%d" % e[i] if e[i] > 1 else "") for i in range(12) if e[i])) for e in f) for k, f in Z2.items() if f and not k.startswith('J3')}
    print("   c1 = c2 = 0: non-zero quadrics / cubics and their monomials: %s" % nzq)
    check("c1 = c2 = 0: the only non-zero condition besides the genus-3 rows is q_7 = kappa' w^2 (so c3_2 = 0, hence c3 = 0 by the genus-2 rows)", list(nzq) == ['q_7'] and nzq['q_7'] == [('w^2',)])
    Z3 = {k: msubs(f, {4: 0}, p) for k, f in Z2.items()}
    check("low part zero: every quadric and cubic vanishes identically in (c4, c5, c6); only the two genus-3 rows remain (pure-high solutions = their kernel)",
          all(not f for k, f in Z3.items() if not k.startswith('J3')) and all(Z3[k] for k in ('J3_0', 'J3_1')))
    R = MRing(12, p); v = Y.low(R, R.zero, R.zero, R.zero, R.var(3), R.var(4))
    print("   the low coordinates on c1 = 0 as functions of (t, w): c2_2 %s, c3_0 %s, c3_1 %s" % tuple(sorted(VN[i] for i in msupport(v[j])) for j in (4, 5, 6)))
    print("REF11 failures:", FAIL)


if __name__ == "__main__":
    main()
