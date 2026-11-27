# -*- coding: utf-8 -*-
"""lg_wkb_check.py -- the identity row of the conifold symbol to HIGHER order, straight from the WKB recursion of the formal-quintic
Picard-Fuchs operator (fq_rows.WKB, read-only import), against the closed recursion of the LG model (lg_gm.symbol_rows).

psi_{i,n} = L phi_n(L), phi_n a Laurent polynomial in L; the symbol of the identity row is f_n = [L^{4n}] phi_n at colour 0 (top degree),
and zeta^{-i n} f_n at colour i.  Statement tested for n <= Nz (default 14; the stored symbol stops at n = 6):
     (i)  deg_L phi_n = 4n exactly (the degree bound of LG_MODEL.md s.4 with non-vanishing top coefficient),
     (ii) [L^{4n}] phi_n = zeta^{-i n} g^0_n   with g^0 = J_0(-4 zhat) from the cyclic recursion (MC).
Run:  py -I lg_wkb_check.py [Nz]
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import WKB
from lg_gm import symbol_rows


def main():
    Nz = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    rows = symbol_rows(Nz); allok = True
    for tag, p in (("89", P2), ("61", P1)):
        R = FpRing(p)
        for i in (0, 1):
            t0 = time.time(); w = WKB(R, i, Nz)
            assert w.failed is None
            degs = [max(w.phi[n]) for n in range(Nz + 1)]
            top = [w.phi[n][4 * n] % p for n in range(Nz + 1)]
            pred = [R.of(rows[0][n]) * R.pow(R.zeta, (-i * n) % 5) % p for n in range(Nz + 1)]
            ok_deg = degs == [4 * n for n in range(Nz + 1)]
            ok_top = top == pred
            first_bad = next((n for n in range(Nz + 1) if top[n] != pred[n]), None)
            print("%s p = 2^%s-1, colour %d, Nz = %d: deg_L phi_n = 4n: %s;  [L^4n] phi_n = zeta^(-i n) g^0_n for all n <= %d: %s%s   (%.1f s)"
                  % ("PASS" if ok_deg and ok_top else "FAIL", tag, i, Nz, ok_deg, Nz, ok_top, "" if ok_top else " first failure n = %s" % first_bad, time.time() - t0), flush=True)
            allok = allok and ok_deg and ok_top
    print("ALL PASS" if allok else "FAILURES")


if __name__ == "__main__":
    main()
