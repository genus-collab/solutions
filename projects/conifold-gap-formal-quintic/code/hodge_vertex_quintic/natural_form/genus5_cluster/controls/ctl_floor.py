# -*- coding: utf-8 -*-
"""ctl_floor.py -- THE FLOOR-ORDER RESPONSE IS INDEPENDENT OF k.

   With tau' = tau + c u^k,  E_M'(a,b) - E_M(a,b) = [(tau + c u^k)^m - tau^m] Eh_0(a,b),  m = a+b+1, so
       [u^k] (E_M' - E_M) = m c Eh_0(a,b)     because [u^0] tau^(m-1) = 1,
   a number that does NOT depend on k.  At the floor order u^(k+E-1) of an outer graph with E edges every other
   factor of the row must be taken at u^0, so the floor-order response of each row is the same number for every k.
   Consequence: a single nat_tk control that shows a non-zero floor-order response proves the same non-zero
   response for every other k in the window -- i.e. the regularity test is sensitive to tau at every order u^1..u^(2g-3).

   This script checks that at genus 3 (k = 1, 2, 3, c = 1, 2): the floor-order coefficients of the vertex row (u^k),
   of the two one-edge rows (u^k) and of the five/six two-edge rows (u^(k+1)) must agree across k.

   py -I -B ctl_floor.py [prime]
"""
import os, sys
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctl_lin import table, G, K


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"
    print("=" * 100)
    print("genus %d, 2^%s-1: floor-order response of every row to tau_k -> tau_k + c" % (G, prime))
    print("=" * 100)
    sp0, base = table(prime, 0, 0)
    p = sp0.p
    E = [g['E'] for g in sp0.inner_graphs]
    res = {}
    for c in (1, 2):
        for k in (1, 2, 3):
            sp, tab = table(prime, k, c)
            d = {i: [(tab.get(i, [0] * K)[j] - base.get(i, [0] * K)[j]) % p for j in range(K)] for i in range(len(E))}
            res[(k, c)] = d
    print("row -> floor-order coefficient [u^(k+max(E,1)-1)] of (row(nat_tk) - row(nat)) / c")
    print("%-5s %-3s %-12s %s" % ("row", "E", "floor order", "k=1, k=2, k=3 (c=1) | k=1, k=2, k=3 (c=2, divided by 2)"))
    agree = True
    for i in range(len(E)):
        e = max(E[i], 1)
        vals = []
        for c in (1, 2):
            inv = pow(c, p - 2, p)
            for k in (1, 2, 3):
                j = k + e - 1
                vals.append(res[(k, c)][i][j] * inv % p if j < K else None)
        shown = [v for v in vals if v is not None]
        same = len(set(shown)) <= 1
        agree = agree and same
        if E[i] <= 2:
            print("%-5d %-3d u^(k+%d)      %s   %s" % (i, E[i], e - 1, vals, "AGREE" if same else "DIFFER"))
    print()
    print("floor-order response independent of k and linear in c, for every row with E <= 2: %s" % ("CONFIRMED" if agree else "NOT CONFIRMED"))
    print()
    print("(rows with E >= 3 only have a floor order inside the window for k = 1, so nothing to compare there)")


if __name__ == "__main__":
    main()
