# -*- coding: utf-8 -*-
"""ms_singlist.py -- task 2(b): EVERY computed correlator <tau_k1(H^a1) .. tau_kn(H^an)>_h with negative v-exponent, with its remainder
      r = <..>_h - F_h^(n)(v) prod PiC(tau_ki H^ai)          (the (ML) prediction on K; h = 0: F_0 = -(5/2) v^2 log v, n >= 3)
and its D-content: the insertions that are NOT in K, i.e. have Pi_m != 0 for some m >= 1 (a component along a descendant tau_m(H),
m >= 1, of the vanishing direction).  (ML) predicts r = 0 whenever every insertion is in K; the script checks that and lists the rest.
Output: summary per (h, n) on stdout; the full list in ms_singlist_full.log.  Run: py -I -B ms_singlist.py
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc
import ms_engine as me
import ms_model as mm


def F0der(n):
    c = Fr(-5)
    for j in range(n - 3): c *= -(j + 1)
    return c


def in_K(ins):
    """tau_k(H^b) is in K iff Pi_m(tau_k H^b) = PiC(tau_(k-m) H^b) = 0 for 1 <= m <= k"""
    b, k = ins
    return all(mm.pic(b, k - m)[0] == 0 for m in range(1, k + 1))


def main():
    plan = [(0, 3), (0, 4), (0, 5), (1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2)]
    full = open(os.path.join(HERE, "ms_singlist_full.log"), "w", encoding="ascii")
    full.write("h | insertions | v-exponent | value at v=1 | F_h^(n) prod PiC | remainder | insertions not in K\n")
    allK_ok = True
    for h, n in plan:
        if not os.path.exists(me._cache_path(h, n)): print("skip (%d,%d): no amplitude" % (h, n)); continue
        nsing = nK = nKbad = nD = nDnz = 0
        for ins in mc.all_insertions(h, n):
            e = mc.vexp(h, ins)
            if e >= 0: continue
            val = mc.corr(h, ins)
            pred = (F0der(n) if h == 0 else mc.Fder(h, n))
            for x in ins: pred *= mc.pic(x)
            r = val - pred; outK = [x for x in ins if not in_K(x)]; nsing += 1
            if not outK:
                nK += 1
                if r != 0: nKbad += 1; allK_ok = False
            else:
                nD += 1; nDnz += (r != 0)
            if val != 0 or pred != 0:
                full.write("%d | %s | %d | %s | %s | %s | %s\n" % (h, mc.fmt_ins(ins), e, val, pred, r, mc.fmt_ins(outK) if outK else "-"))
        print("(%d,%d): %d singular correlators; %d with all insertions in K (remainder != 0 for %d); %d with D-content (remainder != 0 for %d)"
              % (h, n, nsing, nK, nKbad, nD, nDnz), flush=True)
    full.close()
    print("ALL-K singular correlators have zero remainder: %s" % allK_ok)


if __name__ == "__main__":
    main()
