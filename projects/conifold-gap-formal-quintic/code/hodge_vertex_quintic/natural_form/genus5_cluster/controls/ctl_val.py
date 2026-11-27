# -*- coding: utf-8 -*-
"""ctl_val.py -- the u-valuations that drive the prediction: val(E_out) for the natural M, for the bare symbol M_0
   and for the tau_k controls, plus the coefficients of tau and the valuation of t_c(u).

   py -I -B ctl_val.py [genus ...]     (default 4 5)

   For a variant with tau' = tau + d u^k one has E_M' - E_M = (a+b+1) tau^(a+b) d u^k + O(u^(2k)), so
   E_out' - E_out = -(that), of valuation exactly k wherever Eh_0(a, b) != 0.  This script measures that directly.
"""
import os, sys, collections
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)


def report(g, prime="89"):
    from g5_engine import make_split
    from fq4_split import tau_pred
    from fractions import Fraction as Fr
    from fq_fp import P1, P2
    p = P2 if prime == "89" else P1
    K = 2 * g - 2
    print("=" * 100)
    print("genus %d, 2^%s-1, K = %d (u^0..u^%d), zhat-order KZ = %d" % (g, prime, K, K - 1, 3 * g - 3))
    print("   tau = %s" % [str(x) for x in tau_pred(K)])
    out = {}
    for variant in ("nat", "M0") + tuple("nat_t%d" % k for k in range(1, K)):
        sp = make_split(p, g, (0, 0), K, variant)
        h = collections.Counter(sp.valO.values())
        nz = sum(1 for key in sp.EO for i in range(5) for j in range(5) if sp.EO[key][i][j] is not None)
        print("   %-8s min val(E_out) over all %d keys (a,b) = %d;  histogram of the per-key valuation: %s"
              % (variant, len(sp.valO), min(sp.valO.values()), sorted(h.items())))
        out[variant] = sp
        if variant == "nat":
            # per-entry valuation, not just per key
            hh = collections.Counter()
            for key in sp.EO:
                for i in range(5):
                    for j in range(5):
                        hh[sp.val(sp.EO[key][i][j])] += 1
            print("            per-entry valuation histogram (val = K means identically zero mod u^K): %s" % sorted(hh.items()))
            print("            non-zero entries of E_out: %d of %d" % (nz, 25 * len(sp.EO)))
    # difference of E_out between the controls and nat: valuation must be exactly k (M0 behaves as k = 1)
    nat = out["nat"]
    print("   valuation of E_out(control) - E_out(nat), minimised over all keys and entries:")
    for variant in ("M0",) + tuple("nat_t%d" % k for k in range(1, K)):
        sp = out[variant]
        vals = [nat.val(nat.sub(sp.EO[key][i][j], nat.EO[key][i][j])) for key in nat.EO for i in range(5) for j in range(5)]
        h = collections.Counter(vals)
        print("      %-8s min = %d   histogram %s" % (variant, min(vals), sorted(h.items())))
    # t_c(u): the change of variables that turns [u^k] into [t_c^(k-(2g-2))]
    from fq_hat import tc_series
    tc = tc_series(nat.C, nat.R, K + 4)
    print("   t_c(u): valuation %s, leading coefficient %s (non-zero), so u^k <-> t_c^k up to triangular mixing" % (tc.v, tc.coef(tc.v)))


if __name__ == "__main__":
    gs = [int(x) for x in sys.argv[1:]] or [4, 5]
    for g in gs:
        report(g)
