# -*- coding: utf-8 -*-
"""ms_delta4.py -- why the canonical normalisation matters: first-order effect of P -> P exp(c z^(2l-1)) on the genus-g free energy.

For a scalar symplectic factor c(z) = exp(c z^m), m odd, the Givental action changes the leaves by c z^(m+1) 1 and the edges by
c (w^m + z^m)/(w + z) eta^-1 = c Sum_{j=0}^{m-1} (-1)^j w^(m-1-j) z^j eta^-1.  At the leading order (frozen theory, R' = 1) the change
of the genus-g free energy of the model, per unit c, is
   Delta_g(m) = <tau_(m+1)(1)>_g + (1/2) Sum_j (-1)^j [ <tau_(m-1-j)(e_a) tau_j(e^a)>_(g-1)
                                                       + Sum_{h1 + h2 = g, h_i >= 1} <tau_(m-1-j)(e_a)>_h1 <tau_j(e^a)>_h2 ],
eta^-1 = (1/5) Sum_a H^a (x) H^(3-a)  (N_a N_(3-a) = 1).  This is (a multiple of) Int ch_m(E) Omega^LG (Mumford/GRR).
For m = 5 (the z^5 normalisation of the gauge, hatted u^4) and g = 4 this is the predicted coefficient of the change of the t_c^(-2)
term of the genus-4 graph sum when the gauge is NOT canonically normalised; Delta_4(5) != 0 means the genus-4 test sees the normalisation.
A control: Delta_g(m) for the frozen theory must equal d/dc of F_g along the flow, which by homogeneity is c-independent; we also print
the pieces.  Run: py -I -B ms_delta4.py [g] [m]
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_corr as mc
import ms_engine as me


def corr_or_none(h, ins):
    n = len(ins)
    if 2 * h - 2 + n <= 0: return None
    if not os.path.exists(me._cache_path(h, n)): raise RuntimeError("amplitude (%d,%d) missing" % (h, n))
    if sum(k for a, k in ins) > 3 * h - 3 + n: return Fr(0)
    return mc.corr(h, tuple(sorted(ins)))


def delta(g, m):
    leaf = corr_or_none(g, ((0, m + 1),))
    loop = Fr(0); bridge = Fr(0)
    for j in range(m):
        s = (-1) ** j
        for a in range(5):
            b = (3 - a) % 5
            x = corr_or_none(g - 1, ((a, m - 1 - j), (b, j)))
            if x: loop += s * x / 5
            for h1 in range(1, g):
                h2 = g - h1
                y1 = corr_or_none(h1, ((a, m - 1 - j),)); y2 = corr_or_none(h2, ((b, j),))
                if y1 and y2: bridge += s * y1 * y2 / 5
    return leaf, loop / 2, bridge / 2


if __name__ == "__main__":
    g = int(sys.argv[1]) if len(sys.argv) > 1 else 4; m = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    lf, lp, br = delta(g, m)
    print("g = %d, m = %d:  leaf %s  loop %s  bridges %s  ->  Delta = %s   (v-exponent of the change: %s)"
          % (g, m, lf, lp, br, lf + lp + br, me.vexp(g, ((0, m + 1),))))
