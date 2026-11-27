# -*- coding: utf-8 -*-
"""ur_c5.py -- the fitted constant c_5 of ur_gauge_run.log at 2^61-1 against the rational value found at 2^89-1 and against
-p_{5,1}(0) = -p_{5,2}(0) of the gauge normalised by p_{5,0}(0) = 0 (ur_chain_check.log): the true gauge has P_11(0,z) = P_22(0,z) = 1 + O(z^10)."""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from fq_fp import FpRing, P1, P2
c5 = Fr(43397, 807975000)
r61 = 1868382097900760522
ok = FpRing(P1).of(c5) == r61
print("c_5 = %s; residue at 2^61-1 equals the fitted residue %d: %s" % (c5, r61, ok))
print("denominator factorisation: 807975000 = 2^3 * 3^5 * 5^5 * 7 * 19 :", 2**3 * 3**5 * 5**5 * 7 * 19 == 807975000)
