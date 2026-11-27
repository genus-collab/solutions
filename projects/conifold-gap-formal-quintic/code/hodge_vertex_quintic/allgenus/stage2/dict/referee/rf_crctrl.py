"""Referee: detection power of the critic's falsifier (critic/cr_dictionary.py, imported READ-ONLY, nothing written).
The H^2 leg v/p is replaced IN MEMORY by v^2/p^6 (same weight 4, outside Q[p,q]); Test A and Test B must then FAIL.
Usage: py -I -B rf_crctrl.py [KMAX] [WMAX]      (default 12 48; runs in under a second)
"""
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CRITIC = os.path.normpath(os.path.join(HERE, "..", "..", "..", "critic"))
KMAX = sys.argv[1] if len(sys.argv) > 1 else "12"
WMAX = sys.argv[2] if len(sys.argv) > 2 else "48"
sys.argv = [sys.argv[0], KMAX, WMAX]          # cr_dictionary reads KMAX, WMAX from argv at import time
sys.path.insert(0, CRITIC)
import cr_dictionary as CD  # noqa: E402

if __name__ == "__main__":
    print("rf_crctrl.py: critic falsifier with the H^2 leg v/p replaced by v^2/p^6; KMAX=%s WMAX=%s" % (KMAX, WMAX))
    CD.LEG0[2] = {(-6, 0, 2): Fr(1)}            # critic keys (p-exponent, log power, v-exponent)
    T = CD.leg_tables()
    CD.test_a(T)
    CD.test_b(T)
    n = len(CD.FAILS)
    nb = sum(1 for f in CD.FAILS if "NOT in L" in f)
    print("RESULT: %d FAIL lines, %d of them 'log-free combination NOT in L' (expected > 0: the falsifier detects the bad leg)"
          % (n, nb))
