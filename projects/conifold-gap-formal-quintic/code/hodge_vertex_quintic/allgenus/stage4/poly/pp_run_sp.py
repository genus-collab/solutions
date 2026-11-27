# -*- coding: utf-8 -*-
"""pp_run_sp.py -- Stage 4 / poly: compare the formal stationary phase Phi^SP of the Mellin-Barnes summand (pp_sp) with the
QRR-normalised hatted WKB polynomials ph_n (pp_core), mod 2^61-1 and 2^89-1.
Checks per n <= N: (a) [w^n] Phi^SP mod u^K equals ph_n (all coefficients u^0..u^(K-1)); (b) [w^n] Phi^SP vanishes above
u^floor(4n/5) (polynomiality, for the orders inside the window).
Run: py -I -B pp_run_sp.py N K     (log: logs/pp_run_sp_N<N>_K<K>.log)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pp_core import WKB
from pp_sp import phi_sp


def main(N, K):
    t0 = time.time(); W = WKB(N); allok = True
    print("pp_run_sp: N = %d, K = %d (u^0..u^%d), primes 2^61-1, 2^89-1" % (N, K, K - 1))
    for p in (2 ** 61 - 1, 2 ** 89 - 1):
        sp = phi_sp(p, N, K); nsame = 0; npoly = 0; ncoef = 0
        for n in range(N + 1):
            poly = [Fr(c).numerator % p * pow(Fr(c).denominator % p, p - 2, p) % p for c in W.ph[n]] + [0] * K
            same = all(sp[n][k] == poly[k] for k in range(K)); nsame += same
            dk = 4 * n // 5; pol = all(sp[n][k] == 0 for k in range(dk + 1, K)); npoly += pol; ncoef += K
            if not (same and pol):
                allok = False; print("  p=%d n=%d MISMATCH same=%s poly=%s" % (p, n, same, pol))
        print("p = %d: SP == QRR for %d/%d orders w^0..w^%d (%d u-coefficients each); polynomial (deg <= floor(4n/5)) %d/%d  [%.1f s]"
              % (p, nsame, N + 1, N, K, npoly, N + 1, time.time() - t0))
    print("RESULT Phi^SP == Phi^QRR (incl. the free coefficients f_q, q <= %d) mod u^%d: %s" % (N // 5, K, "PASS" if allok else "FAIL"))


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
