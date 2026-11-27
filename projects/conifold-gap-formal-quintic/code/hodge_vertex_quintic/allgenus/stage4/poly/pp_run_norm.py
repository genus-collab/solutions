# -*- coding: utf-8 -*-
"""pp_run_norm.py -- Stage 4 / poly: the normalisation c^act(z) = P^act_11(0, z) of the gauge of the QRR-normalised WKB function,
from the scalar hatted polynomials alone (pp_core) and the first-row extraction (pp_gauge), exact over Q, frame (0, 0).

Prints, for q = 0..Q: [z^5q] c^act (must be 0 for q >= 1 iff D = 1 through z^5q), p^act_(5q,0)(0), the free coefficients f_q of the
QRR solution, the h_auto checks (H) for 5 not | n, and two controls: the zero-free-coefficient solution (c = 1/F_inf) and the stored
normalisation (p_(5q,0)(0) = 0), whose c(z) must start 1 + m5 z^5, m5 = -43397/807975000.
Run: py -I -B pp_run_norm.py Q   (log: logs/pp_run_norm_Q<Q>.log)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pp_core import WKB
from pp_series import periods, tau_series
from pp_gauge import first_row, cnorm

M5 = Fr(-43397, 807975000)


def main(Q):
    t0 = time.time(); NK = 5 * Q; K = 4 * Q + 3
    f0, f2 = periods(K); tau = tau_series(K)
    print("pp_run_norm: Q = %d (z^0..z^%d), WKB to w^%d, u-series mod u^%d, frame (0,0), exact over Q" % (Q, 5 * Q, NK, K))
    W = WKB(NK, free='qrr')
    print("(H) h_auto: ph_n(1) = h_n automatically for all n <= %d, 5 not | n: %s (%d checks)"
          % (NK, all(x for _, x in W.h_auto), len(W.h_auto)))
    print("degrees deg ph_n <= floor(4n/5) for n <= %d: %s" % (NK, all(len(p) - 1 <= 4 * n // 5 for n, p in enumerate(W.ph))))
    p, g = first_row(W.ph, NK, K, C0=f0, tau=tau)
    c = cnorm(p, Q)
    print("first-row extraction: all divisibility checks passed (regularity of P_0.)  [%.1f s]" % (time.time() - t0))
    ok = True
    for q in range(Q + 1):
        print("q=%d  [z^%d] c^act = %s   p^act_(5q,0)(0) = %s" % (q, 5 * q, c[q], p[5 * q][0]))
        if q >= 1 and c[q] != 0: ok = False
        if q == 0 and c[q] != 1: ok = False
    for q in range(1, Q + 1): print("   f_%d (free coefficient [u^%d] ph_%d, QRR) = %s" % (q, 4 * q, 5 * q, W.free[q]))
    print("RESULT (N) c^act = 1 through z^%d: %s" % (5 * Q, "PASS" if ok else "FAIL"))
    # control 1: zero free coefficients -> c = 1/F_inf^QRR
    Z = WKB(NK, free='zero'); pz, _ = first_row(Z.ph, NK, K, C0=f0, tau=tau); cz = cnorm(pz, Q)
    F = [Fr(0)] * (Q + 1); F[0] = Fr(1)
    for q in range(1, Q + 1): F[q] = W.free[q]
    inv = [Fr(1)] + [Fr(0)] * Q
    for n in range(1, Q + 1): inv[n] = -sum((F[k] * inv[n - k] for k in range(1, n + 1)), Fr(0))
    print("CONTROL zero-free solution: c != 1 at z^5: %s; c = 1/F_inf^QRR through z^%d: %s"
          % (cz[1] != 0, 5 * Q, cz == inv))
    # control 2: stored normalisation p_(5q,0)(0) = 0  <=> free coefficients shifted so that p_(5q,0)(0) = 0 (solve q by q)
    free = {}
    for q in range(1, Q + 1):
        free[q] = Fr(0)
        S = WKB(5 * q, free=free); ps, _ = first_row(S.ph, 5 * q, 4 * q + 3, C0=f0[:4 * q + 3], tau=tau[:4 * q + 3])
        free[q] = -ps[5 * q][0]
    S = WKB(NK, free=free); ps, _ = first_row(S.ph, NK, K, C0=f0, tau=tau); cs = cnorm(ps, Q)
    print("CONTROL stored normalisation: p_(5q,0)(0) = %s; [z^5] c = %s (= m5: %s); c != 1: %s"
          % ([str(ps[5 * q][0]) for q in range(Q + 1)], cs[1], cs[1] == M5, any(cs[1:])))
    print("time %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
