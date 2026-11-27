# -*- coding: utf-8 -*-
"""ur_chain_check.py -- UNIFORM_REDUCTION.md, check (C'): the residue relations of the regular gauge at u = 0 (exact, over Q).

With P = Sum_k z^k P_k(u) the regular gauge of ur_gauge.gauge (entry (a, b) = (b + k, b) of P_k is p_{k,b}):
 (R1) P(0, z) commutes with E_{12} (the residue of both connections at u = 0):
        p_{k,1}(0) = 0 for k != 0 mod 5 (column 1),  p_{k,2-k}(0) = 0 for k != 0 mod 5 (row 2),  p_{k,1}(0) = p_{k,2}(0) for k = 0 mod 5.
 (R2) the H (x) H entry (2, 1) vanishes to second order: p_{k,1} = O(u^2) for k = 1 mod 5.
 (R3) the block rows (2, 3) x columns (0, 1) is, to leading order in u, the rank-one matrix pi(z) (u/z, 1)^T (1, -u/z), pi(z) = P_30(0, z):
        p_{k-1,0}'(0) = p_{k,0}(0),   p_{k-1,1}'(0) = -p_{k,0}(0),   p_{k-2,1}''(0)/2 = -p_{k,0}(0)      for k = 3 mod 5.
 (X1) the order-u propagator: p_{1,0}(0) = p_{1,2}(0) = eps_1/2 and p_{1,3}(0) = p_{1,4}(0), from the first-order data of the periods.
Run: py -I ur_chain_check.py [K] [NU]     (defaults 14, 44)
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ur_gauge import gauge


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    NU = int(sys.argv[2]) if len(sys.argv) > 2 else 44
    OUT = open(os.path.join(HERE, "ur_chain_check.log"), "w", encoding="utf-8"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)

    res = gauge(K, NU); p = res['p']
    c0 = lambda k, b, n=0: p[k][b].coef(n)
    r1a = all(c0(k, 1) == 0 for k in range(1, K + 1) if k % 5)
    r1b = all(c0(k, (2 - k) % 5) == 0 for k in range(1, K + 1) if k % 5)
    r1c = all(c0(k, 1) == c0(k, 2) for k in range(0, K + 1, 5))
    check("(R1) column 1 of P(0, z) vanishes off the diagonal (k <= %d)" % K, r1a)
    check("(R1) row 2 of P(0, z) vanishes off the diagonal (k <= %d)" % K, r1b)
    check("(R1) P_11(0, z) = P_22(0, z)", r1c, ["k=%d: %s" % (k, c0(k, 1)) for k in range(0, K + 1, 5)])
    check("(R2) p_{k,1} = O(u^2) for k = 1 mod 5", all(c0(k, 1) == 0 and c0(k, 1, 1) == 0 for k in range(1, K + 1, 5)))
    for k in range(3, K + 1, 5):
        pi = c0(k, 0)
        ok = c0(k - 1, 0, 1) == pi and c0(k - 1, 1, 1) == -pi and c0(k - 2, 1, 2) == -pi
        check("(R3) k = %d: pi_k = P_30,k(0) = %s;  P_20 = u pi/z, P_31 = -u pi/z, P_21 = -u^2 pi/z^2 to leading order" % (k, pi), ok,
              "" if ok else (c0(k - 1, 0, 1), c0(k - 1, 1, 1), c0(k - 2, 1, 2)))
    # first-order data
    bq, bt = res['bq'], res['bt']
    eps = []
    for a in range(5):
        r = bt[a].mul(bq[a].inv())
        eps.append(r.coef(1))
    log("first-order ratios bt_a/b_a = 1 + eps_a u: eps = %s (sum %s)" % ([str(x) for x in eps], sum(eps)))
    log("p_{1,b}(0), b = 0..4: %s" % [str(c0(1, b)) for b in range(5)])
    check("(X1) p_{1,0}(0) = p_{1,2}(0) = eps_1/2 and p_{1,3}(0) = p_{1,4}(0) = eps_3 - eps_1/2 ... (value check)",
          c0(1, 0) == c0(1, 2) == eps[1] / 2 and c0(1, 3) == c0(1, 4), [str(eps[1] / 2), str(c0(1, 3))])
    log("diagonal entries p_{5,b}(0) in the normalisation p_{5,0}(0) = 0: %s" % [str(c0(5, b)) for b in range(5)])
    log("UR_CHAIN_CHECK " + ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); OUT.close()


if __name__ == "__main__":
    main()
