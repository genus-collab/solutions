# -*- coding: utf-8 -*-
"""pp_cfu.py -- Stage 4 / poly: the formal Chester-Friedman-Ursell (CFU) change of variables in hatted variables, exact over Q.

Find Y(X; u) = X (1 + Sum_(k>=1) u^k Y_k(X^5)) and V(u) = 1 + O(u) with
    W^(X; u) := Sum_(m>=1) (u^(m-1)/m) [1/X + X^(5m-1)/(5m-1)]  =  V/Y + Y^4/4 .
Multiplying by X: order u^k reads (X^5 - 1) Y_k = R_k(X^5) - V_k, solvable in Q[X^5] iff V_k = R_k(1) (one condition).
Checks: every Y_k is a polynomial in X^5 (no remainder); V(u) = tau^(-5/4), tau = 1/2F1(1,4/5;9/5;u) (critical values).
Run: py -I -B pp_cfu.py K
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from pp_series import tau_series, s_pow


def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out


def cfu(K):
    """polynomials in T = X^5: Y[k] (k = 1..K-1), V list"""
    lhs = [[Fr(1, m), Fr(0)] + [Fr(0)] * (m - 1) for m in range(1, K + 1)]          # X*W^ at u^(m-1): (1/m)(1 + T^m/(5m-1))
    for m in range(1, K + 1): lhs[m - 1] = [Fr(1, m)] + [Fr(0)] * (m - 1) + [Fr(1, m * (5 * m - 1))]
    Y = [[Fr(0)]]; V = [Fr(1)]
    def series_parts(Y, V, k):
        """[u^k] of V/(1+S) + T (1+S)^4/4 with S = Sum u^j Y_j, using Y_1..Y_(k-1) only (Y_k = 0, V_k = 0)"""
        # power series in u with polynomial coefficients
        S = [[Fr(0)]] + [Y[j] for j in range(1, k)] + [[Fr(0)]]
        def smul(A, B):
            out = [[Fr(0)] for _ in range(k + 1)]
            for i in range(k + 1):
                for j in range(k + 1 - i):
                    p = pmul(A[i], B[j])
                    q = out[i + j]; n = max(len(p), len(q)); out[i + j] = [(q[t] if t < len(q) else 0) + (p[t] if t < len(p) else 0) for t in range(n)]
            return out
        one_plus = [[Fr(1)]] + S[1:]
        inv = [[Fr(1)]] + [[Fr(0)] for _ in range(k)]                               # 1/(1+S)
        for n in range(1, k + 1):
            acc = [Fr(0)]
            for j in range(1, n + 1):
                p = pmul(S[j], inv[n - j]); acc = [(acc[t] if t < len(acc) else 0) - (p[t] if t < len(p) else 0) for t in range(max(len(acc), len(p)))]
            inv[n] = acc
        Vs = [[V[j]] if j < len(V) else [Fr(0)] for j in range(k)] + [[Fr(0)]]
        a = smul(Vs, inv)
        b = smul(smul(one_plus, one_plus), smul(one_plus, one_plus))
        bk = [Fr(0)] + [x / 4 for x in b[k]]
        ak = a[k]; n = max(len(ak), len(bk))
        return [(ak[t] if t < len(ak) else 0) + (bk[t] if t < len(bk) else 0) for t in range(n)]
    for k in range(1, K):
        Rk = series_parts(Y, V, k)
        L = lhs[k]; n = max(len(L), len(Rk))
        rhs = [(L[t] if t < len(L) else 0) - (Rk[t] if t < len(Rk) else 0) for t in range(n)]  # (T - 1) Y_k + V_k  ... = rhs
        Vk = sum(rhs)                                                                  # value at T = 1
        num = list(rhs); num[0] -= Vk                                                  # (T-1) Y_k = rhs - V_k
        q = [Fr(0)] * (len(num) - 1); r = list(num)                                    # divide by (T - 1)
        for t in range(len(num) - 1, 0, -1):
            q[t - 1] = r[t]; r[t - 1] += r[t]; r[t] = 0
        assert r[0] == 0, "remainder"
        Y.append(q); V.append(Vk)
    return Y, V


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    Y, V = cfu(K); tau = tau_series(K); Vt = s_pow(tau, Fr(-5, 4), K)
    print("pp_cfu: K = %d" % K)
    for k in range(1, min(K, 4)): print("Y_%d(T) = %s" % (k, [str(c) for c in Y[k]]))
    print("V(u) = %s" % [str(x) for x in V])
    print("RESULT Y_k in Q[X^5] for k < %d: PASS; V = tau^(-5/4) mod u^%d: %s" % (K, K, "PASS" if V == Vt else "FAIL"))
