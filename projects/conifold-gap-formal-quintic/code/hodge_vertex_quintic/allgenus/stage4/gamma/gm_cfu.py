# -*- coding: utf-8 -*-
"""gm_cfu.py -- the hatted Chester-Friedman-Ursell map of GAMMA.md s.3 (Lemma CFU), exact over Q.

Hatted phase of the full-Stirling 1-form (x = s xi, z = s^4 zhat, u = s^5):
   What(xi; u) = a(u) xi^-1 + Sum_l g_l u^l xi^(5l+4),  a = log(1-u)/u = -Sum_n u^n/(n+1),  g_l = -1/((l+1)(5l+4)).
Model phase at v(u) = u tau^(-5/4) in the variable Xi (X = v^(1/5) Xi):  M(Xi; u) = -(Xi^4/4 + Xi^-1)/tau(u),  1/tau = Sum_n 4u^n/(5n+4).
Lemma CFU: there is a unique xi = Xi psi(Xi^5; u), psi = 1 + Sum_{k>=1} u^k phi_k(y), phi_k in Q[y], with What(xi; u) = M(Xi; u).
Order u^k (y = Xi^5): (1 - y) phi_k(y) = -[u^k] F(psi_<k),  F(psi) := a/psi + Sum_l g_l u^l y^(l+1) psi^(5l+4) + (1 + y/4)/tau;
solvable iff the right side vanishes at y = 1, which is the critical-value identity What(1; u) = -(5/4)/tau(u).
Checks: (K1) What(1; u) = -(5/4)/tau(u) mod u^K; (K2) the right side vanishes at y = 1 for k <= K (phi_k polynomial);
(K3) phi_k(1) = 0 (critical point fixed); (K4) the identity What(xi(Xi); u) = M(Xi; u) holds mod u^(K+1) as Laurent polynomials;
reported: deg phi_k (if deg phi_k <= k the flat map x(X; u) is analytic at (X, u) = (0, 0); not needed for the proof).
Control: with v(u) replaced by u (tau -> 1) the obstruction appears at k = 1.
Run: py -I -B gm_cfu.py [K]   log: logs/gm_cfu_<K>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def padd(a, b, c=1):
    n = max(len(a), len(b)); out = [Fr(0)] * n
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i] += c * x
    while out and out[-1] == 0: out.pop()
    return out
def pmul(a, b):
    if not a or not b: return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] += x * y
    return out
def pval(a, x): return sum(c * x ** i for i, c in enumerate(a))


def umul(A, B, K):
    """product of u-series with polynomial coefficients (lists over u of polys in y), mod u^(K+1)"""
    out = [[] for _ in range(K + 1)]
    for i in range(K + 1):
        for j in range(K + 1 - i):
            if A[i] and B[j]: out[i + j] = padd(out[i + j], pmul(A[i], B[j]))
    return out
def uinv(A, K):
    """1/A, A[0] = [1]"""
    out = [[] for _ in range(K + 1)]; out[0] = [Fr(1)]
    for k in range(1, K + 1):
        s = []
        for j in range(1, k + 1):
            if A[j] and out[k - j]: s = padd(s, pmul(A[j], out[k - j]))
        out[k] = [-x for x in s]
    return out
def upow(A, e, K):
    out = [[Fr(1)]] + [[] for _ in range(K)]; base = A
    while e:
        if e & 1: out = umul(out, base, K)
        base = umul(base, base, K); e >>= 1
    return out


def F_of(psi, K, tau_inv):
    a = [Fr(-1, n + 1) for n in range(K + 1)]
    pinv = uinv(psi, K)
    out = [[] for _ in range(K + 1)]
    for i in range(K + 1):
        for j in range(K + 1 - i):
            if pinv[j]: out[i + j] = padd(out[i + j], [a[i] * x for x in pinv[j]])
    for l in range(K + 1):
        g = Fr(-1, (l + 1) * (5 * l + 4)); P = upow(psi, 5 * l + 4, K - l)
        for j in range(K + 1 - l):
            if P[j]: out[l + j] = padd(out[l + j], [Fr(0)] * (l + 1) + [g * x for x in P[j]])
    for k in range(K + 1):
        out[k] = padd(out[k], [tau_inv[k], tau_inv[k] / 4])
    return out


def solve(K, tau_inv):
    psi = [[Fr(1)]] + [[] for _ in range(K)]; info = []
    for k in range(1, K + 1):
        F = F_of(psi, k, tau_inv); r = [-x for x in F[k]]
        r1 = pval(r, Fr(1)) if r else Fr(0)
        if r1 != 0: info.append((k, False, None, None)); return psi, info
        rem = list(r); n = len(rem)                                     # divide r by (1 - y)
        qq = [Fr(0)] * max(n - 1, 0); carry = Fr(0)
        for i in range(n - 1): carry += rem[i]; qq[i] = carry           # r/(1-y): q_i = r_0 + .. + r_i
        while qq and qq[-1] == 0: qq.pop()
        psi[k] = qq; info.append((k, True, len(qq) - 1, pval(qq, Fr(1)) if qq else Fr(0)))
    return psi, info


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "gm_cfu_%d.log" % K), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + (": " + str(det) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("gm_cfu: hatted CFU map mod u^%d, exact over Q" % (K + 1))
    tau_inv = [Fr(4, 5 * n + 4) for n in range(K + 2)]
    W1 = [Fr(-1, n + 1) - Fr(1, (n + 1) * (5 * n + 4)) for n in range(K + 1)]
    check("(K1) critical value: What(1; u) = -(5/4)/tau(u) mod u^%d" % (K + 1), W1 == [-Fr(5, 4) * t for t in tau_inv[:K + 1]])
    psi, info = solve(K, tau_inv)
    check("(K2) solvability at every order k <= %d (right side vanishes at y = 1)" % K, all(ok for (_, ok, _, _) in info) and len(info) == K)
    check("(K3) phi_k(1) = 0 for k <= %d (the critical point Xi = 1 is fixed)" % K, all(v == 0 for (_, _, _, v) in info))
    log("   deg phi_k:", [d for (_, _, d, _) in info], "  phi_1 =", [str(x) for x in psi[1]], "  phi_2 =", [str(x) for x in psi[2]])
    F = F_of(psi, K, tau_inv)
    check("(K4) What(xi(Xi); u) = M(Xi; u) identically mod u^%d" % (K + 1), all(not f for f in F))
    _, info0 = solve(min(K, 3), [Fr(1)] + [Fr(0)] * (K + 1))
    check("(C) control: with tau = 1 (v = u) the order-1 equation is obstructed", info0 and not info0[0][1])
    log("GM_CFU %d: %s   (%.1f s)" % (K, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
