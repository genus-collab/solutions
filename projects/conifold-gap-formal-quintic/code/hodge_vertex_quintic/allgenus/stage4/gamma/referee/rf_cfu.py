# -*- coding: utf-8 -*-
"""rf_cfu.py -- REFEREE check of Lemma CFU (hatted Chester-Friedman-Ursell map), own code, exact over Q.

What(xi; u) = (l/u) xi^-1 - Sum_(m>=0) u^m xi^(5m+4)/((m+1)(5m+4)),  l = log(1-u);   M(Xi; u) = -(Xi^4/4 + Xi^-1)/tau,
1/tau = 2F1(1, 4/5; 9/5; u) = Sum_n 4 u^n/(5n+4).  Lemma CFU: xi = Xi psi(Xi^5; u), psi = 1 + Sum_k u^k phi_k(y), phi_k in Q[y],
with What(xi; u) = M(Xi; u).  Solved order by order: [u^k] of Xi (What - M) = r_k(y) + (1 - y) phi_k(y) (y = Xi^5).
Independent verification of the result: substitute numerical Xi (rational) and compare u-series directly.
Run: py -I -B rf_cfu.py K     log: logs/rf_cfu_<K>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_ring import FQ, mul, inv, add, scal

Q = FQ()


def padd(a, b, c=1):
    n = max(len(a), len(b)); a = a + [Fr(0)] * (n - len(a)); b = b + [Fr(0)] * (n - len(b))
    out = [x + c * y for x, y in zip(a, b)]
    while out and out[-1] == 0: out.pop()
    return out


def pmul(a, b):
    if not a or not b: return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out


def bmul(A, B, K):
    """bivariate: lists over u of polynomials in y, mod u^K"""
    out = [[] for _ in range(K)]
    for i, a in enumerate(A[:K]):
        for j, b in enumerate(B[:K - i]): out[i + j] = padd(out[i + j], pmul(a, b))
    return out


def binv(A, K):
    """1/A, A[0] = [1]"""
    assert A[0] == [Fr(1)]
    out = [[Fr(1)]] + [[] for _ in range(K - 1)]
    for k in range(1, K):
        acc = []
        for j in range(1, k + 1): acc = padd(acc, pmul(A[j], out[k - j]))
        out[k] = padd([], acc, -1)
    return out


def lhs(psi, K, tau_one=False):
    """Xi (What(Xi psi) - M(Xi)) as bivariate series mod u^K"""
    lu = [[Fr(-1, n + 1)] for n in range(K)]                              # l/u
    it = [[Fr(4, 5 * n + 4)] if not tau_one else ([Fr(1)] if n == 0 else []) for n in range(K)]   # 1/tau
    out = bmul(lu, binv(psi, K), K)
    p4 = bmul(bmul(psi, psi, K), bmul(psi, psi, K), K); p5 = bmul(p4, psi, K); pw = p4
    for m in range(K):
        term = [[] for _ in range(K)]
        for i, a in enumerate(pw[:K - m]): term[i + m] = padd([], pmul([Fr(0)] * (m + 1) + [Fr(1)], a), Fr(-1, (m + 1) * (5 * m + 4)))
        out = [padd(x, y) for x, y in zip(out, term)]
        pw = bmul(pw, p5, K)
    tail = bmul(it, [[Fr(1), Fr(1, 4)]] + [[] for _ in range(K - 1)], K)      # + (1 + y/4)/tau
    return [padd(x, y) for x, y in zip(out, tail)]


def solve(K, tau_one=False):
    psi = [[Fr(1)]] + [[] for _ in range(K)]; info = dict(r1=[], phi1=[])
    for k in range(1, K + 1):
        r = lhs(psi, k + 1, tau_one)[k]
        info['r1'].append(sum(r))
        if sum(r) != 0: return psi, info, k
        q = []; rem = list(r)                                               # phi_k = -r/(1 - y) = r/(y - 1)
        while len(rem) > 1:
            c = rem[-1]; q.insert(0, c); rem[-2] += c; rem.pop()
        assert rem == [] or rem[0] == 0
        psi[k] = q; info['phi1'].append(sum(q))
    return psi, info, None


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "rf_cfu_%d.log" % K), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("rf_cfu: Lemma CFU mod u^%d, exact over Q" % (K + 1))
    wh1 = [Fr(-1, n + 1) - Fr(1, (n + 1) * (5 * n + 4)) for n in range(K + 1)]
    check("(K1) What(1; u) = -(5/4)/tau(u) = -Sum 5u^n/(5n+4) mod u^%d" % (K + 1), wh1 == [Fr(-5, 5 * n + 4) for n in range(K + 1)])
    psi, info, obst = solve(K)
    check("(K2) solvable at every order k <= %d: r_k(1) = 0" % K, obst is None and all(x == 0 for x in info['r1']))
    check("(K3) phi_k(1) = 0 for k <= %d (the five saddles Xi^5 = 1 are fixed)" % K, all(x == 0 for x in info['phi1']))
    degs = [len(psi[k]) - 1 for k in range(1, K + 1)]
    check("(K4) deg phi_k = k for k <= %d" % K, degs == list(range(1, K + 1)), degs)
    check("(K4) phi_1 = (1 - y)/18, phi_2 = 13/567 - (101/4536) y - y^2/1512",
          psi[1] == [Fr(1, 18), Fr(-1, 18)] and psi[2] == [Fr(13, 567), Fr(-101, 4536), Fr(-1, 1512)])
    bad = []
    for X in (Fr(2), Fr(3, 2), Fr(-1, 3), Fr(5, 7)):                     # direct substitution, u-series only
        y = X ** 5; ps = [sum(c * y ** i for i, c in enumerate(psi[k])) for k in range(K + 1)]
        xi = scal(Q, ps, X); n = K + 1
        wh = scal(Q, mul(Q, [Fr(-1, m + 1) for m in range(n)], inv(Q, xi, n), n), 1)
        x5 = mul(Q, mul(Q, mul(Q, xi, xi, n), mul(Q, xi, xi, n), n), xi, n); pw = mul(Q, mul(Q, xi, xi, n), mul(Q, xi, xi, n), n)
        for m in range(n):
            wh = add(Q, wh, [Fr(0)] * m + scal(Q, pw, Fr(-1, (m + 1) * (5 * m + 4)))[:n - m]); pw = mul(Q, pw, x5, n)
        Mv = scal(Q, [Fr(4, 5 * m + 4) for m in range(n)], -(X ** 4 / 4 + 1 / X))
        if wh[:n] != Mv: bad.append(X)
    check("(K5) What(Xi psi(Xi^5; u); u) = M(Xi; u) mod u^%d by direct substitution at Xi = 2, 3/2, -1/3, 5/7" % (K + 1), not bad, bad)
    _, info_c, obst_c = solve(3, tau_one=True)
    check("(C) control: tau -> 1 (v = u) is obstructed at order 1 (r_1(1) != 0)", obst_c == 1, info_c['r1'])
    log("     phi_3 = %s" % [str(c) for c in psi[3]])
    log("RF_CFU %d: %s   (%.1f s)" % (K, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
