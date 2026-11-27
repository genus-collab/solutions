# -*- coding: utf-8 -*-
"""g1_cyclic_check.py -- evidence (not proof) for prop:set-cyclicLR (gap cluster G1).

Exact rational arithmetic (fractions.Fraction), Python stdlib only.  Run:  py -I -B g1_cyclic_check.py [N] [K]
(N = q-order, default 30; K = w-order of the starting series, default 12).  Writes g1_cyclic_check.log next to itself
(relative path only is printed).

Set-up (GAPS_G1.md s.2).  At the zeta point (lambda_i = zeta^i, H^5 = 1) the I-function is
    I/z = q^w Fhat(w, q),   Fhat(w, q) = Sum_d q^d prod_{k=1}^{5d} (5w + k) / prod_{k=1}^{d} ((w + k)^5 - w^5),   w <-> H/z,
where every monomial w^k stands for H^(k mod 5) z^(-k)  (since (H + kz)^5 - 1 = z^5 ((w + k)^5 - w^5) when H^5 = 1).
Birkhoff recursion (the lowered QDE z d_T s_p = kappa_p s_(p+1), s_p = H^p q^w g_p(w)):
    g_0 = Fhat / Fhat(0),   I_(p+1) := 1 + D [w^1] g_p,   g_(p+1) := (1 + D/w) g_p / I_(p+1),   D = q d/dq,
so that C_1 kappa_p = I_(p+1) (p = 0..4).  Claims checked (all as q-series to O(q^(N+1))):
    I_1 = C_1,  I_2 C_0^2 C_1^2 = L^5,  I_3 = C_1,  I_4 = C_0,  I_5 = C_0,  g_5 = g_0  (closing),
i.e. kappa = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1), L^5 = 1/(1 - 5^5 q).  Also, in Q[[q]][ell] (ell = formal log q,
D ell = 1), with c_k = Sum_j ell^j/j! [w^(k-j)] Fhat:  PF4 c_k = 0 (k <= 3), PF4 c_4 = 1, det(D^j c_k)_(j,k<=3) = L^10 (Lemmas D, E).
Controls that MUST fail:  (c1) the non-equivariant series ((w + k)^5 in place of (w + k)^5 - w^5): the closing g_5 = g_0 fails;
(c2) a wrong twist (O(4) or O(6) in place of O(5), same denominators): I_2 C_0^2 C_1^2 = L^5 fails;
(c3) the wrong guess kappa_3 = C_1/C_0: fails.
"""
import os, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
LOGNAME = "g1_cyclic_check.log"


def qmul(a, b, N):
    out = [Fr(0)] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                y = b[j]
                if y:
                    out[i + j] += x * y
    return out


def qinv(a, N):
    assert a[0] != 0
    out = [Fr(0)] * (N + 1)
    out[0] = 1 / a[0]
    for n in range(1, N + 1):
        s = Fr(0)
        for k in range(1, n + 1):
            if a[k]:
                s += a[k] * out[n - k]
        out[n] = -s / a[0]
    return out


def qD(a):
    return [n * x for n, x in enumerate(a)]


def qadd(a, b):
    return [x + y for x, y in zip(a, b)]


def qscal(a, c):
    return [c * x for x in a]


def wpoly_mul(p, r, K):
    out = [Fr(0)] * (K + 1)
    for i, x in enumerate(p):
        if x and i <= K:
            for j, y in enumerate(r):
                if i + j > K:
                    break
                out[i + j] += x * y
    return out


def wpoly_inv(p, K):
    out = [Fr(0)] * (K + 1)
    out[0] = 1 / p[0]
    for n in range(1, K + 1):
        s = Fr(0)
        for k in range(1, min(n, len(p) - 1) + 1):
            s += p[k] * out[n - k]
        out[n] = -s / p[0]
    return out


def fhat(N, K, deg=5, equivariant=True):
    """Fhat as a list over w-order m (0..K) of q-series (length N+1)."""
    F = [[Fr(0)] * (N + 1) for _ in range(K + 1)]
    for d in range(N + 1):
        num = [Fr(1)]
        for k in range(1, deg * d + 1):
            num = wpoly_mul(num, [Fr(k), Fr(deg)], K)
        den = [Fr(1)]
        for k in range(1, d + 1):
            pk = [Fr(k) ** 5, 5 * Fr(k) ** 4, 10 * Fr(k) ** 3, 10 * Fr(k) ** 2, 5 * Fr(k)]
            if not equivariant:
                pk = pk + [Fr(1)]          # (w + k)^5 including w^5
            den = wpoly_mul(den, pk, K)
        A = wpoly_mul(num, wpoly_inv(den, K), K)
        for m in range(K + 1):
            F[m][d] += A[m]
    return F


def birkhoff(F, N, P=5):
    """returns (Is, gs): Is[p] = I_p for p = 1..P, gs[p] = g_p as list over w-order of q-series."""
    K = len(F) - 1
    C0inv = qinv(F[0], N)
    g = [qmul(F[m], C0inv, N) for m in range(K + 1)]
    gs = {0: g}
    Is = {}
    for p in range(P):
        Kp = len(g) - 1
        if Kp < 2:
            break
        h = [qadd(g[m], qD(g[m + 1])) for m in range(Kp)]       # (1 + D/w) g  (D g_0 = 0)
        Ip = h[0]
        Is[p + 1] = Ip
        Iinv = qinv(Ip, N)
        g = [qmul(x, Iinv, N) for x in h]
        gs[p + 1] = g
    return Is, gs


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    lines = []
    def log(s):
        print(s); lines.append(s)
    fails = 0
    def check(name, ok, must_fail=False):
        nonlocal fails
        tag = ("PASS" if ok else "FAIL") if not must_fail else ("PASS (control fails as required)" if not ok else "FAIL (control did not fail)")
        if (not must_fail and not ok) or (must_fail and ok):
            fails += 1
        log("%s  %s" % (tag, name))
    log("g1_cyclic_check: N = %d (q-order), K = %d (w-order)" % (N, K))
    F = fhat(N, K)
    C0 = F[0]
    Is, gs = birkhoff(F, N)
    It1 = F[1]                                   # [w^1] Fhat = Itilde_1
    C1 = qadd([Fr(1)] + [Fr(0)] * N, qD(qmul(It1, qinv(C0, N), N)))      # C_1 = D(log q + Itilde_1/C_0)
    L5 = [Fr(3125) ** n for n in range(N + 1)]   # 1/(1 - 3125 q)
    check("C_0 = Sum (5d)!/(d!)^5 q^d", all(C0[d] == Fr(__import__('math').factorial(5 * d), __import__('math').factorial(d) ** 5) for d in range(N + 1)))
    check("I_1 = C_1  (kappa_0 = 1)", Is[1] == C1)
    check("I_2 C_0^2 C_1^2 = L^5  (kappa_1 = L^5/(C_0^2 C_1^3))", qmul(qmul(Is[2], qmul(C0, C0, N), N), qmul(C1, C1, N), N) == L5)
    check("I_3 = C_1  (kappa_2 = 1)", Is[3] == C1)
    check("I_4 = C_0  (kappa_3 = C_0/C_1)", Is[4] == C0)
    check("I_5 = C_0  (kappa_4 = C_0/C_1)", Is[5] == C0)
    M = min(len(gs[5]), len(gs[0])) - 1
    check("closing g_5 = g_0 through w^%d" % M, all(gs[5][m] == gs[0][m] for m in range(M + 1)))
    check("product C_0 I_1 I_2 I_3 I_4 = L^5", qmul(qmul(qmul(C0, Is[1], N), qmul(Is[2], Is[3], N), N), Is[4], N) == L5)
    # Lemma D / Lemma E ingredients in Q[[q]][ell] (ell = formal log q, D ell = 1): PF4 c_k (k = 0..4) and W(c_0..c_3) = L^10
    from math import factorial
    def lD(f):                      # f: dict ell-power -> q-series
        out = {}
        for j, a in f.items():
            out[j] = qadd(out.get(j, [Fr(0)] * (N + 1)), qD(a))
            if j:
                out[j - 1] = qadd(out.get(j - 1, [Fr(0)] * (N + 1)), qscal(a, Fr(j)))
        return out
    def lcomb(terms):               # terms: list of (q-series coefficient, dict)
        out = {}
        for c, f in terms:
            for j, a in f.items():
                out[j] = qadd(out.get(j, [Fr(0)] * (N + 1)), qmul(c, a, N))
        return out
    def lmul(f, g):
        out = {}
        for i, a in f.items():
            for j, b in g.items():
                out[i + j] = qadd(out.get(i + j, [Fr(0)] * (N + 1)), qmul(a, b, N))
        return out
    one = [Fr(1)] + [Fr(0)] * N
    qser = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)
    c = [{j: qscal(F[k - j], Fr(1, factorial(j))) for j in range(k + 1)} for k in range(5)]
    def PF4(f):
        d1 = lD(f); d2 = lD(d1); d3 = lD(d2); d4 = lD(d3)
        return lcomb([(qadd(one, qscal(qser, Fr(-3125))), d4), (qscal(qser, Fr(-6250)), d3), (qscal(qser, Fr(-4375)), d2),
                      (qscal(qser, Fr(-1250)), d1), (qscal(qser, Fr(-120)), f)])
    def is_const(f, val):
        return all((all(x == 0 for x in a) if j else (a[0] == val and all(x == 0 for x in a[1:]))) for j, a in f.items())
    check("Lemma D: PF4 c_k = 0 (k = 0..3), PF4 c_4 = 1 (in Q[[q]][ell])",
          all(is_const(PF4(c[k]), Fr(0)) for k in range(4)) and is_const(PF4(c[4]), Fr(1)))
    rows = [c[:4]]
    for j in range(3):
        rows.append([lD(f) for f in rows[-1]])
    def ldet(M):                     # Leibniz expansion, 4 x 4
        import itertools
        tot = {}
        for perm in itertools.permutations(range(4)):
            sgn = 1
            for x in range(4):
                for y in range(x + 1, 4):
                    if perm[x] > perm[y]: sgn = -sgn
            prod = {0: one}
            for r in range(4):
                prod = lmul(prod, M[r][perm[r]])
            for j, a in prod.items():
                tot[j] = qadd(tot.get(j, [Fr(0)] * (N + 1)), qscal(a, Fr(sgn)))
        return tot
    W = ldet(rows)
    L10 = [Fr((n + 1) * 3125 ** n) for n in range(N + 1)]          # (1 - 3125 q)^-2
    check("Lemma E: W(c_0..c_3) = L^10 (ell-free)", all(all(x == 0 for x in a) for j, a in W.items() if j) and W[0] == L10)
    # controls
    Fn = fhat(N, K, equivariant=False)
    Isn, gsn = birkhoff(Fn, N)
    check("control c1 (non-equivariant): I_1..I_4 agree with the zeta point", all(Isn[p] == Is[p] for p in range(1, 5)))
    Mn = min(len(gsn[5]), len(gsn[0])) - 1
    check("control c1 (non-equivariant): closing g_5 = g_0", all(gsn[5][m] == gsn[0][m] for m in range(Mn + 1)), must_fail=True)
    for dg in (4, 6):
        Fw = fhat(N, K, deg=dg)
        Isw, _ = birkhoff(Fw, N)
        C0w = Fw[0]; C1w = qadd([Fr(1)] + [Fr(0)] * N, qD(qmul(Fw[1], qinv(C0w, N), N)))
        check("control c2 (twist O(%d)): I_2 C_0^2 C_1^2 = L^5" % dg,
              qmul(qmul(Isw[2], qmul(C0w, C0w, N), N), qmul(C1w, C1w, N), N) == L5, must_fail=True)
    check("control c3 (wrong guess kappa_3 = C_1/C_0, i.e. I_4 = C_1^2/C_0)", Is[4] == qmul(qmul(C1, C1, N), qinv(C0, N), N), must_fail=True)
    log("first coefficients: C_1 = %s ; I_2 = %s" % ([str(x) for x in C1[:4]], [str(x) for x in Is[2][:4]]))
    log("TOTAL: %d unexpected result(s)" % fails)
    with open(os.path.join(HERE, LOGNAME), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("log written to " + LOGNAME)


if __name__ == "__main__":
    main()
