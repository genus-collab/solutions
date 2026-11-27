# -*- coding: utf-8 -*-
"""r3_uadic.py -- referee check (independent code) of the u-adic statements of lem:tra-frames (b)-(d) and lem:tra-conv (a):
the decomposition q = Sum_(k, a != 1) q_(k,a) frame(k,a) + sigma(q)(z) H for INFINITE q in HH^u_+ = H (x) Q[[u, z]].

v(u) = u tau^(-5/4), 1/tau = 2F1(1, 4/5; 9/5; u): the hypergeometric coefficients are taken from the Pochhammer definition and
F^(5/4) from the ODE F g' = (5/4) F' g (not from a binomial series).  Frames: graded frames at v = 1 by Gaussian elimination
(as in r3_frames.py), substituted v -> v(u) with the exponent (4(k-j)+nu_a+1)/5 (part (d) of the prover's statement).
Tests, all mod u^N:
 (U1) finite q (z-degree K, random coefficients in Q[[u]]): the top-down solution s of Pi_i(q - s(z)H) = 0, 0 <= i <= K, equals
      sigma_j = q_(j,1) - Sum_(k > j, a != 1) q_(k,a) f_j^(k,a)(v(u)) (frame coordinates = flat components with a != 1);
 (U2) INFINITE q (q_(k,b) = (1 + k u)(b + 1) for all k): sigma_j is computed from the frame formula with the k-sum cut where the
      u-order of f_j^(k,a) reaches N; then Pi_i(q) - (Theta sigma)_i = 0 mod u^N for 0 <= i <= I, where Pi_i(q) uses the infinite
      q (only finitely many terms are non-zero mod u^N) and (Theta sigma)_i = Sum_(j >= i) sigma_j Pi_0(z^(j-i) H);
 (U3) the cut used in (U2) is sharp: the k-th frame changes [u^m] sigma_j only for k - j <= (5m - 1)/4;
 (C)  CONTROL v = 1 + u (v(0) != 0): the partial sums over k of the [u^0] part of sigma_0 for the same infinite q change at
      infinitely many k (no u-adic limit), and Pi_0 of the infinite q has no u-adic value.
Usage (working directory: the scratch folder):  py -I -B r3_uadic.py N I   -- log on stdout, nothing written, no path printed.
"""
import sys
from fractions import Fraction as Fr
from math import factorial

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)
AS = (0, 2, 3, 4)


# ---- truncated series mod u^N -------------------------------------------------------------------------------------
def smul(a, b):
    N = len(a)
    out = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]:
                    out[i + j] += x * b[j]
    return out


def sadd(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]


def sconst(c, N):
    return [Fr(c)] + [Fr(0)] * (N - 1)


def ordu(a):
    return next((i for i, x in enumerate(a) if x), None)


def model_v(N):
    """v(u) = u * F^(5/4), F = 2F1(1,4/5;9/5;u) from the Pochhammer quotient; returns (v, check of F_n = 4/(5n+4))."""
    F = [Fr(1)]
    for n in range(N):
        F.append(F[-1] * (Fr(4, 5) + n) * (1 + n) / ((Fr(9, 5) + n) * (n + 1)))
    okF = all(F[n] == Fr(4, 5 * n + 4) for n in range(N + 1))
    g = [Fr(1)]
    for n in range(N - 1):
        rhs = Fr(5, 4) * sum((i + 1) * F[i + 1] * g[n - i] for i in range(n + 1))
        rhs -= sum(F[i] * (n - i + 1) * g[n - i + 1] for i in range(1, n + 1))
        g.append(rhs / (n + 1))
    return [Fr(0)] + g[:N - 1], okF


class VP:
    def __init__(self, v):
        self.v, self.N, self.cache = v, len(v), {0: sconst(1, len(v))}

    def __call__(self, e):
        if e not in self.cache:
            self.cache[e] = smul(self(e - 1), self.v)
        return self.cache[e]


# ---- graded data at v = 1 (independent of the prover's code) -------------------------------------------------------
def yc1(b, n):
    """[z^-n] Y^C_b at v = 1 from the residue double sum; v-power returned as well."""
    m, tot, vp = NU[b], Fr(0), None
    for i in range(n + 1):
        l = m + 4 * i + 1
        if l >= 0 and i + l == n:
            tot += Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** l, factorial(l))
            vp = l
    return EPS[b] * tot, vp


def pi0_1(m, b):
    """Pi_0(z^m H^b) = (-1)^(m+1) [z^-m] Y^C_b: (value at v = 1, v-power)."""
    c, vp = yc1(b, m)
    return (-1) ** (m + 1) * c, vp


def gauss(M, rhs):
    n = len(M)
    A = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fac = A[r][c] / A[c][c]
                A[r] = [x - fac * y for x, y in zip(A[r], A[c])]
    return [A[r][n] / A[r][r] for r in range(n)]


def frame1(k, a):
    if k == 0:
        return {}
    M = [[pi0_1(j - i, 1)[0] if j >= i else Fr(0) for j in range(k)] for i in range(k)]
    f = gauss(M, [-pi0_1(k - i, a)[0] for i in range(k)])
    return {j: c for j, c in enumerate(f) if c != 0}


# ---- u-adic objects ------------------------------------------------------------------------------------------------
def pi_u(i, qcoef, kmax, vp, N):
    """Pi_i(q) mod u^N for q given by qcoef(k, b), summed over i <= k <= kmax."""
    out = [Fr(0)] * N
    for k in range(i, kmax + 1):
        for b in range(5):
            c, e = pi0_1(k - i, b)
            if c:
                out = sadd(out, smul(qcoef(k, b), vp(e)), c)
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    I = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    print("r3_uadic.py  N = %d (mod u^N), I = %d; v(u) = u tau^(-5/4); exact over Q" % (N, I))
    v, okF = model_v(N + 2)
    v = v[:N]
    vp = VP(v)
    print("[V] 2F1 coefficients from Pochhammer = 4/(5n+4): %s;  v = %s u + %s u^2 + %s u^3 + .." % (okF, v[1], v[2], v[3]))
    ok = okF and v[1] == 1 and v[2] == Fr(5, 9)
    JMAX = I + (5 * N) // 4 + 2                              # sigma_j needed for j - i < 5N/4
    KMAX = JMAX + (5 * N) // 4 + 2                           # frames needed: k - j <= (5N - 1)/4
    fr = {(k, a): frame1(k, a) for k in range(KMAX + 1) for a in AS}
    fu = {}
    for (k, a), f in fr.items():
        for j, c in f.items():
            e = (4 * (k - j) + NU[a] + 1) // 5
            fu[(k, a, j)] = [c * x for x in vp(e)]
    # (U1) finite q: top-down solve versus the frame formula
    K = 25
    q1 = {(k, b): [Fr((7 * k + 3 * b + i) % 11 - 5) for i in range(N)] for k in range(K + 1) for b in range(5)}
    t = [pi_u(i, lambda kk, bb: q1[(kk, bb)], K, vp, N) for i in range(K + 1)]
    s = [None] * (K + 1)
    for i in range(K, -1, -1):
        acc = t[i][:]
        for qq in range(1, (K - i) // 5 + 1):
            c, e = pi0_1(5 * qq, 1)
            acc = sadd(acc, smul(vp(e), s[i + 5 * qq]), -c)
        s[i] = acc
    sig = []
    for j in range(K + 1):
        acc = q1[(j, 1)][:]
        for k in range(j + 1, K + 1):
            for a in AS:
                if (k, a, j) in fu:
                    acc = sadd(acc, smul(q1[(k, a)], fu[(k, a, j)]), -1)
        sig.append(acc)
    u1 = (s == sig)
    print("[U1] finite q (z-degree %d): top-down solution of Pi_i(q - s H) = 0 equals the frame formula on all %d entries: %s"
          % (K, K + 1, u1))
    ok &= u1

    # (U2) infinite q
    def qinf(k, b):
        return [Fr(b + 1), Fr((b + 1) * k)] + [Fr(0)] * (N - 2)
    cnt = [0]

    def sigma_list(cut):
        """sigma_j, j <= JMAX, from the frame formula with the k-sum restricted to k - j <= cut."""
        out = []
        for j in range(JMAX + 1):
            acc = qinf(j, 1)
            for k in range(j + 1, min(KMAX, j + cut) + 1):
                for a in AS:
                    if (k, a, j) in fu:
                        term = smul(qinf(k, a), fu[(k, a, j)])
                        if ordu(term) is not None and 4 * (k - j) > 5 * N - 1:
                            cnt[0] += 1                       # a frame beyond the predicted cut still visible mod u^N
                        acc = sadd(acc, term, -1)
            out.append(acc)
        return out

    def theta_fail(sg):
        bad = 0
        for i in range(I + 1):
            lhs = pi_u(i, qinf, i + (5 * N) // 4 + 2, vp, N)
            rhs = [Fr(0)] * N
            for j in range(i, JMAX + 1):
                c, e = pi0_1(j - i, 1)
                if c:
                    rhs = sadd(rhs, smul(sg[j], vp(e)), c)
            bad += lhs != rhs
        return bad
    bad2 = theta_fail(sigma_list(KMAX))
    print("[U2] infinite q: Pi_i(q) = (Theta sigma)_i mod u^%d for 0 <= i <= %d: failures %d;  frames beyond the cut visible: %d"
          % (N, I, bad2, cnt[0]))
    ok &= bad2 == 0 and cnt[0] == 0
    cut = max(1, (5 * N - 1) // 4 - 2)                       # (U2-C) failing control: the k-sum cut too early
    badc = theta_fail(sigma_list(cut))
    print("[U2-C] CONTROL, frame sum cut at k - j <= %d: Pi_i(q) != (Theta sigma)_i for %d of %d values of i (expected > 0)"
          % (cut, badc, I + 1))
    ok &= badc > 0
    # (U3) sharpness of the cut: the last k that changes [u^m] sigma_j
    rows = []
    for (j, m) in ((0, 1), (0, N - 1), (3, N - 1), (7, N - 2)):
        last = j
        for k in range(j + 1, KMAX + 1):
            for a in AS:
                if (k, a, j) in fu and smul(qinf(k, a), fu[(k, a, j)])[m] != 0:
                    last = max(last, k)
        rows.append((j, m, last, j + (5 * m - 1) // 4))
    u3 = all(last <= bound for (_, _, last, bound) in rows)
    print("[U3] last frame index k changing [u^m] sigma_j, (j, m, last, bound j + floor((5m-1)/4)): %s: %s" % (rows, u3))
    ok &= u3
    # (C) control v = 1 + u: at u^0 every v-power is 1, so the u^0 parts are plain infinite sums of rationals
    KC = 60
    frc = {(k, a): (fr[(k, a)] if (k, a) in fr else frame1(k, a)) for k in range(KC + 1) for a in AS}
    changes, ps = 0, Fr(0)
    for k in range(1, KC + 1):
        for a in AS:
            c0 = frc[(k, a)].get(0, Fr(0))
            if c0:
                ps -= qinf(k, a)[0] * c0
                changes += 1
    pc, pchanges = Fr(0), 0
    for m in range(0, KC + 1):
        for b in range(5):
            c, e = pi0_1(m, b)
            if c:
                pc += qinf(m, b)[0] * c
                pchanges += 1
    print("[C] CONTROL v = 1 + u: the u^0 part of sigma_0 (frame formula, k <= %d) changes at %d values of (k, a); the u^0 part of "
          "Pi_0(q) changes at %d values of m <= %d (expected: both grow with the cut: no u-adic limit)" % (KC, changes, pchanges, KC))
    ok &= changes >= 10 and pchanges >= 10
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


if __name__ == "__main__":
    main()
