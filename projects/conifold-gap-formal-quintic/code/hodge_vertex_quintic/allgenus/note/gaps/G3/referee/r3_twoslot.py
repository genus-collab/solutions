# -*- coding: utf-8 -*-
"""r3_twoslot.py -- referee check (independent code) of G3-2g: a bivector E in H (x) H (x) Q[[u, psi, psi']] annihilated by every
Pi_i in each slot equals Sum E^(alpha beta) frame_alpha(psi) (x) frame_beta(psi'), E^(alpha beta) its flat components with a, b != 1.

E is built WITHOUT frames: E := (pr (x) pr)(B) for a random bivector B of bidegree <= K, where pr(q) = q - s(z) H and s solves
Pi_i(q - s(z) H) = 0 (0 <= i <= K) from the top (finite in z); v = v(u) = u tau^(-5/4), everything mod u^N.  Frames: Gaussian
elimination at v = 1, substituted v -> v(u) with the exponent (4(k-j)+nu_a+1)/5.
Checks: (A) (Pi_i (x) 1) E = 0 = (1 (x) Pi_i) E, 0 <= i <= K + 3; (R) E = the frame double sum of its flat components with a, b != 1;
(S) E is symmetric when B is.  CONTROL: E' := (pr (x) 1)(B) (second slot not projected): (1 (x) Pi_i) E' != 0 and E' differs from the
frame double sum of its own flat components.
Usage (working directory: the scratch folder):  py -I -B r3_twoslot.py N K   -- log on stdout, nothing written, no path printed.
"""
import sys
import random
from fractions import Fraction as Fr
from math import factorial

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)
AS = (0, 2, 3, 4)


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


def model_v(N):
    F = [Fr(4, 5 * n + 4) for n in range(N + 1)]
    g = [Fr(1)]
    for n in range(N - 1):
        rhs = Fr(5, 4) * sum((i + 1) * F[i + 1] * g[n - i] for i in range(n + 1))
        rhs -= sum(F[i] * (n - i + 1) * g[n - i + 1] for i in range(1, n + 1))
        g.append(rhs / (n + 1))
    return [Fr(0)] + g[:N - 1]


def pi0_1(m, b):
    """Pi_0(z^m H^b) = (-1)^(m+1) [z^-m] Y^C_b at v = 1 (residue double sum) and its v-power."""
    if m < 0:
        return Fr(0), None
    nu, tot, vp = NU[b], Fr(0), None
    for i in range(m + 1):
        l = nu + 4 * i + 1
        if l >= 0 and i + l == m:
            tot += Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** l, factorial(l))
            vp = l
    return (-1) ** (m + 1) * EPS[b] * tot, vp


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    rng = random.Random(4242)
    print("r3_twoslot.py  N = %d (mod u^N), bidegree K = %d; v(u) = u tau^(-5/4); exact over Q" % (N, K))
    v = model_v(N)
    pw = [[Fr(1)] + [Fr(0)] * (N - 1)]
    for _ in range(K + 2):
        pw.append(smul(pw[-1], v))

    def Pi(i, q):                      # q: {(k, b): series}, finite
        out = [Fr(0)] * N
        for (k, b), c in q.items():
            p, e = pi0_1(k - i, b)
            if p:
                out = sadd(out, smul(c, pw[e]), p)
        return out

    def pr(q):
        s = {}
        for i in range(K, -1, -1):
            acc = Pi(i, q)
            for j, sj in s.items():
                p, e = pi0_1(j - i, 1)
                if p:
                    acc = sadd(acc, smul(sj, pw[e]), -p)
            s[i] = acc
        out = {kb: c[:] for kb, c in q.items()}
        for j, sj in s.items():
            out[(j, 1)] = sadd(out.get((j, 1), [Fr(0)] * N), sj, -1)
        return out

    def frame(k, a):                   # Gaussian elimination at v = 1, then v -> v(u)
        M = [[pi0_1(j - i, 1)[0] for j in range(k)] for i in range(k)]
        rhs = [-pi0_1(k - i, a)[0] for i in range(k)]
        A = [row[:] + [r] for row, r in zip(M, rhs)]
        for c in range(k):
            p = next(r for r in range(c, k) if A[r][c] != 0)
            A[c], A[p] = A[p], A[c]
            for r in range(k):
                if r != c and A[r][c] != 0:
                    fac = A[r][c] / A[c][c]
                    A[r] = [x - fac * y for x, y in zip(A[r], A[c])]
        out = {(k, a): [Fr(1)] + [Fr(0)] * (N - 1)}
        for j in range(k):
            f = A[j][k] / A[j][j]
            if f:
                out[(j, 1)] = [f * x for x in pw[(4 * (k - j) + NU[a] + 1) // 5]]
        return out
    fr = {(k, a): frame(k, a) for k in range(K + 1) for a in AS}
    keys = [(k, b) for k in range(K + 1) for b in range(5)]
    Bc = {}
    for i, x in enumerate(keys):
        for y in keys[i:]:
            Bc[(x, y)] = [Fr(rng.randint(-3, 3)) for _ in range(N)]
            Bc[(y, x)] = Bc[(x, y)]

    def slot1(E, op):                  # apply op to the psi-slot for every fixed psi'-index
        out = {}
        for y in keys:
            col = {x: E[(x, y)] for x in keys if (x, y) in E}
            for x, c in op(col).items():
                out[(x, y)] = c
        return out

    def swap(E):
        return {(y, x): c for (x, y), c in E.items()}
    E = swap(slot1(swap(slot1(Bc, pr)), pr))
    E1 = slot1(Bc, pr)                 # control: second slot not projected

    def ann(E, first):
        bad = 0
        for i in range(K + 4):
            for y in keys:
                col = {x: E[(x, y)] for x in keys if (x, y) in E} if first else {x: E[(y, x)] for x in keys if (y, x) in E}
                bad += any(Pi(i, col))
        return bad

    def recon(E):
        out = {}
        for (al, be), c in E.items():
            if al[1] == 1 or be[1] == 1 or not any(c):
                continue
            for x, fa in fr[al].items():
                for y, fb in fr[be].items():
                    out[(x, y)] = sadd(out.get((x, y), [Fr(0)] * N), smul(c, smul(fa, fb)))
        return out

    def same(E, F):
        ks = set(E) | set(F)
        return all(E.get(k, [Fr(0)] * N) == F.get(k, [Fr(0)] * N) for k in ks)
    a1, a2 = ann(E, True), ann(E, False)
    print("[A] (Pi_i (x) 1) E and (1 (x) Pi_i) E, 0 <= i <= %d: non-zero columns %d, %d" % (K + 3, a1, a2))
    rc = same(E, recon(E))
    sym = same(E, swap(E))
    print("[R] E = Sum E^(alpha beta) frame_alpha (x) frame_beta (flat components with a, b != 1): %s;  [S] E symmetric: %s" % (rc, sym))
    c1, c2 = ann(E1, False), same(E1, recon(E1))
    print("[C] CONTROL E' = (pr (x) 1) B: (1 (x) Pi_i) E' non-zero columns %d (expected > 0); E' = its frame double sum: %s "
          "(expected False)" % (c1, c2))
    ok = a1 == 0 and a2 == 0 and rc and sym and c1 > 0 and not c2
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


if __name__ == "__main__":
    main()
