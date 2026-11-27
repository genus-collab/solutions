# -*- coding: utf-8 -*-
"""lit_vertex_check.py -- literature consistency check of the twist-weight sign (Stage 5, task lit).

Computes the degree-0 formal-quintic invariants N_{g,0}, g = 2, 3, directly from the vertex class of
Lho-Pandharipande arXiv:1803.01409 eq. (26),
    H_g^{p_i} = e(E^v (x) T_{p_i}) / e(T_{p_i})  *  (5 l_i) / e(E^v (x) (5 l_i)),
with tangent weights l_i - l_j (j != i) (their s.1.1) and the specialisation l_i = zeta^i (their s.0.2),
N_{g,0} = Sum_i Int_{Mbar_g} H_g^{p_i}.  Hodge integrals: Mumford's relation c(E)c(E^v) = 1 reduces every
top-degree lambda monomial to a multiple of lambda_1^(3g-3); Int_{Mbar_g} lambda_{g-1}^3 from
Faber-Pandharipande arXiv:math/9810173 Theorem 4 (g = 2: 1/2880, g = 3: 1/725760).
Checks:
  (A) H_g^{p_i} = (W_i)^(g-1) * Prod_{a in w_i} Lambda^v(a)/a^g  as classes (all coefficients), with
      W_i = Prod_{j!=i}(l_i - l_j)/(5 l_i) and w_i = {l_i - l_j} U {-5 l_i}   (the project's convention);
  (B) N_{2,0} = -19/9600 and N_{3,0} = 1/224000 (the project's localisation values, FORMAL_QUINTIC_genus2.md s.2);
  (C) control: the twist weight +5 l_i in place of -5 l_i changes N_{2,0} or N_{3,0}.
Exact arithmetic in Q(zeta_5).  Usage: py -I -B lit_vertex_check.py
"""
from fractions import Fraction as Fr
from itertools import product


# ---------------- Q(zeta_5): coefficient lists on 1, z, z^2, z^3 with z^4 = -1 - z - z^2 - z^3
def zn(c): return [Fr(c), Fr(0), Fr(0), Fr(0)]
def zadd(a, b): return [x + y for x, y in zip(a, b)]
def zneg(a): return [-x for x in a]
def zmul(a, b):
    t = [Fr(0)] * 7
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): t[i + j] += x * y
    for k in (6, 5, 4):                                  # z^k = z^(k-4) z^4 = -z^(k-4)(1 + z + z^2 + z^3)
        c = t[k]; t[k] = Fr(0)
        if c:
            for m in range(4): t[k - 4 + m] -= c
    return t[:4]
def zpow_root(e):                                        # zeta^e
    e %= 5
    if e == 4: return [Fr(-1)] * 4
    v = [Fr(0)] * 4; v[e] = Fr(1); return v
def zconj(a, k):                                         # Galois automorphism zeta -> zeta^k
    out = zn(0)
    for i, x in enumerate(a):
        if x: out = zadd(out, [x * y for y in zpow_root(i * k)])
    return out
def zinv(a):
    c = zmul(zmul(zconj(a, 2), zconj(a, 3)), zconj(a, 4))
    nrm = zmul(a, c)
    assert all(x == 0 for x in nrm[1:]) and nrm[0] != 0
    return [x / nrm[0] for x in c]
def zscal(a, f): return [x * f for x in a]
def zrat(a):
    assert all(x == 0 for x in a[1:]), a
    return a[0]


# ---------------- graded polynomials in lambda_1..lambda_g over Q(zeta), truncated at degree 3g-3
class Hodge:
    def __init__(self, g): self.g = g; self.top = 3 * g - 3
    def deg(self, m): return sum((k + 1) * e for k, e in enumerate(m))
    def one(self): return {tuple([0] * self.g): zn(1)}
    def const(self, c): return {tuple([0] * self.g): c}
    def lam(self, k, c):                                 # c * lambda_k (lambda_0 = 1)
        if k == 0: return self.const(c)
        m = [0] * self.g; m[k - 1] = 1; return {tuple(m): c}
    def add(self, P, Q):
        R = dict(P)
        for m, c in Q.items(): R[m] = zadd(R[m], c) if m in R else c
        return R
    def mul(self, P, Q):
        R = {}
        for (m1, c1), (m2, c2) in product(P.items(), Q.items()):
            m = tuple(x + y for x, y in zip(m1, m2))
            if self.deg(m) > self.top: continue
            c = zmul(c1, c2); R[m] = zadd(R[m], c) if m in R else c
        return R
    def lamv(self, a):                                   # Lambda^v(a) = Sum_k (-1)^k lambda_k a^(g-k)
        P = {}
        for k in range(self.g + 1):
            c = zn(1)
            for _ in range(self.g - k): c = zmul(c, a)
            P = self.add(P, self.lam(k, zscal(c, (-1) ** k)))
        return P
    def inv(self, P):                                    # P = c0 (1 + u), u nilpotent
        c0 = P[tuple([0] * self.g)]; ic = zinv(c0)
        u = {m: zmul(c, ic) for m, c in P.items() if self.deg(m) > 0}
        S = self.one(); T = self.one()
        for _ in range(self.top):
            T = self.mul(T, {m: zneg(c) for m, c in u.items()}); S = self.add(S, T)
        return {m: zmul(c, ic) for m, c in S.items()}
    def scale(self, P, c): return {m: zmul(v, c) for m, v in P.items()}
    def normal(self, P):
        """normal form modulo Mumford's relations c(E)c(E^v) = 1: g = 2: lambda_2 = lambda_1^2/2;
        g = 3: lambda_2 = lambda_1^2/2, lambda_1 lambda_3 = lambda_1^4/8, lambda_3^2 = 0 (basis lambda_1^a, lambda_3)"""
        R = {}
        for m, c in P.items():
            e1 = m[0] + 2 * m[1]; f = Fr(1, 2) ** m[1]; key = ('l1', e1)
            if self.g == 3 and m[2] >= 2: continue
            if self.g == 3 and m[2] == 1:
                if e1 == 0: key = ('l3', 0)
                else: key = ('l1', e1 + 3); f = f * Fr(1, 8)
            R[key] = zadd(R[key], zscal(c, f)) if key in R else zscal(c, f)
        return R
    def integral(self, P):
        """Int over Mbar_g of the top-degree part; lambda monomials reduced by Mumford's relation"""
        g = self.g
        if g == 2:   red = lambda m: Fr(1, 2) ** m[1]                       # lambda_2 = lambda_1^2/2
        elif g == 3: red = lambda m: Fr(0) if m[2] >= 2 else Fr(1, 2) ** m[1] * Fr(1, 8) ** m[2]
        else: raise ValueError
        # g = 3: lambda_2 = lambda_1^2/2, lambda_1 lambda_3 = lambda_1^4/8, lambda_3^2 = 0 (every monomial with a lambda_3
        # and total degree 6 contains lambda_1 after substituting lambda_2, except lambda_3^2)
        top = {2: Fr(1, 2880), 3: Fr(1, 725760) * 8}[g]                     # Int lambda_1^(3g-3) from FP Thm 4
        tot = zn(0)
        for m, c in P.items():
            if self.deg(m) == self.top: tot = zadd(tot, zscal(c, red(m) * top))
        return tot


def vertex_LP(H, i, twist_sign=+1):
    """LP (26) at p_i, l_j = zeta^j; twist_sign = -1 replaces 5 l_i/e(E^v (x) 5 l_i) by the control class"""
    li = zpow_root(i); P = H.one()
    for j in range(5):
        if j == i: continue
        a = zadd(li, zneg(zpow_root(j)))
        P = H.mul(P, H.scale(H.lamv(a), zinv(a)))
    w = zscal(li, Fr(5))
    if twist_sign > 0: T = H.scale(H.inv(H.lamv(w)), w)                    # (5 l_i) / Lambda^v(5 l_i)
    else:                                                                # control: w^(1-g) Lambda^v(w)/w^g
        wg = zn(1)
        for _ in range(2 * H.g - 1): wg = zmul(wg, w)
        T = H.scale(H.lamv(w), zinv(wg))
    return H.mul(P, T)


def vertex_project(H, i, wsign=-1):
    """(W_i)^(g-1) Prod_{a in w_i} Lambda^v(a)/a^g, w_i = {l_i - l_j} U {wsign * 5 l_i}"""
    li = zpow_root(i); ws = [zadd(li, zneg(zpow_root(j))) for j in range(5) if j != i] + [zscal(li, Fr(5 * wsign))]
    W = zn(1)
    for a in ws[:4]: W = zmul(W, a)
    W = zmul(W, zinv(zscal(li, Fr(5))))
    P = H.one()
    for a in ws:
        ag = zn(1)
        for _ in range(H.g): ag = zmul(ag, a)
        P = H.mul(P, H.scale(H.lamv(a), zinv(ag)))
    Wp = zn(1)
    for _ in range(H.g - 1): Wp = zmul(Wp, W)
    return H.scale(P, Wp)


def main():
    target = {2: Fr(-19, 9600), 3: Fr(1, 224000)}
    ok = True
    for g in (2, 3):
        H = Hodge(g)
        def eq(P, Q): return all(zadd(P.get(m, zn(0)), zneg(Q.get(m, zn(0)))) == zn(0) for m in set(P) | set(Q))
        raw = all(eq(vertex_LP(H, i), vertex_project(H, i)) for i in range(5))
        same = all(eq(H.normal(vertex_LP(H, i)), H.normal(vertex_project(H, i))) for i in range(5))
        print("g = %d: (A0) the two classes as polynomials in FREE lambda's: %s (expected: differ; they agree only modulo"
              " Mumford's relation)" % (g, "equal" if raw else "differ"))
        N = zrat([sum(x) for x in zip(*[H.integral(vertex_LP(H, i)) for i in range(5)])])
        Nc = zrat([sum(x) for x in zip(*[H.integral(vertex_LP(H, i, -1)) for i in range(5)])])
        Np = zrat([sum(x) for x in zip(*[H.integral(vertex_project(H, i, +1)) for i in range(5)])])
        print("g = %d: (A) LP (26) class = project class (W_i)^(g-1) Hodge_(w_i) modulo Mumford, all degrees, i = 0..4: %s"
              % (g, "PASS" if same else "FAIL"))
        print("g = %d: (B) N_(g,0) from LP (26) = %s; project value %s: %s" % (g, N, target[g], "PASS" if N == target[g] else "FAIL"))
        print("g = %d: (C) control, twist class w^(1-g) Lambda^v(w)/w^g: %s; project formula with weight +5 l_i: %s"
              " (control must differ from %s: %s)" % (g, Nc, Np, target[g], "FAIL as expected" if Nc != target[g] else "DID NOT FAIL"))
        ok = ok and same and N == target[g]
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")


if __name__ == '__main__':
    main()
