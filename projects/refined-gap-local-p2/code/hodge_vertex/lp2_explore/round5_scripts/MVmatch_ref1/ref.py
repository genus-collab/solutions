# Independent referee check: Schur side only; series normalised by x^{|nu|}.
from fractions import Fraction as Fr
from math import factorial as f, comb
from functools import lru_cache
from itertools import product
M = 34


def partitions(n, m=None):
    m = m or n
    if n == 0:
        yield ()
        return
    for k in range(min(n, m), 0, -1):
        for r in partitions(n - k, k):
            yield (k,) + r


def bern(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BN = bern(M + 2)


def mul(a, b):
    c = [Fr(0)] * (M + 1)
    for i, u in enumerate(a):
        if u:
            for j, v in enumerate(b):
                if i + j > M:
                    break
                if v:
                    c[i + j] += u * v
    return c


def hookser(h):
    # h x/(2 sinh(h x/2)) = sum B_{2n}(1/2) (hx)^{2n}/(2n)!, B_n(1/2)=(2^{1-n}-1)B_n
    s = [Fr(0)] * (M + 1)
    for n in range(0, M // 2 + 1):
        s[2 * n] = (Fr(2) ** (1 - 2 * n) - 1) * BN[2 * n] * Fr(h) ** (2 * n) / f(2 * n)
    return s


hs = {}


def hooks(nu):
    conj = [sum(1 for p in nu if p > i) for i in range(nu[0])] if nu else []
    return [nu[i] - j - 1 + conj[j] - i for i in range(len(nu)) for j in range(nu[i])]


def kappa(nu):
    return sum(p * (p - 2 * (i + 1) + 1) for i, p in enumerate(nu))


def W(nu, tau):
    r = [Fr(0)] * (M + 1)
    r[0] = Fr(1)
    prodh = 1
    for h in hooks(nu):
        if h not in hs:
            hs[h] = hookser(h)
        r = mul(r, hs[h])
        prodh *= h
    c = (tau + Fr(1, 2)) * kappa(nu) / 2
    e = [c ** n / f(n) for n in range(M + 1)]
    r = mul(r, e)
    return [v / prodh for v in r]


@lru_cache(None)
def chi(nu, mu):
    if not mu:
        return 1 if not nu else 0
    k = mu[0]
    n = len(nu)
    if n == 0:
        return 0
    L = n + k
    beta = [(nu[i] if i < n else 0) + L - 1 - i for i in range(L)]
    s = 0
    S = set(beta)
    for i, b in enumerate(beta):
        if b - k >= 0 and b - k not in S:
            ht = sum(1 for c in beta if b - k < c < b)
            nb = sorted([c for c in beta if c != b] + [b - k], reverse=True)
            nn = tuple(p for p in (nb[j] - (L - 1 - j) for j in range(L)) if p > 0)
            s += (-1) ** ht * chi(nn, mu[1:])
    return s


def zmu(mu):
    z = 1
    for m in set(mu):
        c = mu.count(m)
        z *= f(c) * m ** c
    return z


class Sys:
    def __init__(s, tau):
        s.tau = tau
        s.wc = {}
        s.Z = {}
        s.C = {}

    def w(s, nu):
        if nu not in s.wc:
            s.wc[nu] = W(nu, s.tau)
        return s.wc[nu]

    def Zc(s, mu):
        mu = tuple(sorted(mu, reverse=True))
        if mu not in s.Z:
            n = sum(mu)
            r = [Fr(0)] * (M + 1)
            if n == 0:
                r[0] = Fr(1)
            for nu in partitions(n):
                if n == 0:
                    break
                c = chi(nu, mu)
                if c:
                    wn = s.w(nu)
                    r = [a + Fr(c, zmu(mu)) * b for a, b in zip(r, wn)]
            s.Z[mu] = r
        return s.Z[mu]

    def Cc(s, mu):
        mu = tuple(sorted(mu, reverse=True))
        if mu in s.C:
            return s.C[mu]
        k = mu[0]
        mk = mu.count(k)
        vals = sorted(set(mu))
        cnt = [mu.count(v) for v in vals]
        r = [Fr(0)] * (M + 1)
        for e in product(*[range(c + 1) for c in cnt]):
            sub = []
            for v, c in zip(vals, e):
                sub += [v] * c
            sub = tuple(sorted(sub, reverse=True))
            mk1 = sub.count(k)
            if mk1 == 0 or sub == mu:
                continue
            comp = list(mu)
            for v in sub:
                comp.remove(v)
            t = mul(s.Cc(sub), s.Zc(comp))
            r = [a + mk1 * b for a, b in zip(r, t)]
        zz = s.Zc(mu)
        res = [(a * mk - b) / mk for a, b in zip(zz, r)]
        s.C[mu] = res
        return res

    def coef(s, mu, n):
        i = n + sum(mu)
        if i < 0:
            return Fr(0)
        return s.Cc(mu)[i]


def lag(n, al, X):
    return sum((-1) ** i * comb(n + al, n - i) * X ** i / Fr(f(i)) for i in range(n + 1))


def Pi(m, tau):
    r = Fr(1)
    for a in range(1, m):
        r *= (m * tau + a)
    return r / f(m - 1)


def Maut(R):
    z = 1
    for m in set(R):
        z *= f(R.count(m))
    return z


def Xi(S, h, R, tau):
    sg = tau * (1 + tau)
    r = len(R)
    J = sum(R)
    N = 3 * h - 3 + r
    pre = Fr(Maut(R)) / (sg ** (h + r - 1))
    for j in R:
        pre /= Pi(j, tau)
    tot = Fr(0)
    for a in range(0, N + 1):
        mu = tuple(R) + (1,) * a
        tot += (-1) ** a * sg ** (-a) * lag(N - a, 2 * h - 2 + r + a, Fr(-J - a)) * S.coef(mu, 2 * h - 2 + r + a)
    return pre * tot


tau = Fr(2, 7)
sg = tau * (1 + tau)
S = Sys(tau)
print("[p1^2] x^0:", S.coef((1, 1), 0), "sigma/4=", sg / 4)
print("[p2]x^-1:", S.coef((2,), -1), "(2t+1)/4=", (2 * tau + 1) / 4)
print("forbidden powers (1,1):", [S.coef((1, 1), n) for n in (-3, -2, -1, 1)], " (2,2):", [S.coef((2, 2), n) for n in (-3, -2, -1, 1)])


def wk1(ks):
    ks = sorted(ks)
    n = len(ks)
    if sum(ks) != n:
        return Fr(0)
    if ks == [1]:
        return Fr(1, 24)
    if ks[0] == 0:
        rest = ks[1:]
        t = Fr(0)
        for i in range(len(rest)):
            if rest[i] > 0:
                r = rest[:]
                r[i] -= 1
                t += wk1(r)
        return t
    return (n - 1) * wk1(ks[1:])


def comps(tot, n):
    if n == 1:
        yield (tot,)
        return
    for i in range(tot + 1):
        for r in comps(tot - i, n - 1):
            yield (i,) + r


def M1(mu):
    l = len(mu)
    t = Fr(0)
    for k in comps(l, l):
        c = Fr(1)
        for m, e in zip(mu, k):
            c *= m ** e
        t += c * wk1(list(k))
    u = Fr(0)
    for k in comps(l - 1, l):
        c = Fr(1)
        for m, e in zip(mu, k):
            c *= m ** e
        mn = f(l - 1)
        for e in k:
            mn //= f(e)
        u += c * Fr(mn, 24)
    return t - (1 + 1 / sg) * u


ok = True
for mu in [(1,), (2,), (3,), (1, 1), (2, 1), (2, 2), (3, 1), (1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 2, 2), (1, 1, 1, 1)]:
    l = len(mu)
    pm = Fr(1)
    for j in mu:
        pm *= Pi(j, tau)
    rhs = sg ** (1 + l - 1) * pm * M1(mu) / Maut(mu)
    lhs = S.coef(mu, l)
    ok &= lhs == rhs
    print(mu, lhs == rhs)
print("g1 (A) all:", ok)

A2 = Fr(1, 5760)
Ma = {a: S.coef((1,) * a, 2 + a) * f(a) / sg ** (1 + a) for a in range(1, 9)}


def sco(m, M0):
    t = lag(m, 1, Fr(0)) * M0 * (-1) ** m
    for a in range(1, m + 1):
        t += Ma[a] / f(a) * (-1) ** (m - a) * lag(m - a, 1 + a, Fr(-a))
    return t


for m in (4, 5, 6, 7):
    v0 = sco(m, Fr(0))
    v1 = sco(m, Fr(1))
    print("m", m, "M0 solved", -v0 / (v1 - v0), "pred", -A2 / sg)
M0 = -A2 / sg
Qc = [sco(m, M0) for m in range(0, 4)]
print("Q_2(s;empty) coeffs", Qc, "value at -1", sum((-1) ** i * c for i, c in enumerate(Qc)))
print("(1+s)^2 QQ_2(s;empty) at s=-1 =", -M0)
print("Xi_2(2)", Xi(S, 2, [2], tau), " Xi_2(3)", Xi(S, 2, [3], tau))
print("Xi_2(2,2)", Xi(S, 2, [2, 2], tau))
print("Xi_3(2)", Xi(S, 3, [2], tau), " Xi_3(3)", Xi(S, 3, [3], tau), " Xi_3(4)", Xi(S, 3, [4], tau), " Xi_3(5)", Xi(S, 3, [5], tau))
print("Xi_3(2,2)", Xi(S, 3, [2, 2], tau), "Xi_3(3,2)", Xi(S, 3, [3, 2], tau))
