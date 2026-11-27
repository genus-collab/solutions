# Exact check (evidence only) of the formulas exactly as written in LEMMA_A_PROOF.md:
# Prop 2.1 on the Schur side -> M_g(mu); Lemma 3.2 (Laguerre form) -> X_g(U); (5.1) = Theorem M_k; target coefficient.
# Fresh code: imports nothing from the project.
import sys
from collections import Counter
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial, comb
from itertools import product

tau = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 7)
sigma = tau * (1 + tau)
PREC = 44


def Pi(k):
    r = Fr(1)
    for a in range(1, k):
        r *= (k * tau + a)
    return r / factorial(k - 1)


def S_series(n):
    s = [Fr(0)] * (n + 1)
    for i in range(0, n + 1, 2):
        s[i] = Fr(1, factorial(i + 1) * 2 ** i)
    inv = [Fr(0)] * (n + 1)
    inv[0] = Fr(1)
    for m in range(1, n + 1):
        inv[m] = -sum(s[i] * inv[m - i] for i in range(1, m + 1))
    return inv


SS = S_series(PREC)


def mul(a, b):
    r = [Fr(0)] * (PREC + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j > PREC:
                break
            if y:
                r[i + j] += x * y
    return r


def parts(n, m=None):
    m = n if m is None else m
    if n == 0:
        yield ()
        return
    for k in range(min(n, m), 0, -1):
        for p in parts(n - k, k):
            yield (k,) + p


@lru_cache(None)
def chi(nu, mu):
    if not mu:
        return 1 if not nu else 0
    k = mu[0]
    rest = mu[1:]
    l = len(nu)
    b = [nu[i] + l - 1 - i for i in range(l)]
    tot = 0
    S = set(b)
    for bi in b:
        c = bi - k
        if c < 0 or c in S:
            continue
        ht = sum(1 for x in b if c < x < bi)
        nb = sorted([x for x in b if x != bi] + [c], reverse=True)
        n2 = tuple(v for v in (nb[j] - (l - 1 - j) for j in range(l)) if v > 0)
        tot += (-1) ** ht * chi(n2, rest)
    return tot


def hooks(nu):
    conj = [sum(1 for x in nu if x > j) for j in range(nu[0])] if nu else []
    return [nu[i] - j + conj[j] - i - 1 for i in range(len(nu)) for j in range(nu[i])]


def kappa(nu):
    return sum(nu[i] * (nu[i] - 2 * (i + 1) + 1) for i in range(len(nu)))


@lru_cache(None)
def Znu(nu):
    # x^{|nu|} e^{(tau+1/2) kappa x/2} prod 1/(2 sinh(h x/2)) as a power series
    c = (tau + Fr(1, 2)) * kappa(nu) / 2
    r = [c ** i / factorial(i) for i in range(PREC + 1)]
    for h in hooks(nu):
        sh = [SS[i] * Fr(h) ** i / h for i in range(PREC + 1)]
        r = mul(r, sh)
    return tuple(r)


def norm(mu):
    return tuple(sorted(mu, reverse=True))


@lru_cache(None)
def zz(mu):
    # x^{|mu|} * |Aut mu| [p_mu] Z
    n = sum(mu)
    r = [Fr(0)] * (PREC + 1)
    den = 1
    for k in mu:
        den *= k
    for nu in parts(n):
        c = chi(nu, mu)
        if c:
            z = Znu(nu)
            for i in range(PREC + 1):
                r[i] += Fr(c, den) * z[i]
    return tuple(r)


def submultisets(mu):
    C = Counter(mu)
    keys = sorted(C)
    for choice in product(*[range(C[k] + 1) for k in keys]):
        w = 1
        for k, c in zip(keys, choice):
            w *= comb(C[k], c)
        sub = []
        rest = []
        for k, c in zip(keys, choice):
            sub += [k] * c
            rest += [k] * (C[k] - c)
        yield norm(sub), norm(rest), w


@lru_cache(None)
def cc(mu):
    # x^{|mu|} * |Aut mu| [p_mu] F, labelled-part cumulant
    if not mu:
        return tuple([Fr(0)] * (PREC + 1))
    s0 = mu[0]
    rest = mu[1:]
    r = list(zz(mu))
    for sub, oth, w in submultisets(rest):
        if not oth:
            continue
        B = norm((s0,) + sub)
        t = mul(cc(B), zz(oth))
        for i in range(PREC + 1):
            r[i] -= w * t[i]
    return tuple(r)


def M(g, mu):
    # Prop 2.1: [x^{2g-2+l}] hat c = sigma^{g+l-1} Pi_mu M_g(mu)
    mu = norm(mu)
    l = len(mu)
    n = sum(mu)
    idx = 2 * g - 2 + l + n
    assert idx <= PREC
    v = cc(mu)[idx]
    p = sigma ** (g + l - 1)
    for k in mu:
        p *= Pi(k)
    return v / p


def lag(n, al, x):
    return sum((-1) ** i * comb(n + al, n - i) * Fr(x) ** i / factorial(i) for i in range(n + 1))


@lru_cache(None)
def X(g, U):
    U = norm(U)
    r = len(U)
    J = sum(U)
    assert r >= 1
    if g == 0:
        v = Fr(1)
        if r >= 3:
            for m in range(1, r - 2):
                v *= J - m
            return v
        if r == 2:
            return Fr(1, J)
        return Fr(1, J * (J + 1))
    core = tuple(u for u in U if u >= 2)
    v = r - len(core)
    if v and core:
        J0 = sum(core)
        f = Fr(1)
        for i in range(v):
            f *= J0 + i
        return f * X(g, core)
    if v and g >= 2:
        return Fr(0)
    if v and g == 1 and v > 1:
        return factorial(v - 1) * X(1, (1,))
    N = 3 * g - 3 + r
    return sum(Fr((-1) ** a, factorial(a)) * lag(N - a, 2 * g - 2 + r + a, -(J + a)) * M(g, U + (1,) * a)
               for a in range(N + 1))


def ff(a, k):
    r = 1
    for i in range(k):
        r *= a - i
    return r


def set_partitions(L):
    if not L:
        yield []
        return
    h = L[0]
    for sp in set_partitions(L[1:]):
        yield [[h]] + sp
        for i in range(len(sp)):
            yield sp[:i] + [[h] + sp[i]] + sp[i + 1:]


def calpha(alpha, k0):
    C = Counter(alpha)
    v = Fr(1)
    for i, m in C.items():
        v *= Fr(-1, i) ** m / factorial(m)
    if alpha == (k0,):
        v += 1
    return v


OPT = {'codd': -1, 'tau': tau}


def block_phi(Bparts, TB, g, k, bump=None):
    b = len(Bparts)
    if b % 2 == 1:
        c = 0 if k == 0 else OPT['codd']
    else:
        c = 2 if k == 0 else 1
    if c == 0:
        return Fr(0)
    U = norm(tuple(TB) + tuple(Bparts))
    d = sum(U) - len(U) - 2 * g + 2
    fa = ff(d, k)
    if fa == 0:
        return Fr(0)
    pref = Fr(c)
    for i in Bparts:
        pref *= i * Pi(i)
    pref *= sigma ** (b + g - 1) * OPT['tau'] ** k / factorial(k) * fa
    x = X(g, U)
    if bump and (g, U) == bump:
        x += 1
    return pref * x


def relation(k0, R, n, bump=None):
    tot = Fr(0)
    for alpha in parts(k0):
        ca = calpha(alpha, k0)
        for sp in set_partitions(list(range(len(alpha)))):
            blocks = [[alpha[i] for i in B] for B in sp]
            nb = len(blocks)
            for asg in product(range(nb), repeat=len(R)):
                TBs = [[R[i] for i in range(len(R)) if asg[i] == b] for b in range(nb)]

                def rec(b, rem):
                    if b == nb:
                        return Fr(1) if rem == 0 else Fr(0)
                    s = Fr(0)
                    B = blocks[b]
                    sz = len(B)
                    for g in range(0, rem // 2 + 1):
                        kmin = 1 if sz % 2 else 0
                        for k in range(kmin, rem - 2 * g + 2 - sz + 1):
                            nB = 2 * g - 2 + sz + k
                            if nB < 0 or nB > rem:
                                continue
                            p = block_phi(B, TBs[b], g, k, bump)
                            if p == 0:
                                continue
                            s += p * rec(b + 1, rem - nB)
                    return s

                tot += ca * rec(0, n)
    return tot


if __name__ == '__main__':
    print('tau', tau)
    print('M_1((1))', M(1, (1,)), 'expected', -Fr(1, 24) / sigma)
    print('X_1(1)', X(1, (1,)), 'expected', -(1 + 1 / sigma) / 24)
    print('X_2(2)', X(2, (2,)), ' X_2(3)', X(2, (3,)), ' X_2(2,2)', X(2, (2, 2)))
    print('X_3(2),(3),(4)', X(3, (2,)), X(3, (3,)), X(3, (4,)), ' X_3(2,2)', X(3, (2, 2)), ' X_3(5)', X(3, (5,)))
    fails = tot = 0
    cases = [(2, k0, R) for k0 in (2, 3, 4) for R in ((), (2,), (3,), (2, 2), (4,))] + \
            [(3, k0, R) for k0 in (2, 3, 4) for R in ((), (2,), (3,))]
    for G, k0, R in cases:
        for n in (2 * G, 2 * G + 1):
            v = relation(k0, R, n)
            tot += 1
            if v != 0:
                fails += 1
                print('FAIL', G, k0, R, n, v)
    print('relations', tot, 'fails', fails)
    ok = bad = 0
    for G, k0, R in cases:
        T = norm(R + (k0,))
        D = sum(T) - len(T)
        v0 = relation(k0, R, 2 * G)
        v1 = relation(k0, R, 2 * G, bump=(G, T))
        d = D - 2 * G + 2
        ex = -(k0 - 1) * Pi(k0) * sigma ** G * tau * d
        if v1 - v0 == ex:
            ok += 1
        else:
            bad += 1
            print('resp', G, k0, R, v1 - v0, ex)
    print('target coefficient ok', ok, 'bad', bad)
    OPT['codd'] = 1
    c1 = sum(1 for G, k0, R in cases if relation(k0, R, 2 * G) != 0)
    OPT['codd'] = -1
    print('control flip c(odd,k): nonzero', c1, 'of', len(cases))
    OPT['tau'] = -tau
    c2 = sum(1 for G, k0, R in cases if relation(k0, R, 2 * G) != 0)
    OPT['tau'] = tau
    print('control beta->-beta: nonzero', c2, 'of', len(cases))
