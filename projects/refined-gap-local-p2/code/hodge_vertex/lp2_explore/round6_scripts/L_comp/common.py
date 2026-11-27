# L_comp: own combinatorics; no reuse of authors' scripts
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial, comb
from collections import Counter

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n-k, k):
            yield (k,)+rest

def zee(mu):
    z = 1
    for k, m in Counter(mu).items():
        z *= k**m * factorial(m)
    return z

def aut(mu):
    a = 1
    for k, m in Counter(mu).items():
        a *= factorial(m)
    return a

def conj(nu):
    if not nu: return ()
    return tuple(sum(1 for x in nu if x > i) for i in range(nu[0]))

def contents_sum(nu):
    return sum(j-i for i, r in enumerate(nu) for j in range(r))

def hooks(nu):
    c = conj(nu)
    return [r-j-1 + c[j]-i-1 + 1 for i, r in enumerate(nu) for j in range(r)]

def beta_set(nu, L):
    nu = list(nu)+[0]*(L-len(nu))
    return tuple(nu[i]+L-1-i for i in range(L))

def from_beta(bs):
    bs = sorted(bs, reverse=True)
    L = len(bs)
    nu = [bs[i]-(L-1-i) for i in range(L)]
    return tuple(x for x in nu if x > 0)

@lru_cache(maxsize=None)
def chi(nu, mu):
    """Murnaghan-Nakayama by removal of rim hooks using beta numbers. mu tuple, order irrelevant"""
    if not mu:
        return 1 if not nu else 0
    if sum(nu) != sum(mu):
        return 0
    k = mu[0]; rest = mu[1:]
    # strip trailing ones by hook length count for speed if all parts ones
    if k == 1 and all(x == 1 for x in rest):
        n = sum(nu); h = 1
        for x in hooks(nu): h *= x
        return factorial(n)//h
    L = len(nu)+k
    bs = list(beta_set(nu, L))
    S = set(bs)
    tot = 0
    for b in bs:
        nb = b-k
        if nb >= 0 and nb not in S:
            height = sum(1 for c in bs if nb < c < b)
            new = [x for x in bs if x != b]+[nb]
            tot += (-1)**height * chi(from_beta(new), rest)
    return tot

def norm_mu(mu):
    """order parts descending so ones are stripped last"""
    return tuple(sorted(mu, reverse=True))

def schur_in_p(nu):
    """dict mu->coefficient of p_mu in s_nu"""
    n = sum(nu)
    return {mu: Fr(chi(nu, mu), zee(mu)) for mu in partitions(n)}

def bernoulli_table(N):
    B = [Fr(0)]*(N+1)
    B[0] = Fr(1)
    for m in range(1, N+1):
        s = Fr(0)
        for k in range(m):
            s += comb(m+1, k)*B[k]
        B[m] = -s/(m+1)
    return B
BERN = bernoulli_table(80)

def Pi(m, tau):
    r = Fr(1)
    for a in range(1, m):
        r *= (m*tau+a)
    return r/factorial(m-1)

def falling(a, k):
    r = 1
    for i in range(k):
        r = r*(a-i)
    return r

def rising(a, k):
    r = 1
    for i in range(k):
        r = r*(a+i)
    return r
