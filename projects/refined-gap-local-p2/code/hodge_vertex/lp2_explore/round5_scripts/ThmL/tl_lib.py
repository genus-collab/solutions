# ThmL round 5: independent exact library (stdlib only, written from scratch).
# Partitions are descending tuples. Polys in p are dicts {partition: Fraction}.
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial
import itertools

def partitions(n, m=None):
    if m is None: m = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, m), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def f2(nu):
    return sum(j - i for i, r in enumerate(nu) for j in range(r))

def zmu(mu):
    z = 1
    for k in set(mu):
        m = mu.count(k)
        z *= k**m * factorial(m)
    return z

def aut(mu):
    a = 1
    for k in set(mu):
        a *= factorial(mu.count(k))
    return a

@lru_cache(maxsize=None)
def chi(nu, mu):
    """Murnaghan-Nakayama on beta sets. chi^nu at cycle type mu."""
    if not mu:
        return 1 if not nu else 0
    k = mu[0]
    rest = mu[1:]
    L = len(nu)
    beta = [nu[i] + L - 1 - i for i in range(L)]
    bs = set(beta)
    tot = 0
    for i, b in enumerate(beta):
        c = b - k
        if c < 0 or c in bs:
            continue
        # leg length = number of beta strictly between c and b
        leg = sum(1 for x in beta if c < x < b)
        nb = sorted([x for x in beta if x != b] + [c], reverse=True)
        nn = tuple(y for y in (nb[j] - (L - 1 - j) for j in range(L)) if y > 0)
        tot += (-1)**leg * chi(nn, rest)
    return tot

def schur_in_p(nu):
    n = sum(nu)
    return {mu: Fr(chi(nu, mu), zmu(mu)) for mu in partitions(n) if chi(nu, mu) != 0}

def eval_p(mu, pv):
    r = Fr(1)
    for k in mu:
        r *= pv[k]
    return r

def poly_mul(a, b, maxdeg):
    r = {}
    for m1, c1 in a.items():
        d1 = sum(m1)
        for m2, c2 in b.items():
            if d1 + sum(m2) > maxdeg:
                continue
            m = tuple(sorted(m1 + m2, reverse=True))
            r[m] = r.get(m, 0) + c1 * c2
    return r

def d_t(a, k):
    """d/dt_k = k d/dp_k."""
    r = {}
    for m, c in a.items():
        e = m.count(k)
        if e:
            l = list(m)
            l.remove(k)
            t = tuple(l)
            r[t] = r.get(t, 0) + c * e * k
    return r

def d_multi(a, idx):
    for k in idx:
        a = d_t(a, k)
    return a

def hirota(f, g, alpha, maxdeg):
    """D^alpha f.g = sum over subsets S of alpha (-1)^{|alpha - S|} d^S f d^(alpha-S) g, truncated."""
    r = {}
    n = len(alpha)
    for bits in itertools.product([0, 1], repeat=n):
        S = [alpha[i] for i in range(n) if bits[i]]
        C = [alpha[i] for i in range(n) if not bits[i]]
        sgn = (-1)**len(C)
        pr = poly_mul(d_multi(f, S), d_multi(g, C), maxdeg)
        for m, c in pr.items():
            r[m] = r.get(m, 0) + sgn * c
    return r

def member_Pk(k):
    """P_k = D_k + p_k(-Dtilde) as {alpha (partition of k): coefficient}."""
    out = {}
    for lam in partitions(k):
        c = Fr(1)
        for i in set(lam):
            m = lam.count(i)
            c *= Fr(-1, i)**m / factorial(m)
        out[lam] = c
    out[(k,)] = out.get((k,), 0) + 1
    return {a: c for a, c in out.items() if c != 0}

def set_partitions(items):
    """All set partitions of the list `items` (positions labelled)."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for sp in set_partitions(rest):
        yield [[first]] + sp
        for i in range(len(sp)):
            yield sp[:i] + [[first] + sp[i]] + sp[i+1:]

def Pi(k, tau):
    r = Fr(1)
    for a in range(1, k):
        r *= (k * tau + a)
    return r / factorial(k - 1)

def bernoulli_list(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    from math import comb
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1)
    return B

# ---- truncated series in one variable: lists of Fractions ----
def ser_mul(a, b, L):
    r = [Fr(0)] * (L + 1)
    for i, x in enumerate(a):
        if x == 0 or i > L: continue
        for j, y in enumerate(b):
            if i + j > L: break
            if y != 0:
                r[i + j] += x * y
    return r

def ser_exp(a, L):
    """exp of series a with a[0]=0."""
    e = [Fr(0)] * (L + 1)
    e[0] = Fr(1)
    for n in range(1, L + 1):
        s = Fr(0)
        for k in range(1, n + 1):
            if k < len(a) and a[k] != 0:
                s += k * a[k] * e[n - k]
        e[n] = s / n
    return e

def hooks(nu):
    conj = [sum(1 for r in nu if r > j) for j in range(nu[0])] if nu else []
    return [nu[i] - j + conj[j] - i - 1 for i in range(len(nu)) for j in range(nu[i])]
