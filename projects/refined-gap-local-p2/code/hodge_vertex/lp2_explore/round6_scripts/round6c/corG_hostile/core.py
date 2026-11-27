# Fresh engine, imports no project script. Exact arithmetic only.
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import permutations, product
from math import factorial, comb
import sys
sys.setrecursionlimit(100000)

@lru_cache(None)
def bern(n):
    if n == 0:
        return Fr(1)
    return -sum(comb(n + 1, k) * bern(k) for k in range(n)) / (n + 1)

def df(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r

# Witten-Kontsevich by DVV
@lru_cache(None)
def wk(g, d):
    d = tuple(sorted(d))
    n = len(d)
    if g < 0 or 2 * g - 2 + n <= 0 or sum(d) != 3 * g - 3 + n or n == 0:
        return Fr(0)
    if g == 0 and n == 3:
        return Fr(1)
    if g == 1 and d == (1,):
        return Fr(1, 24)
    if d[0] == 0:
        rest = d[1:]
        s = Fr(0)
        for i in range(len(rest)):
            if rest[i] > 0:
                s += wk(g, rest[:i] + (rest[i] - 1,) + rest[i + 1:])
        return s
    if 1 in d:
        i = d.index(1)
        rest = d[:i] + d[i + 1:]
        return (2 * g - 2 + len(rest)) * wk(g, rest)
    k = d[-1] - 1
    rest = d[:-1]
    s = Fr(0)
    for j in range(len(rest)):
        s += Fr(df(2 * k + 2 * rest[j] + 1), df(2 * rest[j] - 1)) * wk(g, rest[:j] + (rest[j] + k,) + rest[j + 1:])
    for r in range(k):
        t = k - 1 - r
        c = df(2 * r + 1) * df(2 * t + 1)
        x = wk(g - 1, (r, t) + rest)
        m = len(rest)
        for mask in range(1 << m):
            A = tuple(rest[i] for i in range(m) if mask >> i & 1)
            B = tuple(rest[i] for i in range(m) if not mask >> i & 1)
            da = r + sum(A)
            na = len(A) + 1
            if (da - na + 3) % 3:
                continue
            g1 = (da - na + 3) // 3
            g2 = g - g1
            if g1 < 0 or g2 < 0:
                continue
            x += wk(g1, (r,) + A) * wk(g2, (t,) + B)
        s += Fr(c, 2) * x
    return s / df(2 * k + 3)

# Mumford (Faber-Pandharipande s.1.1 form): integrals of psi and ch_odd
@lru_cache(None)
def mum(g, psi, ch):
    psi = tuple(sorted(psi))
    ch = tuple(sorted(ch))
    n = len(psi)
    if g < 0 or 2 * g - 2 + n <= 0:
        return Fr(0)
    if sum(psi) + sum(ch) != 3 * g - 3 + n:
        return Fr(0)
    if not ch:
        return wk(g, psi) if n > 0 else Fr(0)
    if g == 0:
        return Fr(0)
    k = ch[-1]
    rest = ch[:-1]
    if k % 2 == 0:
        return Fr(0)
    l = (k + 1) // 2
    pref = bern(2 * l) / factorial(2 * l)
    s = mum(g, psi + (2 * l,), rest)
    for i in range(n):
        s -= mum(g, psi[:i] + (psi[i] + k,) + psi[i + 1:], rest)
    bd = Fr(0)
    for i in range(2 * l - 1):
        j = 2 * l - 2 - i
        sg = (-1) ** i
        bd += sg * mum(g - 1, psi + (j, i), rest)
        m = len(rest)
        for mask in range(1 << m):
            ra = tuple(rest[x] for x in range(m) if mask >> x & 1)
            rb = tuple(rest[x] for x in range(m) if not mask >> x & 1)
            for pm in range(1 << n):
                A = tuple(psi[x] for x in range(n) if pm >> x & 1)
                B = tuple(psi[x] for x in range(n) if not pm >> x & 1)
                da = sum(A) + sum(ra) + j
                na = len(A) + 1
                if (da - na + 3) % 3:
                    continue
                g1 = (da - na + 3) // 3
                g2 = g - g1
                if g1 < 0 or g2 < 0:
                    continue
                bd += sg * mum(g1, A + (j,), ra) * mum(g2, B + (i,), rb)
    s += bd / 2
    return pref * s

def odd_parts(total, maxp=None):
    if total == 0:
        yield ()
        return
    if maxp is None:
        maxp = total
    if maxp % 2 == 0:
        maxp -= 1
    for p in range(maxp, 0, -2):
        if p <= total:
            for rest in odd_parts(total - p, p):
                yield (p,) + rest

@lru_cache(None)
def omega_terms(g, psi):
    n = len(psi)
    dim = 3 * g - 3 + n
    deg = dim - sum(psi)
    if deg < 0 or 2 * g - 2 + n <= 0:
        return ()
    out = []
    for part in odd_parts(deg):
        v = mum(g, tuple(psi), part)
        if v == 0:
            continue
        c = v
        mult = {}
        for k in part:
            mult[k] = mult.get(k, 0) + 1
        for k, m in mult.items():
            c *= Fr((-factorial(k - 1)) ** m, factorial(m))
        out.append((tuple(sorted(mult.items())), c))
    return tuple(out)

def stable_graphs(g):
    res = {}
    for V in range(1, 2 * g - 1):
        pairs = [(i, j) for i in range(V) for j in range(i, V)]
        for gs in product(range(g + 1), repeat=V):
            if sum(gs) > g:
                continue
            Ecount = g - sum(gs) + V - 1
            if Ecount < 0 or list(gs) != sorted(gs, reverse=True):
                continue
            def rec(idx, left, m):
                if idx == len(pairs):
                    if left == 0:
                        yield dict(m)
                    return
                for c in range(left + 1):
                    m[pairs[idx]] = c
                    yield from rec(idx + 1, left - c, m)
                m[pairs[idx]] = 0
            for m in rec(0, Ecount, {}):
                val = [0] * V
                for (i, j), c in m.items():
                    val[i] += c
                    val[j] += c
                if any(2 * gs[i] - 2 + val[i] <= 0 for i in range(V)):
                    continue
                seen = {0}
                st = [0]
                while st:
                    u = st.pop()
                    for (i, j), c in m.items():
                        if c and (i == u or j == u):
                            w = j if i == u else i
                            if w not in seen:
                                seen.add(w)
                                st.append(w)
                if len(seen) < V:
                    continue
                best = None
                aut = 0
                for P in permutations(range(V)):
                    gp = tuple(gs[P[i]] for i in range(V))
                    mat = tuple(m.get((min(P[i], P[j]), max(P[i], P[j])), 0) for (i, j) in pairs)
                    k = (gp, mat)
                    if best is None or k < best:
                        best = k
                    if gp == gs and all(m.get((min(P[i], P[j]), max(P[i], P[j])), 0) == m[(i, j)] for (i, j) in pairs):
                        aut += 1
                if best in res:
                    continue
                A = aut
                for (i, j), c in m.items():
                    A *= factorial(c)
                    if i == j:
                        A *= 2 ** c
                res[best] = (gs, dict(m), A)
    return list(res.values())

def beta(n, e1, e2):
    g = (n + 2) // 2
    def s(m, e):
        if m % 2:
            return Fr(0)
        return (Fr(2) ** (1 - m) - 1) * bern(m) * e ** m / factorial(m)
    return factorial(n - 1) * sum(s(a, e1) * s(2 * g - a, e2) for a in range(0, 2 * g + 1)) / (e1 * e2)

def omega_num(g, psi, e1, e2):
    E = e1 + e2
    s = Fr(0)
    for items, c in omega_terms(g, tuple(sorted(psi))):
        t = c
        for k, m in items:
            t *= (Fr(1) / E ** k + (-1) ** k * (Fr(1) / e1 ** k + Fr(1) / e2 ** k)) ** m
        s += t
    return s

def excess_multisets(left, minv=2):
    yield ()
    for d in range(minv, left + 2):
        for rest in excess_multisets(left - (d - 1), d):
            yield (d,) + rest

def Xi(g, R, e1, e2):
    """Xi_g(R)= Sum_m (1/m!) Int Omega_w prod lambda_j(psi_i) prod T(psi_k), T=-psi^2/(E-psi)."""
    E = e1 + e2
    s = Fr(0)
    r = len(R)
    dim = 3 * g - 3 + r
    def lamcoef(j, a):
        return Fr(j) ** a / E ** (a + 1)
    from itertools import product as P
    for ds in excess_multisets(dim):
        mult = {}
        for d in ds:
            mult[d] = mult.get(d, 0) + 1
        w = Fr(1)
        for d, c in mult.items():
            w /= factorial(c)
        for d in ds:
            w *= -(E ** (1 - d))
        budget = dim - sum(d - 1 for d in ds)
        for aa in P(range(budget + 1), repeat=r):
            if sum(aa) > budget:
                continue
            c = w
            for j, a in zip(R, aa):
                c *= lamcoef(j, a)
            s += c * omega_num(g, tuple(aa) + ds, e1, e2)
    return s

if __name__ == "__main__":
    for g in (2, 3, 4):
        print('graphs', g, len(stable_graphs(g)))
    print('lambda1 M11 (1/24):', mum(1, (0,), (1,)))
    print('lambda1^3 M2 (1/2880):', mum(2, (), (1, 1, 1)))
    print('lambda1^6 M3 (1/90720):', mum(3, (), (1, 1, 1, 1, 1, 1)))
    print('3!ch3 M2 (-1/5760):', 6 * mum(2, (), (3,)))
    for (e1, e2) in [(Fr(3), Fr(7)), (Fr(1), Fr(1)), (Fr(5, 2), Fr(-3, 7))]:
        for g in (2, 3, 4):
            E = e1 + e2
            print('Thm VI', g, e1, e2, (-e1 * e2 * E ** 3) ** (g - 1) * Xi(g, (), e1, e2) == beta(2 * g - 2, e1, e2))
    # M_irr sample: Xi_h(R)=0, R parts >=2, D(R)<=2h-3; and first non-vanishing just outside the range
    e1, e2 = Fr(3), Fr(7)
    print('Xi_2({2})', Xi(2, (2,), e1, e2), ' Xi_2({3})', Xi(2, (3,), e1, e2), ' Xi_2({2,2})', Xi(2, (2, 2), e1, e2))
    for R in [(2,), (3,), (4,), (2, 2), (2, 3), (2, 2, 2)]:
        print('Xi_3', R, Xi(3, R, e1, e2))
    print('Xi_3 (5,) D=4 outside range', Xi(3, (5,), e1, e2))
