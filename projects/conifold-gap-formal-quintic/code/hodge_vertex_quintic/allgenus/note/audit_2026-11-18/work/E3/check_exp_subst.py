# Independent sanity check (not a justification) for lem:fsp-res (iii):
#   T_chi(E) := sum_j E^{(j)} chi^j / j!  equals  E * exp(-(f0(X+chi) - f0(X))/zh),
# with E = exp(-f0/zh), f0 = 1/X + X^4/4, chi in u Q[X,1/X][[u]].
# Elements: dict (a, j, b) -> Fraction, meaning u^a X^j w^b, w = 1/zh.
# Truncation: u^a with a <= N, w^b with b <= K.  All operations used are w-linear
# and do not lower the w-degree, so truncation mod w^(K+1) is exact.
from fractions import Fraction as F
from math import factorial
N, K = 4, 7

def add(p, q, c=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + c * v
        if r[k] == 0:
            del r[k]
    return r

def mul(p, q):
    r = {}
    for (a1, j1, b1), v1 in p.items():
        for (a2, j2, b2), v2 in q.items():
            a, b = a1 + a2, b1 + b2
            if a > N or b > K:
                continue
            k = (a, j1 + j2, b)
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}

def scal(p, c):
    return {k: c * v for k, v in p.items() if c * v != 0}

def dX(p):
    r = {}
    for (a, j, b), v in p.items():
        if j != 0:
            r[(a, j - 1, b)] = r.get((a, j - 1, b), 0) + j * v
    return r

ONE = {(0, 0, 0): F(1)}

def pw(p, n):
    r = ONE
    for _ in range(n):
        r = mul(r, p)
    return r

def exp_small(p, maxn):
    # p has positive u-order or positive w-order; sum p^n/n!
    r, t = ONE, ONE
    for n in range(1, maxn + 1):
        t = scal(mul(t, p), F(1, n))
        r = add(r, t)
    return r

# f0 = 1/X + X^4/4 ; g = -w f0
f0 = {(0, -1, 0): F(1), (0, 4, 0): F(1, 4)}
g = scal(mul({(0, 0, 1): F(1)}, f0), -1)
E = exp_small(g, K)  # exact mod w^(K+1)

chi = {(1, 1, 0): F(1), (1, -1, 0): F(2), (2, 3, 0): F(1, 3), (3, -2, 0): F(-5, 7)}

# LHS: T_chi(E)
lhs, Ej = {}, E
for j in range(0, N + 1):
    lhs = add(lhs, scal(mul(Ej, pw(chi, j)), F(1, factorial(j))))
    Ej = dX(Ej)

# Laurent substitution of f0: (X+chi)^{-1} = sum (-1)^j X^{-1-j} chi^j ; (X+chi)^4/4
inv = {}
for j in range(0, N + 1):
    inv = add(inv, scal(mul({(0, -1 - j, 0): F(1)}, pw(chi, j)), F((-1) ** j)))
Xc = add({(0, 1, 0): F(1)}, chi)
f0sub = add(inv, scal(pw(Xc, 4), F(1, 4)))
h = add(f0sub, f0, -1)  # in u Q[X^+-1][[u]]
assert all(a >= 1 for (a, j, b) in h)
rhs = mul(E, exp_small(scal(mul({(0, 0, 1): F(1)}, h), -1), N))

diff = add(lhs, rhs, -1)
print("number of monomials in LHS:", len(lhs), " RHS:", len(rhs))
print("LHS == RHS (mod u^%d, w^%d):" % (N + 1, K + 1), len(diff) == 0)

# residue invariance check: Res_X[Omega(X+chi)(1+chi')] = Res[Omega] for Omega = E*beta
beta = {(0, 2, 0): F(3), (0, -3, 0): F(1), (1, 1, 0): F(-2), (2, 0, 0): F(1, 5)}
Om = mul(E, beta)
T, Oj = {}, Om
for j in range(0, N + 1):
    T = add(T, scal(mul(Oj, pw(chi, j)), F(1, factorial(j))))
    Oj = dX(Oj)
lhsR = mul(T, add(ONE, dX(chi)))
res = lambda p: {(a, b): v for (a, j, b), v in p.items() if j == -1}
print("Res invariance:", res(lhsR) == res(Om))
