"""Independent low-order checks for appendix A (large-radius cyclic structure) in q.
Elements of Q[[q]][l] (l = formal log q) are lists [a_0, a_1, ...] of truncated series (coefficient of l^j).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
from ser import S

N = 9   # precision in q
WN = 6  # powers of w kept: w^0..w^5

def zero():
    return S(0, [0] * N)
def one():
    return S(0, [1] + [0] * (N - 1))

# polynomial arithmetic in w (truncated at w^WN) with Fraction coefficients
def pmul(a, b):
    c = [Fr(0)] * WN
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j < WN:
                c[i + j] += x * y
    return c
def pinv(a):
    b = [Fr(0)] * WN
    b[0] = 1 / a[0]
    for k in range(1, WN):
        b[k] = -sum(a[j] * b[k - j] for j in range(1, k + 1)) / a[0]
    return b

# P(w,q) = sum_d q^d prod_{k=1}^{5d}(5w+k) / prod_{k=1}^d((w+k)^5 - w^5)
ptil = [[Fr(0)] * N for _ in range(WN)]   # ptil[m][d]
for d in range(N):
    num = [Fr(1)] + [Fr(0)] * (WN - 1)
    for k in range(1, 5 * d + 1):
        num = pmul(num, [Fr(k), Fr(5)] + [Fr(0)] * (WN - 2))
    den = [Fr(1)] + [Fr(0)] * (WN - 1)
    for k in range(1, d + 1):
        # (w+k)^5 - w^5 = 5k w^4 + 10k^2 w^3 + 10 k^3 w^2 + 5 k^4 w + k^5
        fac = [Fr(k ** 5), Fr(5 * k ** 4), Fr(10 * k ** 3), Fr(10 * k ** 2), Fr(5 * k)] + [Fr(0)] * (WN - 5)
        den = pmul(den, fac[:WN])
    term = pmul(num, pinv(den))
    for m in range(WN):
        ptil[m][d] = term[m]
pt = [S(0, ptil[m]) for m in range(WN)]

C0 = pt[0]
print("C0 =", C0.c[:4])
assert C0.c[:4] == [1, 120, 113400, 168168000]
# tilde I_1 = 5 sum (5d)!/(d!)^5 (sum_{m=d+1}^{5d} 1/m) q^d
tI1 = S(0, [5 * Fr(factorial(5 * d), factorial(d) ** 5) * sum(Fr(1, m) for m in range(d + 1, 5 * d + 1)) for d in range(N)])
print("p~_1 == tilde I_1:", (pt[1] - tI1).iszero())

# log-polynomials
def Dq(f):
    # q d/dq on a truncated series
    return S(f.v, [(f.v + k) * f.c[k] for k in range(len(f.c))])
def LD(x):
    n = len(x)
    out = []
    for j in range(n):
        t = Dq(x[j])
        if j + 1 < n:
            t = t + (j + 1) * x[j + 1]
        out.append(t)
    return out
def Ladd(x, y):
    n = max(len(x), len(y))
    x = x + [zero()] * (n - len(x)); y = y + [zero()] * (n - len(y))
    return [x[j] + y[j] for j in range(n)]
def Lscal(s, x):
    return [s * t for t in x]
def Lmul(x, y):
    out = [zero() for _ in range(len(x) + len(y) - 1)]
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            out[i + j] = out[i + j] + a * b
    return out
def LDk(x, k):
    for _ in range(k):
        x = LD(x)
    return x
def Lshow(x):
    return [t.c[:4] for t in x]
def Liszero(x):
    return all(t.iszero() for t in x)

# p_k = sum_j l^j/j! ptil_{k-j}
def pk(k):
    return [pt[k - j] * Fr(1, factorial(j)) for j in range(k + 1)]
p = [pk(k) for k in range(5)]

q = S(1, [1] + [0] * (N - 1))
def PF4(x):
    # D^4 - 5 q prod_{k=1}^4 (5D+k)
    t = x
    for k in range(1, 5):
        t = Ladd(Lscal(Fr(5), LD(t)), Lscal(Fr(k), t))
    return Ladd(LDk(x, 4), Lscal(Fr(-5), [q * s for s in t]))
for k in range(5):
    r = PF4(p[k])
    const = [t.c[0] for t in r]
    print("PF4 p_%d: l-coefficients constant terms %s, higher q-orders zero: %s" % (k, const, all(all(c == 0 for c in t.c[1:]) for t in r)))

# L^5 = 1/(1-5^5 q)
L5 = S(0, [Fr(3125) ** d for d in range(N)])
L10 = L5 * L5

# Wronskian det(D^j p_k), j,k=0..3, in Q[[q]][l]
def det(mat):
    n = len(mat)
    if n == 1:
        return mat[0][0]
    tot = [zero()]
    for j in range(n):
        minor = [[mat[i][k] for k in range(n) if k != j] for i in range(1, n)]
        t = Lmul(mat[0][j], det(minor))
        if j % 2:
            t = Lscal(Fr(-1), t)
        tot = Ladd(tot, t)
    return tot
Wr = det([[LDk(p[k], j) for k in range(4)] for j in range(4)])
print("Wronskian: l-free:", all(t.iszero() for t in Wr[1:]), "; equals L^10:", (Wr[0] - L10).iszero())

# mirror map and kappa's
T = [tI1 / C0, one()]          # T = l + tI1/C0
C1 = LD(T)
assert Liszero(C1[1:])
C1 = C1[0]
print("C1 =", C1.c[:3])
def dT(x):
    return [t / C1 for t in LD(x)]
def divs(x, s):
    return [t / s for t in x]
P0 = [divs(p[k], C0) for k in range(5)]   # p^(0)_k
print("p^(0)_1 == T:", Liszero(Ladd(P0[1], Lscal(Fr(-1), T))))
k1 = dT(dT(P0[2]))
print("kappa_1 l-free:", Liszero(k1[1:]), "; equals L^5/(C0^2 C1^3):", (k1[0] - L5 / (C0 * C0 * C1 * C1 * C1)).iszero(), "; kappa_1 =", k1[0].c[:3])
k2 = dT(divs(dT(dT(P0[3])), k1[0]))
print("kappa_2 l-free:", Liszero(k2[1:]), "; equals 1:", (k2[0] - one()).iszero())
k3 = dT(dT(divs(dT(dT(P0[4])), k1[0])))
print("kappa_3 l-free:", Liszero(k3[1:]), "; equals C0/C1:", (k3[0] - C0 / C1).iszero(), "; kappa_3 =", k3[0].c[:3])
# product of the five kappa^LR = (L/C1)^5
prod = k1[0] * k3[0] * k3[0]
print("prod kappa^LR == L^5/C1^5:", (prod - L5 / (C1 ** 5)).iszero())
# Yukawa check: 5*kappa_1 in the variable Q = e^T should be 5 + 2875 Q + 4876875 Q^2 (quintic 3-point function)
eT = S(1, [1] + [0] * (N - 1))
# exp(tI1/C0)
x = tI1 / C0
ex = one(); term = one()
for n in range(1, N + 1):
    term = term * x * Fr(1, n)
    ex = ex + term
Qv = q * ex
k1Q = (k1[0] - one() - 575 * Qv - Fr(4876875, 5) * Qv * Qv)
print("5*kappa_1 = 5 + 2875 Q + 4876875 Q^2 + O(Q^3):", k1Q.norm().v >= 3, "; Q =", Qv.c[:3])
