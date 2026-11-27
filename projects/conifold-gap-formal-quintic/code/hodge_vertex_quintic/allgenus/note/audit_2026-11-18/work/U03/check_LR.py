# Independent sanity checks for Appendix A (large-radius cyclic structure) -- exact arithmetic.
# Elements of Q[[q]][l] truncated at q^N : dict {(n,k): c} = c q^n l^k.
from fractions import Fraction as F
from math import factorial, comb

N = 7

class S:
    __slots__ = ("d",)
    def __init__(self, d=None):
        self.d = {}
        if d:
            for key, c in d.items():
                if c != 0 and key[0] <= N:
                    self.d[key] = F(c)
    @staticmethod
    def const(c): return S({(0, 0): F(c)})
    def __add__(a, b):
        if not isinstance(b, S): b = S.const(b)
        d = dict(a.d)
        for k, c in b.d.items(): d[k] = d.get(k, 0) + c
        return S(d)
    __radd__ = __add__
    def __neg__(a): return S({k: -c for k, c in a.d.items()})
    def __sub__(a, b):
        if not isinstance(b, S): b = S.const(b)
        return a + (-b)
    def __rsub__(a, b): return S.const(b) - a
    def __mul__(a, b):
        if not isinstance(b, S): b = S.const(b)
        d = {}
        for (n1, k1), c1 in a.d.items():
            for (n2, k2), c2 in b.d.items():
                n = n1 + n2
                if n > N: continue
                key = (n, k1 + k2); d[key] = d.get(key, 0) + c1*c2
        return S(d)
    __rmul__ = __mul__
    def D(a):
        d = {}
        for (n, k), c in a.d.items():
            if n: d[(n, k)] = d.get((n, k), 0) + n*c
            if k: d[(n, k-1)] = d.get((n, k-1), 0) + k*c
        return S(d)
    def Dk(a, k):
        for _ in range(k): a = a.D()
        return a
    def inv(a):
        assert all(k == 0 for (_, k) in a.d) and (0, 0) in a.d
        c0 = a.d[(0, 0)]
        r = S({key: c/c0 for key, c in a.d.items() if key != (0, 0)})
        res = S.const(1); term = S.const(1)
        for _ in range(N + 1):
            term = term * (-r)
            if not term.d: break
            res = res + term
        return res * (1/c0)
    def __repr__(a):
        return " + ".join(f"({c})q^{n}" + (f" l^{k}" if k else "") for (n, k), c in sorted(a.d.items())) or "0"
    def iszero(a): return not a.d

q = S({(1, 0): 1}); l = S({(0, 1): 1})

# --- P(w,q) = sum_d q^d prod_{k=1}^{5d}(5w+k) / prod_{k=1}^{d}((w+k)^5 - w^5), power series in w up to w^W
W = 5
def pmul(a, b):
    r = [F(0)]*(W+1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j > W: break
            r[i+j] += x*y
    return r
def pinv(a):
    r = [F(0)]*(W+1); r[0] = 1/a[0]
    for n in range(1, W+1):
        r[n] = -sum(a[k]*r[n-k] for k in range(1, n+1))/a[0]
    return r
ptilde = [S() for _ in range(W+1)]
for d in range(0, N+1):
    num = [F(1)] + [F(0)]*W
    for k in range(1, 5*d+1):
        num = pmul(num, [F(k), F(5)] + [F(0)]*(W-1))
    den = [F(1)] + [F(0)]*W
    for k in range(1, d+1):
        fac = [F(comb(5, j)*k**(5-j)) for j in range(0, 5)] + [F(0)]*(W-4)   # (w+k)^5 - w^5
        den = pmul(den, fac[:W+1])
    ratio = pmul(num, pinv(den))
    for m in range(W+1):
        ptilde[m] = ptilde[m] + S({(d, 0): ratio[m]})
p = []
for k in range(W+1):
    tot = S()
    lj = S.const(1)
    for j in range(0, k+1):
        tot = tot + lj * ptilde[k-j] * F(1, factorial(j))
        lj = lj * l
    p.append(tot)

C0 = S({(d, 0): F(factorial(5*d), factorial(d)**5) for d in range(N+1)})
I1t = S({(d, 0): 5*F(factorial(5*d), factorial(d)**5)*sum(F(1, m) for m in range(d+1, 5*d+1)) for d in range(1, N+1)})
assert (p[0] - C0).iszero() and (p[1] - (C0*l + I1t)).iszero()
print("p_0 = C_0 and p_1 = I_1: ok")

def PF4(f):
    g = f
    for k in range(1, 5): g = 5*g.D() + k*g
    return f.Dk(4) - 5*q*g
for k in range(0, 5):
    print("PF4 p_%d =" % k, PF4(p[k]))
print("D PF4 p_5 - p_0 =", PF4(p[5]).D() - p[0], " (coefficient of z^0 of PF5(I/z)=0)")

# L = (1-5^5 q)^(-1/5)
def binom_series(x, alpha):
    # (1+x)^alpha, x in q Q[[q]]
    res = S.const(1); term = S.const(1)
    for n in range(1, N+1):
        term = term * x * ((alpha - (n-1))/n)
        res = res + term
    return res
L = binom_series(-(5**5)*q, F(-1, 5))
L5 = L*L*L*L*L
assert (L5*(1 - 5**5*q) - 1).iszero()
assert (L.D() - (L5*L - L)*F(1, 5)).iszero()

# Wronskian
def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = S()
    for j in range(n):
        if M[0][j].iszero(): continue
        tot = tot + (-1)**j * M[0][j] * det([row[:j] + row[j+1:] for row in M[1:]])
    return tot
Wr = det([[p[k].Dk(j) for k in range(4)] for j in range(4)])
print("Wr(p_0..p_3) - L^10 =", Wr - L5*L5)

# mirror map and Birkhoff recursion
C0i = C0.inv()
T = l + I1t*C0i
C1 = T.D(); C1i = C1.inv()
dT = lambda f: C1i * f.D()
P0 = [p[k]*C0i for k in range(W+1)]
assert (P0[1] - T).iszero()
k1 = dT(dT(P0[2]))
print("kappa_1 - L^5/(C0^2 C1^3) =", k1 - L5*C0i*C0i*C1i*C1i*C1i)
k1i = k1.inv()
k2 = dT(k1i*dT(dT(P0[3])))
print("kappa_2 - 1 =", k2 - 1)
k3 = dT(dT(k1i*dT(dT(P0[4]))))
print("kappa_3 - C0/C1 =", k3 - C0*C1i)
# kappa_4 from the next step: J_4 = sum p^(4)_k H^k z^(4-k), z dT J_4 = kappa_4 J_0  => kappa_4 = dT p^(4)_5 (coefficient of H^0 z^0)
p4_5 = dT(dT(k1i*dT(dT(P0[5])))) * k3.inv()
k4 = dT(p4_5)
print("kappa_4 - C0/C1 =", k4 - C0*C1i, " (independent check; the note gets kappa_4 = kappa_3 from the Frobenius property)")
print("kappa_3 =", k3)
