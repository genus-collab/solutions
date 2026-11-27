# Independent sanity check (exact arithmetic) of the statement of prop:set-cyclicLR (i),
# of eq:set-McycLR, of cor:set-cyclicPF5 and of the first WKB coefficient.
# Elements of Q[[q]][l] (l = formal log q), truncated at q^N: dict j -> list of N+1 Fractions.
from fractions import Fraction as Fr
from math import factorial, comb

N = 6  # q-order

def zero():
    return {}

def const(c):
    return {0: [Fr(c)] + [Fr(0)] * N}

def clean(a):
    return {j: v for j, v in a.items() if any(x != 0 for x in v)}

def add(a, b, s=1):
    r = {j: list(v) for j, v in a.items()}
    for j, v in b.items():
        if j not in r:
            r[j] = [Fr(0)] * (N + 1)
        for n in range(N + 1):
            r[j][n] += s * v[n]
    return clean(r)

def mul(a, b):
    r = {}
    for j1, v1 in a.items():
        for j2, v2 in b.items():
            j = j1 + j2
            if j not in r:
                r[j] = [Fr(0)] * (N + 1)
            for n1 in range(N + 1):
                if v1[n1] == 0:
                    continue
                for n2 in range(N + 1 - n1):
                    r[j][n1 + n2] += v1[n1] * v2[n2]
    return clean(r)

def scal(c, a):
    return clean({j: [Fr(c) * x for x in v] for j, v in a.items()})

def D(a):
    r = {}
    for j, v in a.items():
        r.setdefault(j, [Fr(0)] * (N + 1))
        for n in range(N + 1):
            r[j][n] += n * v[n]
        if j >= 1:
            r.setdefault(j - 1, [Fr(0)] * (N + 1))
            for n in range(N + 1):
                r[j - 1][n] += j * v[n]
    return clean(r)

def inv(a):
    # inverse of an l-free unit
    assert set(a.keys()) <= {0}
    v = a[0]
    assert v[0] != 0
    w = [Fr(0)] * (N + 1)
    w[0] = 1 / v[0]
    for n in range(1, N + 1):
        w[n] = -sum(v[k] * w[n - k] for k in range(1, n + 1)) / v[0]
    return {0: w}

def powser(a, e):
    # a = 1 + O(q), l-free; rational exponent e: binomial series
    v = a[0]
    assert v[0] == 1
    x = {0: [Fr(0)] + v[1:]}
    r = const(1)
    term = const(1)
    for k in range(1, N + 1):
        term = scal(Fr(e - (k - 1)) / k, mul(term, x))
        r = add(r, term)
    return r

def eq(a, b):
    return clean(add(a, b, -1)) == {}

# ---- the series P(w,q) = sum_d q^d prod_{k<=5d}(5w+k)/prod_{k<=d}((w+k)^5-w^5), in w up to w^4
W = 5
def wmul(a, b):
    r = [Fr(0)] * W
    for i in range(W):
        for j in range(W - i):
            r[i + j] += a[i] * b[j]
    return r
def winv(a):
    r = [Fr(0)] * W
    r[0] = 1 / a[0]
    for n in range(1, W):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r

ptilde = [[Fr(0)] * (N + 1) for _ in range(W)]  # ptilde[m][d]
for d in range(N + 1):
    num = [Fr(1)] + [Fr(0)] * (W - 1)
    for k in range(1, 5 * d + 1):
        num = wmul(num, [Fr(k), Fr(5)] + [Fr(0)] * (W - 2))
    den = [Fr(1)] + [Fr(0)] * (W - 1)
    for k in range(1, d + 1):
        # (w+k)^5 - w^5 = sum_{i<5} C(5,i) k^(5-i) w^i
        f = [Fr(comb(5, i) * k ** (5 - i)) for i in range(5)]
        den = wmul(den, f[:W])
    c = wmul(num, winv(den))
    for m in range(W):
        ptilde[m][d] = c[m]

def lpow(j):
    return {j: [Fr(1)] + [Fr(0)] * N}

p = []
for k in range(W):
    acc = zero()
    for j in range(k + 1):
        acc = add(acc, scal(Fr(1, factorial(j)), mul(lpow(j), {0: list(ptilde[k - j])})))
    p.append(acc)

C0 = {0: [Fr(factorial(5 * d), factorial(d) ** 5) for d in range(N + 1)]}
I1t = {0: [Fr(0)] + [5 * Fr(factorial(5 * d), factorial(d) ** 5) * sum(Fr(1, m) for m in range(d + 1, 5 * d + 1)) for d in range(1, N + 1)]}
assert eq(p[0], C0), "p0 != C0"
I1 = add(mul(C0, lpow(1)), I1t)
assert eq(p[1], I1), "p1 != I1"
C0i = inv(C0)
T = add(lpow(1), mul(I1t, C0i))
C1 = D(T)
C1i = inv(C1)
qser = {0: [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)}
u = add(const(1), scal(-5 ** 5, qser))
L = powser(u, Fr(-1, 5))
L5 = inv(u)

def dT(a):
    return mul(C1i, D(a))

def PF4(a):
    # D^4 - 5 q prod_{k=1..4}(5D+k)
    d1 = D(a); d2 = D(d1); d3 = D(d2); d4 = D(d3)
    b = add(add(add(add(scal(3125, d4), scal(6250, d3)), scal(4375, d2)), scal(1250, d1)), scal(120, a))
    return add(d4, mul(qser, b), -1)

print("PF4 p_k (k=0..4):", [PF4(pk) for pk in p])

P0 = [mul(pk, C0i) for pk in p]
k1 = dT(dT(P0[2]))
k1_claim = mul(L5, mul(mul(C0i, C0i), mul(C1i, mul(C1i, C1i))))
print("kappa_1 = L^5/(C0^2 C1^3):", eq(k1, k1_claim))
k1i = inv(k1)
k2 = dT(mul(k1i, dT(dT(P0[3]))))
print("kappa_2 = 1:", eq(k2, const(1)))
k3 = dT(dT(mul(k1i, dT(dT(P0[4])))))
print("kappa_3 = C0/C1:", eq(k3, mul(C0, C1i)))
# product of kappas = (L/C1)^5
prod = mul(k1_claim, mul(mul(C0, C1i), mul(C0, C1i)))
print("prod kappa = (L/C1)^5:", eq(prod, mul(L5, powser(C1, -5))))
# Delta_0^{-1} = eta(e_0,e_0) = kappa_0 kappa_1 kappa_2 mu^{-3} = L^2/C0^2
mu = mul(L, C1i)
lhs = mul(k1_claim, powser(mul(mu, const(1)) if False else mu, -3)) if mu[0][0] == 1 else None
print("eta(e_0,e_0) = L^2/C0^2:", eq(lhs, mul(mul(L, L), mul(C0i, C0i))))
# Wronskian L^10 (lem:set-wronskian)
rows = [[pk for pk in p[:4]]]
for j in range(1, 4):
    rows.append([D(x) for x in rows[-1]])
def det(m):
    n = len(m)
    if n == 1:
        return m[0][0]
    acc = zero()
    for c in range(n):
        minor = [r[:c] + r[c + 1:] for r in m[1:]]
        acc = add(acc, mul(m[0][c], det(minor)), 1 if c % 2 == 0 else -1)
    return acc
print("Wr = L^10:", eq(det(rows), powser(u, -2)))
# D L = (L^6 - L)/5
print("D L = (L^6-L)/5:", eq(D(L), scal(Fr(1, 5), add(mul(L5, L), L, -1))))
# colour-0 WKB amplitude: A = L - (3/20) L^5 z + ... : check z^2 coefficient of PF5^{(0)} A vanishes.
# PF5^{(0)} = (zD+L)^5 - 1 - q prod_{k=1..5}(5(zD+L)+k z); represent series in z as list of coefficient elements.
ZN = 3
def zadd(a, b, s=1):
    return [add(a[i], b[i], s) for i in range(ZN)]
def Mop(a):  # (zD + L) a
    r = [mul(L, a[i]) for i in range(ZN)]
    for i in range(1, ZN):
        r[i] = add(r[i], D(a[i - 1]))
    return r
def Nop(a, k):  # (5(zD+L) + k z) a
    m = Mop(a)
    r = [scal(5, m[i]) for i in range(ZN)]
    for i in range(1, ZN):
        r[i] = add(r[i], scal(k, a[i - 1]))
    return r
def PF5col(a):
    x = a
    for _ in range(5):
        x = Mop(x)
    y = a
    for k in range(1, 6):
        y = Nop(y, k)
    y = [mul(qser, t) for t in y]
    return zadd(zadd(x, a, -1), y, -1)
A = [L, scal(Fr(-3, 20), mul(L5, mul(L, const(1)))) if False else scal(Fr(-3, 20), mul(L5, const(1))), zero()]
# note: A_1 = -(3/20) L^5
res = PF5col(A)
print("PF5^(0)(L - 3/20 L^5 z) z^0,z^1,z^2 coefficients vanish:", [clean(r) == {} for r in res])
A2 = [L, add(scal(Fr(-3, 20), L5), scal(Fr(1, 7), L)), zero()]
res2 = PF5col(A2)
print("same with A_1 + (1/7) L (kernel direction):", [clean(r) == {} for r in res2])
A3 = [L, scal(Fr(3, 20), L5), zero()]
res3 = PF5col(A3)
print("wrong sign +3/20 L^5:", [clean(r) == {} for r in res3])
