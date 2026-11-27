"""Prop 2.12: exact y-series through high order directly from the I-function (eq:Ifun) as printed in BS, with no shortcut;
Lemma 1.2(c): float check of BS line 5145 against t_c.  Evidence only; controls included."""
from fractions import Fraction as Fr
from math import factorial, log, pi, sqrt
import cmath

N = 60

def mul(a, b, n=N):
    out = [Fr(0)] * (n + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j > n: break
            out[i + j] += x * y
    return out

def inv(a, n=N):
    out = [Fr(0)] * (n + 1)
    out[0] = 1 / a[0]
    for k in range(1, n + 1):
        out[k] = -sum(a[i] * out[k - i] for i in range(1, k + 1)) * out[0]
    return out

def theta(a):
    return [i * x for i, x in enumerate(a)]

# I-function coefficients from eq:Ifun: treat H, s=1/z as commuting nilpotents to the needed order.
# summand_d = y^d prod_{m=0}^{3d-1} (-m - 3 H s) / prod_j prod_{m=1}^d (m + (H+alpha_j)s)
# expand as a polynomial in (H, s) keeping the total s-degree <= 2  (z y^{H/z} multiplies by z) with alpha generic numeric
alpha = [Fr(2, 7), Fr(-3, 5), Fr(11, 13)]
e1 = sum(alpha)

def poly_mul(p, q):
    out = {}
    for (a, b), x in p.items():
        for (c, d), y in q.items():
            if b + d <= 2:
                out[(a + c, b + d)] = out.get((a + c, b + d), Fr(0)) + x * y
    return out

def poly_inv_linear(c0, cH, cs):
    """1/(c0 + cH H s + cs s), series in s to order 2.  Terms: H^a s^b"""
    # =1/c0 * 1/(1+u), u=(cH H + cs) s /c0
    u1 = {(1, 1): cH / c0, (0, 1): cs / c0}
    out = {(0, 0): Fr(1)}
    u = {(0, 0): Fr(1)}
    sign = 1
    for k in range(1, 3):
        u = poly_mul(u, u1)
        sign = -sign
        for key, val in u.items():
            out[key] = out.get(key, Fr(0)) + sign * val
    return {k: v / c0 for k, v in out.items()}

def summand(d):
    p = {(0, 0): Fr(1)}
    for m in range(0, 3 * d):
        p = poly_mul(p, {(0, 0): Fr(-m), (1, 1): Fr(-3)})
    for a in alpha:
        for m in range(1, d + 1):
            p = poly_mul(p, poly_inv_linear(Fr(m), Fr(1), a))
    return p

# I = z y^{H/z} Sum_d y^d summand_d.
# Write I/z y^{-H/z} = 1 + Sum_d y^d (c1_d H s + s^2 (c2_d H + c3_d H^2) + ...)
c1 = [Fr(0)] * (N + 1)
c2 = [Fr(0)] * (N + 1)
c3 = [Fr(0)] * (N + 1)
for d in range(1, N + 1):
    p = summand(d)
    # no H^0 terms
    assert p.get((0, 0), 0) == 0 and p.get((0, 1), 0) == 0 and p.get((0, 2), 0) == 0
    c1[d] = p.get((1, 1), Fr(0))
    c2[d] = p.get((1, 2), Fr(0))
    c3[d] = p.get((2, 2), Fr(0))
# I = z(1 + H s L + H^2 s^2 L^2/2 + ...)(1 + H s C1 + s^2 (H C2 + H^2 C3)) with L=log y
# coefficient z^0 H^1: L + C1 =: I_0^[1];  z^-1: H^1: C2 ; H^2: L^2/2 + L C1 + C3
# theta I_0^[1] = 1 + theta C1 = A
A = [a for a in theta(c1)]; A[0] = Fr(1)
# theta I_1^[1]/e1 = G
G = [x / e1 for x in theta(c2)]
# theta I_1^[2] = L A + (C1 + theta C3)  [theta(L^2/2)=L; theta(L C1)=C1+L theta C1]
h = [x + y for x, y in zip(c1, theta(c3))]
Ainv = inv(A)
# D = theta(P_I/A)= theta(L + h/A) = 1 + theta(h/A)
D = theta(mul(h, Ainv)); D[0] += 1
Ey = theta(mul(G, Ainv))
X = [Fr((-27) ** k) for k in range(N + 1)]
lhs1 = mul(mul(A, A), D)
print("A^2 D = X through y^%d:" % N, lhs1 == X)
print("2E(y) = D - X through y^%d:" % N, all(2 * Ey[k] == D[k] - X[k] for k in range(N + 1)))
print("A = Sum a_d (-y)^d:", all(A[d] == Fr(factorial(3 * d), factorial(d) ** 3) * (-1) ** d for d in range(N + 1)))
Hn = lambda n: sum(Fr(1, k) for k in range(1, n + 1))
print("G = -Sum a_d H_d (-y)^d:", all(G[d] == -Fr(factorial(3 * d), factorial(d) ** 3) * Hn(d) * (-1) ** d for d in range(N + 1)))
print("h = Sum 3 a_d (H_3d-H_d)(-y)^d:", all(h[d] == 3 * Fr(factorial(3 * d), factorial(d) ** 3) * (Hn(3 * d) - Hn(d)) * (-1) ** d for d in range(N + 1)))
# (K): Sum a_d(3H_3d-H_d)(-y)^d = -A log(1+27y)
logD = [Fr(0)] + [Fr((-1) ** (k + 1) * 27 ** k, k) for k in range(1, N + 1)]
rhs = [-x for x in mul(A, logD)]
print("(K) through y^%d:" % N, all(rhs[d] == Fr(factorial(3 * d), factorial(d) ** 3) * (3 * Hn(3 * d) - Hn(d)) * (-1) ** d for d in range(N + 1)))
# dictionary: C1i=1/A, d2=E(C1i^2-1)/(2 Dc) = e1 E(y)
# control: perturb h
h2 = list(h); h2[3] += 1
D2 = theta(mul(h2, Ainv)); D2[0] += 1
print("control (h perturbed at y^3): A^2 D = X?", mul(mul(A, A), D2) == X)

# Lemma 1.2(c): floats.  BS 5145 at real y in (-1/27,0) with log y = log|y| + i pi s
def tcf(y, s=+1, nmax=4000):
    L = complex(log(abs(y)), s * pi)
    a = 1.0
    I01 = L
    I12 = L * L / 2
    term = 1.0
    # (3d-1)!/(d!)^3 = a_d/(3d)
    ad = 1.0
    for d in range(1, nmax):
        ad = ad * 3 * (3 * d - 1) * (3 * d - 2) / d ** 2
        H3dm1 = None
        term = ad / (3 * d) * (-y) ** d
        I01 += 3 * term
        I12 += 3 * L * term + 9 * term * (hs[3 * d - 1] - hs[d])
    return I01, I12

hs = [0.0]
for k in range(1, 3 * 4000 + 2):
    hs.append(hs[-1] + 1.0 / k)

# t_c from exact series in x = 1+27 y
M = 200
wc = [1.0]
for n in range(M):
    wc.append((n * (n - 1) + 2 * n + 2 / 9) * wc[n] / ((n + 1) ** 2))
# t_c = int wc/(1-x)
wc_over = []
s = 0.0
for n in range(M + 1):
    s += wc[n]
    wc_over.append(s)
def tc_of(x):
    return sum(wc_over[n] * x ** (n + 1) / (n + 1) for n in range(M))

for y in [-1 / 54, -1 / 40, -1 / 100, -1 / 33]:
    x = 1 + 27 * y
    for s_ in (+1, -1):
        I01, I12 = tcf(y, s_)
        val = sqrt(3) / (2 * pi) * (I12 - 1j * pi * I01 - pi ** 2)
        print(f"y={y:.5f} branch s={s_:+d}: BS t_cf = {val:.8f}; t_c(x) = {tc_of(x):.8f}")
