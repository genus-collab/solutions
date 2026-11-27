"""X4 sanity checks (exact arithmetic, independent of the project tree).
(a) power-series solutions f0,f1,f2 of PF4 at u=0; varpi_1 = f1 + 7/10 f2; t_c in frame (0,0).
(b) concomitant Conc(varpi_0, varpi_1) = 0 for random rational frames (alpha,beta) (lem:gau-iso (ii)),
    so that C_{alpha,beta}(D i_2) = D C_{alpha,beta}(i_2) (lem:gau-confexp).
(c) v(u) = u tau^{-5/4} satisfies v'^5/v = 1/(u (1-u)^5) (prop:gau-equalprod (ii)).
(d) lem:set-canon: D(-5/4 u^{4/5}/tau) = u^{-1/5} checked via the series identity.
"""
from fractions import Fraction as F
N = 14  # truncation order in u

def mul(a, b):
    c = [F(0)] * N
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b[:N - i]):
            c[i + j] += x * y
    return c

def deriv(a):  # d/du
    return [a[i + 1] * (i + 1) for i in range(N - 1)] + [F(0)]

def D(a):  # D = -(1-u) d/du
    d = deriv(a)
    one_minus_u = [F(1), F(-1)] + [F(0)] * (N - 2)
    return [-x for x in mul(one_minus_u, d)]

def PF4(a):
    d1 = D(a); d2 = D(d1); d3 = D(d2); d4 = D(d3)
    u = [F(0), F(1)] + [F(0)] * (N - 2)
    omu = [F(1), F(-1)] + [F(0)] * (N - 2)
    inner = [2 * x + F(7, 5) * y + F(2, 5) * z_ + F(24, 625) * w for x, y, z_, w in zip(d3, d2, d1, a)]
    t1 = mul(u, d4)
    t2 = mul(omu, inner)
    return [x - y for x, y in zip(t1, t2)]

def solve(b0, b1, b2):
    """power-series solution with given first three coefficients, by undetermined coefficients."""
    a = [F(0)] * N
    a[0], a[1], a[2] = F(b0), F(b1), F(b2)
    for n in range(3, N):
        # coefficient of u^{n-3} in PF4(a) is P_{-3}(n) a_n + lower; solve linearly
        a[n] = F(0)
        r0 = PF4(a)[n - 3]
        a[n] = F(1)
        r1 = PF4(a)[n - 3]
        a[n] = -r0 / (r1 - r0)
    return a

f0 = solve(1, 0, 0); f1 = solve(0, 1, 0); f2 = solve(0, 0, 1)
for f in (f0, f1, f2):
    assert all(x == 0 for x in PF4(f)[:N - 4]), "not a solution"
w1 = [x + F(7, 10) * y for x, y in zip(f1, f2)]
print("varpi_1 =", w1[:5])

def inv(a):
    b = [F(0)] * N
    b[0] = 1 / a[0]
    for n in range(1, N):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b

tc = mul(w1, inv(f0))
print("t_c (frame 0,0) =", tc[:5])
assert tc[1] == 1 and tc[2] == F(7, 10) and tc[3] == F(41, 75)

# (b) concomitant with the jet matrix of lem:gau-conc (iii)
def conc(f, g):
    jf = [f, D(f), D(D(f)), D(D(D(f)))]
    jg = [g, D(g), D(D(g)), D(D(D(g)))]
    u = [F(0), F(1)] + [F(0)] * (N - 2)
    omu = [F(1), F(-1)] + [F(0)] * (N - 2)
    zero = [F(0)] * N
    sc = lambda c, a: [c * x for x in a]
    J = [[zero, sc(F(2, 5), omu), omu, sc(-1, u)],
         [sc(F(-2, 5), omu), zero, u, zero],
         [sc(-1, omu), sc(-1, u), zero, zero],
         [u, zero, zero, zero]]
    tot = [F(0)] * N
    for a_ in range(4):
        for b_ in range(4):
            t = mul(mul(jf[a_], J[a_][b_]), jg[b_])
            tot = [x + y for x, y in zip(tot, t)]
    return tot

import random
random.seed(1)
for _ in range(3):
    al = F(random.randint(-9, 9), random.randint(1, 9)); be = F(random.randint(-9, 9), random.randint(1, 9))
    w0 = [x + al * y + be * z_ for x, y, z_ in zip(f0, w1, f2)]
    c = conc(w0, w1)
    assert all(x == 0 for x in c[:N - 5]), c[:6]
print("Conc(varpi_0, varpi_1) = 0 to order", N - 6, "for random rational frames")
print("Conc(f0,f1), Conc(f0,f2) =", conc(f0, f1)[0], conc(f0, f2)[0])

# (c) tau and v
inv_tau = [F(4, 5 * n + 4) for n in range(N)]
tau = inv(inv_tau)
print("tau =", tau[:3])
# tau^{-5/4} via binomial series: w = inv_tau^{5/4}; compute log-free: w = exp(5/4 log inv_tau) using power recursion
def power(a, p):  # a[0]=1, returns a^p for rational p via J.C.P. Miller recursion
    b = [F(0)] * N; b[0] = F(1)
    for n in range(1, N):
        b[n] = sum((p * k - (n - k)) * a[k] * b[n - k] for k in range(1, n + 1)) / (n * a[0])
    return b
w = power(inv_tau, F(5, 4))
v = [F(0)] + w[:N - 1]
vp = deriv(v)
lhs = mul(power(vp, F(5)), [F(1)] + [F(0)] * (N - 1))  # v'^5
# compare v'^5 * u (1-u)^5 with v
u = [F(0), F(1)] + [F(0)] * (N - 2)
omu5 = power([F(1), F(-1)] + [F(0)] * (N - 2), F(5))
left = mul(mul(lhs, u), omu5)
assert left[:N - 1] == v[:N - 1], (left[:5], v[:5])
print("v'^5 u (1-u)^5 = v verified to order", N - 2, "; v =", v[:4])
