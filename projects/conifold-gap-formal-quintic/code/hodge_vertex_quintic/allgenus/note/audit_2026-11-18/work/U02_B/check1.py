# Independent sanity checks for unit U02 (lens B).  Exact arithmetic with fractions.
from fractions import Fraction as Fr

# ---------- Laurent polynomials / series in u as dicts {exp: coeff} ----------
def add(a, b, cb=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + cb * v
        if r[k] == 0:
            del r[k]
    return r

def scal(a, c):
    return {k: v * c for k, v in a.items() if v * c != 0}

def mulu(a, n):
    return {k + n: v for k, v in a.items()}

# an element of Q((u))[ell] is a dict {ellpower: laurent}; D = -(1-u) d/du, D ell = 1 - 1/u
def Dser(a):
    # -(1-u) d/du
    r = {}
    for k, v in a.items():
        if k != 0:
            r[k - 1] = r.get(k - 1, 0) - k * v
            r[k] = r.get(k, 0) + k * v
    return {k: v for k, v in r.items() if v != 0}

def Dfull(x):
    # x = {p: series}; D(s ell^p) = D(s) ell^p + p s (1-1/u) ell^{p-1}
    r = {}
    for p, s in x.items():
        r[p] = add(r.get(p, {}), Dser(s))
        if p >= 1:
            t = add(s, mulu(s, -1), -1)
            r[p - 1] = add(r.get(p - 1, {}), scal(t, p))
    return {p: s for p, s in r.items() if s}

def addf(x, y, c=1):
    r = dict(x)
    for p, s in y.items():
        r[p] = add(r.get(p, {}), s, c)
    return {p: s for p, s in r.items() if s}

def scalf(x, c):
    return {p: scal(s, c) for p, s in x.items()}

def muluf(x, n):
    return {p: mulu(s, n) for p, s in x.items()}

def PF4(x):
    d1 = Dfull(x); d2 = Dfull(d1); d3 = Dfull(d2); d4 = Dfull(d3)
    inner = addf(addf(addf(scalf(d3, 2), scalf(d2, Fr(7, 5))), scalf(d1, Fr(2, 5))), scalf(x, Fr(24, 625)))
    one_minus_u = addf(inner, muluf(inner, 1), -1)
    return addf(muluf(d4, 1), one_minus_u, -1)

# P_j(n): coefficient of u^{n+j} in PF4 u^n
def P(j, n):
    r = PF4({0: {n: Fr(1)}})
    return r.get(0, {}).get(n + j, Fr(0))

print("support of PF4 u^7:", sorted(PF4({0: {7: Fr(1)}})[0].keys()))
for x in range(-3, 8):
    x_ = Fr(x)
    assert P(-3, x) == x_ * (x_ - 1) ** 2 * (x_ - 2), x
    assert P(-2, x) == -x_ * (x_ - 1) * (4 * x_ * x_ - 6 * x_ + Fr(17, 5)), x
print("P_{-3}, P_{-2} confirmed on x=-3..7")
print("P_{-1}(0) =", P(-1, 0))
# P_{-1} polynomial of degree <= 4: interpolate derivative at 0 via PF4(ell)
r = PF4({1: {0: Fr(1)}})
print("PF4(ell): ell^0 part:", r.get(0), " ell^1 part:", r.get(1))
# derivative P_j'(n) = coefficient of u^{n+j} ell^0 in PF4(u^n ell)
def Pp(j, n):
    r = PF4({1: {n: Fr(1)}})
    return r.get(0, {}).get(n + j, Fr(0))
print("P_{-3}'(0), P_{-3}'(1), P_{-3}'(2):", Pp(-3, 0), Pp(-3, 1), Pp(-3, 2))
print("P_{-2}'(0), P_{-2}'(1):", Pp(-2, 0), Pp(-2, 1))
print("P_{-1}'(0):", Pp(-1, 0))

# ---------- power series solutions f_0,f_1,f_2 mod u^N ----------
N = 9
def solve(b0, b1, b2):
    b = [Fr(b0), Fr(b1), Fr(b2)]
    # coefficient of u^M in PF4 g: sum_j P_j(M-j) b_{M-j}, j=-3..1
    for M in range(0, N - 3):
        s = Fr(0)
        for j in (-2, -1, 0, 1):
            n = M - j
            if 0 <= n < len(b):
                s += P(j, n) * b[n]
        b.append(-s / P(-3, M + 3))
    return b
f0 = solve(1, 0, 0); f1 = solve(0, 1, 0); f2 = solve(0, 0, 1)
print("f0:", f0[:6]); print("f1:", f1[:6]); print("f2:", f2[:6])
# check the coefficients of negative powers vanish as claimed
for f in (f0, f1, f2):
    r = PF4({0: {k: c for k, c in enumerate(f) if c != 0}})
    bad = {k: v for k, v in r.get(0, {}).items() if k < N - 4}
    assert not bad, bad
w1 = [a + Fr(7, 10) * c for a, c in zip(f1, f2)]
print("varpi_1:", w1[:6])
def ser_mul(a, b):
    r = [Fr(0)] * N
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < N:
                r[i + j] += x * y
    return r
def ser_inv(a):
    r = [Fr(0)] * N
    r[0] = 1 / a[0]
    for n in range(1, N):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r
tc = ser_mul(w1, ser_inv(f0))
print("t_c (frame 0,0):", tc[:5])
# logarithmic solution: obstruction of g = b0 f0 + b1 f1 + b2 f2
def obstruction(g):
    r = PF4({1: {k: c for k, c in enumerate(g) if c != 0}})
    return [r.get(0, {}).get(m, Fr(0)) for m in (-3, -2, -1)], r.get(1, {})
for name, f in (("f0", f0), ("f1", f1), ("f2", f2), ("varpi_1", w1)):
    ob, l1 = obstruction(f)
    print("obstruction", name, ob, " ell-part low orders:", {k: v for k, v in l1.items() if k < N - 4})

# ---------- tau ----------
F21 = [Fr(4, 5 * n + 4) for n in range(N)]
tau = ser_inv(F21)
print("tau:", tau[:4])

# ---------- T-system, colour 0 ----------
M = 8
def T(k, F):
    # (-1 + t/5 (k - 4 t d/dt)) F
    r = [Fr(0)] * (M + 1)
    for n, c in enumerate(F):
        if n <= M:
            r[n] += -c
        if n + 1 <= M:
            r[n + 1] += Fr(k - 4 * n, 5) * c
    return r
def neg(F):
    return [-c for c in F]
def chain(F0):
    F1 = T(-1, F0); F2 = neg(T(-2, F1)); F3 = T(2, F2); F4 = neg(T(1, F3)); F0b = neg(T(0, F4))
    return F1, F2, F3, F4, F0b
# solve closing relation order by order
F0 = [Fr(1)] + [Fr(0)] * M
for n in range(1, M):
    # choose F0[n] so that coefficient of t^{n+1} in F0b - F0 vanishes
    F0[n] = Fr(0)
    _, _, _, _, F0b = chain(F0)
    resid0 = F0b[n + 1] - F0[n + 1]
    F0[n] = Fr(1)
    _, _, _, _, F0b = chain(F0)
    resid1 = F0b[n + 1] - F0[n + 1]
    F0[n] = -resid0 / (resid1 - resid0)
    assert resid1 - resid0 == 4 * n, (n, resid1 - resid0)
F1, F2, F3, F4, F0b = chain(F0)
print("F_0:", F0[:5])
print("F_b(0):", F0[0], F1[0], F2[0], F3[0], F4[0])
print("closing residual:", [F0b[k] - F0[k] for k in range(M)])
print("F_1:", F1[:4]); print("F_2:", F2[:4]); print("F_3:", F3[:4]); print("F_4:", F4[:4])

# h = exp(-sum b_k z^{2k-1} sum a^{1-2k}) to low order: first coefficient -3/20, second 9/800
import json
json.dump({"F": [[str(c) for c in F[:M]] for F in (F0, F1, F2, F3, F4)]}, open("Fs.json", "w"))
