# Check thm:con-Pa(ii)/(iv) at a rational L: c_n(L) = L^{4n} ph_n(L^{-5}) for n = 1, 2, 3
# using the Taylor data of def:gam-c (first Mellin-Barnes summand), Wick contraction (eq:fsp-univ).
from fractions import Fraction as Fr
from math import comb, factorial
import sys

L = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2)
K = 3
TD = 6 * K + 4
# power sums p_k(L) = sum_j (zeta^j - L)^{-k}: generating function 5(L+t)^4/(1-(L+t)^5) = sum_k p_k t^{k-1}
def series_inv(a, n):
    r = [Fr(0)] * n
    r[0] = 1 / a[0]
    for i in range(1, n):
        s = sum(a[j] * r[i - j] for j in range(1, min(i, len(a) - 1) + 1))
        r[i] = -s / a[0]
    return r
N = TD + 4
num = [5 * comb(4, i) * L ** (4 - i) for i in range(5)] + [Fr(0)] * N
den = [Fr(1) - L ** 5] + [-comb(5, i) * L ** (5 - i) for i in range(1, 6)] + [Fr(0)] * N
inv = series_inv(den[:N], N)
gen = [sum(num[j] * inv[i - j] for j in range(i + 1)) for i in range(N)]
p = {k: gen[k - 1] for k in range(1, N)}
B2k = {1: Fr(1, 6), 2: Fr(-1, 30), 3: Fr(1, 42), 4: Fr(-1, 30)}
bk = {k: B2k[k] / (2 * k * (2 * k - 1)) for k in B2k}
# phase Taylor coefficients f_d (d>=2)
f = {}
for k in range(1, TD):
    f[k + 1] = (5 * Fr((-1) ** (k + 1)) * L ** (-k) + p[k]) / (k * (k + 1))
C = -1 / (2 * f[2])
assert C == -L * (1 - L ** 5) / 5
F3 = {d: f[d] for d in range(3, TD + 1)}
# amplitude log: S1(L+t)-S1(L) + sum_k z^{2k-1} S_{2k}(L+t)
def binom_gen(alpha, d):
    r = Fr(1)
    for i in range(d):
        r *= (alpha - i)
    return r / factorial(d)
# polys in (z, t): dict
def mul(a, b, tmax, zmax):
    r = {}
    for (za, ta), ca in a.items():
        for (zb, tb), cb in b.items():
            if za + zb > zmax or ta + tb > tmax: continue
            k = (za + zb, ta + tb)
            r[k] = r.get(k, 0) + ca * cb
    return {k: v for k, v in r.items() if v != 0}
loga = {}
for k in range(1, TD + 1):
    loga[(0, k)] = (Fr((-1) ** (k + 1)) * L ** (-k) + p[k]) / (2 * k)
for k in range(1, 3):
    for r in range(0, TD + 1):
        c = bk[k] * binom_gen(1 - 2 * k, r) * (Fr(5) ** (1 - 2 * k) * L ** (1 - 2 * k - r) + Fr((-1) ** r) * p[2 * k - 1 + r])
        if c != 0:
            loga[(2 * k - 1, r)] = loga.get((2 * k - 1, r), 0) + c
amp = {(0, 0): Fr(1)}
term = {(0, 0): Fr(1)}
for n in range(1, TD + 2):
    term = {k: v / n for k, v in mul(term, loga, TD, K).items()}
    if not term: break
    for k, v in term.items():
        amp[k] = amp.get(k, 0) + v
def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r
F3p = [{0: Fr(1)}]
for r in range(1, 2 * K + 1):
    prev = F3p[-1]; new = {}
    for d1, c1 in prev.items():
        for d2, c2 in F3.items():
            if d1 + d2 <= TD:
                new[d1 + d2] = new.get(d1 + d2, 0) + c1 * c2
    F3p.append(new)
c = []
for k in range(K + 1):
    acc = Fr(0)
    for m in range(k + 1):
        am = {t: v for (z, t), v in amp.items() if z == m}
        for r in range(0, 2 * (k - m) + 1):
            deg = 2 * (k - m + r)
            s = Fr(0)
            for d1, c1 in F3p[r].items():
                d2 = deg - d1
                if d2 in am:
                    s += c1 * am[d2]
            acc += s * dfact(deg - 1) * C ** (deg // 2) / factorial(r)
    c.append(acc)
u = L ** (-5)
ph = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000) - Fr(126, 625) * u + Fr(67, 1250) * u ** 2]
print("L =", L)
for n in range(K + 1):
    print(" n=%d  c_n(L) = %s   L^{4n} ph_n(L^-5) = %s   equal: %s" % (n, c[n], L ** (4 * n) * ph[n], c[n] == L ** (4 * n) * ph[n]))
