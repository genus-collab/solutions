# Checks around def:set-confexp: the concomitant relation Conc(C_0, I_1) = 0 at q = 0 side,
# and Conc(varpi_0, varpi_1) = 0 at the conifold side; the image of W; D-compatibility of the image of i_3.
from fractions import Fraction as Fr
from math import factorial

N = 10
# ---------- large-radius side: series in q, with a formal log (dict {power of log: list}) ----------
def smul(a, b):
    r = [Fr(0)] * N
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j < N:
                r[i + j] += x * y
    return r
def sadd(a, b, c=1):
    return [x + c * y for x, y in zip(a, b)]
def Dq(a):  # q d/dq
    return [n * x for n, x in enumerate(a)]
C0 = [Fr(factorial(5 * d), factorial(d) ** 5) for d in range(N)]
It = [Fr(0)] + [5 * Fr(factorial(5 * d), factorial(d) ** 5) * sum(Fr(1, m) for m in range(d + 1, 5 * d + 1)) for d in range(1, N)]
# PF4 C0 = 0 ?
u = [Fr(1), Fr(-5 ** 5)] + [Fr(0)] * (N - 2)
one_minus_u = [Fr(0), Fr(5 ** 5)] + [Fr(0)] * (N - 2)
def jets(a):
    j = [a]
    for _ in range(4):
        j.append(Dq(j[-1]))
    return j
def PF4(a):
    j = jets(a)
    inner = sadd(sadd(sadd([2 * x for x in j[3]], [Fr(7, 5) * x for x in j[2]]), [Fr(2, 5) * x for x in j[1]]), [Fr(24, 625) * x for x in j[0]])
    return sadd(smul(u, j[4]), smul(one_minus_u, inner), -1)
print("PF4 C0 =", PF4(C0)[:6])
# I_1 = C0 log + It ; D^j I_1 = log * D^j C0 + j D^{j-1} C0 + D^j It
jc = jets(C0); ji = jets(It)
zero = [Fr(0)] * N
I1_log = [jc[j] for j in range(4)]
I1_reg = [sadd(ji[j], ([j * x for x in jc[j - 1]] if j >= 1 else zero)) for j in range(4)]
# PF4 I_1 = 0: log part = PF4 C0, regular part:
D4reg = sadd(ji[4], [4 * x for x in jc[3]])
inner = sadd(sadd(sadd([2 * x for x in I1_reg[3]], [Fr(7, 5) * x for x in I1_reg[2]]), [Fr(2, 5) * x for x in I1_reg[1]]), [Fr(24, 625) * x for x in I1_reg[0]])
print("PF4 I1 (regular part) =", sadd(smul(u, D4reg), smul(one_minus_u, inner), -1)[:6])
# jet matrix
def J(a, b):
    m = {(0, 1): [Fr(2, 5) * x for x in one_minus_u], (0, 2): one_minus_u, (0, 3): [-x for x in u],
         (1, 0): [-Fr(2, 5) * x for x in one_minus_u], (1, 2): u, (2, 0): [-x for x in one_minus_u],
         (2, 1): [-x for x in u], (3, 0): u}
    return m.get((a, b), zero)
conc_log = zero; conc_reg = zero
for a in range(4):
    for b in range(4):
        conc_log = sadd(conc_log, smul(smul(jc[a], J(a, b)), I1_log[b]))
        conc_reg = sadd(conc_reg, smul(smul(jc[a], J(a, b)), I1_reg[b]))
print("Conc(C0,I1): log part", conc_log[:6], " regular part", conc_reg[:6])
# W = C0^2 + C0 D It - D C0 It ; C1 = D T
W = sadd(sadd(smul(C0, C0), smul(C0, ji[1])), smul(jc[1], It), -1)
print("ev(W) low orders:", W[:3])

# ---------- conifold side: power series in u ----------
def Du(a):  # -(1-u) d/du on power series (list index = power of u)
    r = [Fr(0)] * N
    for k, x in enumerate(a):
        if k >= 1:
            r[k - 1] -= k * x
            r[k] += k * x
    return r
def P4coef(j, n):
    x = Fr(n)
    # use operator directly on monomial
    a = [Fr(0)] * (N + 6); a[n + 3] = Fr(1)  # shift by 3 to keep negative powers
    return None
# solve PF4 f = 0 by direct recursion using the operator on truncated series with enough extra terms
import itertools
def PF4u(a):
    j = [a]
    for _ in range(4):
        j.append(Du(j[-1]))
    uu = [Fr(0), Fr(1)] + [Fr(0)] * (N - 2)
    omu = [Fr(1), Fr(-1)] + [Fr(0)] * (N - 2)
    inner = sadd(sadd(sadd([2 * x for x in j[3]], [Fr(7, 5) * x for x in j[2]]), [Fr(2, 5) * x for x in j[1]]), [Fr(24, 625) * x for x in j[0]])
    return sadd(smul(uu, j[4]), smul(omu, inner), -1)
def solve(b0, b1, b2):
    f = [Fr(b0), Fr(b1), Fr(b2)] + [Fr(0)] * (N - 3)
    for n in range(3, N):
        r0 = PF4u(f)[n - 3]
        f[n] = Fr(1)
        r1 = PF4u(f)[n - 3]
        f[n] = -r0 / (r1 - r0)
    return f
f0 = solve(1, 0, 0); f1 = solve(0, 1, 0); f2 = solve(0, 0, 1)
w1 = sadd(f1, f2, Fr(7, 10))
print("f0:", f0[:5]); print("varpi1:", w1[:5])
uu = [Fr(0), Fr(1)] + [Fr(0)] * (N - 2)
omu = [Fr(1), Fr(-1)] + [Fr(0)] * (N - 2)
def Ju(a, b):
    m = {(0, 1): [Fr(2, 5) * x for x in omu], (0, 2): omu, (0, 3): [-x for x in uu],
         (1, 0): [-Fr(2, 5) * x for x in omu], (1, 2): uu, (2, 0): [-x for x in omu],
         (2, 1): [-x for x in uu], (3, 0): uu}
    return m.get((a, b), zero)
def jetsu(a):
    j = [a]
    for _ in range(3):
        j.append(Du(j[-1]))
    return j
def conc(f, g):
    jf = jetsu(f); jg = jetsu(g)
    acc = zero
    for a in range(4):
        for b in range(4):
            acc = sadd(acc, smul(smul(jf[a], Ju(a, b)), jg[b]))
    return acc[:N - 4]
print("Conc(f0,f1) =", conc(f0, f1))
print("Conc(f0,f2) =", conc(f0, f2))
print("Conc(f1,f2) =", conc(f1, f2))
print("Conc(f0,varpi1) =", conc(f0, w1))
print("Conc(f2,varpi1) =", conc(f2, w1))
# frame (alpha,beta) = (1/3, -2/7): t_c and D t_c(0)
al, be = Fr(1, 3), Fr(-2, 7)
w0 = sadd(sadd(f0, w1, al), f2, be)
def sinv(a):
    r = [Fr(0)] * N
    r[0] = 1 / a[0]
    for n in range(1, N):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) / a[0]
    return r
tc = smul(w1, sinv(w0))
print("t_c in frame (1/3,-2/7):", tc[:4], " expected u^2 coeff:", Fr(7, 10) - al)
print("Conc(varpi0,varpi1) =", conc(w0, w1))
