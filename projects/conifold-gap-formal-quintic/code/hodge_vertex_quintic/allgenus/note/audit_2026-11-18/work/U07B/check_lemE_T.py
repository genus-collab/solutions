# Independent check of lem:tra-E, cor:tra-edges, lem:tra-T on a random graded symplectic dressing.
# Truncated arithmetic: z in [-ZN, Z], u in [0, K).  v = v(u) = u + u^2/3 (v'(0)=1).
from fractions import Fraction as F
import random, itertools
random.seed(7)
K = 7          # u-truncation
Z = 14         # positive z truncation
ZN = 12        # negative z truncation
n5 = range(5)

def padd(a, b, s=1):
    r = dict(a)
    for k, c in b.items():
        r[k] = r.get(k, 0) + s * c
        if r[k] == 0: del r[k]
    return r
def pmul(a, b):
    r = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            i, j = i1 + i2, j1 + j2
            if j >= K or i > Z or i < -ZN: continue
            r[(i, j)] = r.get((i, j), 0) + c1 * c2
    return {k: c for k, c in r.items() if c != 0}
def pscal(c, a): return {k: c * x for k, x in a.items() if c * x != 0}
def zsub(a):  # z -> -z
    return {(i, j): (c if i % 2 == 0 else -c) for (i, j), c in a.items()}
def mzero(): return [[{} for _ in n5] for _ in n5]
def mid():
    m = mzero()
    for a in n5: m[a][a] = {(0, 0): F(1)}
    return m
def mmul(A, B):
    C = mzero()
    for a in n5:
        for b in n5:
            s = {}
            for c in n5:
                if A[a][c] and B[c][b]: s = padd(s, pmul(A[a][c], B[c][b]))
            C[a][b] = s
    return C
def madd(A, B, s=1): return [[padd(A[a][b], B[a][b], s) for b in n5] for a in n5]
def mT(A): return [[A[b][a] for b in n5] for a in n5]
def mzsub(A): return [[zsub(A[a][b]) for b in n5] for a in n5]
def mconst(M): return [[({(0, 0): F(M[a][b])} if M[a][b] != 0 else {}) for b in n5] for a in n5]
def mvec(A, x): return [ (lambda a: (lambda s: s)(sum_p([pmul(A[a][c], x[c]) for c in n5])))(a) for a in n5]
def sum_p(lst):
    s = {}
    for t in lst: s = padd(s, t)
    return s
eta = [[F(5) if (a + b) % 5 == 3 else F(0) for b in n5] for a in n5]
etai = [[eta[a][b] / 25 for b in n5] for a in n5]
ETA, ETAI = mconst(eta), mconst(etai)
# v(u)
vser = {(0, 1): F(1), (0, 2): F(1, 3)}
def vpow(e):
    r = {(0, 0): F(1)}
    for _ in range(e): r = pmul(r, vser)
    return r
# graded symplectic factor exp(a z^k), a of shift k, a^* = (-1)^{k+1} a  (a^* = eta^-1 a^T eta)
def shiftmat(k, coeffs):
    m = [[F(0)] * 5 for _ in n5]
    for b in n5: m[(b + k) % 5][b] = coeffs[b]
    return m
def adj(m):
    return [[sum(etai[a][c] * m[d][c] * eta[d][b] for c in n5 for d in n5) for b in n5] for a in n5]
def symm_gen(k):
    sgn = 1 if k % 2 == 1 else -1
    while True:
        c = [F(random.randint(-3, 3), random.randint(1, 3)) for _ in n5]
        m = shiftmat(k, c)
        ma = adj(m)
        g = [[(m[a][b] + sgn * ma[a][b]) / 2 for b in n5] for a in n5]
        if any(g[a][b] != 0 for a in n5 for b in n5): return g
def mexp_z(g, k, upow=0):
    # exp(g z^k u^upow), truncated
    X = [[({(k, upow): g[a][b]} if g[a][b] != 0 else {}) for b in n5] for a in n5]
    R = mid(); T = mid(); n = 1
    while True:
        T = mmul(T, X); T = [[pscal(F(1, n), T[a][b]) for b in n5] for a in n5]
        if all(not T[a][b] for a in n5 for b in n5): break
        R = madd(R, T); n += 1
    return R
# P0: diagonal symplectic, p_a p_{3-a} = 1, p_4 = +-1; with u-dependence
p0 = [{(0, 0): F(2), (0, 1): F(1)}, None, None, None, {(0, 0): F(1)}]
p0[1] = {(0, 0): F(3), (0, 2): F(-1, 2)}
# inverses mod u^K
def pinv_u(a):
    c0 = a[(0, 0)]; r = {(0, 0): 1 / c0}
    rest = padd(a, {(0, 0): c0}, -1)
    term = {(0, 0): 1 / c0}
    x = pscal(-1 / c0, rest)
    for _ in range(K):
        term = pmul(term, x); r = padd(r, term)
    return r
p0[3] = pinv_u(p0[0]); p0[2] = pinv_u(p0[1])
P0 = mzero()
for a in n5: P0[a][a] = p0[a]
# R = product of graded symplectic factors, some with u
Rm = mid()
for (k, up) in [(1, 0), (2, 1), (3, 0), (1, 2), (5, 1), (4, 0), (6, 1), (7, 0)]:
    Rm = mmul(Rm, mexp_z(symm_gen(k), k, up))
P = mmul(P0, Rm)   # dressing in lowered components (acts on lowered comps), graded
# check (Sp): P(-z)^T eta P(z) = eta
S = mmul(mmul(mT(mzsub(P)), ETA), P)
print("Sp ok:", all(padd(S[a][b], ETA[a][b], -1) == {} for a in n5 for b in n5))
# vanishing solution YC at v(u)
from math import factorial
nu = [0, -1, 3, 2, 1]; eps = [1, -1, -1, 1, 1]
YC = [{} for _ in n5]
for b in n5:
    for n in range(0, ZN + 1):
        if (n - nu[b] - 1) % 5 != 0: continue
        m = n // 5
        c = F(eps[b] * (-1) ** n, 4 ** m * factorial(m) * factorial(n - m))
        t = pscal(c, vpow(n - m))
        YC[b] = padd(YC[b], {(-n, j): x for (i, j), x in t.items()})
def ptrunc_ok(): pass
# Pi_i(q) = -[z^i] <YC(-z), q(z)>
def Pi(i, q):
    s = {}
    for a in n5: s = padd(s, pmul(zsub(YC[a]), q[a]))
    return {j: -c for (ii, j), c in s.items() if ii == i}
def neg(d): return {k: -c for k, c in d.items()}
# Y' = P YC, check grading claim and T_P
Yp = mvec(P, YC)
TP = {j: -c for (i, j), c in Yp[0].items() if i == -1}
print("T_P mod u^K:", TP)
# vector form
Pvec = mmul(mmul(ETAI, P), ETA)
# inverse of Pvec: Pvec = P0vec Rvec; invert Rvec = 1 + N via geometric series
P0vec = mmul(mmul(ETAI, P0), ETA)
Rvec = mmul(mmul(ETAI, Rm), ETA)
def minv_unip(Rv):
    N = madd(Rv, mid(), -1); Inv = mid(); T = mid()
    for _ in range(Z + K + 2):
        T = mmul(T, N); T = [[neg(T[a][b]) for b in n5] for a in n5]
        if all(not T[a][b] for a in n5 for b in n5): break
        Inv = madd(Inv, T)
    return Inv
Rvi = minv_unip(Rvec)
P0vi = mzero()
for a in n5:
    P0vi[a][a] = pinv_u(P0vec[a][a])
Pvi = mmul(Rvi, P0vi)
mz1 = [{(1, 0): F(-1)}, {}, {}, {}, {}]
qT = mvec(Pvi, mz1)
d = P0[3][3]
print("d == P0vec[0][0]:", padd(d, P0vec[0][0], -1) == {})
# Lemma T(b): Pi_i(q_T) = (-1)^{i+1} [z^{i-1}] Y'_0  (valid for i small enough given truncation)
okb = True
for i in range(0, 4):
    lhs = Pi(i, qT)
    rhs = {j: ((-1) ** (i + 1)) * c for (ii, j), c in Yp[0].items() if ii == i - 1}
    lhs = {j: c for j, c in lhs.items() if j < K}
    if lhs != rhs: okb = False; print("Lemma T(b) mismatch at i=", i, lhs, rhs)
print("Lemma T(b) identity (without (A)) ok:", okb)
# Lemma T(d): T_P - p00 v in u^2
diff = padd(TP, {j: c for (i, j), c in pmul(P0[0][0], vser).items()}, -1)
print("T_P - p00 v lowest u-order:", min(diff) if diff else None)
# Lemma E: contraction of Omega_s(f,.) into first slot of E_R vs -R^{-1}(psi')(Rf)_{>=0}(psi')
# with R = Rvec (vector form), f = z^{-1} J^C,  J^C = eta^{-1} YC
JC = mvec(ETAI, YC)
f = [{(i - 1, j): c for (i, j), c in JC[a].items()} for a in n5]
# coefficients of R^{-1}(psi) in z
def zcoef_mat(M, k): return [[{(0, j): c for (i, j), c in M[a][b].items() if i == k} for b in n5] for a in n5]
L = 6  # psi' degree to test
KMAX = Z
Rinv_k = [zcoef_mat(Rvi, k) for k in range(KMAX + 1)]
def Ncoef(k, l):
    T = mmul(mmul(Rinv_k[k], ETAI), mT(Rinv_k[l]))
    T = [[neg(T[a][b]) for b in n5] for a in n5]
    if k == 0 and l == 0: T = madd(T, ETAI)
    return T
def Ecoef(k, l):
    S = mzero()
    for j in range(0, l + 1):
        if k + 1 + j > KMAX: raise Exception("trunc")
        Nc = Ncoef(k + 1 + j, l - j)
        S = madd(S, Nc, (-1) ** j)
    return S
# covector: Omega_s(f, psi^m e) = ([z^{-m-1}] f(-z))^T eta e
fm = [zsub(f[a]) for a in n5]
def cov(m): return [{(0, j): c for (i, j), c in fm[a].items() if i == -m - 1} for a in n5]
Mmax = ZN  # f has z^{-1}..z^{-ZN-1}
okE = True
Rf = mvec(Rvec, f)
RfPos = [{(i, j): c for (i, j), c in Rf[a].items() if i >= 0} for a in n5]
for l in range(0, L):
    C = [{} for _ in n5]
    for m in range(0, 9):
        c = cov(m)
        ce = [sum_p([pmul(c[a], ETA[a][b]) for a in n5]) for b in n5]
        if all(not x for x in ce): continue
        Em = Ecoef(m, l)
        for b in n5:
            C[b] = padd(C[b], sum_p([pmul(ce[a], Em[a][b]) for a in n5]))
    # RHS coefficient of psi'^l
    Rhs = mvec(Rvi, RfPos)
    rl = [{(0, j): -c for (i, j), c in Rhs[b].items() if i == l} for b in n5]
    for b in n5:
        lhs = {k: c for k, c in C[b].items() if k[1] < K - 1}
        rhs = {k: c for k, c in rl[b].items() if k[1] < K - 1}
        if lhs != rhs:
            okE = False; print("Lemma E mismatch l=", l, "b=", b, lhs, rhs)
print("Lemma E formula ok (psi'-degree < %d, mod u^%d):" % (L, K - 1), okE)
# Since this P generally violates (A), the contraction should be nonzero iff Y'_{>=1} != 0
Ypos = any(i >= 1 for a in n5 for (i, j) in Yp[a])
print("Y' has positive z-powers (A fails):", Ypos)
