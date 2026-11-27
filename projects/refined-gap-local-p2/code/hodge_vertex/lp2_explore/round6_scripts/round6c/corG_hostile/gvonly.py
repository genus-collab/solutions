# Independent of rows and graph sums: uses only GV numbers (test_lr tables), the mirror map and hypergeometric series.
# (1) genus 0: I_1^[2] = -3 F_0'(t); Yukawa = -(1/3) X/A^3.
# (2) genus 1: F_1 = -(1/2) log A - (1/12) log(Delta_c y^7) as Q-series.
# (3) genus 2: fit an ansatz in F, X from N_{2,d}, d<=10; substitute F_cf; read the polar part.
# (4) eq:F2ref substitution directly.
# (5) t_cf numerically, both branches.
from fractions import Fraction as Fr
from math import factorial, comb, pi, sqrt, log
import cmath
GV0 = [3, -6, 27, -192, 1695, -17064, 188454, -2228160, 27748899, -360012150]
GV1 = [0, 0, -10, 231, -4452, 80948, -1438086, 25301295, -443384578, 7760515332]
GV2 = [0, 0, 0, -102, 5430, -194022, 5784837, -155322234, 3894455457, -93050366010]
D = 10
N = D + 1
def mul(a, b, n=N):
    r = [Fr(0)] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]):
                r[i + j] += x * y
    return r
def inv(a, n=N):
    r = [Fr(0)] * n
    r[0] = 1 / Fr(a[0])
    for k in range(1, n):
        r[k] = -sum(a[i] * r[k - i] for i in range(1, k + 1)) * r[0]
    return r
def exps(a, n=N):
    r = [Fr(0)] * n
    r[0] = Fr(1)
    for m in range(1, n):
        r[m] = sum(k * a[k] * r[m - k] for k in range(1, m + 1)) / m
    return r
def logs(a, n=N):
    # a[0]=1
    r = [Fr(0)] * n
    for m in range(1, n):
        r[m] = (m * a[m] - sum(k * r[k] * a[m - k] for k in range(1, m))) / m
    return r
def H(n):
    return sum(Fr(1, k) for k in range(1, n + 1))
# y-series
A = [Fr(factorial(3 * d), factorial(d) ** 3) * (-1) ** d for d in range(N)]
g1 = [Fr(0)] + [3 * Fr(factorial(3 * d - 1), factorial(d) ** 3) * (-1) ** d for d in range(1, N)]
g2 = [Fr(0)] + [9 * Fr(factorial(3 * d - 1), factorial(d) ** 3) * (H(3 * d - 1) - H(d)) * (-1) ** d for d in range(1, N)]
Qy = exps(g1)
def powQ(d):
    # Q^d/y^d
    r = [Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(d):
        r = mul(r, Qy)
    return r
def q_to_y(coeffs):
    """coeffs[d] of Q^d, d>=0 -> series in y"""
    r = [Fr(0)] * N
    for d, c in enumerate(coeffs):
        if d < N and c:
            p = powQ(d)
            for i in range(N - d):
                r[d + i] += c * p[i]
    return r
def Ngd(table, h, d, g):
    # contribution of n^h to N_{g,d} is computed in N_all below
    pass
def sinc_pow(p, n):
    s = [Fr(0)] * n
    for m in range(0, (n + 1) // 2):
        s[2 * m] = Fr((-1) ** m, factorial(2 * m + 1) * 4 ** m)
    if p == 0:
        return [Fr(1)] + [Fr(0)] * (n - 1)
    base = s if p > 0 else inv(s, n)
    r = [Fr(1)] + [Fr(0)] * (n - 1)
    for _ in range(abs(p)):
        r = mul(r, base, n)
    return r
TAB = {0: GV0, 1: GV1, 2: GV2}
def N_all(g, d):
    t = Fr(0)
    for k in range(1, d + 1):
        if d % k == 0:
            dd = d // k
            for h in range(0, g + 1):
                n = TAB[h][dd - 1]
                if n:
                    t += n * Fr(1, k) * Fr(k) ** (2 * g - 2) * sinc_pow(2 * h - 2, 2 * g - 2 * h + 1)[2 * g - 2 * h]
    return t
# genus 0: N_{0,d}=Sum_{k|d} n_{d/k}/k^3
N0 = [Fr(0)] + [sum(Fr(GV0[d // k - 1], k ** 3) for k in range(1, d + 1) if d % k == 0) for d in range(1, N)]
print('N_0,d', [str(x) for x in N0[1:5]])
# (1) g2 - g^2/2 = -3 Sum d N Q^d
lhs = [a - b / 2 for a, b in zip(g2, mul(g1, g1))]
rhs = q_to_y([-3 * d * N0[d] for d in range(N)])
print('(1a) I_1^[2] = -3 F_0prime through y^%d:' % D, lhs == rhs)
X = inv([Fr(1), Fr(27)] + [Fr(0)] * (N - 2))
Ai = inv(A)
yuk = [-x / 3 for x in mul(X, mul(Ai, mul(Ai, Ai)))]
rhs2 = q_to_y([Fr(-1, 3)] + [d ** 3 * N0[d] for d in range(1, N)])
print('(1b) Yukawa -(1/3) X/A^3 = -1/3 + Sum d^3 N_0d Q^d:', yuk == rhs2)
# (2) genus 1: F_1 = Sum N_{1,d} Q^d, N_{1,d}=Sum_{k|d} (n^0_{d/k}/12 + n^1_{d/k})/k
N1 = [Fr(0)] + [sum((Fr(GV0[d // k - 1], 12) + GV1[d // k - 1]) / k for k in range(1, d + 1) if d % k == 0) for d in range(1, N)]
# -1/2 log A -1/12 log(Delta_c) -7/12 log y , with log y= t - g1 : F_1 + t*7/12 ...
la = logs(A)
ld = logs([Fr(1), Fr(27)] + [Fr(0)] * (N - 2))
ser = [-x / 2 - w / 12 + Fr(1, 12) * h for x, w, h in zip(la, ld, g1)]
print('(2) F_1 + t/12 = -(1/2)log A -(1/12)log(Delta_c) +(1/12) g, y-series equals Sum N_1d Q^d:', ser == q_to_y(N1))
# ---- conifold series
M = 16
W = []
x = Fr(1)
for n in range(M):
    W.append(x)
    x = x * (Fr(1, 3) + n) * (Fr(2, 3) + n) / (n + 1) ** 2
def mulM(a, b):
    return mul(a, b, M)
dW = [(n + 1) * W[n + 1] for n in range(M - 1)] + [Fr(0)]
thW = [(dW[n - 1] if n >= 1 else 0) - dW[n] for n in range(M)]  # (xi-1) W'
Fcf = mulM(thW, inv(W, M))
pc = []
s = Fr(0)
for n in range(M):
    s += W[n]
    pc.append(s)
tc = [Fr(0)] + [pc[n] / (n + 1) for n in range(M - 1)]
print('t_c =', [str(c) for c in tc[:5]], ' F_cf(0)=', Fcf[0])
# (3) ansatz: F^k X^j, -1 <= j <= 2-k, k<=3
terms = [(k, j) for k in range(4) for j in range(-1, 3 - k)]
Fy = mul([d * A[d] for d in range(N)], Ai)
def xpow(j):
    if j >= 0:
        r = [Fr(1)] + [Fr(0)] * (N - 1)
        for _ in range(j):
            r = mul(r, X)
        return r
    r = [Fr(1)] + [Fr(0)] * (N - 1)
    base = [Fr(1), Fr(27)] + [Fr(0)] * (N - 2)
    for _ in range(-j):
        r = mul(r, base)
    return r
def fpow(k):
    r = [Fr(1)] + [Fr(0)] * (N - 1)
    for _ in range(k):
        r = mul(r, Fy)
    return r
basis = [mul(fpow(k), xpow(j)) for (k, j) in terms]
target = q_to_y([Fr(0)] + [N_all(2, d) for d in range(1, N)])
# solve for y^1..y^10 ; constant term unknown, so add unknown const "1" is (k=0,j=0), set by gauge
rows = [[basis[i][d] for i in range(len(terms))] + [target[d]] for d in range(1, N)]
nv = len(terms)
piv = []
r = 0
for c in range(nv):
    p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
    if p is None:
        continue
    rows[r], rows[p] = rows[p], rows[r]
    f = rows[r][c]
    rows[r] = [x / f for x in rows[r]]
    for i in range(len(rows)):
        if i != r and rows[i][c] != 0:
            f = rows[i][c]
            rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
    piv.append(c)
    r += 1
print('(3) ansatz terms', len(terms), 'rank', r, 'of', len(rows), 'equations; inconsistent rows:', sum(1 for i in range(r, len(rows)) if rows[i][-1] != 0))
free = [terms[c] for c in range(nv) if c not in piv]
print('    free directions:', free)
sol = {terms[c]: rows[i][-1] for i, c in enumerate(piv)}
print('    fitted coefficients (free set to 0):', {k: str(v) for k, v in sol.items()})
# Laurent in xi: X=1/xi
def cf_eval(coefs):
    # returns dict power->coef for xi^-3.. up to 5
    out = {}
    for (k, j), c in coefs.items():
        p = [Fr(1)] + [Fr(0)] * (M - 1)
        for _ in range(k):
            p = mulM(p, Fcf)
        for i, x in enumerate(p):
            n = i - j
            out[n] = out.get(n, Fr(0)) + c * x
    return out
L = cf_eval(sol)
print('    F_{0,2}|F->F_cf Laurent in xi:', {n: str(L[n]) for n in sorted(L) if n <= 1})
# polar in t_c: -1/80 /t_c^2 ; 1/t_c^2 = xi^-2 (1 - 11/9 xi + ...)
tcr = tc[1:]
itc = inv(mulM(tcr + [Fr(0)], tcr + [Fr(0)]), M - 1)
print('    1/t_c^2 =xi^-2*', [str(x) for x in itc[:3]])
print('    coefficient xi^-2', L.get(-2), ' expected -1/80; xi^-1 minus (-1/80)(-11/9):', L.get(-1, 0) + Fr(1, 80) * itc[1])
# (4) eq:F2ref printed
F02 = {(3, -1): Fr(5, 8), (2, -1): Fr(5, 8), (2, 0): Fr(-1, 2), (1, 1): Fr(13, 96), (1, -1): Fr(20, 96), (1, 0): Fr(-32, 96),
       (0, 2): Fr(-27, 2160), (0, 1): Fr(98, 2160), (0, 0): Fr(-121, 2160), (0, -1): Fr(50, 2160)}
L2 = cf_eval(F02)
print('(4) eq:F2ref F_{0,2}: xi^-2', L2.get(-2), ' xi^-1 residual', L2.get(-1, 0) + Fr(1, 80) * itc[1])
same = all(sol.get(k, 0) == v or k == (0, 0) for k, v in F02.items())
print('    fit vs eq:F2ref coefficients agree except the constant direction:', same)
# (5) t_cf
def tcf_num(y, branch):
    l = log(-y)
    ly = complex(l, branch * pi)
    NN = 400
    # terms via ratios
    I0 = ly
    I12 = -ly ** 2 / 2 + ly * 0
    g = 0.0
    g2v = 0.0
    from math import lgamma
    for d in range(1, NN):
        c = 3 * (lgamma(3 * d) - 3 * lgamma(d + 1))
    g = 0.0
    term = 0.0
    a = 1.0  # (3d-1)!/(d!)^3 (-y)^d
    Hd = 0.0
    H3 = 0.0
    for d in range(1, NN):
        lg = lgamma(3 * d) - 3 * lgamma(d + 1) + d * log(-y)
        t = (2.718281828459045 ** lg) * ((-1) ** d) * ((-1) ** d)  # -y>0 and (-y)^d positive
        H1 = sum(1.0 / k for k in range(d + 1, 3 * d))
        g += 3 * t
        g2v += 9 * t * H1
    I0 = ly + g
    I2 = -ly ** 2 / 2 + ly * I0 + g2v
    return sqrt(3) / (2 * pi) * (I2 - 1j * pi * I0 - pi ** 2)
for y in (-1 / 54, -1 / 81):
    xi = 1 + 27 * y
    val = sum(float(c) * xi ** i for i, c in enumerate(tc))
    print('(5) y', y, 't_c(trunc 16)', val, 'BS expr log y=l+i pi', tcf_num(y, 1), ' l-i pi', tcf_num(y, -1))
