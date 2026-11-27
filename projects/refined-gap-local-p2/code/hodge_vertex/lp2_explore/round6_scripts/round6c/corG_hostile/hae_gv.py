# Comparison only (no rows, no graph sum, no gap input): BS's unrefined HAE fixes the F-dependence of F_{0,g};
# the X-only ambiguity is fitted to GV numbers N_{g,d}, d<=10, and then F -> F_cf gives the conifold Laurent expansion.
from fractions import Fraction as Fr
from math import factorial
from test_lr import N_gd
D = 10
N = D + 1
def padd(a, b, f=Fr(1)):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, Fr(0)) + f * v
    return {k: v for k, v in r.items() if v != 0}
def pmul(a, b):
    r = {}
    for (k1, j1), v1 in a.items():
        for (k2, j2), v2 in b.items():
            key = (k1 + k2, j1 + j2)
            r[key] = r.get(key, Fr(0)) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}
ThF = {(2, 0): Fr(-1), (1, 0): Fr(-1), (1, 1): Fr(1), (0, 0): Fr(-2, 9), (0, 1): Fr(2, 9)}
def theta(p):
    r = {}
    for (k, j), v in p.items():
        if k:
            r = padd(r, pmul({(k - 1, j): v * k}, ThF))
        if j:
            r = padd(r, {(k, j): -v * j, (k, j + 1): v * j})
    return r
th1 = {(1, 0): Fr(-1, 2), (0, 0): Fr(-1, 6), (0, 1): Fr(1, 12)}   # theta F_{0,1}
F2 = {(3, -1): Fr(5, 8), (2, -1): Fr(5, 8), (2, 0): Fr(-1, 2), (1, 1): Fr(13, 96), (1, -1): Fr(20, 96), (1, 0): Fr(-32, 96)}
# series
A = [Fr(factorial(3 * d), factorial(d) ** 3) * (-1) ** d for d in range(N)]
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
Fy = mul([d * A[d] for d in range(N)], inv(A))
X = inv([Fr(1), Fr(27)] + [Fr(0)] * (N - 2))
Xi = [Fr(1), Fr(27)] + [Fr(0)] * (N - 2)
def ev_y(p):
    r = [Fr(0)] * N
    for (k, j), v in p.items():
        s = [Fr(1)] + [Fr(0)] * (N - 1)
        for _ in range(k):
            s = mul(s, Fy)
        for _ in range(abs(j)):
            s = mul(s, X if j > 0 else Xi)
        r = [a + v * b for a, b in zip(r, s)]
    return r
g1 = [Fr(0)] + [3 * Fr(factorial(3 * d - 1), factorial(d) ** 3) * (-1) ** d for d in range(1, N)]
def exps(a, n=N):
    r = [Fr(0)] * n
    r[0] = Fr(1)
    for m in range(1, n):
        r[m] = sum(k * a[k] * r[m - k] for k in range(1, m + 1)) / m
    return r
Qy = exps(g1)
def q_to_y(c):
    r = [Fr(0)] * N
    pw = [Fr(1)] + [Fr(0)] * (N - 1)
    for d in range(N):
        for i in range(N - d):
            r[d + i] += c[d] * pw[i]
        pw = mul(pw, Qy)
    return r
# conifold series
M = 22
W = []
x = Fr(1)
for n in range(M):
    W.append(x)
    x = x * (Fr(1, 3) + n) * (Fr(2, 3) + n) / (n + 1) ** 2
dW = [(n + 1) * W[n + 1] for n in range(M - 1)] + [Fr(0)]
thW = [(dW[n - 1] if n >= 1 else Fr(0)) - dW[n] for n in range(M)]
Fcf = mul(thW, inv(W, M), M)
pc = []
s = Fr(0)
for n in range(M):
    s += W[n]
    pc.append(s)
tcr = [pc[n] / (n + 1) for n in range(M - 1)]  # t_c/xi
def cf_laurent(p, nmax=1):
    out = {}
    for (k, j), c in p.items():
        s = [Fr(1)] + [Fr(0)] * (M - 1)
        for _ in range(k):
            s = mul(s, Fcf, M)
        for i, v in enumerate(s):
            n = i - j
            out[n] = out.get(n, Fr(0)) + c * v
    return out
def tc_neg(power):
    r = [Fr(1)] + [Fr(0)] * (M - 2)
    it = inv(tcr, M - 1)
    for _ in range(power):
        r = mul(r, it, M - 1)
    return r
def fit(g, Fdep, label):
    target = q_to_y([Fr(0)] + [N_gd(g, d) for d in range(1, N)])
    base = ev_y(Fdep)
    js = list(range(-(g - 1), 2 * g - 1))
    bas = [ev_y({(0, j): Fr(1)}) for j in js]
    rows = [[bas[i][d] for i in range(len(js))] + [target[d] - base[d]] for d in range(0, N)]
    # y^0: gauge value N_{g,0}=0, the T_F value of F_{0,g}
    n = len(js)
    r = 0
    piv = []
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        f = rows[r][c]
        rows[r] = [x / f for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
    bad = sum(1 for i in range(r, len(rows)) if rows[i][-1] != 0)
    print(label, 'unknowns', n, 'rank', r, 'equations', len(rows), 'inconsistent', bad)
    amb = {(0, js[c]): rows[i][-1] for i, c in enumerate(piv)}
    print('  ambiguity:', {k[1]: str(v) for k, v in amb.items()})
    P = padd(Fdep, amb)
    return P
def show(g, P):
    L = cf_laurent(P)
    tn = tc_neg(2 * g - 2)
    b = Fr(3) ** (g - 1) * [Fr(1, 6), Fr(-1, 30), Fr(1, 42), Fr(-1, 30)][g - 1] / (2 * g * (2 * g - 2))
    print('  xi-Laurent:', {n: str(L.get(n, 0)) for n in range(-(2 * g - 2) - 1, 1)})
    print('  3^(g-1)B_2g/(2g(2g-2)) =', b)
    res = {n: L.get(n, Fr(0)) - b * tn[n + 2 * g - 2] for n in range(-(2 * g - 2), 0)}
    print('  polar part of F^cf_{0,%d} - b/t_c^%d :' % (g, 2 * g - 2), {k: str(v) for k, v in res.items()})
    print('  xi^0 coefficient of F^cf_{0,%d}:' % g, L.get(0))
def integrate_F(p):
    return {(k + 1, j): v / (k + 1) for (k, j), v in p.items()}
# genus 2 control
rhs2 = padd(padd(theta(th1), pmul({(1, 0): Fr(-1)}, th1)), pmul(th1, th1))
dF2 = pmul({(0, -1): Fr(3, 2)}, rhs2)
print('genus 2 HAE reproduces eq:F2ref F-dependence:', integrate_F(dF2) == {k: v for k, v in F2.items()})
P2 = fit(2, integrate_F(dF2), 'g=2')
show(2, P2)
# genus 3
thF2 = theta(P2)
rhs3 = padd(padd(theta(thF2), pmul({(1, 0): Fr(-1)}, thF2)), pmul(th1, thF2), Fr(2))
P3 = fit(3, integrate_F(pmul({(0, -1): Fr(3, 2)}, rhs3)), 'g=3')
show(3, P3)
thF3 = theta(P3)
rhs4 = padd(padd(padd(theta(thF3), pmul({(1, 0): Fr(-1)}, thF3)), pmul(th1, thF3), Fr(2)), pmul(thF2, thF2))
P4 = fit(4, integrate_F(pmul({(0, -1): Fr(3, 2)}, rhs4)), 'g=4')
show(4, P4)
