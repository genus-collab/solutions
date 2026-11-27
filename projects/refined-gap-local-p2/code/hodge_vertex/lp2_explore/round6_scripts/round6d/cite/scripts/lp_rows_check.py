# Evidence only. Compare BS's I|_{P_j} at the Lho-Pandharipande weights alpha_i = -zeta^i (lambda_0 = 1,
# e1 = e2 = 0) with LP eq. (27): I|_{H=1} = exp(mu/z) (1 + R1 z + R2 z^2 + ...), 1 + D mu = L,
# R1 = (1-L^2)/18, R2 = (1-24L-2L^2+25L^4)/648, L = (1+27 y)^{-1/3}.
# Exact Fractions: series in y through y^N with coefficients truncated Laurent in z.
from fractions import Fraction as F
N = 8
ZMAX = 4   # keep z^k up to k <= ZMAX

def lmul(a, b):
    r = {}
    for i, x in a.items():
        for j, y in b.items():
            if i+j <= ZMAX:
                r[i+j] = r.get(i+j, F(0))+x*y
    return {k: v for k, v in r.items() if v != 0}

def linv_poly_lead(p):
    """inverse of Laurent series p, with lowest term at z^m, m=min(p), truncated"""
    m = min(p)
    c = p[m]
    # q = p/(c z^m) = 1 + ...
    q = {k-m: v/c for k, v in p.items()}
    prec = ZMAX+m+N*3+20
    inv = {0: F(1)}
    cur = {0: F(1)}
    r = {k: -v for k, v in q.items() if k != 0}
    res = [F(0)]*(prec+1)
    res[0] = F(1)
    for n in range(1, prec+1):
        s = F(0)
        for k, v in r.items():
            if k <= n:
                s += v*res[n-k]
        res[n] = s
    return {n-m: res[n]/c for n in range(prec+1) if res[n] != 0 and n-m <= ZMAX}

# y-coefficients of I/z at H = 1, alpha_i = -zeta^i: denominator prod_{m=1}^d ((1+m z)^3 - 1), numerator prod_{m=0}^{3d-1}(-3-m z)
coef = [{0: F(1)}]
ser = {0: F(1)}
for d in range(1, N+1):
    num = {0: F(1)}
    # wide internal ZMAX through a local large truncation
    ZSAVE = ZMAX
    ZMAX = ZSAVE+3*d+6
    for m in range(0, 3*d):
        num = lmul(num, {0: F(-3), 1: F(-m)} if m else {0: F(-3)})
    den = {0: F(1)}
    for m in range(1, d+1):
        den = lmul(den, {1: F(3*m), 2: F(3*m*m), 3: F(m*m*m)})
    val = lmul(num, linv_poly_lead(den))
    ZMAX = ZSAVE
    coef.append({k: v for k, v in val.items() if k <= ZMAX+3*N})

# log of series in y with Laurent coefficients: l = log(1 + sum c_d y^d), d*l_d = d*c_d - sum_{k<d} k l_k c_{d-k}
ZMAX = 2
def lmul2(a, b):
    r = {}
    for i, x in a.items():
        for j, y in b.items():
            r[i+j] = r.get(i+j, F(0))+x*y
    return r
l = [dict() for _ in range(N+1)]
for d in range(1, N+1):
    tot = {k: d*v for k, v in coef[d].items()}
    for k in range(1, d):
        p = lmul2(l[k], coef[d-k])
        for e, v in p.items():
            tot[e] = tot.get(e, F(0))-k*v
    l[d] = {e: v/d for e, v in tot.items() if v != 0}

def lser(d, e):  # coefficient y^d z^e of log(I/z) , trusted for e <= 2
    return l[d].get(e, F(0))

# L = (1+27y)^{-1/3} and powers
def powseries(alpha):  # (1+27 y)^alpha through y^N
    c = [F(1)]
    for d in range(1, N+1):
        c.append(c[-1]*(alpha-(d-1))/d*27)
    return c
def L(p): return powseries(F(-p, 3))

L1, L2, L4 = L(1), L(2), L(4)
# mu_d from D mu = L-1: mu_d = L1[d]/d
ok = True
print("d  z^-1 log I  vs mu_d | z^0 | z^1 vs R1 | z^2 vs R2 - R1^2/2")
R1 = [(F(1) if d == 0 else F(0))/18-L2[d]/18 for d in range(N+1)]
R2 = [((F(1)if d == 0 else F(0))-24*L1[d]-2*L2[d]+25*L4[d])/648 for d in range(N+1)]
R1sq = [sum(R1[i]*R1[d-i] for i in range(d+1)) for d in range(N+1)]
for d in range(1, N+1):
    a = lser(d, -1) == L1[d]/d
    b = lser(d, 0) == 0
    c = lser(d, 1) == R1[d]
    e = lser(d, 2) == R2[d]-R1sq[d]/2
    ok = ok and a and b and c and e
    low = min(l[d]) if l[d] else None
    print(d, a, b, c, e, "lowest z power in log coefficient:", low)
print("all equal through y^%d:" % N, ok)
# Negative control: a wrong R2 must fail
print("control R2 with 25 -> 24 fails:", all(lser(d, 2) != ((F(1)if d == 0 else F(0))-24*L1[d]-2*L2[d]+24*L4[d])/648-R1sq[d]/2 for d in range(1, 4)))
