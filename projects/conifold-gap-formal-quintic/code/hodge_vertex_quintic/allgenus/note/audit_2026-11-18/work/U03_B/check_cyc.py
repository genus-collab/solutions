"""Independent low-order checks for appendix B (cyclic form of PF4) in the variable u.
D = -(1-u) d/du.  PF4 = u D^4 - (1-u)(2 D^3 + 7/5 D^2 + 2/5 D + 24/625).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from ser import S

N = 26  # absolute precision in u

def U(N=N):
    return S(1, [1] + [0] * (N - 2))
def ONE(N=N):
    return S(0, [1] + [0] * (N - 1))

def D(f):
    one_minus_u = S(0, [1, -1] + [0] * (len(f.c)))
    return -(one_minus_u * f.d())

def Dk(f, k):
    for _ in range(k):
        f = D(f)
    return f

def PF4(f):
    u = S(1, [1] + [0] * len(f.c))
    omu = S(0, [1, -1] + [0] * len(f.c))
    return u * Dk(f, 4) - omu * (2 * Dk(f, 3) + Fr(7, 5) * Dk(f, 2) + Fr(2, 5) * D(f) + Fr(24, 625) * f)

# --- power series solutions f_k = u^k + O(u^3), solved coefficient by coefficient
def solve_series(b0, b1, b2, N=N):
    c = [Fr(b0), Fr(b1), Fr(b2)]
    # unknown b_n for n>=3 determined from coefficient of u^(n-3) of PF4 f = 0
    for n in range(3, N):
        # trial with b_n = 0 and b_n = 1 : linear
        f0 = S(0, c + [Fr(0)] + [Fr(0)] * 4)
        f1 = S(0, c + [Fr(1)] + [Fr(0)] * 4)
        r0 = PF4(f0).coef(n - 3)
        r1 = PF4(f1).coef(n - 3)
        assert r1 != r0
        c.append(-r0 / (r1 - r0))
    return S(0, c)

f0 = solve_series(1, 0, 0)
f1 = solve_series(0, 1, 0)
f2 = solve_series(0, 0, 1)
for name, f in (("f0", f0), ("f1", f1), ("f2", f2)):
    r = PF4(f)
    # PF4 lowers precision by 4 derivatives; check the available part vanishes
    print(name, f.c[:6], "PF4 residual zero up to prec", r.prec(), ":", r.iszero())

gam = Fr(7, 10)
w1 = f1 + gam * f2
print("varpi_1 =", w1.c[:5])

# --- indicial data P_{-3}, P_{-2}, P_{-1}'(0) by action on u^n and on u^n*log
# PF4 u^n = sum_j P_j(n) u^{n+j}
def PF4_mono(n):
    f = S(n, [1] + [0] * 10)
    return PF4(f)
for n in range(0, 6):
    r = PF4_mono(n)
    Pm3 = r.coef(n - 3); Pm2 = r.coef(n - 2)
    assert Pm3 == n * (n - 1) ** 2 * (n - 2), (n, Pm3)
    assert Pm2 == -n * (n - 1) * (4 * n * n - 6 * n + Fr(17, 5)), (n, Pm2)
print("P_{-3}, P_{-2} as stated: OK")

# --- logarithmic solution y3 = w1*ell + rho3 : represent a + b*ell with D ell = 1 - 1/u
class LG:
    # a + b*ell
    def __init__(self, a, b):
        self.a = a; self.b = b
def Dlog(x):
    dl = S(-1, [-1, 1] + [0] * (len(x.b.c)))  # 1 - 1/u = -u^-1 + 1
    return LG(D(x.a) + x.b * dl, D(x.b))
def PF4log(x):
    d = [x]
    for _ in range(4):
        d.append(Dlog(d[-1]))
    def comb(sel):
        u = S(1, [1] + [0] * 40)
        omu = S(0, [1, -1] + [0] * 40)
        return u * sel(d[4]) - omu * (2 * sel(d[3]) + Fr(7, 5) * sel(d[2]) + Fr(2, 5) * sel(d[1]) + Fr(24, 625) * sel(d[0]))
    return comb(lambda t: t.a), comb(lambda t: t.b)

# obstruction r(g) = PF4(g*ell) - ell*PF4(g), coefficients of u^-3,u^-2,u^-1
def obstruction(g):
    ra, rb = PF4log(LG(S(0, [0] * len(g.c)), g))
    return [ra.coef(-3), ra.coef(-2), ra.coef(-1)]
print("obstruction columns f0,f1,f2:", obstruction(f0), obstruction(f1), obstruction(f2))
print("obstruction of varpi_1:", obstruction(w1))

# solve rho3 in u^3 Q[[u]] with PF4(w1*ell + rho3) = 0
ra, rb = PF4log(LG(S(0, [0] * len(w1.c)), w1))
# ra = r(w1) power series; need PF4 rho3 = -ra
c = [Fr(0), Fr(0), Fr(0)]
M = N - 6
for n in range(3, M):
    t0 = S(0, c + [Fr(0)] + [Fr(0)] * 4)
    t1 = S(0, c + [Fr(1)] + [Fr(0)] * 4)
    r0 = PF4(t0).coef(n - 3) + ra.coef(n - 3)
    r1 = PF4(t1).coef(n - 3) + ra.coef(n - 3)
    c.append(-r0 / (r1 - r0))
rho3 = S(0, c)
y3 = LG(rho3, S(0, w1.c[:len(rho3.c)]))
ra2, rb2 = PF4log(y3)
print("PF4(y3): ell-part zero:", rb2.iszero(), " series part zero:", ra2.iszero(), "prec", ra2.prec())
print("rho3 =", rho3.c[:6])

# --- Jet matrix
def jet(M=12):
    u = S(1, [1] + [0] * M); omu = S(0, [1, -1] + [0] * M)
    z = S(0, [0] * (M + 1))
    return [[z, Fr(2, 5) * omu, omu, -u], [-Fr(2, 5) * omu, z, u, z], [-omu, -u, z, z], [u, z, z, z]]
J = jet(N)

def conc_series(f, g):
    tot = None
    for a in range(4):
        for b in range(4):
            t = Dk(f, a) * J[a][b] * Dk(g, b)
            tot = t if tot is None else tot + t
    return tot

def conc_def(f, g):
    # definition: sum_{k=1}^4 sum_{j=0}^{k-1} (-1)^j D^j(A_k g) D^{k-1-j} f
    u = S(1, [1] + [0] * N); omu = S(0, [1, -1] + [0] * N)
    A = {4: u, 3: -2 * omu, 2: -Fr(7, 5) * omu, 1: -Fr(2, 5) * omu}
    tot = None
    for k in range(1, 5):
        for j in range(k):
            t = ((-1) ** j) * Dk(A[k] * g, j) * Dk(f, k - 1 - j)
            tot = t if tot is None else tot + t
    return tot

for (n1, a), (n2, b) in [(("f0", f0), ("f1", f1)), (("f0", f0), ("f2", f2)), (("f1", f1), ("f2", f2)), (("f0", f0), ("w1", w1)), (("f1", f1), ("f0", f0))]:
    c1 = conc_series(a, b); c2 = conc_def(a, b)
    d = c1 - c2
    print("Conc(%s,%s) = %s (const term), nonconstant part zero: %s ; jet form == definition: %s" % (n1, n2, c1.coef(0), all(x == 0 for x in c1.c[1:]) and c1.v == 0 or c1.norm().v == 0 and all(x == 0 for x in c1.norm().c[1:]) or c1.iszero(), d.iszero()))

# Conc(f, y3) for f power series: ell-part and series part
def conc_log(f, y):
    d = [y]
    for _ in range(3):
        d.append(Dlog(d[-1]))
    ta = None; tb = None
    for a in range(4):
        for b in range(4):
            fa = Dk(f, a) * J[a][b]
            xa = fa * d[b].a; xb = fa * d[b].b
            ta = xa if ta is None else ta + xa
            tb = xb if tb is None else tb + xb
    return ta, tb
for n1, a in (("f0", f0), ("f1", f1), ("f2", f2)):
    ta, tb = conc_log(a, y3)
    tn = ta.norm()
    print("Conc(%s,y3): ell-coefficient zero: %s ; series part = %s" % (n1, tb.iszero(), ta))

# --- Wronskian of standard basis: u^2 det Wr_f
# det of 4x4 with last column a+b*ell : determinant is linear in last column
def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
def det4(m):
    tot = None
    for j in range(4):
        minor = [[m[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        t = ((-1) ** j) * m[0][j] * det3(minor)
        tot = t if tot is None else tot + t
    return tot
cols = [f0, f1, f2]
dy = [y3]
for _ in range(3):
    dy.append(Dlog(dy[-1]))
Wa = [[Dk(cols[i], a) for i in range(3)] + [dy[a].a] for a in range(4)]
Wb = [[Dk(cols[i], a) for i in range(3)] + [dy[a].b] for a in range(4)]
da = det4(Wa); db = det4(Wb)
u2 = S(2, [1] + [0] * 30)
print("det Wr_f: ell part zero:", db.iszero(), "; u^2 det Wr_f =", (u2 * da))

# --- the cyclic operator
def Mcyc(y0, y1, f):
    t = y1 / y0
    td = D(t)
    u = S(1, [1] + [0] * 60)
    g = f / y0
    g = D(g)
    g = g / td
    g = D(g)
    g = g * (u * y0 * y0 * td * td)
    g = D(g)
    g = g / td
    g = D(g)
    return g / y0

def test_cyc(y0, y1, label):
    ok = True
    worst = None
    for n in range(-2, 7):
        f = S(n, [1] + [0] * (N - 1))
        r = Mcyc(y0, y1, f) - PF4(f)
        if not r.iszero():
            ok = False
        p = r.prec() - n
        worst = p if worst is None else min(worst, p)
    print("PF4 == Mcyc(%s): %s (relative precision checked: u^%d)" % (label, ok, worst))

for al, be in ((0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(-5, 2), Fr(11, 3))):
    w0 = f0 + al * w1 + be * f2
    test_cyc(w0, w1, "varpi_0[%s,%s], varpi_1" % (al, be))
test_cyc(f0 + f1, 3 * w1, "f0+f1, 3 varpi_1")
test_cyc(f2, w1, "f2, varpi_1 (y0(0)=0)")
test_cyc(w1, f0, "varpi_1, f0 (y0 = varpi_1)")
# non-isotropic planes must fail
test_cyc(f0, f1, "f0, f1 (NOT isotropic: expected False)")
test_cyc(f0, f2, "f0, f2 (NOT isotropic: expected False)")

# frame data
w0 = f0
tc = w1 / w0
print("t_c (frame 0,0) =", tc.c[:5] if tc.v == 0 else (tc.v, tc.c[:4]))
