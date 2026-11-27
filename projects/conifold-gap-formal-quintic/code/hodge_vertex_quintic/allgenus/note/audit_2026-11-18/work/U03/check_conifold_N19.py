# Independent sanity checks for Appendix B (cyclic form of PF4) -- exact arithmetic.
# Laurent series in u with coefficients polynomial in ell (formal log u), truncated.
from fractions import Fraction as F
import itertools

N = 19          # keep exponents <= N
MINE = -12      # lowest exponent tracked

class S:
    """element of Q((u))[ell]: dict {(n, k): coeff} meaning coeff * u^n * ell^k, truncated at u^N"""
    __slots__ = ("d",)
    def __init__(self, d=None):
        self.d = {}
        if d:
            for key, c in d.items():
                if c != 0 and key[0] <= N:
                    self.d[key] = F(c)
    @staticmethod
    def const(c): return S({(0, 0): F(c)})
    @staticmethod
    def u(n=1): return S({(n, 0): F(1)})
    @staticmethod
    def ell(): return S({(0, 1): F(1)})
    def __add__(a, b):
        if not isinstance(b, S): b = S.const(b)
        d = dict(a.d)
        for k, c in b.d.items():
            d[k] = d.get(k, 0) + c
        return S(d)
    __radd__ = __add__
    def __neg__(a): return S({k: -c for k, c in a.d.items()})
    def __sub__(a, b):
        if not isinstance(b, S): b = S.const(b)
        return a + (-b)
    def __rsub__(a, b): return S.const(b) - a
    def __mul__(a, b):
        if not isinstance(b, S): b = S.const(b)
        d = {}
        for (n1, k1), c1 in a.d.items():
            for (n2, k2), c2 in b.d.items():
                n = n1 + n2
                if n > N: continue
                key = (n, k1 + k2)
                d[key] = d.get(key, 0) + c1 * c2
        return S(d)
    __rmul__ = __mul__
    def D(a):
        # D = -(1-u) d/du ; D ell = 1 - 1/u ; i.e. D = (1 - 1/u) theta, theta u^n ell^k = n u^n ell^k + k u^n ell^(k-1)
        d = {}
        for (n, k), c in a.d.items():
            for (nn, cc) in ((n, 1), (n - 1, -1)):
                if n != 0:
                    key = (nn, k); d[key] = d.get(key, 0) + cc * c * n
                if k != 0:
                    key = (nn, k - 1); d[key] = d.get(key, 0) + cc * c * k
        return S(d)
    def Dk(a, k):
        for _ in range(k): a = a.D()
        return a
    def val(a):
        return min(n for (n, k) in a.d) if a.d else None
    def inv(a):
        # inverse for an ell-free series
        assert all(k == 0 for (_, k) in a.d)
        v = a.val()
        c0 = a.d[(v, 0)]
        # a = c0 u^v (1 + r)
        r = S({(n - v, 0): c / c0 for (n, _), c in a.d.items() if n != v})
        # geometric series, truncated at relative order N - (-v) ... keep generous
        res = S.const(1); term = S.const(1)
        for _ in range(N - MINE + 2):
            term = term * (-r)
            if not term.d: break
            res = res + term
        return res * S({(-v, 0): 1 / c0})
    def trunc(a, M):
        return S({k: c for k, c in a.d.items() if k[0] <= M})
    def coef(a, n, k=0): return a.d.get((n, k), F(0))
    def __repr__(a):
        return " + ".join(f"({c})u^{n}" + (f" l^{k}" if k else "") for (n, k), c in sorted(a.d.items())) or "0"

u = S.u(); one = S.const(1); ell = S.ell()
A = [ -F(24,625)*(1-u), -F(2,5)*(1-u), -F(7,5)*(1-u), -2*(1-u), u ]
def PF4(f): return sum((A[k]*f.Dk(k) for k in range(5)), S())

# 0. u-form of PF4 against the q-form D^4 - 5 q prod(5D+k), q=(1-u)/5^5
def PF4q(f):
    g = f
    for k in range(1, 5):
        g = 5*g.D() + k*g
    q = (1-u)*F(1, 5**5)
    return f.Dk(4) - 5*q*g
for n in range(0, 5):
    assert (PF4(S.u(n)) - PF4q(S.u(n))).d == {}, n
print("u-form of PF4 agrees with q-form on u^0..u^4")

# 1. indicial polynomials: PF4 u^n = sum_j P_j(n) u^{n+j}
def Pm3(x): return x*(x-1)**2*(x-2)
def Pm2(x): return -x*(x-1)*(4*x*x-6*x+F(17,5))
for n in range(0, 9):
    r = PF4(S.u(n))
    assert r.coef(n-3) == Pm3(F(n)) and r.coef(n-2) == Pm2(F(n)), n
    assert r.val() is None or r.val() >= n-3
print("P_{-3}, P_{-2} as stated")
r = PF4(ell)
print("PF4(ell) =", r, " -> coefficient of u^-1 =", r.coef(-1), "(claimed P_{-1}'(0) = -7/5)")

# 2. power series solutions f_k = u^k + O(u^3) by solving coefficient equations
def solve_ps(b0, b1, b2, rhs=None):
    # find g = sum b_n u^n with PF4 g = rhs (power series), given b0,b1,b2
    b = [F(b0), F(b1), F(b2)]
    for M in range(0, N - 3):
        # coefficient of u^M in PF4 g involves b_{M+3} with factor Pm3(M+3)
        g = S({(n, 0): c for n, c in enumerate(b)})
        r = PF4(g)
        target = rhs.coef(M) if rhs is not None else F(0)
        b.append((target - r.coef(M)) / Pm3(F(M+3)))
    return S({(n, 0): c for n, c in enumerate(b)})
f0 = solve_ps(1, 0, 0); f1 = solve_ps(0, 1, 0); f2 = solve_ps(0, 0, 1)
for f in (f0, f1, f2):
    res = PF4(f)
    assert all(n > N - 5 for (n, k) in res.d), res
print("f0 =", f0.trunc(5)); print("f1 =", f1.trunc(5)); print("f2 =", f2.trunc(5))
gam = F(7, 10)
w1 = f1 + gam*f2
print("varpi1 =", w1.trunc(5))
# obstruction: PF4(g*ell) = ell*PF4(g) + r(g); check negative-power coefficients
for name, g in (("f0", f0), ("f1", f1), ("f2", f2), ("varpi1", w1)):
    r = PF4(g*ell)
    print("  obstruction of", name, ":", [r.coef(m) for m in (-3, -2, -1)], " ell-part low orders:", [r.coef(m, 1) for m in range(-3, 4)])
# rho_3: PF4(w1*ell + rho3) = 0 with rho3 in u^3 Q[[u]]
rr = PF4(w1*ell)
rhs = -S({(n, 0): c for (n, k), c in rr.d.items() if k == 0})
rho3 = solve_ps(0, 0, 0, rhs=rhs)
y3 = w1*ell + rho3
res = PF4(y3)
assert all(n > N - 5 for (n, k) in res.d), res
print("rho3 =", rho3.trunc(6))

# 3. Jet matrix from the definition of the concomitant
def conc(f, g):
    tot = S()
    for k in range(1, 5):
        for j in range(0, k):
            tot = tot + (-1)**j * (A[k]*g).Dk(j) * f.Dk(k-1-j)
    return tot
Jet = [[S() for _ in range(4)] for _ in range(4)]
from math import comb
for a in range(4):
    for b in range(4):
        e = S()
        for k in range(a+b+1, 5):
            e = e + (-1)**(k-1-a) * comb(k-1-a, b) * A[k].Dk(k-1-a-b)
        Jet[a][b] = e
claimed = [[S(), F(2,5)*(1-u), 1-u, -u], [-F(2,5)*(1-u), S(), u, S()], [-(1-u), -u, S(), S()], [u, S(), S(), S()]]
for a in range(4):
    for b in range(4):
        assert (Jet[a][b] - claimed[a][b]).d == {}, (a, b, Jet[a][b])
print("Jet matrix as stated")
def conc_jet(f, g):
    return sum((f.Dk(a)*Jet[a][b]*g.Dk(b) for a in range(4) for b in range(4)), S())
basis = [f0, f1, f2, y3]; names = ["f0", "f1", "f2", "y3"]
print("Gram matrix (entries shown up to u^6; must be constants):")
for i in range(4):
    row = []
    for j in range(4):
        c1 = conc(basis[i], basis[j]).trunc(6); c2 = conc_jet(basis[i], basis[j]).trunc(6)
        assert (c1 - c2).d == {}
        row.append(str(c1))
    print("  ", row)

# 4. Wronskian
def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = S()
    for j in range(n):
        if not M[0][j].d: continue
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        tot = tot + (-1)**j * M[0][j] * det(minor)
    return tot
Wr = [[basis[i].Dk(a) for i in range(4)] for a in range(4)]
dW = det(Wr)
print("u^2 det Wr_f =", (S.u(2)*dW).trunc(5))

# 5. Mcyc(y0,y1) as an operator: compare with PF4 on test functions
def Mcyc_apply(y0, y1, f):
    t = y1 * y0.inv(); td = t.D(); tdi = td.inv(); y0i = y0.inv()
    g = y0i * f
    g = g.D(); g = tdi * g
    g = g.D(); g = (u * y0 * y0 * td * td) * g
    g = g.D(); g = tdi * g
    g = g.D(); g = y0i * g
    return g
def compare(y0, y1, label):
    worst = None
    for test in [S.u(n) for n in range(-2, 7)] + [ell, u*ell, S.u(3)*ell*ell]:
        d = Mcyc_apply(y0, y1, test) - PF4(test)
        v = d.val()
        if v is not None and (worst is None or v < worst): worst = v
    print("  Mcyc - PF4 on test functions,", label, ": lowest surviving order =", worst, "(truncation artefacts expected near u^%d)" % (N-6))
compare(f0, w1, "(f0, varpi1)")
compare(f0 + F(3,7)*w1 - F(5,3)*f2, w1, "frame alpha=3/7, beta=-5/3")
compare(f0 + 2*f2, 5*w1 + 3*(f0 + 2*f2), "(y0, 5 varpi1 + 3 y0)")
compare(f0, f1, "NON-isotropic plane (f0, f1): should fail")
compare(f2 , w1, "(f2, varpi1) (y0(0)=0, allowed by the theorem)")

# 6. tc in the frame (0,0)
tc = w1 * f0.inv()
print("tc (frame 0,0) =", tc.trunc(4))
