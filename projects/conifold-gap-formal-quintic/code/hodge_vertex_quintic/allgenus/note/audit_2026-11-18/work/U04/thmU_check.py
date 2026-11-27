# Independent sanity check (referee unit U04) of:
#  (1) prop:gau-equalprod (ii): v = u*tau^(-5/4) satisfies v'^5/v = 1/(u(1-u)^5);
#  (2) the construction in the proof of Theorem U (a): commutant projection (lem:gau-comm (e)),
#      lifting (lem:gau-lift), scalar equation (lem:gau-scal (a),(b)), for a random pair of
#      one-pole cyclic systems with equal products; integrality at each step and the gauge
#      equation are verified a posteriori.
# Exact arithmetic, truncated Laurent series with precision tracking.
from fractions import Fraction as Fr
import random

class S:
    """Laurent series in u known modulo u^prec."""
    def __init__(self, d, prec):
        self.prec = prec
        self.d = {k: Fr(v) for k, v in d.items() if v != 0 and k < prec}
    def val(self):
        return min(self.d) if self.d else self.prec
    def __add__(a, b):
        p = min(a.prec, b.prec); d = dict(a.d)
        for k, v in b.d.items(): d[k] = d.get(k, 0) + v
        return S(d, p)
    def __neg__(a): return S({k: -v for k, v in a.d.items()}, a.prec)
    def __sub__(a, b): return a + (-b)
    def __mul__(a, b):
        if not isinstance(b, S): return S({k: v * Fr(b) for k, v in a.d.items()}, a.prec)
        p = min(a.prec + b.val(), b.prec + a.val()); d = {}
        for k, v in a.d.items():
            for l, w in b.d.items():
                if k + l < p: d[k + l] = d.get(k + l, 0) + v * w
        return S(d, p)
    __rmul__ = __mul__
    def der(a): return S({k - 1: k * v for k, v in a.d.items()}, a.prec - 1)
    def inv(a):
        v = a.val(); assert a.d, "inverse of zero"
        n = a.prec - v  # unit part known mod u^n
        f = [a.d.get(v + i, Fr(0)) for i in range(n)]
        g = [Fr(0)] * n; g[0] = 1 / f[0]
        for i in range(1, n):
            g[i] = -sum(f[k] * g[i - k] for k in range(1, i + 1)) / f[0]
        return S({i - v: g[i] for i in range(n)}, n - v)
    def __truediv__(a, b): return a * b.inv()
    def integ(a):
        assert a.d.get(-1, 0) == 0, "log term"
        return S({k + 1: v / (k + 1) for k, v in a.d.items()}, a.prec + 1)
    def powq(a, al):
        """binomial power of a series 1 + O(u)"""
        assert a.val() == 0 and a.d[0] == 1
        n = a.prec; f = [a.d.get(i, Fr(0)) for i in range(n)]; g = [Fr(0)] * n; g[0] = Fr(1)
        for m in range(1, n):
            g[m] = sum((al * k - (m - k)) * f[k] * g[m - k] for k in range(1, m + 1)) / m
        return S({i: g[i] for i in range(n)}, n)
    def c(a, k): return a.d.get(k, Fr(0))
    def iszero(a, upto):
        assert a.prec >= upto, ("insufficient precision", a.prec, upto)
        return all(k >= upto for k in a.d)

N = 34
def const(x): return S({0: x}, N)
u = S({1: 1}, N)
ZERO = S({}, N)

# ---------------- (1) equal products for v(u) = u tau^(-5/4)
tauinv = S({n: Fr(4, 5 * n + 4) for n in range(N)}, N)      # 1/tau = sum 4u^n/(5n+4)
tau = tauinv.inv()
v = u * tau.powq(Fr(-5, 4))
one_minus_u = S({0: 1, 1: -1}, N)
lhs = v.der() * v.der() * v.der() * v.der() * v.der() / v
rhs = (u * one_minus_u * one_minus_u * one_minus_u * one_minus_u * one_minus_u).inv()
print("(1) v'^5/v - 1/(u(1-u)^5) vanishes mod u^20:", (lhs - rhs).iszero(20),
      "| tau = 1 - 4/9 u + ...:", tau.c(1) == Fr(-4, 9))

# ---------------- (2) Theorem U (a) for a random pair with equal products
random.seed(20261117)
def runit(c0=None):
    d = {0: Fr(c0 if c0 is not None else random.choice([1, -1, 2, -3, Fr(1, 2)]))}
    for k in range(1, 4): d[k] = Fr(random.randint(-3, 3), random.randint(1, 3))
    return S(d, N)
uinv = S({-1: 1}, N)
bt = [runit(), runit() * uinv, runit(), runit(), runit()]            # tilde b
b = [runit(), runit() * uinv, runit(), runit(), None]
Bt = bt[0] * bt[1] * bt[2] * bt[3] * bt[4]
b[4] = Bt / (b[0] * b[1] * b[2] * b[3])                               # (EP)
assert b[4].val() == 0
I5 = range(5)
def Bm(m, a):                                                          # B^{(m)}_a
    r = const(1)
    for i in range(m): r = r * bt[(a + i) % 5]
    return r
dU = [const(1)]
for a in range(4): dU.append(dU[a] * bt[a] / b[a])
print("    cyclic closure of d^U:", (dU[4] * bt[4] / b[4] - dU[0]).iszero(20))
det = dU[0] * dU[1] * dU[2] * dU[3] * dU[4]
dd = (det * (1 / det.c(0))).powq(Fr(-1, 5))                            # fraktur d
dlog = [dU[a].der() / dU[a] for a in I5]
trdlog = dlog[0] + dlog[1] + dlog[2] + dlog[3] + dlog[4]
print("    d' = -(1/5) tr(dlog Theta) d:", (dd.der() + Fr(1, 5) * trdlog * dd).iszero(20))

def mat0(): return [[S({}, N) for _ in I5] for _ in I5]
def mder(M): return [[M[a][c].der() for c in I5] for a in I5]
def madd(M, Q): return [[M[a][c] + Q[a][c] for c in I5] for a in I5]
def dlogmul(M): return [[dlog[a] * M[a][c] for c in I5] for a in I5]
Ycal = lambda j: [[(Bm(j, c).inv() if (c + j) % 5 == a else S({}, N)) for c in I5] for a in I5]  # Y^j: entries (c+j,c)
def integral(M, upto): return all(M[a][c].val() >= 0 or M[a][c].iszero(upto) for a in I5 for c in I5)

Bcirc = u * Bt
K = 3
Pc = []
P0 = mat0()
for a in I5: P0[a][a] = dd
Pc.append(P0)
ok = True
for k in range(K):
    X = madd(mder(Pc[k]), dlogmul(Pc[k]))
    Zl = mat0()
    for mp in range(1, 6):                 # shift m+1 = mp
        m = mp - 1
        xi = [X[a][(a + mp) % 5] for a in I5]
        sig = [xi[a] / Bm(mp, a) for a in I5]
        tot = sig[0] + sig[1] + sig[2] + sig[3] + sig[4]
        ok &= tot.iszero(12)                                         # (S_k): pr(X)=0
        ok &= all(s.val() >= 0 for s in sig)
        rho = [None] * 5; rho[1] = S({}, N)
        for a in (1, 2, 3, 4): rho[(a + 1) % 5] = rho[a] + sig[a]
        for a in I5: Zl[a][(a + m) % 5] = rho[a] * Bm(m, a)
    ok &= integral(Zl, 12)
    Zz = madd(mder(Zl), dlogmul(Zl))
    w = mat0()
    for j in range(5):
        mp = 5 - j
        phi = S({}, N)
        for a in I5: phi = phi + Zz[a][(a + mp) % 5] * Bm(5 - mp, (a + mp) % 5)
        phi = Fr(1, 5) * phi
        if j == 0:
            wj = dd * (const(Fr(random.randint(-2, 2))) - (phi / dd).integ())
        else:
            assert phi.val() >= -1
            chic = Fr(-j, 5) * Bcirc.der() / Bcirc + Fr(1, 5) * trdlog
            uphi = u * phi
            nmax = min(uphi.prec, chic.prec) - 1
            co = {}
            for n in range(nmax):
                rhs_ = -uphi.c(n) - sum(chic.c(n - 1 - l) * co[l] for l in range(n))
                co[n] = rhs_ / (n + Fr(j, 5))
            wj = S(co, nmax)
        Yj = Ycal(j)
        for a in I5:
            for c in I5:
                if (c + j) % 5 == a: w[a][c] = w[a][c] + wj * Yj[a][c]
    Pn = madd(Zl, w)
    ok &= integral(Pn, 12)
    Pc.append(Pn)
print("(2) (S_k), integrality of sigma, Z_k, check-P_k for k<=%d:" % K, ok)

# gauge equation for P = Theta_U * sum z^k checkP_k
P = [[[dU[a] * Pc[k][a][c] for c in I5] for a in I5] for k in range(K + 1)]
def AP_minus_PAt(M):
    R = mat0()
    for a in I5:
        for c in I5:
            R[a][c] = b[a] * M[(a + 1) % 5][c] - M[a][(c - 1) % 5] * bt[(c - 1) % 5]
    return R
res0 = AP_minus_PAt(P[0])
good = all(res0[a][c].iszero(10) for a in I5 for c in I5)
for k in range(K):
    R = AP_minus_PAt(P[k + 1]); D = mder(P[k])
    good &= all((R[a][c] - D[a][c]).iszero(10) for a in I5 for c in I5)
print("    gauge equation z P' = A P - P At holds at z^0..z^%d mod u^10:" % K, good)
print("    P_0 = d*Theta_U diagonal with unit entries:", all(P[0][a][a].val() == 0 for a in I5))
# non-triviality: off-diagonal entries of P_1 are non-zero
print("    P_1 has non-zero entries:", sum(1 for a in I5 for c in I5 if P[1][a][c].d), "of 25")
