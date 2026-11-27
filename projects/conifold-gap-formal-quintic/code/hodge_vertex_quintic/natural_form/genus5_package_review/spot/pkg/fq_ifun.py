# -*- coding: utf-8 -*-
"""fq_ifun.py -- the equivariant I-function of the formal quintic at the five fixed points, its two expansions,
and the Picard-Fuchs operator.

Conventions (loc_check.py / PLAN Key facts): H|_{p_i} = lambda_i, O(5)|_{p_i} = 5 lambda_i, tangent weights
lambda_i - lambda_j.  With s := 1/z,
  I_i(q,z) = z q^{lambda_i/z} Sum_d q^d c_d^{(i)}(z),
  c_d^{(i)}(z) = prod_{k=1}^{5d}(5 lambda_i + k z) / prod_j prod_{k=1}^d (lambda_i - lambda_j + k z)
              = prod_{k=1}^{5d}(5 lambda_i s + k) / prod_j prod_{k=1}^d ((lambda_i - lambda_j) s + k)   (degree 0 in z).
Picard-Fuchs:  [prod_j (M - lambda_j) - q prod_{k=1}^5 (5M + kz)] I = 0,  M = z q d/dq;  on the q-series
F_i = Sum_d q^d c_d(z) (the factor q^{lambda_i/z} stripped) M acts as M + lambda_i.

Classes:
  IFun(R, lams, D)     c_inf(d,i,N): c_d^{(i)} as Ser in s to O(s^N); c_zero(d,i,N): Laurent Ser in z to O(z^N);
                       pf_check(i, Nz): residual of the PF operator on z F_i (Laurent-in-z q-series), must vanish.
  SL(R)                series in s and l := log q with q-series coefficients: dict {(n, m): Ser(q)} for s^n l^m.
                       add, mul, scale, theta (q d/dq, with theta l = 1), s_shift(k), part(n, m), trunc_s(Ns).
  vandermonde(R, lams) matrices V[i][k] = lambda_i^k and its inverse: f = Sum_k f_k H^k <-> f(p_i).
"""
import math, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_series import Ser


def lin(R, a, b, N):
    """a + b x + O(x^N)"""
    return Ser.from_list(R, [a, b], 0, N)


class IFun(object):
    def __init__(self, R, lams, D):
        self.R = R; self.lams = [R.of(l) for l in lams]; self.D = D
        self.n = len(lams)
        self._ci = {}; self._cz = {}

    def c_inf(self, d, i, N):
        """c_d^{(i)}(s) = prod_{k=1}^{5d}(5 lam_i s + k) / prod_j prod_{k=1}^d ((lam_i - lam_j) s + k), O(s^N)"""
        key = (d, i, N)
        if key in self._ci: return self._ci[key]
        R = self.R; l = self.lams
        num = Ser.const(R, 1, N)
        for k in range(1, 5 * d + 1):
            num = num.mul(lin(R, k, R.mul(R.of(5), l[i]), N))
        den = Ser.const(R, 1, N)
        for j in range(self.n):
            for k in range(1, d + 1):
                den = den.mul(lin(R, k, R.sub(l[i], l[j]), N))
        r = num.mul(den.inv())
        self._ci[key] = r
        return r

    def c_zero(self, d, i, N):
        """c_d^{(i)}(z) as a Laurent series in z (pole order d), to O(z^N)"""
        key = (d, i, N)
        if key in self._cz: return self._cz[key]
        R = self.R; l = self.lams
        M = N + d + 1
        num = Ser.const(R, 1, M)
        for k in range(1, 5 * d + 1):
            num = num.mul(lin(R, R.mul(R.of(5), l[i]), k, M))
        den = Ser.const(R, 1, M)
        for j in range(self.n):
            if j == i: continue
            for k in range(1, d + 1):
                den = den.mul(lin(R, R.sub(l[i], l[j]), k, M))
        r = num.mul(den.inv()).scale(Fr(1, math.factorial(d))).shift(-d).truncate(N)
        self._cz[key] = r
        return r

    def F_zero(self, i, Nz):
        """list over d of c_d^{(i)}(z) Laurent Ser, i.e. F_i = Sum_d q^d c_d(z) = I_i/(z q^{lam_i/z})"""
        return [self.c_zero(d, i, Nz) for d in range(self.D + 1)]

    def pf_check(self, i, Nz):
        """apply [prod_j (M + lam_i - lam_j) - q prod_{k=1}^5 (5(M + lam_i) + k z)] to z F_i; returns list over d of
        Laurent Ser residuals (each should be zero to its precision)"""
        R = self.R; l = self.lams; D = self.D
        F = [c.shift(1) for c in self.F_zero(i, Nz + 6)]        # z F_i
        out = []
        for d in range(D + 1):
            # first part on q^d: M -> d z
            t1 = F[d]
            for j in range(self.n):
                t1 = t1.mul(lin(R, R.sub(l[i], l[j]), d, t1.prec + 2))
            if d >= 1:
                t2 = F[d - 1]
                for k in range(1, 6):
                    t2 = t2.mul(lin(R, R.mul(R.of(5), l[i]), R.of(5 * (d - 1) + k), t2.prec + 2))
                t1 = t1.sub(t2)
            out.append(t1.truncate(Nz))
        return out


class SL(object):
    """series in s and l = log q with q-series coefficients: {(n, m): Ser(q)}"""
    __slots__ = ('R', 'd')

    def __init__(self, R, d=None):
        self.R = R; self.d = {} if d is None else {k: v for k, v in d.items() if not v.iszero()}

    def add(self, o):
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d[k].add(v) if k in d else v
        return SL(self.R, d)

    def sub(self, o): return self.add(o.scale(-1))

    def scale(self, c):
        return SL(self.R, {k: v.scale(c) for k, v in self.d.items()})

    def mulq(self, qs):
        """multiply by a q-series"""
        return SL(self.R, {k: v.mul(qs) for k, v in self.d.items()})

    def mul(self, o, Ns=None):
        d = {}
        for (n1, m1), v1 in self.d.items():
            for (n2, m2), v2 in o.d.items():
                k = (n1 + n2, m1 + m2)
                if Ns is not None and k[0] > Ns: continue
                p = v1.mul(v2)
                d[k] = d[k].add(p) if k in d else p
        return SL(self.R, d)

    def theta(self):
        """q d/dq with theta(l) = 1"""
        R = self.R; d = {}
        for (n, m), v in self.d.items():
            t = v.theta()
            d[(n, m)] = d[(n, m)].add(t) if (n, m) in d else t
            if m >= 1:
                t2 = v.scale(m); k = (n, m - 1)
                d[k] = d[k].add(t2) if k in d else t2
        return SL(R, d)

    def s_shift(self, k): return SL(self.R, {(n + k, m): v for (n, m), v in self.d.items()})
    def trunc_s(self, Ns): return SL(self.R, {k: v for k, v in self.d.items() if k[0] <= Ns})
    def part(self, n, m, N):
        return self.d.get((n, m), Ser.zero(self.R, N))
    def s_part(self, n):
        """dict m -> Ser for the s^n coefficient"""
        return {m: v for (nn, m), v in self.d.items() if nn == n}
    def max_m(self): return max((m for (n, m) in self.d), default=-1)

    @staticmethod
    def exp_ls(R, a, Ns, N):
        """exp(a s l) = Sum_n (a l s)^n / n!, to s^Ns, q-precision N"""
        d = {}
        for n in range(Ns + 1):
            d[(n, n)] = Ser.const(R, R.mul(R.pow(R.of(a), n), R.inv(R.of(math.factorial(n)))), N)
        return SL(R, d)

    @staticmethod
    def from_sser(R, ss, N):
        """s-series with constant coefficients (Ser in s) -> SL with q-constant entries"""
        d = {}
        for n in range(ss.v, ss.prec):
            c = ss.coef(n)
            if not R.iszero(c): d[(n, 0)] = Ser.const(R, c, N)
        return SL(R, d)


def vandermonde(R, lams):
    """V[i][k] = lams[i]^k and Vinv (field matrix inverse by Gauss-Jordan)"""
    n = len(lams)
    V = [[R.pow(R.of(lams[i]), k) for k in range(n)] for i in range(n)]
    A = [row[:] + [R.one if r == c else R.zero for c in range(n)] for r, row in enumerate(V)]
    for col in range(n):
        piv = next(r for r in range(col, n) if not R.iszero(A[r][col]))
        A[col], A[piv] = A[piv], A[col]
        pinv = R.inv(A[col][col])
        A[col] = [R.mul(x, pinv) for x in A[col]]
        for r in range(n):
            if r != col and not R.iszero(A[r][col]):
                f = A[r][col]
                A[r] = [R.sub(A[r][j], R.mul(f, A[col][j])) for j in range(2 * n)]
    Vinv = [row[n:] for row in A]
    return V, Vinv


def components(R, Vinv, vec):
    """f_k = Sum_i Vinv[k][i] f(p_i) for a vector of ring elements or of Ser (needs .add/.scale)"""
    n = len(Vinv)
    if isinstance(vec[0], Ser):
        out = []
        for k in range(n):
            acc = None
            for i in range(n):
                t = vec[i].scale(Vinv[k][i])
                acc = t if acc is None else acc.add(t)
            out.append(acc)
        return out
    return [sum((R.mul(Vinv[k][i], vec[i]) for i in range(n)), R.zero) % R.p if R.is_fp else
            sum((R.mul(Vinv[k][i], vec[i]) for i in range(n)), R.zero) for k in range(n)]


if __name__ == "__main__":
    from fq_fp import FpRing
    R = FpRing()
    lams = [R.z(i) for i in range(5)]
    IF = IFun(R, lams, 6)
    # c_d(s) constant term = (5d)!/(d!)^5
    for d in range(1, 5):
        assert IF.c_inf(d, 0, 4).coef(0) == R.of(Fr(math.factorial(5 * d), math.factorial(d) ** 5))
    # zeta covariance: I_i(q,z) = zeta^i I_0(q, zeta^-i z):  c_d^{(i)}(z) = c_d^{(0)}(zeta^-i z)
    for i in range(5):
        for d in range(1, 4):
            a = IF.c_zero(d, i, 5); b = IF.c_zero(d, 0, 5)
            for k in range(-d, 5):
                assert R.eq(a.coef(k), R.mul(b.coef(k), R.z(-i * k)))
    for i in range(5):
        res = IF.pf_check(i, 6)
        assert all(r.iszero() for r in res), "PF fails at colour %d" % i
    print("fq_ifun ok")
