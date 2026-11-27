# -*- coding: utf-8 -*-
"""fq_series.py -- truncated Laurent series in one variable over a coefficient ring (fq_fp.FpRing or QRing)
with tracked precision.  Ser(R, v, c) represents  sum_i c[i] x^(v+i) + O(x^(v+len(c))).
  .prec     exclusive precision (first unknown exponent)
  .coef(k)  coefficient of x^k (raises if k >= prec)
  add/sub/mul/scale/neg/inv/pow/shift(k)=x^k*  /deriv (d/dx, loses one order) /theta (x d/dx)
  exp (needs valuation >= 1) / log (needs constant term 1) / root(n) / sqrt1 (constant term 1)
  compose(g) (g with valuation >= 1) / reverse() (f = x + O(x^2) -> compositional inverse)
  truncate(N) lowers the precision to O(x^N); eq(other) tests equality within the common precision.
solve_lin(R, M, b): Gaussian elimination over Ser for a square system with unit pivots.
Pure Python.
"""
from fractions import Fraction as Fr


class Ser(object):
    __slots__ = ('R', 'v', 'c')

    def __init__(self, R, v, c):
        self.R = R
        # normalise: strip leading zeros
        i = 0; n = len(c)
        while i < n and R.iszero(c[i]): i += 1
        if i:
            v += i; c = c[i:]
        self.v = v; self.c = list(c)

    # ---- constructors
    @staticmethod
    def zero(R, N): return Ser(R, N, [])
    @staticmethod
    def const(R, a, N):
        a = R.of(a)
        if R.iszero(a) or N <= 0: return Ser(R, N, [])
        return Ser(R, 0, [a] + [R.zero] * (N - 1))
    @staticmethod
    def mono(R, k, a, N):
        """a x^k + O(x^N)"""
        a = R.of(a)
        if N <= k or R.iszero(a): return Ser(R, N, [])
        return Ser(R, k, [a] + [R.zero] * (N - k - 1))
    @staticmethod
    def x(R, N): return Ser.mono(R, 1, 1, N)
    @staticmethod
    def from_list(R, lst, v=0, N=None):
        """sum lst[i] x^(v+i), exact up to O(x^(v+len)) (or to O(x^N) if N given, padding with zeros)"""
        c = [R.of(a) for a in lst]
        if N is not None:
            if N > v + len(c): c = c + [R.zero] * (N - v - len(c))
            else: c = c[:max(N - v, 0)]
        return Ser(R, v, c)

    # ---- access
    @property
    def prec(self): return self.v + len(self.c)
    def coef(self, k):
        if k >= self.prec: raise ValueError("coefficient x^%d beyond precision O(x^%d)" % (k, self.prec))
        if k < self.v: return self.R.zero
        return self.c[k - self.v]
    def coefs(self, k0, k1):
        return [self.coef(k) for k in range(k0, k1)]
    def iszero(self): return len(self.c) == 0
    def val(self): return self.v
    def truncate(self, N):
        if N >= self.prec: return self
        return Ser(self.R, self.v, self.c[:max(N - self.v, 0)])
    def copy(self): return Ser(self.R, self.v, self.c)

    # ---- arithmetic
    def add(self, o):
        R = self.R
        prec = min(self.prec, o.prec); v = min(self.v, o.v)
        out = [R.zero] * max(prec - v, 0)
        for i, x in enumerate(self.c):
            j = self.v + i - v
            if j < len(out): out[j] = x
        for i, x in enumerate(o.c):
            j = o.v + i - v
            if j < len(out): out[j] = R.add(out[j], x)
        return Ser(R, v, out)
    def neg(self): return Ser(self.R, self.v, [self.R.neg(x) for x in self.c])
    def sub(self, o): return self.add(o.neg())
    def scale(self, k):
        k = self.R.of(k)
        if self.R.iszero(k): return Ser(self.R, self.prec, [])
        return Ser(self.R, self.v, [self.R.mul(k, x) for x in self.c])
    def shift(self, k): return Ser(self.R, self.v + k, self.c)
    def mul(self, o, N=None):
        R = self.R
        prec = min(self.prec + o.v, o.prec + self.v)
        if N is not None: prec = min(prec, N)
        v = self.v + o.v
        n = prec - v
        if n <= 0: return Ser(R, prec, [])
        out = [R.zero] * n
        a = self.c; b = o.c
        if R.is_fp:
            p = R.p
            for i, x in enumerate(a):
                if i >= n: break
                if x == 0: continue
                for j in range(min(len(b), n - i)):
                    y = b[j]
                    if y: out[i + j] = (out[i + j] + x * y) % p
        else:
            for i, x in enumerate(a):
                if i >= n: break
                if x == 0: continue
                for j in range(min(len(b), n - i)):
                    out[i + j] += x * b[j]
        return Ser(R, v, out)
    def inv(self):
        R = self.R
        if not self.c: raise ZeroDivisionError("inverse of a series that is zero to its precision")
        n = len(self.c); a = self.c
        out = [R.zero] * n; out[0] = R.inv(a[0])
        for k in range(1, n):
            s = R.zero
            for j in range(1, k + 1):
                if not R.iszero(a[j]): s = R.add(s, R.mul(a[j], out[k - j]))
            out[k] = R.neg(R.mul(s, out[0]))
        return Ser(R, -self.v, out)
    def div(self, o): return self.mul(o.inv())
    def pow(self, n):
        if n < 0: return self.inv().pow(-n)
        if n == 0:
            return Ser.const(self.R, 1, self.prec)  # precision convention: O(x^prec)
        r = self
        for _ in range(n - 1): r = r.mul(self)
        return r
    def deriv(self):
        R = self.R
        return Ser(R, self.v - 1, [R.mul(R.of(self.v + i), x) for i, x in enumerate(self.c)])
    def theta(self):
        R = self.R
        return Ser(R, self.v, [R.mul(R.of(self.v + i), x) for i, x in enumerate(self.c)])
    def integrate_theta(self, c0=0):
        """f with theta f = self (requires zero constant term); constant term c0"""
        R = self.R
        out = [R.zero] * (self.prec)
        for i, x in enumerate(self.c):
            k = self.v + i
            if k == 0:
                if not R.iszero(x): raise ValueError("theta-integral of a constant term")
                continue
            if k < 0: raise ValueError("negative exponent in integrate_theta")
            out[k] = R.mul(x, R.inv(R.of(k)))
        out[0] = R.of(c0)
        return Ser(R, 0, out)
    def exp(self):
        R = self.R
        if self.c and self.v < 1: raise ValueError("exp needs valuation >= 1")
        N = self.prec
        f = [self.coef(k) if k < N else R.zero for k in range(N)]
        out = [R.zero] * N; out[0] = R.one
        for n in range(1, N):
            t = R.zero
            for k in range(1, n + 1):
                if not R.iszero(f[k]): t = R.add(t, R.mul(R.mul(R.of(k), f[k]), out[n - k]))
            out[n] = R.mul(t, R.inv(R.of(n)))
        return Ser(R, 0, out)
    def log(self):
        """log of a series with constant term 1"""
        R = self.R
        if self.v != 0 or not R.eq(self.c[0], R.one): raise ValueError("log needs constant term 1")
        return self.theta().mul(self.inv()).integrate_theta(0)
    def root(self, n):
        return self.log().scale(Fr(1, n)).exp()
    def sqrt1(self): return self.root(2)
    def compose(self, g):
        """self(g), g with valuation >= 1 (Horner over the power-series part; self must have v >= 0)"""
        R = self.R
        if self.v < 0: raise ValueError("compose needs a power series")
        if g.c and g.v < 1: raise ValueError("compose: inner series needs valuation >= 1")
        N = min(self.prec * max(g.v, 1), g.prec) if g.c else self.prec
        res = Ser.zero(R, N)
        for k in range(self.prec - 1, -1, -1):
            res = res.mul(g, N).add(Ser.const(R, self.coef(k), N))
        return res
    def reverse(self):
        """compositional inverse of f = x + O(x^2)"""
        R = self.R
        if self.v != 1 or not R.eq(self.c[0], R.one): raise ValueError("reverse needs f = x + O(x^2)")
        N = self.prec
        g = [R.zero, R.one] + [R.zero] * (N - 2)
        for n in range(2, N):
            gs = Ser(R, 0, g[:n]) if any(not R.iszero(t) for t in g[:n]) else Ser.zero(R, n)
            gs = Ser.from_list(R, g, 0, N)
            fg = self.compose(gs)
            g[n] = R.neg(R.sub(fg.coef(n), g[n]))
        return Ser.from_list(R, g, 0, N)
    def eq(self, o, N=None):
        d = self.sub(o)
        if N is not None: d = d.truncate(N)
        return d.iszero()
    def fmt(self, kmax=None):
        R = self.R
        terms = []
        for i, x in enumerate(self.c):
            if R.iszero(x): continue
            k = self.v + i
            if kmax is not None and k > kmax: break
            terms.append("(%s)x^%d" % (R.fmt(x), k))
        return (" + ".join(terms) if terms else "0") + " + O(x^%d)" % self.prec
    def __repr__(self): return self.fmt()


def solve_lin(R, M, b):
    """solve M x = b for square M of Ser with unit pivots (Gaussian elimination); returns list of Ser"""
    n = len(M)
    A = [row[:] + [b[i]] for i, row in enumerate(M)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            s = A[r][col]
            if s.c and s.v == 0:
                piv = r; break
        if piv is None: raise ValueError("no unit pivot in column %d" % col)
        A[col], A[piv] = A[piv], A[col]
        pinv = A[col][col].inv()
        A[col] = [x.mul(pinv) for x in A[col]]
        for r in range(n):
            if r == col: continue
            f = A[r][col]
            if f.iszero(): continue
            A[r] = [A[r][j].sub(f.mul(A[col][j])) for j in range(n + 1)]
    return [A[i][n] for i in range(n)]


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from fq_fp import FpRing, QRing
    for R in (QRing(), FpRing()):
        N = 12
        x = Ser.x(R, N)
        f = Ser.const(R, 1, N).add(x)           # 1 + x
        g = f.inv()
        assert g.coef(5) == R.of(-1) and g.prec == N
        assert f.mul(g).eq(Ser.const(R, 1, N))
        e = x.exp(); l = e.log()
        assert l.eq(x)
        assert f.log().coef(3) == R.of(Fr(1, 3))
        r = f.root(5)
        assert r.pow(5).eq(f)
        h = x.add(x.mul(x).scale(3))        # x + 3x^2
        hi = h.reverse()
        assert h.compose(hi).eq(x) and hi.compose(h).eq(x)
        L = Ser(R, -2, [R.of(1), R.of(2)])  # x^-2 + 2 x^-1 + O(1)
        Li = L.inv()
        assert Li.v == 2 and Li.prec == 4 and L.mul(Li).coef(0) == R.one
        assert x.deriv().prec == N - 1 and f.theta().eq(x)
        M = [[f, x], [x, f]]; bb = [Ser.const(R, 1, N), x]
        sol = solve_lin(R, M, bb)
        assert f.mul(sol[0]).add(x.mul(sol[1])).eq(bb[0]) and x.mul(sol[0]).add(f.mul(sol[1])).eq(bb[1])
        print("Ser ok over", type(R).__name__)
