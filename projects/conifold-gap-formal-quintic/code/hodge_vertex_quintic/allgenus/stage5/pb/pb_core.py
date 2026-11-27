# -*- coding: utf-8 -*-
"""pb_core.py -- Stage-5 'pb' audit, own core (stdlib only; imports no project or worker module).

Field: Q (MOD = 0, Fractions) or F_p (MOD = p).  LS = truncated Laurent series in u with ABSOLUTE precision:
LS(v, c, prec) means sum_i c[i] u^(v+i) + O(u^prec), len(c) = prec - v.  Every operation propagates precision,
so precision loss in the gauge recursion is tracked, not guessed.
Also: Bernoulli numbers, the quintic PF4 at the conifold point u = 0 (Frobenius basis f0, f1, f2, the log
obstruction that fixes the vanishing period g = f1 + gamma f2 + c0 f0), conifold frames, tau, v(u).
"""
from fractions import Fraction as Fr

MOD = 0


def set_field(p):
    global MOD
    MOD = p


def red(x):
    if MOD:
        if isinstance(x, Fr):
            return x.numerator % MOD * pow(x.denominator % MOD, MOD - 2, MOD) % MOD
        return x % MOD
    return x if isinstance(x, Fr) else Fr(x)


def finv(x):
    if MOD:
        x %= MOD
        if x == 0:
            raise ZeroDivisionError('zero mod p')
        return pow(x, MOD - 2, MOD)
    return Fr(1) / x


def iszero(x):
    return (x % MOD == 0) if MOD else (x == 0)


class LS(object):
    __slots__ = ('v', 'c', 'prec')

    def __init__(self, v, c, prec):
        assert len(c) == prec - v, (v, len(c), prec)
        self.v = v; self.c = c; self.prec = prec

    @staticmethod
    def from_list(lst, v=0, prec=None):
        if prec is None:
            prec = v + len(lst)
        c = [red(x) for x in lst[:prec - v]]
        c += [red(0)] * (prec - v - len(c))
        return LS(v, c, prec)

    @staticmethod
    def const(a, prec):
        return LS.from_list([a], 0, prec)

    def coef(self, n):
        if n >= self.prec:
            raise ValueError('coefficient u^%d beyond precision %d' % (n, self.prec))
        return self.c[n - self.v] if n >= self.v else red(0)

    def norm(self):
        i = 0
        while i < len(self.c) and iszero(self.c[i]):
            i += 1
        return LS(self.v + i, self.c[i:], self.prec)

    def trunc(self, prec):
        prec = min(prec, self.prec)
        if prec <= self.v:
            return LS(prec, [], prec)
        return LS(self.v, self.c[:prec - self.v], prec)

    def _binop(self, o, sgn):
        prec = min(self.prec, o.prec); v = min(self.v, o.v)
        if prec <= v:
            return LS(prec, [], prec)
        c = [0] * (prec - v)
        for i, x in enumerate(self.c[:prec - self.v]):
            c[self.v - v + i] += x
        for i, x in enumerate(o.c[:prec - o.v]):
            c[o.v - v + i] += sgn * x
        return LS(v, [red(x) for x in c], prec)

    def add(self, o): return self._binop(o, 1)
    def sub(self, o): return self._binop(o, -1)
    def scale(self, k):
        k = red(k)
        return LS(self.v, [red(k * x) for x in self.c], self.prec)
    def neg(self): return self.scale(-1)
    def shift(self, m): return LS(self.v + m, list(self.c), self.prec + m)

    def mul(self, o):
        a = self.norm(); b = o.norm()
        v = a.v + b.v
        prec = min(a.v + b.prec, b.v + a.prec)
        n = prec - v
        if n <= 0:
            return LS(prec, [], prec)
        res = [0] * n
        bc = b.c
        for i, x in enumerate(a.c[:n]):
            if x:
                m = min(len(bc), n - i)
                for j in range(m):
                    res[i + j] += x * bc[j]
        return LS(v, [red(x) for x in res], prec)

    def deriv(self):
        return LS(self.v - 1, [red((self.v + i) * x) for i, x in enumerate(self.c)], self.prec - 1)

    def integ(self):
        """int_0^u; the integrand must have no negative powers of u (asserted)."""
        s = self.norm()
        if not s.c:
            return LS(s.prec + 1, [], s.prec + 1)
        assert s.v >= 0, 'negative powers in integrand'
        return LS(s.v + 1, [red(x * finv(s.v + 1 + i)) for i, x in enumerate(s.c)], s.prec + 1)

    def inv(self):
        a = self.norm()
        if not a.c:
            raise ZeroDivisionError('series is O(u^prec)')
        w = a.v; r = a.prec - w
        i0 = finv(a.c[0]); b = [i0]
        for n in range(1, r):
            s = 0
            for k in range(1, n + 1):
                s += a.c[k] * b[n - k]
            b.append(red(-s * i0))
        return LS(-w, b, -w + r)

    def div(self, o): return self.mul(o.inv())

    def power(self, r):
        """self^r for self = 1 + O(u) (v = 0, c[0] = 1), r rational."""
        a = self.norm(); assert a.v == 0 and a.c[0] == red(1)
        r = red(Fr(r)); n = a.prec; b = [red(1)]
        for m in range(1, n):
            s = 0
            for k in range(1, m + 1):
                s += (r * k - (m - k)) * a.c[k] * b[m - k]
            b.append(red(s * finv(m)))
        return LS(0, b, n)

    def ok_zero(self, upto=None):
        top = self.prec if upto is None else min(upto, self.prec)
        return all(iszero(self.coef(n)) for n in range(self.v, top))


def bernoulli(n):
    """B_0..B_n exact (B_1 = -1/2)."""
    B = [Fr(1)]
    for m in range(1, n + 1):
        s = Fr(0); binom = 1
        for k in range(m):
            s += binom * B[k]
            binom = binom * (m + 1 - k) // (k + 1)
        B.append(-s / (m + 1))
    return B


# ---------- the quintic PF4 at the conifold point ----------
PK = {4: (0, 1), 3: (-2, 2), 2: (Fr(-7, 5), Fr(7, 5)), 1: (Fr(-2, 5), Fr(2, 5)), 0: (Fr(-24, 625), Fr(24, 625))}
# p_k(u) = PK[k][0] + PK[k][1] u ;  P = sum_k p_k D^k,  D = q d/dq = -(1-u) d/du


def D_pair(A, B):
    """D on A + B log u (A, B in LS): D(log u) = -(1-u)/u."""
    one_m = LS.from_list([1, -1], 0, A.prec + 2)
    dA = one_m.mul(A.deriv()).neg()
    dB = one_m.mul(B.deriv()).neg()
    corr = B.mul(one_m).shift(-1).neg()
    return dA.add(corr), dB


def PF4_pair(A, B):
    terms = []; cur = (A, B)
    for k in range(5):
        pk = LS.from_list([PK[k][0], PK[k][1]], 0, cur[0].prec + 5)
        terms.append((pk.mul(cur[0]), pk.mul(cur[1])))
        cur = D_pair(*cur)
    rA = terms[0][0]; rB = terms[0][1]
    for x, y in terms[1:]:
        rA = rA.add(x); rB = rB.add(y)
    return rA, rB


def e_coeffs(n, prec_pad=8):
    """P(u^n) = sum_j e_j u^(n+j), j = -3..1 (exact)."""
    mono = LS.from_list([1], n, n + prec_pad)
    zero = LS(n - 6, [red(0)] * (prec_pad + 6), n + prec_pad)
    r, _ = PF4_pair(mono, zero)
    return {j: r.coef(n + j) for j in range(-3, 2)}


class Conifold(object):
    """Frobenius basis of PF4 at u = 0, exact or mod p, to O(u^N)."""

    def __init__(self, N, alpha=0, beta=0):
        self.N = N
        self.E = [e_coeffs(n) for n in range(N + 8)]
        assert all(iszero(self.E[n][-3] - n * (n - 1) ** 2 * (n - 2)) for n in range(N + 8)), 'indicial'
        self.f = [self.analytic(k) for k in range(3)]
        obs = [self.log_obstruction(self.f[k]) for k in range(3)]
        self.obs = obs
        # g = c0 f0 + f1 + gamma f2 with all three obstructions zero: kernel of the 3x3 matrix (columns obs_k)
        M = [[obs[k][i] for k in range(3)] for i in range(3)]
        sol = kernel_f1_normalised(M)
        if sol is None:
            raise RuntimeError('log obstruction: kernel not one-dimensional with f1-coefficient 1')
        self.c0, one, self.gamma = sol
        self.g = self.f[0].scale(self.c0).add(self.f[1]).add(self.f[2].scale(self.gamma))
        self.set_frame(alpha, beta)

    def analytic(self, k):
        N = self.N; a = [red(0)] * (N + 4); a[k] = red(1)
        for n in range(3, N + 4):
            s = 0
            # coefficient of u^(n-3): sum_j e_j(n-3-j) a_(n-3-j) = 0
            for j in (-2, -1, 0, 1):
                idx = n - 3 - j
                if 0 <= idx < n:
                    s += self.E[idx][j] * a[idx]
            a[n] = red(-s * finv(self.E[n][-3]))
        # consistency at n = 0, 1, 2 (resonances): coefficients of u^-3, u^-2, u^-1 vanish
        for m in (-3, -2, -1):
            s = 0
            for j in (-3, -2, -1, 0, 1):
                idx = m - j
                if 0 <= idx < len(a):
                    s += self.E[idx][j] * a[idx]
            assert iszero(s), ('resonance obstructed', k, m)
        return LS.from_list(a, 0, N + 4)

    def log_obstruction(self, g):
        zero = LS(0, [red(0)] * g.prec, g.prec)
        A, B = PF4_pair(zero, g)
        assert B.trunc(g.prec - 4).ok_zero(), 'g is not a solution'
        return [A.coef(m) for m in (-3, -2, -1)]

    def set_frame(self, alpha, beta):
        self.alpha = red(Fr(alpha)); self.beta = red(Fr(beta))
        self.varpi0 = self.f[0].add(self.g.scale(self.alpha)).add(self.f[2].scale(self.beta))
        self.tc = self.g.div(self.varpi0)


def kernel_f1_normalised(M):
    """Solve M x = 0 (3x3) with x[1] = 1; return x if the solution exists and is unique, else None."""
    # columns 0 and 2 are the unknowns; rhs = -column 1
    rows = [[M[i][0], M[i][2], red(-M[i][1])] for i in range(3)]
    # Gaussian elimination on a 3x2 system
    piv = []
    r = 0
    for col in range(2):
        pr = None
        for i in range(r, 3):
            if not iszero(rows[i][col]):
                pr = i; break
        if pr is None:
            return None
        rows[r], rows[pr] = rows[pr], rows[r]
        iv = finv(rows[r][col])
        rows[r] = [red(x * iv) for x in rows[r]]
        for i in range(3):
            if i != r and not iszero(rows[i][col]):
                f = rows[i][col]
                rows[i] = [red(a - f * b) for a, b in zip(rows[i], rows[r])]
        piv.append(r); r += 1
    if not iszero(rows[2][2]):
        return None
    return rows[0][2], red(1), rows[1][2]


def tau_inv(N):
    """1/tau = 2F1(1,4/5;9/5;u) = sum 4/(4+5n) u^n."""
    return LS.from_list([Fr(4, 4 + 5 * n) for n in range(N)], 0, N)


def model_v(N):
    """V = tau^(-5/4), v = u V."""
    V = tau_inv(N).power(Fr(5, 4))
    return V, V.shift(1)
