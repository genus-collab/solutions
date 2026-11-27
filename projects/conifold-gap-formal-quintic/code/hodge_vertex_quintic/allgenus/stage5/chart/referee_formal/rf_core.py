# -*- coding: utf-8 -*-
"""rf_core.py -- Stage 5 / chart / referee_formal: shared exact arithmetic (stdlib only; written for this referee report,
shares no code with the worker's ch_*.py or with any Stage-4 module).

Rings (all with +, -, *, / and Fraction scalars):
  Fraction (Q);  QZ = Q(zeta_5) exact, basis 1, zeta, zeta^2, zeta^3;  FP = F_p with a primitive 5th root of unity;
  LS = truncated Laurent series in eps over Q with absolute-precision tracking;  TS = Q[[u]]/(u^n).
One-variable formal stationary phase by the Morse lemma + Lagrange-Buermann inversion (sp1), NOT by Wick contraction:
  N_n = Sum_(m+k=n) (2k-1)!! C^k [t^(2k)] ( A_m(t) r(t)^(-(2k+1)/2) ),  r = (phase(t) - phase(0))/(ph2 t^2),  C = -1/(2 ph2).
"""
from fractions import Fraction as Fr
from math import comb

INF = 10 ** 9


class QZ:
    __slots__ = ('c',)

    def __init__(s, c): s.c = tuple(c)

    @staticmethod
    def of(x):
        return x if isinstance(x, QZ) else QZ((Fr(x), Fr(0), Fr(0), Fr(0)))

    @staticmethod
    def zeta(k=1):
        k %= 5
        if k == 4: return QZ((Fr(-1),) * 4)
        c = [Fr(0)] * 4; c[k] = Fr(1); return QZ(c)

    def __add__(s, o):
        if not isinstance(o, QZ): o = QZ.of(o)
        a, b = s.c, o.c; return QZ((a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]))
    __radd__ = __add__

    def __neg__(s): return QZ(tuple(-a for a in s.c))
    def __sub__(s, o): return s + (-QZ.of(o))
    def __rsub__(s, o): return QZ.of(o) + (-s)

    def __mul__(s, o):
        a0, a1, a2, a3 = s.c
        if not isinstance(o, QZ):
            o = Fr(o); return QZ((a0 * o, a1 * o, a2 * o, a3 * o))
        b0, b1, b2, b3 = o.c
        r0 = a0 * b0 + a2 * b3 + a3 * b2; r1 = a0 * b1 + a1 * b0 + a3 * b3
        r2 = a0 * b2 + a1 * b1 + a2 * b0; r3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
        r4 = a1 * b3 + a2 * b2 + a3 * b1
        return QZ((r0 - r4, r1 - r4, r2 - r4, r3 - r4))
    __rmul__ = __mul__

    def gal(s, k):
        r = [Fr(0)] * 5
        for i, a in enumerate(s.c): r[(i * k) % 5] += a
        return QZ((r[0] - r[4], r[1] - r[4], r[2] - r[4], r[3] - r[4]))

    def inv(s):
        t = s.gal(2) * s.gal(3) * s.gal(4); n = s * t
        assert n.c[1] == 0 and n.c[2] == 0 and n.c[3] == 0 and n.c[0] != 0, 'QZ: not invertible'
        return t * (1 / n.c[0])

    def __truediv__(s, o): return s * (o.inv() if isinstance(o, QZ) else 1 / Fr(o))
    def __rtruediv__(s, o): return QZ.of(o) * s.inv()
    def __eq__(s, o): return s.c == QZ.of(o).c
    def __hash__(s): return hash(s.c)
    def rat(s): return s.c[0] if s.c[1] == s.c[2] == s.c[3] == 0 else None
    def __repr__(s): return 'QZ' + str(tuple(str(a) for a in s.c))


class FP:
    __slots__ = ('v',)
    p = 0; z1 = 0

    def __init__(s, v): s.v = v

    @staticmethod
    def setup(p):
        FP.p = p
        assert (p - 1) % 5 == 0
        for g in range(2, 500):
            z = pow(g, (p - 1) // 5, p)
            if z != 1: FP.z1 = z; return

    @staticmethod
    def of(x):
        if isinstance(x, FP): return x
        if isinstance(x, int): return FP(x % FP.p)
        x = Fr(x); return FP(x.numerator * pow(x.denominator, -1, FP.p) % FP.p)

    @staticmethod
    def zeta(k=1): return FP(pow(FP.z1, k % 5, FP.p))

    def __add__(s, o): return FP((s.v + (o.v if isinstance(o, FP) else FP.of(o).v)) % FP.p)
    __radd__ = __add__
    def __neg__(s): return FP(-s.v % FP.p)
    def __sub__(s, o): return FP((s.v - (o.v if isinstance(o, FP) else FP.of(o).v)) % FP.p)
    def __rsub__(s, o): return FP((FP.of(o).v - s.v) % FP.p)
    def __mul__(s, o): return FP(s.v * (o.v if isinstance(o, FP) else FP.of(o).v) % FP.p)
    __rmul__ = __mul__
    def inv(s): return FP(pow(s.v, -1, FP.p))
    def __truediv__(s, o): return s * FP.of(o).inv()
    def __rtruediv__(s, o): return FP.of(o) * s.inv()
    def __eq__(s, o): return s.v == FP.of(o).v
    def __hash__(s): return hash(s.v)
    def __repr__(s): return 'FP(%d)' % s.v


def is_prime(n):
    """Miller-Rabin with the bases 2..41: deterministic for n < 3.317e24 (Sorenson-Webster)."""
    if n < 2: return False
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    for b in bases:
        if n % b == 0: return n == b
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for b in bases:
        x = pow(b, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


def prime_1mod10(start):
    n = start - start % 10 + 11
    while not is_prime(n): n += 10
    return n


def bernoulli(n):
    B = [Fr(1)]
    for m in range(1, n + 1):
        B.append(-sum(comb(m + 1, j) * B[j] for j in range(m)) / (m + 1))
    return B


def bk_list(K):
    """b_k = B_2k/(2k(2k-1)), k = 0..K (b_0 unused = 0)"""
    B = bernoulli(2 * K)
    return [Fr(0)] + [B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, K + 1)]


def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r


def gbinom(a, r):
    """generalized binomial coefficient binom(a, r), a a Fraction or int"""
    out = Fr(1)
    for i in range(r): out = out * (Fr(a) - i) / (i + 1)
    return out


def ser_mul(a, b, n, zero):
    out = [zero] * (n + 1)
    for i in range(min(n, len(a) - 1) + 1):
        ai = a[i]
        for j in range(min(n - i, len(b) - 1) + 1): out[i + j] = out[i + j] + ai * b[j]
    return out


def ser_exp(a, n, one):
    """exp of a series with a[0] = 0, to t^n"""
    zero = one - one; e = [one] + [zero] * n
    for k in range(1, n + 1):
        s = zero
        for i in range(1, min(k, len(a) - 1) + 1): s = s + a[i] * e[k - i] * i
        e[k] = s * Fr(1, k)
    return e


def ser_pow1(r, alpha, n, one):
    """r^alpha for a series with r[0] = 1 (J.C.P. Miller recurrence), to t^n"""
    zero = one - one; f = [one] + [zero] * n
    for k in range(1, n + 1):
        s = zero
        for i in range(1, min(k, len(r) - 1) + 1): s = s + r[i] * f[k - i] * (alpha * i - k + i)
        f[k] = s * Fr(1, k)
    return f


def sp1(ph, amps, N, one):
    """one-variable formal SP, unnormalised N_0..N_N. ph[k] = [t^k](phase(t) - phase(0)) for k = 2..2N+2 (ph[0], ph[1]
    ignored; ph[1] must vanish); amps[m][j] = [z^m t^j] amplitude, j <= 2(N-m)."""
    zero = one - one; i2 = 1 / ph[2]
    r = [ph[i + 2] * i2 for i in range(2 * N + 1)]
    C = -(i2 * Fr(1, 2)); Cp = [one]
    for k in range(N): Cp.append(Cp[-1] * C)
    fk = [ser_pow1(r, Fr(-(2 * k + 1), 2), 2 * k, one) for k in range(N + 1)]
    res = []
    for n in range(N + 1):
        tot = zero
        for m in range(n + 1):
            k = n - m; a = amps[m] if m < len(amps) else []; f = fk[k]; s = zero
            for j in range(min(2 * k, len(a) - 1) + 1): s = s + a[j] * f[2 * k - j]
            tot = tot + s * Cp[k] * dfact(2 * k - 1)
        res.append(tot)
    return res


def newton_sums(c, M):
    """power sums p_1..p_M of the roots of the monic x^d + c[0] x^(d-1) + ... + c[d-1]"""
    d = len(c); p = [Fr(d)]
    for m in range(1, M + 1):
        s = Fr(0)
        for i in range(1, min(m - 1, d) + 1): s += c[i - 1] * p[m - i]
        if m <= d: s += m * c[m - 1]
        p.append(-s)
    return p


def h_series(N, twist=Fr(-5)):
    """h_(w_0)(z) through z^N, w_0 = {1 - zeta^j} U {twist}; Newton sums of y = 1/(1 - zeta^j): y^4 - 2y^3 + 2y^2 - y + 1/5"""
    K = N // 2 + 1; b = bk_list(K); p = newton_sums([Fr(-2), Fr(2), Fr(-1), Fr(1, 5)], 2 * K)
    lg = [Fr(0)] * (N + 1)
    for k in range(1, K + 1):
        if 2 * k - 1 <= N: lg[2 * k - 1] = -b[k] * (twist ** (1 - 2 * k) + p[2 * k - 1])
    return ser_exp(lg, N, Fr(1))
