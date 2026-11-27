# cylib.py -- fresh exact implementation (F_p, truncated Laurent series in t = log Q) of the Givental-Teleman graph sum for
# the resolved conifold x C^2 on a general weight point, with Hodge vertices in arbitrary (possibly colour-dependent) frames.
# Pure python, no third-party modules.  Written from rc_lib.m/hvc_lib.m conventions, not a translation of stored results.
import sys, os, itertools, json, math, random
from fractions import Fraction as Fr
from functools import lru_cache

P = (1 << 127) - 1
def setprime(q):
    global P
    P = q
    _invcache.clear()
_invcache = {}
def inv(a):
    a %= P
    r = _invcache.get(a)
    if r is None:
        if a == 0: raise ZeroDivisionError("inverse of 0 mod p")
        r = pow(a, P - 2, P)
        if len(_invcache) < 100000: _invcache[a] = r
    return r
def fr2p(x):
    x = Fr(x)
    return x.numerator % P * inv(x.denominator) % P
def ratrec(a):
    """rational reconstruction of a mod P; returns Fraction (or None)"""
    a %= P
    if a == 0: return Fr(0)
    r0, r1 = P, a; s0, s1 = 0, 1
    B = math.isqrt(P // 2)
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1)

# ------------------------------------------------------------------ truncated Laurent series: (v, c): sum c[i] t^(v+i) + O(t^(v+len c))
def s_norm(v, c):
    i = 0; n = len(c)
    while i < n and c[i] == 0: i += 1
    if i: return (v + i, c[i:])
    return (v, c)
def s_const(a, N):
    a %= P
    if a == 0: return (N, [])
    return (0, [a] + [0] * (N - 1))
def s_add(a, b):
    va, ca = a; vb, cb = b
    prec = min(va + len(ca), vb + len(cb))
    v = min(va, vb)
    out = [0] * max(prec - v, 0)
    for i, x in enumerate(ca):
        j = va + i - v
        if j < len(out): out[j] = x
    for i, x in enumerate(cb):
        j = vb + i - v
        if j < len(out): out[j] = (out[j] + x) % P
    return s_norm(v, out)
def s_scale(k, a):
    k %= P
    return s_norm(a[0], [x * k % P for x in a[1]])
def s_neg(a): return s_scale(-1, a)
def s_mul(a, b, L=None):
    va, ca = a; vb, cb = b
    n = min(len(ca), len(cb))
    if L is not None and L < n: n = L
    if n <= 0: return (va + vb + max(min(len(ca), len(cb)), 0) if False else va + vb, [])
    out = [0] * n
    for i in range(n):
        x = ca[i]
        if x:
            for j in range(n - i):
                out[i + j] += x * cb[j]
    return s_norm(va + vb, [y % P for y in out])
def s_inv(a):
    v, c = a
    if not c or c[0] == 0: raise ZeroDivisionError("series inverse")
    n = len(c); i0 = inv(c[0]); out = [i0] + [0] * (n - 1)
    for k in range(1, n):
        s = 0
        for j in range(1, k + 1):
            s += c[j] * out[k - j]
        out[k] = (-s % P) * i0 % P
    return (-v, out)
def s_pow(a, k):
    if k < 0: return s_pow(s_inv(a), -k)
    r = None
    for _ in range(k):
        r = a if r is None else s_mul(r, a)
    if r is None: return (0, [1] + [0] * (len(a[1]) - 1))
    return r
def s_der(a):
    v, c = a
    # derivative: t^(v+i) -> (v+i) t^(v+i-1); precision drops by one in absolute order
    out = [(v + i) * x % P for i, x in enumerate(c)]
    return s_norm(v - 1, out)
def s_coef(a, k):
    v, c = a
    if k >= v + len(c): raise ValueError("coefficient beyond precision: t^%d, known to O(t^%d)" % (k, v + len(c)))
    if k < v: return 0
    return c[k - v]
def s_prec(a): return a[0] + len(a[1])
def s_iszero(a): return all(x == 0 for x in a[1])

# ------------------------------------------------------------------ Witten-Kontsevich and Hodge integrals (exact)
def df(n):
    r = 1
    while n > 0:
        r *= n; n -= 2
    return r
@lru_cache(maxsize=None)
def wk(g, d):
    d = tuple(sorted(d)); return _wk(g, d)
@lru_cache(maxsize=None)
def _wk(g, d):
    n = len(d)
    if g < 0 or n == 0 or d[0] < 0 or sum(d) != 3 * g - 3 + n or 2 * g - 2 + n <= 0: return Fr(0)
    if g == 0 and d == (0, 0, 0): return Fr(1)
    if g == 1 and d == (1,): return Fr(1, 24)
    if d[0] == 0:
        rest = d[1:]; tot = Fr(0)
        for j in range(len(rest)):
            if rest[j] > 0:
                r = list(rest); r[j] -= 1
                tot += wk(g, tuple(r))
        return tot
    k = d[0] - 1; rest = d[1:]; tot = Fr(0)
    for j in range(len(rest)):
        r = list(rest); r[j] += k
        tot += Fr(df(2 * k + 2 * rest[j] + 1), df(2 * rest[j] - 1)) * wk(g, tuple(r))
    for r_ in range(k):
        c = df(2 * r_ + 1) * df(2 * (k - 1 - r_) + 1)
        tot += Fr(c, 2) * wk(g - 1, (r_, k - 1 - r_) + rest)
        s = 0
        idx = range(len(rest))
        for m in range(len(rest) + 1):
            for S in itertools.combinations(idx, m):
                A = tuple(rest[i] for i in S); B = tuple(rest[i] for i in idx if i not in S)
                for g1 in range(g + 1):
                    x = wk(g1, (r_,) + A)
                    if x: s += x * wk(g - g1, (k - 1 - r_,) + B)
        tot += Fr(c, 2) * s
    return tot / df(2 * k + 3)

_bern = {}
def bernoulli(n):
    if n in _bern: return _bern[n]
    B = [Fr(0)] * (n + 1)
    A = [Fr(0)] * (n + 1)
    for m in range(n + 1):
        A[m] = Fr(1, m + 1)
        for j in range(m, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        B[m] = A[0]
    if n >= 1: B[1] = Fr(-1, 2)
    for i, b in enumerate(B): _bern[i] = b
    return _bern[n]

def dmax(g):
    if g <= 0: return 0
    if g == 1: return 1
    return min(g * (g + 1) // 2, 3 * g - 3)
USE_CUTOFF = [True]
def hi(g, d, K):
    return _hi(g, tuple(sorted(d)), tuple(sorted(K)))
@lru_cache(maxsize=None)
def _hi(g, d, K):
    """Int_{Mbar_{g,n}} prod ch_k (k in K) prod psi^d via Mumford's GRR formula"""
    n = len(d)
    if g < 0 or 2 * g - 2 + n <= 0: return Fr(0)
    if sum(d) + sum(K) != 3 * g - 3 + n: return Fr(0)
    if not K: return wk(g, d) if n else Fr(0)
    if g == 0: return Fr(0)
    if USE_CUTOFF[0]:
        if sum(K) > dmax(g): return Fr(0)
        if K[-1] >= 2 * g: return Fr(0)
    k = K[-1]; Kr = K[:-1]; l = (k + 1) // 2
    tot = hi(g, d + (k + 1,), Kr)
    for i in range(n):
        dd = list(d); dd[i] += k
        tot -= hi(g, dd, Kr)
    b = Fr(0)
    for a in range(k):
        sgn = -1 if a % 2 else 1
        bb = hi(g - 1, d + (a, k - 1 - a), Kr)
        for m in range(n + 1):
            for S in itertools.combinations(range(n), m):
                A = tuple(d[i] for i in S); B = tuple(d[i] for i in range(n) if i not in S)
                for h in range(g + 1):
                    for mk in range(len(Kr) + 1):
                        for T in itertools.combinations(range(len(Kr)), mk):
                            K1 = tuple(Kr[i] for i in T); K2 = tuple(Kr[i] for i in range(len(Kr)) if i not in T)
                            x = hi(h, A + (a,), K1)
                            if x: bb += x * hi(g - h, B + (k - 1 - a,), K2)
        b += sgn * bb
    tot += b / 2
    return bernoulli(2 * l) / math.factorial(2 * l) * tot

def partitions_odd(D, maxpart):
    """multisets of odd integers <= maxpart with sum D"""
    def rec(D, mp):
        if D == 0:
            yield ()
            return
        for k in range(min(mp, D), 0, -1):
            if k % 2 == 1:
                for r in rec(D - k, k):
                    yield (k,) + r
    return list(rec(D, maxpart))

def lam_from_ch(g, nmax):
    """lambda_i as polynomials in ch_odd:  c_t(E) = exp(sum (-1)^(k-1) (k-1)! ch_k t^k); returns dict i-> {K(tuple): Fraction}"""
    # exponent
    res = {0: {(): Fr(1)}}
    # use recursion  i lambda_i = sum_{k=1..i} (-1)^(k-1) k! ch_k lambda_{i-k}, ch_even = 0
    for i in range(1, nmax + 1):
        cur = {}
        for k in range(1, i + 1, 2):
            c = Fr(math.factorial(k), i)  # k odd: (-1)^(k-1) = +1
            for K, v in res[i - k].items():
                K2 = tuple(sorted(K + (k,)))
                cur[K2] = cur.get(K2, Fr(0)) + c * v
        res[i] = cur
    return res

def int_lam(g, mono, d):
    """Int prod_{i in mono} lambda_i prod psi^d, via ch"""
    L = lam_from_ch(g, max(mono) if mono else 0)
    poly = {(): Fr(1)}
    for i in mono:
        new = {}
        for K, v in poly.items():
            for K2, w in L[i].items():
                KK = tuple(sorted(K + K2))
                new[KK] = new.get(KK, Fr(0)) + v * w
        poly = new
    return sum((v * hi(g, d, K) for K, v in poly.items()), Fr(0))

# ------------------------------------------------------------------ Hodge vertex
def class_s(w):
    """exponents s_k with prod_a Lambda^vee(w_a)/w_a^g = exp(sum s_k ch_k), k odd;  w = list of residues mod p (or None for the point class).
       Returns function k -> s_k mod p"""
    if w is None: return lambda k: 0
    def s(k):
        tot = 0
        for x in w:
            tot += pow(inv(x), k, P)
        return (-math.factorial(k - 1) * tot) % P
    return s

class VertexTable:
    """integrals of the class with weights w on Mbar_{g,n} with psi^d; cached"""
    def __init__(self, w):
        self.w = None if w is None else [x % P for x in w]
        self.s = class_s(self.w)
        self.cache = {}
        self.sk = {}
    def sval(self, k):
        r = self.sk.get(k)
        if r is None:
            r = self.s(k); self.sk[k] = r
        return r
    def val(self, g, d):
        d = tuple(sorted(d))
        key = (g, d)
        r = self.cache.get(key)
        if r is not None: return r
        n = len(d)
        D = 3 * g - 3 + n - sum(d)
        if D < 0 or 2 * g - 2 + n <= 0:
            r = 0
        elif self.w is None:
            r = fr2p(wk(g, d)) if D == 0 and n else 0
        elif D > dmax(g):
            r = 0
        else:
            r = 0
            for K in partitions_odd(D, max(2 * g - 1, 1)):
                x = hi(g, d, K)
                if x == 0: continue
                c = fr2p(x)
                from collections import Counter
                for k, m in Counter(K).items():
                    c = c * pow(self.sval(k), m, P) % P * inv(math.factorial(m)) % P
                r = (r + c) % P
        self.cache[key] = r
        return r

def h_coeffs(w, zmax, s=-1):
    """coefficients in z of Exp[s * ex], ex = -Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k).  w residues mod p; None -> 1"""
    out = [0] * (zmax + 1); out[0] = 1
    if w is None: return out
    ex = [0] * (zmax + 1)
    k = 1
    while 2 * k - 1 <= zmax:
        c = fr2p(bernoulli(2 * k) / (2 * k * (2 * k - 1)))
        tot = 0
        for x in w: tot += pow(inv(x), 2 * k - 1, P)
        ex[2 * k - 1] = (-c * tot * s) % P
        k += 1
    return series_exp(ex, zmax)
def series_exp(ex, zmax):
    out = [0] * (zmax + 1); out[0] = 1
    for n in range(1, zmax + 1):
        t = 0
        for k in range(1, n + 1):
            t += k * ex[k] * out[n - k]
        out[n] = t % P * inv(n) % P
    return out
def ffc_closed(W, zmax):
    """QRR constant: Exp[-Sum B_2k/(2k(2k-1)) z^(2k-1) (a1^-(2k-1) + (a0+a2)^-(2k-1) - e1^-(2k-1) - e2^-(2k-1))]"""
    ex = [0] * (zmax + 1)
    k = 1
    while 2 * k - 1 <= zmax:
        c = fr2p(bernoulli(2 * k) / (2 * k * (2 * k - 1)))
        m = 2 * k - 1
        S = pow(inv(W['a1']), m, P) + pow(inv((W['a0'] + W['a2']) % P), m, P) - pow(inv(W['e1']), m, P) - pow(inv(W['e2']), m, P)
        ex[m] = (-c * S) % P
        k += 1
    return series_exp(ex, zmax)

def qrr_target(wts, zmax):
    """Exp[-Sum_k B_(2k+2)/((2k+2)(2k+1)) z^(2k+1) Sum_w w^-(2k+1)]"""
    ex = [0] * (zmax + 1)
    k = 0
    while 2 * k + 1 <= zmax:
        c = fr2p(bernoulli(2 * k + 2) / ((2 * k + 2) * (2 * k + 1)))
        tot = 0
        for w in wts: tot += pow(inv(w), 2 * k + 1, P)
        ex[2 * k + 1] = (-c * tot) % P
        k += 1
    return series_exp(ex, zmax)

# ------------------------------------------------------------------ Wick combinatorics
@lru_cache(maxsize=None)
def wick_terms(m):
    """list of (coefficient Fraction, n, tuple of (k, mult)): z^m coefficient of the normalised Laplace expansion"""
    out = []
    if m == 0: return [(Fr(1), 0, ())]
    for n in range(m, 3 * m + 1):
        j = n - m
        if j == 0: continue
        for p in _parts(2 * n, j, 3):
            from collections import Counter
            cnt = Counter(p)
            c = Fr(math.factorial(2 * n) * (-1) ** n, math.factorial(n) * 2 ** n)
            for k, mu in cnt.items():
                c /= math.factorial(k) ** mu * math.factorial(mu)
            out.append((c, n, tuple(sorted(cnt.items()))))
    return out
def _parts(total, nparts, minp):
    if nparts == 0:
        if total == 0: yield ()
        return
    for first in range(minp, total // nparts + 1):
        for r in _parts(total - first, nparts - 1, first):
            yield (first,) + r

def thetapolys(kmax):
    """P_k with theta^k log(1 - cU) = -P_k(A), A = cU/(1-cU): P_1 = A, P_{k+1} = (A + A^2) dP_k/dA.  list of coefficient lists"""
    Ps = {1: [0, 1]}
    for k in range(1, kmax):
        p = Ps[k]
        dp = [i * p[i] for i in range(1, len(p))]
        new = [0] * (len(dp) + 2)
        for i, x in enumerate(dp):
            new[i + 1] += x; new[i + 2] += x
        Ps[k + 1] = new
    return Ps

# ------------------------------------------------------------------ Frobenius data and rows
def weights(a0, a1, a2, e1, e2):
    return {'a0': fr2p(a0), 'a1': fr2p(a1), 'a2': fr2p(a2), 'e1': fr2p(e1), 'e2': fr2p(e2)}
def cyw(a1, a2, e1, e2):
    a1, a2, e1, e2 = Fr(a1), Fr(a2), Fr(e1), Fr(e2)
    return weights(e1 + e2 - a1 - a2, a1, a2, e1, e2)

def exp_series(N):
    c = [1] * N
    for i in range(1, N): c[i] = c[i - 1] * inv(i) % P
    return (0, c)
def eval_poly(pc, x, N):
    """Horner, pc integer coefficients; x Laurent series"""
    r = None
    for c in reversed(pc):
        if r is None:
            r = s_const(c, N) if c % P else (10 ** 6, [])
            if c % P == 0: r = None
            continue
        r = s_mul(r, x)
        if c % P: r = s_add(r, s_scale(c, ('const', None)) if False else s_const(c, N * 3))
    return r

class Rows:
    """R-matrix rows (point frame R_Z incl. QRR constants) at the two critical points, as t-series"""
    def __init__(self, W, zmax, N=24, ffc=None, Qmode='t'):
        self.W = W; self.zmax = zmax; self.N = N
        a0, a1, a2, e1, e2 = W['a0'], W['a1'], W['a2'], W['e1'], W['e2']
        self.E_alpha = (a0 + a1 + a2) % P
        Q = exp_series(N + 4)
        # the two roots of Q(a0+a1)U^2 - (a0+Q(a1-a2))U - a2 = 0 as t-series
        A = s_scale((a0 + a1), Q)
        B = s_add(s_const(a0, N + 4), s_scale(a1 - a2, Q))
        C = s_const(a2, N + 4)
        def leq(U): return s_add(s_add(s_mul(A, s_mul(U, U)), s_neg(s_mul(B, U))), s_neg(C))
        def dleq(U): return s_add(s_scale(2, s_mul(A, U)), s_neg(B))
        roots = {}
        for c, u0 in ((1, 1), (2, (-a2 * inv(a0 + a1)) % P)):
            U = (0, [u0] + [0] * (N + 3))
            for _ in range(8):
                U = s_add(U, s_neg(s_mul(leq(U), s_inv(dleq(U)))))
            roots[c] = U
        self.U = roots
        self.Q = Q
        self.ffc = ffc_closed(W, zmax) if ffc is None else ffc
        Ps = thetapolys(2 * zmax + 3)
        self.L = {}; self.Delta = {}; self.wick = {}
        for c in (1, 2):
            U = roots[c]
            QU = s_mul(Q, U)
            X1 = s_mul(QU, s_inv(s_add(s_const(1, N + 4), s_neg(QU))))
            X2 = s_mul(U, s_inv(s_add(s_const(1, N + 4), s_neg(U))))
            # powers
            def Pk(k, X):
                r = (10 ** 7, [])
                pc = Ps[k]
                res = None
                for cf in reversed(pc):
                    if res is not None: res = s_mul(res, X)
                    if cf % P:
                        t = s_const(cf, N + 4)
                        res = t if res is None else s_add(res, t)
                return res
            Wk = {}
            a02 = (a0 + a2) % P
            for k in range(2, 2 * zmax + 3):
                Wk[k] = s_add(s_scale(-a1, Pk(k, X1)) if a1 else (10 ** 7, []), s_scale(-a02, Pk(k, X2)) if a02 else (10 ** 7, []))
            self_W2inv = s_inv(Wk[2])
            ipow = {0: (0, [1] + [0] * (N + 3))}
            for n in range(1, 3 * zmax + 1): ipow[n] = s_mul(ipow[n - 1], self_W2inv)
            wick = {0: (0, [1] + [0] * (N + 3))}
            for m in range(1, zmax + 1):
                tot = None
                for cf, n, pr in wick_terms(m):
                    t = ipow[n]
                    for k, mu in pr:
                        for _ in range(mu): t = s_mul(t, Wk[k])
                    t = s_scale(fr2p(cf), t)
                    tot = t if tot is None else s_add(tot, t)
                wick[m] = tot
            self.wick[c] = wick
            # eigenvalue of H*: L = -Q((a0+a1)U + a2)/(1-Q) = Q((a0+a1)U+a2)/(Q-1)
            Qm1 = s_add(Q, s_const(-1, N + 4))
            self.L[c] = s_mul(s_mul(Q, s_add(s_scale(a0 + a1, U), s_const(a2, N + 4))), s_inv(Qm1))
        # metric
        d = inv((a0 + a1) * (a0 + a2) % P * e1 % P * e2 % P)
        E = self.E_alpha
        self.eta = [[E * inv(a1 * a2 % P) % P * d % P, d], [d, (-a0) * d % P]]
        det = (self.eta[0][0] * self.eta[1][1] - self.eta[0][1] * self.eta[1][0]) % P
        di = inv(det)
        self.etainv = [[self.eta[1][1] * di % P, -self.eta[0][1] * di % P], [-self.eta[1][0] * di % P, self.eta[0][0] * di % P]]
        for c in (1, 2):
            o = 3 - c
            Li, Lj = self.L[c], self.L[o]
            diff = s_add(Lj, s_neg(Li))
            den = s_add(s_add(s_scale(self.eta[0][0], s_mul(Lj, Lj)), s_scale(-2 * self.eta[0][1], Lj)), s_const(self.eta[1][1], N + 4))
            self.Delta[c] = s_mul(s_mul(diff, diff), s_inv(den))
        # identity row with constants, then the H row
        self.r0 = {}; self.r1 = {}
        for c in (1, 2):
            r0 = []
            for m in range(zmax + 1):
                tot = None
                for k in range(m + 1):
                    if self.ffc[k] % P:
                        t = s_scale(self.ffc[k], self.wick[c][m - k])
                        tot = t if tot is None else s_add(tot, t)
                r0.append(tot)
            dlog = s_mul(s_der(self.Delta[c]), s_inv(s_scale(2, self.Delta[c])))
            r1 = [s_mul(self.L[c], r0[0])]
            for m in range(1, zmax + 1):
                r1.append(s_add(s_mul(self.L[c], r0[m]), s_add(s_der(r0[m - 1]), s_neg(s_mul(r0[m - 1], dlog)))))
            self.r0[c] = r0; self.r1[c] = r1

class Frame:
    """frame: per-colour class weights (mod p) or None.  rows r'_{a,c} = r_{a,c} h_{w_c}^{s}; edges, dilatons, Delta powers"""
    def __init__(self, rows, wcol, s=-1, gmax=4):
        self.rows = rows; self.wcol = wcol; self.zmax = rows.zmax
        zmax = rows.zmax
        self.vt = {c: VertexTable(wcol[c]) for c in (1, 2)}
        self.r = {}
        for c in (1, 2):
            hc = h_coeffs(wcol[c], zmax, s)
            rr = {}
            for a, src in ((0, rows.r0[c]), (1, rows.r1[c])):
                new = []
                for m in range(zmax + 1):
                    tot = None
                    for j in range(m + 1):
                        if hc[m - j] % P:
                            t = s_scale(hc[m - j], src[j])
                            tot = t if tot is None else s_add(tot, t)
                    new.append(tot)
                rr[a] = new
            self.r[c] = rr
        self.dil = {c: {n: s_scale((-1) ** n, self.r[c][0][n - 1]) for n in range(1, zmax + 2)} for c in (1, 2)}
        ei = rows.etainv
        self.ed = {}
        for i in (1, 2):
            for j in (1, 2):
                pre = {}
                for p in range(zmax + 1):
                    for q in range(zmax + 1 - p):
                        tot = None
                        for a in (0, 1):
                            for b in (0, 1):
                                if ei[a][b] % P:
                                    t = s_scale(-(-1) ** (p + q) * ei[a][b], s_mul(self.r[i][a][p], self.r[j][b][q]))
                                    tot = t if tot is None else s_add(tot, t)
                        pre[(p, q)] = tot
                tab = {}
                for m in range(zmax):
                    for n in range(zmax - m):
                        tot = None
                        for k in range(n + 1):
                            t = s_scale((-1) ** k, pre[(m + k + 1, n - k)])
                            tot = t if tot is None else s_add(tot, t)
                        tab[(m, n)] = tot
                self.ed[(i, j)] = tab
        self.dp = {}
        for c in (1, 2):
            D = rows.Delta[c]; Di = s_inv(D)
            self.dp[(c, 0)] = (0, [1] + [0] * (len(D[1]) - 1))
            self.dp[(c, -1)] = Di
            cur = self.dp[(c, 0)]
            for k in range(1, gmax):
                cur = s_mul(cur, D); self.dp[(c, k)] = cur

# ------------------------------------------------------------------ stable graphs, decorations, orbits
DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
def load_graphs(g):
    b = open(os.path.join(DATADIR, "StableGraph-%d.bin" % g), "rb").read()
    raw = [x for x in b[1:-1] if x != 13]
    out = []; i = 0
    while i < len(raw):
        K = raw[i]; gv = raw[i + 1:i + 1 + K]; m = raw[i + 1 + K:i + 1 + 2 * K]; l = raw[i + 1 + 2 * K:i + 1 + 3 * K]
        i2 = i + 1 + 3 * K
        A = [[0] * K for _ in range(K)]
        for j in range(K): A[j][j] = l[j]
        for j in range(K):
            for ii in range(j + 1, K):
                A[ii][j] = A[j][ii] = raw[i2]; i2 += 1
        out.append((tuple(gv), tuple(m), A))
        i = i2
    return out
def edges_of(A):
    K = len(A); ed = []
    for i in range(K):
        ed += [(i, i)] * A[i][i]
    for i in range(K):
        for j in range(i + 1, K):
            ed += [(i, j)] * A[i][j]
    return ed
def autos_of(gv, A):
    V = len(gv); ed = edges_of(A); out = []
    for sig in itertools.permutations(range(V)):
        if any(gv[sig[v]] != gv[v] for v in range(V)): continue
        if any(A[sig[i]][sig[j]] != A[i][j] for i in range(V) for j in range(V)): continue
        img = [tuple(sorted((sig[u], sig[v]))) for u, v in ed]
        classes = {}
        for e, im in enumerate(img): classes.setdefault(im, []).append(e)
        choices = []
        for im, cl in classes.items():
            targets = [e for e, x in enumerate(ed) if x == im]
            choices.append([dict(zip(cl, p)) for p in itertools.permutations(targets)])
        for ch in itertools.product(*choices):
            emap = {}
            for d in ch: emap.update(d)
            flipsets = [(0, 1) if ed[e][0] == ed[e][1] else (0,) for e in range(len(ed))]
            for fl in itertools.product(*flipsets):
                out.append((sig, emap, fl))
    return out
def aut_count(gv, A):
    V = len(gv); n = 0
    for sig in itertools.permutations(range(V)):
        if any(gv[sig[v]] != gv[v] for v in range(V)): continue
        if any(A[sig[i]][sig[j]] != A[i][j] for i in range(V) for j in range(V)): continue
        n += 1
    r = n
    for i in range(V):
        r *= 2 ** A[i][i] * math.factorial(A[i][i])
        for j in range(i + 1, V): r *= math.factorial(A[i][j])
    return r
def compositions(s, nd, minp):
    """ordered tuples of nd parts >= minp with sum s"""
    if nd == 0:
        if s == 0: yield ()
        return
    for f in range(minp, s - minp * (nd - 1) + 1):
        for r in compositions(s - f, nd - 1, minp):
            yield (f,) + r
def tuples_sum_le(n, budget):
    if n == 0:
        yield ()
        return
    for f in range(budget + 1):
        for r in tuples_sum_le(n - 1, budget - f):
            yield (f,) + r

_nonzero_memo = {}
def hodge_poly_nonzero(g, d):
    """is the Hodge correlator (exp class) a non-zero polynomial of the s_k?"""
    d = tuple(sorted(d)); key = (g, d)
    r = _nonzero_memo.get(key)
    if r is None:
        n = len(d); D = 3 * g - 3 + n - sum(d)
        if D < 0 or 2 * g - 2 + n <= 0: r = False
        elif D > dmax(g): r = False
        else:
            r = any(hi(g, d, K) != 0 for K in partitions_odd(D, max(2 * g - 1, 1)))
        _nonzero_memo[key] = r
    return r
def vertex_options(g, hs, hodge):
    """list of (hp tuple, dils tuple, 1/nd!) for a vertex of genus g with hs half-edges; hodge=False: point vertices"""
    res = []
    for nd in range(0, max(0, 3 * g - 3 + hs) + 1):
        dim = 3 * g - 3 + hs + nd
        if dim < 2 * nd: continue
        dilss = [()] if nd == 0 else [c for s in range(2 * nd, dim + 1) for c in compositions(s, nd, 2)]
        for dils in dilss:
            budget = dim - sum(dils)
            for hp in tuples_sum_le(hs, budget):
                if hodge:
                    if hodge_poly_nonzero(g, hp + dils): res.append((hp, dils))
                else:
                    if sum(hp) == budget and wk(g, hp + dils) != 0: res.append((hp, dils))
    return res

class Decorated:
    """orbits of labelled decorations for all stable graphs of genus g"""
    def __init__(self, g, hodge):
        self.g = g; self.hodge = hodge
        self.graphs = []
        for gi, (gv, m, A) in enumerate(load_graphs(g)):
            V = len(gv); ed = edges_of(A)
            H = [(e, s) for e in range(len(ed)) for s in (0, 1)]
            hv = [ed[e][s] for e, s in H]
            hsAt = [[i for i, x in enumerate(hv) if x == v] for v in range(V)]
            opts = [vertex_options(gv[v], len(hsAt[v]), hodge) for v in range(V)]
            auts = autos_of(gv, A)
            assert len(auts) == aut_count(gv, A)
            orbs = {}
            nlab = 0
            for dec in itertools.product(*opts):
                nlab += 1
                psi = [[0, 0] for _ in ed]
                for v in range(V):
                    for r, h in enumerate(hsAt[v]):
                        e, s = H[h]; psi[e][s] = dec[v][0][r]
                dil = [tuple(sorted(dec[v][1])) for v in range(V)]
                best = None
                for sig, emap, fl in auts:
                    psi2 = [[0, 0] for _ in ed]
                    for e in range(len(ed)):
                        e2 = emap[e]
                        tgt = (sig[ed[e][0]], sig[ed[e][1]])
                        ordr = (0, 1) if tgt == ed[e2] else (1, 0)
                        if ed[e2][0] == ed[e2][1] and fl[e] == 1: ordr = (ordr[1], ordr[0])
                        psi2[e2][ordr[0]] = psi[e][0]; psi2[e2][ordr[1]] = psi[e][1]
                    dil2 = [None] * V
                    for v in range(V): dil2[sig[v]] = dil[v]
                    key = (tuple(tuple(x) for x in psi2), tuple(dil2))
                    if best is None or key < best: best = key
                orbs.setdefault(best, []).append(dec)
            self.graphs.append(dict(idx=gi + 1, gv=gv, A=A, ed=ed, H=H, hsAt=hsAt, aut=len(auts), orbs=orbs, nlab=nlab, V=V))
    def counts(self):
        return dict(graphs=len(self.graphs), labelled=sum(G['nlab'] for G in self.graphs), orbits=sum(len(G['orbs']) for G in self.graphs))

def eval_dec(G, dec, col, frame, T):
    """value at the colouring col of a labelled decoration, up to and including t^T; returns series or None if O(t^(T+1))"""
    V = G['V']; ed = G['ed']; H = G['H']; hsAt = G['hsAt']
    scal = 1
    facs = []
    hp = [0] * len(H)
    for v in range(V):
        c = col[v]; hps, dils = dec[v]
        vt = frame.vt[c]
        val = vt.val(G['gv'][v], hps + dils)
        if val == 0: return None
        scal = scal * val % P
        from collections import Counter
        for _, mu in Counter(dils).items():
            pass
        scal = scal * inv(math.factorial(len(dils))) % P
        for r, h in enumerate(hsAt[v]): hp[h] = hps[r]
        facs.append(frame.dp[(c, G['gv'][v] - 1)])
        for d in dils: facs.append(frame.dil[c][d])
    for e in range(len(ed)):
        i, j = col[ed[e][0]], col[ed[e][1]]
        facs.append(frame.ed[(i, j)][(hp[2 * e], hp[2 * e + 1])])
    vt_ = 0
    for f in facs:
        if not f[1]:
            # known zero up to its precision
            vt_ += f[0]
        else:
            vt_ += f[0]
    L = T - vt_ + 1
    if L <= 0: return None
    # require enough precision
    res = None
    for f in facs:
        if len(f[1]) < L:
            if not f[1]:
                if f[0] > T + 40: return None
            raise ValueError("insufficient precision: factor length %d, need %d" % (len(f[1]), L))
        g = (f[0], f[1][:L])
        res = g if res is None else s_mul(res, g, L)
    if res is None: res = (0, [1])
    return s_scale(scal, res)

class GraphSum:
    def __init__(self, deco, frame, T):
        self.deco = deco; self.frame = frame; self.T = T
        g = deco.g
        Kmin = -(2 * g + 2)  # lowest power we record; a guard below checks that nothing is lower
        self.low = -(6 * g)
        self.results = []   # list of dict(graph idx, key, coeffs: dict power->value, percol: dict col-> dict power->value)
    def run(self, percol=False, verbose=False):
        T = self.T
        for G in self.deco.graphs:
            cols = list(itertools.product((1, 2), repeat=G['V']))
            ai = inv(G['aut'])
            for key, decs in sorted(G['orbs'].items()):
                total = {}; pc = {}
                # value is Aut-invariant: use the representative and multiplicity
                rep = decs[0]
                mult = len(decs) % P * ai % P
                for col in cols:
                    x = eval_dec(G, rep, col, self.frame, T)
                    if x is None: continue
                    v, c = x
                    d = {}
                    for i, y in enumerate(c):
                        if y:
                            k = v + i
                            if k <= T:
                                y2 = y * mult % P
                                total[k] = (total.get(k, 0) + y2) % P
                                d[k] = y2
                    if percol and d: pc[col] = d
                total = {k: y for k, y in total.items() if y}
                self.results.append(dict(graph=G['idx'], gv=G['gv'], ed=G['ed'], key=key, c=total, pc=pc, V=G['V'], aut=G['aut'], bound=sum(max(2 * h - 2, 0) for h in G['gv'])))
        return self
    def total(self):
        tot = {}
        for r in self.results:
            for k, y in r['c'].items():
                tot[k] = (tot.get(k, 0) + y) % P
        return {k: y for k, y in tot.items() if y}
    def by_graph(self):
        out = {}
        for r in self.results:
            d = out.setdefault(r['graph'], dict(gv=r['gv'], ed=r['ed'], c={}, nsing=0, n=0, maxpole=0))
            d['n'] += 1
            neg = [k for k, y in r['c'].items() if k < 0]
            if neg:
                d['nsing'] += 1; d['maxpole'] = max(d['maxpole'], -min(neg))
            for k, y in r['c'].items():
                d['c'][k] = (d['c'].get(k, 0) + y) % P
        for d in out.values():
            d['c'] = {k: y for k, y in d['c'].items() if y}
        return out

# ------------------------------------------------------------------ the expected total
def P_g_poly(g, e1, e2):
    """t^(2g) coefficient of e1 e2 t^2/(4 sinh(e1 t/2) sinh(e2 t/2)); e1,e2 mod p"""
    # (x/2)/sinh(x/2) = sum_i c_i x^(2i), c_i = -(1-2^(1-2i)) B_2i/(2i)! (standard)  (i=0: 1)
    def c(i):
        if i == 0: return Fr(1)
        return -(1 - Fr(2) ** (1 - 2 * i)) * bernoulli(2 * i) / math.factorial(2 * i)
    tot = 0
    for i in range(g + 1):
        j = g - i
        tot += fr2p(c(i) * c(j)) * pow(e1, 2 * i, P) % P * pow(e2, 2 * j, P)
    return tot % P
def zeta_neg(n):
    """zeta(-n), n >= 0"""
    if n == 0: return Fr(-1, 2)
    return -bernoulli(n + 1) / (n + 1) * (1 if n >= 0 else 0) if True else None
def li_neg_exp(n, kmax):
    """Laurent coefficients in t of Li_{-n}(e^t), n >= 1: dict power -> Fraction, to t^kmax"""
    out = {-(n + 1): Fr(math.factorial(n) * (-1) ** (n + 1))}
    for k in range(0, kmax + 1):
        # zeta(-n-k) t^k/k!
        out[k] = zeta_neg_correct(n + k) / math.factorial(k)
    return out
def zeta_neg_correct(n):
    if n == 0: return Fr(-1, 2)
    # zeta(-n) = (-1)^n B_{n+1}^+/(n+1) = -B_{n+1}/(n+1) for n >= 1
    return -bernoulli(n + 1) / (n + 1)
def expected_total(g, W, kmax):
    e1, e2 = W['e1'], W['e2']
    Pg = P_g_poly(g, e1, e2) * inv(e1 * e2 % P) % P
    li = li_neg_exp(2 * g - 3, kmax)
    return {k: fr2p(v) * Pg % P for k, v in li.items()}
def const_map(g, W):
    """F_g(0) from the two fixed points: Sum_p (prod w)^(g-1) Int_{Mbar_g} prod_w Lambda^vee(w)/w^g, five weights"""
    a0, a1, a2, e1, e2 = W['a0'], W['a1'], W['a2'], W['e1'], W['e2']
    tot = 0
    for wts in ([a0, a1, a2, -e1 % P, -e2 % P], [-a0 % P, (a0 + a1) % P, (a0 + a2) % P, -e1 % P, -e2 % P]):
        pr = 1
        for x in wts: pr = pr * x % P
        vt = VertexTable(wts)
        tot += pow(pr, g - 1, P) * vt.val(g, ())
    return tot % P

def fmt(y):
    r = ratrec(y)
    return str(r) if r is not None else "F_p:%d" % y

def const_closed_fibrewise(g, W):
    e1, e2 = W['e1'], W['e2']
    c = abs(bernoulli(2 * g) * bernoulli(2 * g - 2)) / (4 * g * (2 * g - 2) * math.factorial(2 * g - 2))
    return (-fr2p(c) * (pow(e1, 2 * g - 2, P) + pow(e2, 2 * g - 2, P))) % P
_const_map0 = const_map
def const_map(g, W):
    return const_closed_fibrewise(g, W) if W['a0'] == 0 else _const_map0(g, W)

_const_map1 = const_map
def const_map(g, W):
    if W['a0'] == 0 and not (W['a1'] in (W['e1'], W['e2']) and W['a2'] in (W['e1'], W['e2'])):
        return None
    return _const_map1(g, W)
