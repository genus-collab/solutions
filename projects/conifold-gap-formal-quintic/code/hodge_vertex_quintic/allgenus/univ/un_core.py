# -*- coding: utf-8 -*-
"""un_core.py -- series arithmetic with a formal first-order (or higher) parameter eps, hatted rows of the formal quintic as plain
lists, and the edge / leaf tables of the graph sum computed from ANY table of hatted rows.

Everything here is read-only with respect to the project: the project modules are imported (sys.path), nothing is written there.

Series.  A series in (u, eps) truncated at u^KU and eps^L is a list of KU*L residues mod p, entry [e*KU + k] = coefficient of
eps^e u^k; None stands for zero (the convention of fq4_split.Split).  L = 1 is the plain u-series of fq4_split.

Rows.  Qh[a][i][n] = [zhat^n] Qhat_a^(i)  (a flat index 0..4, i colour 0..4, n = 0..KZ), the hatted rows of fq3_msum.HatRowsU:
        Qhat_a^(i)(zhat) = eta(H^a, R(zhat) e~_i) s^(1 - c_a),   c = (0,-1,3,2,1),  zhat = z/s^4,  u = s^5.
The whole graph sum is a function of these rows: edges (EdgeLeaf formula of fq_sser), dilaton leaves T(d)_i = (-1)^d Qh[0][i][d-1],
vertex factor (zeta^(4i) xhat)^(2g-2+n) with 1/(zeta^(4i) xhat) = Qh[0][i][0].

API: Arith(p, KU, L);  nat_rows(prime, g, ab, K) -> dict;  tables(ar, Qh, KZ, eta_inv) -> (ER, TR);  xhat_of(ar, Qh, R);
     left_dress(ar, Qh, X) (rows -> (1 + X) rows, X = {(a, b, k): series});  right_scalar(ar, Qh, c) (rows -> rows c(zhat)).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2

CEXP = (0, -1, 3, 2, 1)
PRIMES = {"89": P2, "61": P1}


def mu(b, k):
    """u-exponent of the flat-frame term z^k E_{b+k,b} in hatted variables: zhat^k u^mu E_{b+k,b}"""
    v = 4 * k + CEXP[b % 5] - CEXP[(b + k) % 5]
    assert v % 5 == 0
    return v // 5


class Arith(object):
    """arithmetic of (u, eps)-series as lists of KU*L residues (None = 0)"""
    def __init__(self, p, KU, L=1):
        self.p = p; self.KU = KU; self.L = L; self.N = KU * L

    def zero(self): return None
    def one(self): return [1] + [0] * (self.N - 1)

    def const(self, c):
        c %= self.p
        return ([c] + [0] * (self.N - 1)) if c else None

    def mono(self, c, k, e=0):
        """c u^k eps^e"""
        c %= self.p
        if not c or k >= self.KU or e >= self.L: return None
        out = [0] * self.N; out[e * self.KU + k] = c
        return out

    def from_u(self, coefs, e=0):
        """u-series (list of residues, any length) placed at eps^e"""
        if e >= self.L: return None
        out = [0] * self.N; p = self.p
        for k, c in enumerate(coefs[:self.KU]): out[e * self.KU + k] = c % p
        return out if any(out) else None

    def add(self, a, b):
        if a is None: return b
        if b is None: return a
        p = self.p; out = [(x + y) % p for x, y in zip(a, b)]
        return out if any(out) else None

    def neg(self, a):
        if a is None: return None
        p = self.p; return [(-x) % p for x in a]

    def sub(self, a, b): return self.add(a, self.neg(b))

    def scale(self, a, c):
        c %= self.p
        if a is None or not c: return None
        p = self.p; return [x * c % p for x in a]

    def mul(self, a, b):
        if a is None or b is None: return None
        K = self.KU; L = self.L; p = self.p; out = [0] * self.N
        for ea in range(L):
            oa = ea * K
            if not any(a[oa:oa + K]): continue
            for eb in range(L - ea):
                ob = eb * K; oo = oa + ob
                for i in range(K):
                    x = a[oa + i]
                    if x:
                        for j in range(K - i):
                            y = b[ob + j]
                            if y: out[oo + i + j] += x * y
        out = [c % p for c in out]
        return out if any(out) else None

    def inv(self, a):
        """inverse of a series with invertible constant term (eps^0 u^0)"""
        K = self.KU; p = self.p
        a0 = a[:K]; assert a0[0] % p, "series not invertible"
        y0 = [0] * K; y0[0] = pow(a0[0], p - 2, p)
        for k in range(1, K):
            y0[k] = (-sum(a0[j] * y0[k - j] for j in range(1, k + 1)) * y0[0]) % p
        y = y0 + [0] * (self.N - K); two = self.const(2); n = 1
        while n < self.L:
            y = self.mul(y, self.sub(two, self.mul(a, y))); n *= 2
        return y

    def part(self, a, e):
        """the u-coefficients of eps^e as a list of KU residues"""
        if a is None: return [0] * self.KU
        return list(a[e * self.KU:(e + 1) * self.KU])

    def val(self, a, e=0):
        """u-valuation of the eps^e part (KU if it vanishes)"""
        if a is None: return self.KU
        for k in range(self.KU):
            if a[e * self.KU + k]: return k
        return self.KU


def nat_rows(prime, g, ab=("0", "0"), K=None):
    """the hatted rows of the formal quintic in the conifold frame ab as lists (the first lines of fq4_split.Split.__init__, repeated
    here because Split does not keep the row object).  Returns dict(R, p, Qh, x, eta_inv, Rh, C, KZ, K)."""
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq_symbol import Symbol
    from fq3_msum import HatRowsU
    p = PRIMES[prime]; KZ = 3 * g - 3; K = (2 * g - 2) if K is None else K
    R = FpRing(p)
    rows = Rows(R, Nz=KZ); C = Conifold(R, Nu=K + 4, alpha=Fr(ab[0]), beta=Fr(ab[1])); S = SRows(rows, C, prec_abs=5 * K + 4)
    sym = Symbol(S); assert sym.val_ok
    H = HatRowsU(S, K); assert H.val_ok and H.min_prec >= K
    def lst(s):
        assert s.prec >= K
        return [s.coef(k) % p for k in range(K)]
    Qh = [[[lst(H.Q[a][i][n]) for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    return dict(R=R, p=p, Qh=Qh, x=lst(H.x), eta_inv=dict(H.eta_inv), Rh=sym.Rh, C=C, KZ=KZ, K=K)


def lift_rows(ar, Qh):
    """plain u-lists (length >= KU) -> (u, eps)-series at eps^0"""
    return [[[ar.from_u(q) if q is not None else None for q in col] for col in row] for row in Qh]


def tables(ar, Qh, KZ, eta_inv):
    """(ER, TR): ER[(m, n)][i][j] for m + n + 1 <= KZ, TR[d][i] for d = 2..KZ+1  (fq_sser.EdgeLeaf, on lists)
    G(M, N)[i][j] = delta_ij delta_M0 delta_N0 - (-1)^(M+N) Sum_ab eta^ab Q_a^(i)[M] Q_b^(j)[N];  E(m, n) = Sum_k (-1)^k G(m+1+k, n-k)"""
    p = ar.p; G = {}
    pairs = [(a, b, e % p) for (a, b), e in eta_inv.items()]
    for M_ in range(KZ + 1):
        for N_ in range(KZ + 1 - M_):
            tab = [[None] * 5 for _ in range(5)]
            sg = p - 1 if (M_ + N_) % 2 == 0 else 1
            for i in range(5):
                for j in range(5):
                    acc = None
                    for a, b, e in pairs:
                        acc = ar.add(acc, ar.scale(ar.mul(Qh[a][i][M_], Qh[b][j][N_]), e))
                    acc = ar.scale(acc, sg)
                    if i == j and M_ == 0 and N_ == 0: acc = ar.add(acc, ar.one())
                    tab[i][j] = acc
            G[(M_, N_)] = tab
    ER = {}
    for m in range(KZ):
        for n in range(KZ - m):
            tab = [[None] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    acc = None
                    for k in range(n + 1):
                        t = G[(m + 1 + k, n - k)][i][j]
                        acc = ar.add(acc, ar.neg(t) if k % 2 else t)
                    tab[i][j] = acc
            ER[(m, n)] = tab
    TR = {d: [ar.scale(Qh[0][i][d - 1], 1 if d % 2 == 0 else p - 1) for i in range(5)] for d in range(2, KZ + 2)}
    return ER, TR, G


def symplectic_check(ar, Qh, KZ, eta_inv):
    """failures of Sum_ab Sum_m Q_a^(i)[m] eta^ab (-1)^(n-m) Q_b^(j)[n-m] = delta_ij delta_n0 for n <= KZ"""
    p = ar.p; bad = []
    pairs = [(a, b, e % p) for (a, b), e in eta_inv.items()]
    for n in range(KZ + 1):
        for i in range(5):
            for j in range(5):
                acc = None
                for a, b, e in pairs:
                    for m in range(n + 1):
                        t = ar.scale(ar.mul(Qh[a][i][m], Qh[b][j][n - m]), e if (n - m) % 2 == 0 else p - e)
                        acc = ar.add(acc, t)
                if i == j and n == 0: acc = ar.sub(acc, ar.one())
                if acc is not None: bad.append((n, i, j))
    return bad


def xhat_of(ar, Qh, R):
    """xhat with Qh[0][i][0] = 1/(zeta^(4i) xhat); asserts the Z/5 covariance of the unit row at zhat^0"""
    p = ar.p; x = ar.inv(Qh[0][0][0])
    for i in range(5):
        chk = ar.scale(ar.mul(Qh[0][i][0], x), R.z(4 * i))
        assert chk == ar.one(), "unit row at zhat^0 is not zeta^(-4i)/xhat (colour %d)" % i
    return x


def left_dress(ar, Qh, X, KZ):
    """rows -> (1 + X) rows:  Qh'[a][i][n] = Qh[a][i][n] + Sum_{b, k} X[(a, b, k)] Qh[b][i][n - k]"""
    out = [[[Qh[a][i][n] for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
    for (a, b, k), ser in X.items():
        if ser is None: continue
        for i in range(5):
            for n in range(k, KZ + 1):
                out[a][i][n] = ar.add(out[a][i][n], ar.mul(ser, Qh[b][i][n - k]))
    return out


def right_scalar(ar, Qh, c, KZ):
    """rows -> rows c(zhat), c = {k: series} (the term k = 0 must be given explicitly if wanted)"""
    out = [[[None] * (KZ + 1) for i in range(5)] for a in range(5)]
    for a in range(5):
        for i in range(5):
            for n in range(KZ + 1):
                acc = None
                for k, ser in c.items():
                    if k <= n: acc = ar.add(acc, ar.mul(ser, Qh[a][i][n - k]))
                out[a][i][n] = acc
    return out


def mat_exp_left(ar, X, KZ, order):
    """exp(X) - 1 for X = {(a, b, k): series} with every k >= 1 (or nilpotent enough): truncated at zhat^KZ and at 'order' terms"""
    def mmul(A, B):
        out = {}
        for (a, b, k), s in A.items():
            for (b2, c, k2), t in B.items():
                if b2 != b or k + k2 > KZ: continue
                key = (a, c, k + k2); out[key] = ar.add(out.get(key), ar.mul(s, t))
        return {key: v for key, v in out.items() if v is not None}
    total = dict(X); power = dict(X); fact = 1
    for m in range(2, order + 1):
        power = mmul(power, X); fact *= m
        if not power: break
        inv = pow(fact, ar.p - 2, ar.p)
        for key, v in power.items(): total[key] = ar.add(total.get(key), ar.scale(v, inv))
    return {key: v for key, v in total.items() if v is not None}
