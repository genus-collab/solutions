# Graph sum G_g (Hodge vertices) and Tbar_g (point vertices) from the closed-form rows.
# Coefficients: F_p (p=2^127-1, omega in F_p) or exact Q(omega). Truncated Laurent series with tracked precision.
from fractions import Fraction as Fr
from math import factorial, comb
from core import *
import itertools, sys, time

P = 2 ** 127 - 1

class FpField:
    name = 'Fp'
    def __init__(self, p=P):
        self.p = p
        g = 2
        while True:
            w = pow(g, (p - 1) // 3, p)
            if w != 1:
                break
            g += 1
        self.omega = w
        self.zero = 0
        self.one = 1
    def fr(self, x):
        x = Fr(x)
        return x.numerator % self.p * pow(x.denominator, -1, self.p) % self.p
    def red(self, s):
        return s % self.p
    def inv(self, a):
        return pow(a, -1, self.p)
    def is0(self, a):
        return a % self.p == 0

class QW:
    __slots__ = ('a', 'b')
    def __init__(self, a, b=0):
        self.a = a
        self.b = b
    def __add__(self, o):
        if isinstance(o, QW):
            return QW(self.a + o.a, self.b + o.b)
        return QW(self.a + o, self.b)
    __radd__ = __add__
    def __neg__(self):
        return QW(-self.a, -self.b)
    def __sub__(self, o):
        if isinstance(o, QW):
            return QW(self.a - o.a, self.b - o.b)
        return QW(self.a - o, self.b)
    def __mul__(self, o):
        if isinstance(o, QW):
            # w^2 = -1-w
            ac = self.a * o.a
            bd = self.b * o.b
            return QW(ac - bd, self.a * o.b + self.b * o.a - bd)
        return QW(self.a * o, self.b * o)
    __rmul__ = __mul__

class QWField:
    name = 'Qw'
    def __init__(self):
        self.omega = QW(Fr(0), Fr(1))
        self.zero = QW(Fr(0), Fr(0))
        self.one = QW(Fr(1), Fr(0))
    def fr(self, x):
        return QW(Fr(x), Fr(0))
    def red(self, s):
        return s
    def inv(self, x):
        # norm a^2-ab+b^2; conj(a+bw)=a+b w^2 = a-b -b w
        n = x.a * x.a - x.a * x.b + x.b * x.b
        return QW((x.a - x.b) / n, -x.b / n)
    def is0(self, x):
        return x.a == 0 and x.b == 0

class L:
    """Truncated Laurent series over a field: sum_{i<len(c)} c[i] x^(val+i); exact up to x^(val+len(c))"""
    __slots__ = ('F', 'val', 'c')
    def __init__(self, F, val, c):
        self.F = F
        self.val = val
        self.c = c
        self._norm()
    def _norm(self):
        c = self.c
        F = self.F
        i = 0
        n = len(c)
        while i < n and F.is0(c[i]):
            i += 1
        if i:
            self.c = c[i:]
            self.val += i
    @property
    def prec(self):
        return self.val + len(self.c)
    def __add__(self, o):
        F = self.F
        if not isinstance(o, L):
            o = L.const(F, o, max(self.prec, 1))
        v = min(self.val, o.val)
        pr = min(self.prec, o.prec)
        n = pr - v
        if n <= 0:
            return L(F, pr, [])
        r = [F.zero] * n
        for i, x in enumerate(self.c):
            k = self.val + i - v
            if k < n:
                r[k] = x
        for i, x in enumerate(o.c):
            k = o.val + i - v
            if k < n:
                r[k] = F.red(r[k] + x)
        return L(F, v, r)
    def __neg__(self):
        F = self.F
        return L(F, self.val, [F.red(-x) for x in self.c])
    def __sub__(self, o):
        if not isinstance(o, L):
            o = L.const(self.F, o, max(self.prec, 1))
        return self + (-o)
    def __mul__(self, o):
        F = self.F
        if not isinstance(o, L):
            return L(F, self.val, [F.red(x * o) for x in self.c])
        a = self.c
        b = o.c
        n = min(len(a), len(b))
        v = self.val + o.val
        if n == 0:
            # absolute precision
            if len(a) == 0:
                return L(F, self.val + o.val, [])
            return L(F, self.val + o.val, [])
        red = F.red
        z = F.zero
        r = [z] * n
        for k in range(n):
            s = z
            for i in range(k + 1):
                s = s + a[i] * b[k - i]
            r[k] = red(s)
        return L(F, v, r)
    def inv(self):
        F = self.F
        a = self.c
        n = len(a)
        if n == 0:
            raise ZeroDivisionError('inverse of series with no significant terms')
        i0 = F.inv(a[0])
        r = [F.zero] * n
        r[0] = i0
        for k in range(1, n):
            s = F.zero
            for i in range(1, k + 1):
                s = s + a[i] * r[k - i]
            r[k] = F.red(F.zero - s * i0)
        return L(F, -self.val, r)
    def deriv(self):
        F = self.F
        return L(F, self.val - 1, [F.red(x * F.fr(self.val + i)) for i, x in enumerate(self.c)])
    @staticmethod
    def const(F, x, prec):
        if prec <= 0:
            return L(F, prec, [])
        return L(F, 0, [x] + [F.zero] * (prec - 1))
    def coeff(self, n):
        if n < self.val:
            return self.F.zero
        if n >= self.prec:
            raise ValueError('beyond precision %d >= %d' % (n, self.prec))
        return self.c[n - self.val]

def rec_r(n):
    r = [Fr(1)]
    prev = Fr(0)
    for m in range(n):
        nxt = (3 * (3 * m * m + 3 * m + 1) * r[m] - m ** 3 * (r[m - 1] if m >= 1 else 0)) / (27 * (m + 1))
        r.append(nxt)
    return r

def S_polys(E, K):
    """[z^k] S(v,z) as dict {power:Fraction}, k<=K."""
    r = rec_r(K + 1)
    out = [dict() for _ in range(K + 1)]
    for l in range(K + 1):
        # 1/prod_{j=1}^{l+1} (1+jz/E) series
        ser = [Fr(0)] * (K + 1 - l)
        ser[0] = Fr(1)
        for j in range(1, l + 2):
            a = Fr(j) / E
            for i in range(1, len(ser)):
                ser[i] = ser[i] - a * ser[i - 1]
        for i, x in enumerate(ser):
            if x != 0:
                out[l + i][l] = out[l + i].get(l, Fr(0)) + r[l] * x
    return out

def theta_poly(d, E):
    # theta v = 1/E - v/3
    o = {}
    for n, c in d.items():
        if n != 0:
            o[n - 1] = o.get(n - 1, Fr(0)) + c * n / E
            o[n] = o.get(n, Fr(0)) - c * Fr(n, 3)
    return o

def add_to(o, d, f=Fr(1), shift=0):
    for n, c in d.items():
        o[n + shift] = o.get(n + shift, Fr(0)) + c * f

def T_apply(seq, E, K):
    """T = z theta + (1/v - E/3) + z (1/(E v) - 1/3) acting on z-series of Laurent polys in v."""
    out = [dict() for _ in range(K + 1)]
    for k in range(K + 1):
        add_to(out[k], seq[k], Fr(1), -1)
        add_to(out[k], seq[k], -E / 3, 0)
        if k >= 1:
            add_to(out[k], theta_poly(seq[k - 1], E), Fr(1), 0)
            add_to(out[k], seq[k - 1], Fr(1) / E, -1)
            add_to(out[k], seq[k - 1], Fr(-1, 3), 0)
    return out

def hw_series(e1, e2, K):
    E = e1 + e2
    logh = [Fr(0)] * (K + 1)
    for m in range(1, K // 2 + 2):
        if 2 * m - 1 <= K:
            logh[2 * m - 1] = bern(2 * m) / (2 * m * (2 * m - 1)) * (-E ** (1 - 2 * m) + e1 ** (1 - 2 * m) + e2 ** (1 - 2 * m))
    h = [Fr(0)] * (K + 1)
    h[0] = Fr(1)
    for n in range(1, K + 1):
        h[n] = sum(k * logh[k] * h[n - k] for k in range(1, n + 1)) / n
    return h

class Data:
    """Rows at the three colours. expansions: v[i] (L), Dc (L, Delta_c), c1 (L), FF (L)."""
    def __init__(self, F, e1, e2, v, Dc, c1, FF, K, frame='hodge'):
        self.F = F
        self.e1, self.e2 = e1, e2
        self.E = e1 + e2
        E = self.E
        self.K = K
        self.frame = frame
        Sp = S_polys(E, K)
        Qp = T_apply(Sp, E, K)
        TQp = T_apply(Qp, E, K)
        c1i = c1.inv()
        d2 = (c1i * c1i - F.one) * Dc.inv() * F.fr(E / 2)
        eta00 = e1 * e2 * E ** 3 / 9
        eta11 = -3 * e1 * e2 * E
        eta12 = -3 * e1 * e2
        self.eta = {(0, 0): F.fr(eta00), (1, 1): F.fr(eta11), (1, 2): F.fr(eta12), (2, 1): F.fr(eta12)}
        self.ncol = len(v)
        self.Delta = []
        self.P = []
        kap = F.fr(3 * e1 * e2 * E)
        h = hw_series(e1, e2, K) if frame == 'point' else None
        for i, vi in enumerate(v):
            vinv = vi.inv()
            pw = {0: L.const(F, F.one, len(vi.c) + 5)}
            def power(n):
                if n in pw:
                    return pw[n]
                if n > 0:
                    pw[n] = power(n - 1) * vi
                else:
                    pw[n] = power(n + 1) * vinv
                return pw[n]
            def ev(d):
                s = None
                for n, c in d.items():
                    if c == 0:
                        continue
                    t = power(n) * F.fr(c)
                    s = t if s is None else s + t
                if s is None:
                    return L(F, 10 ** 6, [])
                return s
            self.Delta.append(vinv * vinv * kap)
            S = [ev(Sp[k]) for k in range(K + 1)]
            Q = [ev(Qp[k]) for k in range(K + 1)]
            TQ = [ev(TQp[k]) for k in range(K + 1)]
            P0 = S
            P1 = [Q[k] * c1i for k in range(K + 1)]
            P2 = []
            Dcc1 = Dc * c1
            for k in range(K + 1):
                br = TQ[k] - Q[k] * d2
                if k >= 1:
                    br = br - FF * Q[k - 1]
                P2.append(br * Dcc1)
            rows = [P0, P1, P2]
            if frame == 'point':
                newrows = []
                for row in rows:
                    nr = []
                    for k in range(K + 1):
                        s = row[k]
                        for j in range(1, k + 1):
                            if h[j] != 0:
                                s = s + row[k - j] * F.fr(h[j])
                        nr.append(s)
                    newrows.append(nr)
                rows = newrows
            self.P.append(rows)
        self.ecache = {}
        self.dcache = {}
        self.vcache = {}
        self.pcache = {}

    def pair(self, i, j, M, b):
        key = (i, j, M, b)
        r = self.pcache.get(key)
        if r is None:
            A = self.P[i]
            B = self.P[j]
            et = self.eta
            r = A[0][M] * B[0][b] * et[(0, 0)] + A[1][M] * B[1][b] * et[(1, 1)] + (A[1][M] * B[2][b] + A[2][M] * B[1][b]) * et[(1, 2)]
            self.pcache[key] = r
        return r

    def edge(self, i, j, a, b):
        """Ehat_ij(a,b) = Sum_k (-1)^k N_{a+k+1,b-k}, N_{M,c} = -(-1)^(M+c) <P_M,P_c>"""
        key = (i, j, a, b)
        r = self.ecache.get(key)
        if r is None:
            s = None
            for k in range(b + 1):
                M = a + k + 1
                c = b - k
                sign = (-1) ** k * (-(-1) ** (M + c))
                t = self.pair(i, j, M, c) * self.F.fr(sign)
                s = t if s is None else s + t
            r = s
            self.ecache[key] = r
        return r

    def dil(self, i, d):
        key = (i, d)
        r = self.dcache.get(key)
        if r is None:
            r = self.P[i][0][d - 1] * self.F.fr((-1) ** d)
            self.dcache[key] = r
        return r

    def vint(self, h, psi):
        if self.frame == 'hodge':
            return omega_num(h, tuple(sorted(psi)), self.e1, self.e2)
        return wk(h, tuple(sorted(psi))) if 2 * h - 2 + len(psi) > 0 else Fr(0)

    def vertex(self, h, i, aa):
        """Delta_i^(h-1) Sum_m (1/m!) Sum_d Int ... prod Dil."""
        aa = tuple(sorted(aa))
        key = (h, i, aa)
        r = self.vcache.get(key)
        if r is not None:
            return r
        F = self.F
        n = len(aa)
        dim = 3 * h - 3 + n
        budget = dim - sum(aa)
        s = None
        for ds in excess_multisets(budget):
            if n == 0 and len(ds) == 0 and h < 2:
                continue
            I = self.vint(h, aa + ds)
            if I == 0:
                continue
            mult = {}
            for d in ds:
                mult[d] = mult.get(d, 0) + 1
            w = I
            for d, c in mult.items():
                w = w / factorial(c)
            t = None
            for d in ds:
                x = self.dil(i, d)
                t = x if t is None else t * x
            if t is None:
                t = L.const(F, F.fr(w), 10 ** 3)
            else:
                t = t * F.fr(w)
            s = t if s is None else s + t
        if s is None:
            r = 0
        else:
            D = self.Delta[i]
            if h == 0:
                r = s * D.inv()
            elif h == 1:
                r = s
            else:
                DD = D
                for _ in range(h - 2):
                    DD = DD * D
                r = s * DD
        self.vcache[key] = r
        return r

    def T(self, gs, m, col):
        """decoration sum T(Gamma,c)"""
        V = len(gs)
        edges = []
        for (i, j), c in m.items():
            edges += [(i, j)] * c
        val = [0] * V
        for (i, j) in edges:
            val[i] += 1
            val[j] += 1
        budgets = [3 * gs[i] - 3 + val[i] for i in range(V)]
        ex = [[] for _ in range(V)]
        total = [None]
        F = self.F
        def rec(e, pr):
            if e == len(edges):
                t = pr
                for v in range(V):
                    x = self.vertex(gs[v], col[v], tuple(ex[v]))
                    if x == 0:
                        return
                    t = x if t is None else t * x
                total[0] = t if total[0] is None else total[0] + t
                return
            (i, j) = edges[e]
            bi = budgets[i] - sum(ex[i])
            bj = budgets[j] - sum(ex[j])
            for a in range(bi + 1):
                ex[i].append(a)
                bj2 = bj - a if i == j else bj
                for b in range(bj2 + 1):
                    ex[j].append(b)
                    x = self.edge(col[i], col[j], a, b)
                    rec(e + 1, x if pr is None else pr * x)
                    ex[j].pop()
                ex[i].pop()
        rec(0, None)
        return total[0]

    def total(self, g, per=None):
        tot = None
        for gs, m, A in stable_graphs_cached(g):
            for col in itertools.product(range(self.ncol), repeat=len(gs)):
                t = self.T(gs, m, col)
                if t is None:
                    continue
                t = t * self.F.fr(Fr(1, A))
                if per is not None:
                    per.append((gs, tuple(sorted(m.items())), col, t))
                tot = t if tot is None else tot + t
        return tot

_sg = {}
def stable_graphs_cached(g):
    if g not in _sg:
        _sg[g] = stable_graphs(g)
    return _sg[g]

# ---------- series in Q
def conifold_setup(F, N):
    """Series in xi over F: v_i, Dc, c1=C1cf=-wc, FF=F_cf, t_c."""
    # wc
    w = []
    x = Fr(1)
    for n in range(N + 3):
        w.append(x)
        x = x * (Fr(1, 3) + n) * (Fr(2, 3) + n) / (n + 1) ** 2
    # (1-xi)^{1/3}
    bin_ = []
    x = Fr(1)
    for n in range(N + 3):
        bin_.append(x)
        x = x * (Fr(1, 3) - n) / (n + 1) * (-1)
    t = [-b / 3 for b in bin_]
    # t_c' = wc/(1-xi)
    pc = []
    s = Fr(0)
    for n in range(N + 3):
        s += w[n]
        pc.append(s)
    tc = [Fr(0)] + [pc[n] / (n + 1) for n in range(N + 2)]
    fs = lambda a: L(F, 0, [F.fr(x) for x in a[:N + 2]])
    W = fs(w)
    xi = L(F, 1, [F.one] + [F.zero] * (N + 1))
    theta_W = (xi - F.one) * W.deriv()
    FF = theta_W * W.inv()
    c1 = -W
    return dict(t=fs(t), tc=L(F, 0, [F.fr(x) for x in tc[:N + 2]]), xi=xi, c1=c1, FF=FF, W=W)

def conifold_data(F, e1, e2, N, K, frame='hodge', dF=None):
    su = conifold_setup(F, N)
    E = e1 + e2
    v = []
    om = F.one
    for i in range(3):
        ti = su['t'] * om
        one = L.const(F, F.one, N + 2)
        vi = (ti * F.fr(3) + one) * (ti * F.fr(E)).inv()
        v.append(vi)
        om = F.red(om * F.omega)
    FF = su['FF']
    if dF is not None:
        FF = FF + dF(su)
    return Data(F, e1, e2, v, su['xi'], su['c1'], FF, K, frame), su

def lr_setup(F, Dmax, N):
    """Large radius: local parameter t, y=t^3. A, F, Q."""
    a = []
    for d in range(Dmax + 2):
        a.append(Fr(factorial(3 * d), factorial(d) ** 3 ) * (-1) ** d)
    # series of A, theta A in y
    # F = theta log A
    A = a
    thA = [d * a[d] for d in range(Dmax + 2)]
    invA = [Fr(0)] * (Dmax + 2)
    invA[0] = Fr(1)
    for n in range(1, Dmax + 2):
        invA[n] = -sum(A[k] * invA[n - k] for k in range(1, n + 1))
    Fs = [sum(thA[k] * invA[n - k] for k in range(n + 1)) for n in range(Dmax + 2)]
    # mirror map: log(Q/y)=3 Sum (3d-1)!/(d!)^3 (-y)^d
    lg = [Fr(0)] * (Dmax + 2)
    for d in range(1, Dmax + 2):
        lg[d] = 3 * Fr(factorial(3 * d - 1), factorial(d) ** 3) * (-1) ** d
    ex = [Fr(0)] * (Dmax + 2)
    ex[0] = Fr(1)
    for n in range(1, Dmax + 2):
        ex[n] = sum(k * lg[k] * ex[n - k] for k in range(1, n + 1)) / n
    return dict(A=A, F=Fs, Qy=ex)

def lr_data(F, e1, e2, Dmax, Nt, K, frame='hodge'):
    su = lr_setup(F, Dmax + Nt // 3 + 2, Nt)
    E = e1 + e2
    def in_t(ser):
        c = [F.zero] * Nt
        for d, x in enumerate(ser):
            if 3 * d < Nt:
                c[3 * d] = F.fr(x)
        return L(F, 0, c)
    A = in_t(su['A'])
    FF = in_t(su['F'])
    t = L(F, 1, [F.one] + [F.zero] * (Nt - 1))
    Dc = (t * t * t) * F.fr(27) + L.const(F, F.one, Nt)
    v = []
    om = F.one
    for i in range(3):
        ti = t * om
        one = L.const(F, F.one, Nt)
        v.append((ti * F.fr(3) + one) * (ti * F.fr(E)).inv())
        om = F.red(om * F.omega)
    return Data(F, e1, e2, v, Dc, A, FF, K, frame), su

def needed_K(g):
    return 3 * g - 3 + 2

def rr(x, p=P):
    """rational reconstruction"""
    r0, r1 = p, x % p
    s0, s1 = 0, 1
    B = int((p // 2) ** 0.5)
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 < 0:
        r1, s1 = -r1, -s1
    return Fr(r1, s1)
