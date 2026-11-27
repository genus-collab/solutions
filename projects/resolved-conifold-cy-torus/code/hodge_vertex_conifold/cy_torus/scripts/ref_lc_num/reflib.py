# reflib.py -- independent referee implementation (F_p, p = 2^127 - 1) of the frame-B rows of the
# resolved conifold x C^2 at Q = 1 (t = log Q), built from the superpotential
#     W(U) = a1 log(1 - Q U) - a2 log U + A log(1 - U),  A = a0 + a2,
# by the Wick/Gaussian expansion (D1), the Q = 0 constant h_nw (Lemma 13.3) and the frame-B twist h_w^-1,
# the H-row (L_i + z D_i) r_0, eta from (C4'), Delta_i = 1/eta(e_i, e_i), the pairing eta^-1,
# dilaton factors Dil_i(d) = (-1)^d r_{0,i}[d-1], edge factors EdN (C11).
# Nothing is imported from earlier project code.  Pure python, no sympy.
import math
from fractions import Fraction as Fr

P = (1 << 127) - 1

def inv(a):
    a %= P
    if a == 0:
        raise ZeroDivisionError("inverse of 0 mod p")
    return pow(a, P - 2, P)

def F(x):
    x = Fr(x)
    return x.numerator % P * inv(x.denominator) % P

def ratrec(a):
    """rational reconstruction of a mod P (returns Fraction or None)"""
    a %= P
    if a == 0:
        return Fr(0)
    r0, r1 = P, a
    s0, s1 = 0, 1
    Bd = math.isqrt(P // 2)
    while r1 > Bd:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > Bd:
        return None
    return Fr(r1, s1)

def show(a):
    r = ratrec(a)
    return str(r) if r is not None else ("0x%x" % (a % P))

# ---------------------------------------------------------------- truncated Laurent series in t
# S(v, c): sum_i c[i] t^(v+i), known for exponents v .. v+len(c)-1; prec = v + len(c) (unknown from t^prec on)

class S:
    __slots__ = ('v', 'c')
    def __init__(self, v, c):
        c = [x % P for x in c]
        i = 0
        n = len(c)
        while i < n and c[i] == 0:
            i += 1
        self.v = v + i
        self.c = c[i:]
    @property
    def prec(self):
        return self.v + len(self.c)
    def iszero(self):
        return len(self.c) == 0
    def coef(self, k):
        if k >= self.prec:
            raise ValueError("coefficient t^%d beyond precision %d" % (k, self.prec))
        if k < self.v:
            return 0
        return self.c[k - self.v]
    def __add__(self, o):
        if isinstance(o, int):
            o = const(o, self.prec)
        prec = min(self.prec, o.prec)
        v = min(self.v, o.v)
        out = [0] * max(prec - v, 0)
        for i, x in enumerate(self.c):
            j = self.v + i - v
            if j < len(out):
                out[j] = x
        for i, x in enumerate(o.c):
            j = o.v + i - v
            if j < len(out):
                out[j] = (out[j] + x) % P
        return S(v, out)
    __radd__ = __add__
    def __neg__(self):
        return S(self.v, [(-x) % P for x in self.c])
    def __sub__(self, o):
        if isinstance(o, int):
            o = const(o, self.prec)
        return self + (-o)
    def __rsub__(self, o):
        return (-self) + o
    def scale(self, k):
        k %= P
        return S(self.v, [x * k % P for x in self.c])
    def __mul__(self, o):
        if isinstance(o, int):
            return self.scale(o)
        prec = min(self.prec + o.v, o.prec + self.v)
        v = self.v + o.v
        n = prec - v
        if n <= 0 or not self.c or not o.c:
            return S(prec, [])
        out = [0] * n
        ca, cb = self.c, o.c
        for i, x in enumerate(ca):
            if i >= n:
                break
            if x:
                m = min(len(cb), n - i)
                for j in range(m):
                    out[i + j] += x * cb[j]
        return S(v, [y % P for y in out])
    __rmul__ = __mul__
    def inverse(self):
        if not self.c:
            raise ZeroDivisionError("inverse of zero series")
        n = len(self.c)
        c = self.c
        i0 = inv(c[0])
        out = [i0] + [0] * (n - 1)
        for k in range(1, n):
            s = 0
            for j in range(1, k + 1):
                if j < n:
                    s += c[j] * out[k - j]
            out[k] = (-s * i0) % P
        return S(-self.v, out)
    def __truediv__(self, o):
        if isinstance(o, int):
            return self.scale(inv(o))
        return self * o.inverse()
    def der(self):
        # d/dt; the coefficient t^(prec-1) is known, its derivative is t^(prec-2); so prec drops by one
        out = [(self.v + i) % P * x % P for i, x in enumerate(self.c)]
        return S(self.v - 1, out)
    def trunc(self, prec):
        if prec >= self.prec:
            return self
        return S(self.v, self.c[:max(prec - self.v, 0)])
    def shift(self, k):
        return S(self.v + k, self.c)
    def __repr__(self):
        return "S(v=%d, prec=%d, [%s])" % (self.v, self.prec, ", ".join(show(x) for x in self.c[:6]) + (", ..." if len(self.c) > 6 else ""))

def const(a, prec):
    a %= P
    if a == 0 or prec <= 0:
        return S(prec, [])
    return S(0, [a] + [0] * (prec - 1))

def tvar(prec):
    return S(1, [1] + [0] * (prec - 2))

def s_exp(f):
    """exp of a series with valuation >= 1 (regular, no constant term); result known to f.prec"""
    prec = f.prec
    if f.v <= 0 and f.c:
        raise ValueError("exp needs valuation >= 1")
    n = prec
    fc = [f.coef(k) if k < prec else 0 for k in range(n)]
    # E' = f' E  ->  k E_k = sum_{j=1}^{k} j f_j E_{k-j}
    E = [1] + [0] * (n - 1)
    for k in range(1, n):
        s = 0
        for j in range(1, k + 1):
            if fc[j]:
                s += j * fc[j] % P * E[k - j]
        E[k] = s % P * inv(k) % P
    return S(0, E)

def s_log1(f):
    """log of a unit with constant term 1"""
    prec = f.prec
    if f.v != 0 or f.c[0] != 1:
        raise ValueError("log needs constant term 1")
    n = prec
    g = f - 1  # valuation >= 1
    # log(1+g) : L' = g'/(1+g)
    return s_integrate(g.der() / f, prec)

def s_integrate(f, prec):
    """antiderivative with zero constant term, for f with valuation >= 0"""
    out = [0] * prec
    for k in range(f.v, min(f.prec, prec - 1)):
        out[k + 1] = f.coef(k) * inv(k + 1) % P
    return S(0, out).trunc(min(prec, f.prec + 1))

def s_sqrt1(f):
    """square root of a unit with constant term 1, branch with value 1"""
    if f.v != 0 or f.c[0] != 1:
        raise ValueError("sqrt needs constant term 1")
    return s_exp(s_log1(f).scale(inv(2)))

def s_cmp(a, b):
    """first exponent where a and b differ (within common precision), or None; also returns compared precision"""
    prec = min(a.prec, b.prec)
    d = a - b
    for k in range(d.v, prec):
        if d.coef(k) != 0:
            return k, prec
    return None, prec

# ---------------------------------------------------------------- z-Laurent polynomials with S coefficients: dict k -> S
def zmul(a, b, zcap):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            k = i + j
            if k > zcap:
                continue
            p = x * y
            if k in out:
                out[k] = out[k] + p
            else:
                out[k] = p
    return out

def zadd(a, b):
    out = dict(a)
    for k, y in b.items():
        out[k] = out[k] + y if k in out else y
    return out

def zscale(a, s):
    return {k: x * s for k, x in a.items()}

def zconst(ser):
    return {0: ser}

def zprec(a):
    return min(x.prec for x in a.values())

def zcoef(a, k, prec):
    return a.get(k, S(prec, []))

# ---------------------------------------------------------------- Bernoulli numbers, h_v(z)
def bernoulli(nmax):
    Bn = [Fr(0)] * (nmax + 1)
    Bn[0] = Fr(1)
    for m in range(1, nmax + 1):
        s = Fr(0)
        for k in range(m):
            s += Fr(math.comb(m + 1, k)) * Bn[k]
        Bn[m] = -s / (m + 1)
    return Bn

def h_series(v, zmax, prec):
    """h_v(z) = Exp[-Sum_k bk_k z^(2k-1) Sum_a v_a^(1-2k)], bk_k = B_2k/(2k(2k-1)); v a list of mod-p weights.
    Returns dict k -> S (constant series of precision prec)."""
    Bn = bernoulli(2 * ((zmax + 1) // 2) + 2)
    ex = {}
    for k in range(1, zmax // 2 + 2):
        e = 2 * k - 1
        if e > zmax:
            break
        bk = F(Bn[2 * k] / (2 * k * (2 * k - 1)))
        ps = 0
        for va in v:
            ps = (ps + pow(inv(va), 2 * k - 1, P)) % P
        ex[e] = (-bk * ps) % P
    # exp of a scalar z-series with zero constant term
    E = [1] + [0] * zmax
    ec = [ex.get(e, 0) for e in range(zmax + 1)]
    for k in range(1, zmax + 1):
        s = 0
        for j in range(1, k + 1):
            if ec[j]:
                s += j * ec[j] % P * E[k - j]
        E[k] = s % P * inv(k) % P
    return {k: const(E[k], prec) for k in range(zmax + 1) if E[k]}

def zinv_scalar(a, zmax, prec):
    """inverse of a scalar z-power series given as dict k->S with constant-series values (a[0] != 0)"""
    ac = [0] * (zmax + 1)
    for k, x in a.items():
        if 0 <= k <= zmax and x.c:
            ac[k] = x.coef(0)
    i0 = inv(ac[0])
    out = [i0] + [0] * zmax
    for k in range(1, zmax + 1):
        s = 0
        for j in range(1, k + 1):
            s += ac[j] * out[k - j]
        out[k] = (-s * i0) % P
    return {k: const(out[k], prec) for k in range(zmax + 1) if out[k]}

def geom(E_, j, zmax, prec):
    """1/(E + j z) as a scalar z-power series"""
    iE = inv(E_)
    r = (-j * iE) % P
    out = {}
    cur = iE
    for k in range(zmax + 1):
        if cur:
            out[k] = const(cur, prec)
        cur = cur * r % P
    return out

# ---------------------------------------------------------------- the construction at a weight point
class Point:
    def __init__(self, a0, a1, a2, e1, e2, N, Z, frame=None, use_hnw=True, tag=""):
        """weights as Fractions; N = number of t-coefficients kept (series known mod t^N at the start);
        Z = max z-order of the rows; frame = multiset (list of Fractions) whose h^-1 multiplies the rows
        (default w = (E, -e1, -e2)); use_hnw: include the Q = 0 constant h_nw = h_{(a1, A, -e1, -e2)}."""
        self.tag = tag
        self.fr = dict(a0=Fr(a0), a1=Fr(a1), a2=Fr(a2), e1=Fr(e1), e2=Fr(e2))
        self.fr['A'] = self.fr['a0'] + self.fr['a2']
        self.fr['B'] = self.fr['a0'] + self.fr['a1']
        self.fr['E'] = self.fr['a0'] + self.fr['a1'] + self.fr['a2']
        self.fr['u'] = (self.fr['E'] - 2 * self.fr['a1']) * (self.fr['E'] - 2 * self.fr['a2'])
        for k in ('a1', 'a2', 'A', 'B', 'E', 'e1', 'e2'):
            if self.fr[k] == 0:
                raise ValueError("(ND) violated: %s = 0" % k)
        self.a0, self.a1, self.a2, self.e1, self.e2 = (F(self.fr[k]) for k in ('a0', 'a1', 'a2', 'e1', 'e2'))
        self.A, self.B, self.E, self.u = (F(self.fr[k]) for k in ('A', 'B', 'E', 'u'))
        self.N, self.Z = N, Z
        self.isCY = (self.fr['E'] == self.fr['e1'] + self.fr['e2'])
        a0, a1, a2, A, B, E = self.a0, self.a1, self.a2, self.A, self.B, self.E
        # t, Q, tau
        self.t = tvar(N)
        self.Q = s_exp(self.t)
        self.Qi = s_exp(-self.t)
        self.tau = 1 - self.Qi
        # critical points by Newton on LeqU = B Q U^2 - (a0 + Q (a1 - a2)) U - a2
        def leq(U):
            return B * self.Q * U * U - (a0 + (a1 - a2) * self.Q) * U - a2
        def dleq(U):
            return 2 * B * self.Q * U - (a0 + (a1 - a2) * self.Q)
        self.U = {}
        for i, u0 in ((1, 1), (2, (-a2 * inv(B)) % P)):
            U = const(u0, N)
            for _ in range(8):
                U = U - leq(U) / dleq(U)
            self.U[i] = U
            r = leq(U)
            assert r.iszero() or r.v >= N - 1, ("Newton failed", i, r)
        # L_i = -a1 Q U/(1 - Q U)
        self.L = {}
        for i in (1, 2):
            QU = self.Q * self.U[i]
            self.L[i] = (-a1) * QU / (1 - QU)
        # W_k(U_i): theta^k W, W = a1 log(1 - QU) - a2 log U + A log(1 - U);  theta log(1 - aU) = 1 - p_a, theta p = p^2 - p
        self.kmax = 2 * Z + 2
        pol = [[0, 1]]  # Pol_0(p) = p  (coefficients in p: index = power)
        for n in range(1, self.kmax):
            prev = pol[-1]
            # Pol_n = (p^2 - p) Pol_{n-1}'
            d = [(j * prev[j]) % P for j in range(1, len(prev))]  # derivative coefficients, powers 0..len-2
            new = [0] * (len(prev) + 1)
            for j, x in enumerate(d):
                new[j + 2] = (new[j + 2] + x) % P
                new[j + 1] = (new[j + 1] - x) % P
            pol.append(new)
        self.pol = pol
        self.W = {}
        for i in (1, 2):
            U = self.U[i]
            p1 = (1 - U).inverse()
            pQ = (1 - self.Q * U).inverse()
            def horner(coeffs, p):
                acc = const(0, max(p.prec, N))
                for cf in reversed(coeffs):
                    acc = acc * p + cf
                return acc
            Wk = {}
            Wk[1] = a1 * (1 - pQ) - a2 + A * (1 - p1)   # should vanish at a critical point
            for k in range(2, self.kmax + 1):
                Wk[k] = -(a1 * horner(pol[k - 1], pQ)) - A * horner(pol[k - 1], p1)
            self.W[i] = Wk
        # Wick rows (point frame)
        self.rpt = {i: self.wick_row(self.W[i]) for i in (1, 2)}
        # constants: h_nw (Q = 0 constant) and frame twist h_frame^-1
        prec = zprec(self.rpt[1])
        self.rowprec = prec
        if frame is None:
            frame = [self.fr['E'], -self.fr['e1'], -self.fr['e2']]
        self.frame = [Fr(x) for x in frame]
        nw = [a1, A, (-self.e1) % P, (-self.e2) % P]
        hnw = h_series(nw, Z, prec) if use_hnw else {0: const(1, prec)}
        hfr = h_series([F(x) for x in self.frame], Z, prec)
        hfr_inv = zinv_scalar(hfr, Z, prec)
        self.twist = zmul(hnw, hfr_inv, Z)   # the scalar multiplying every row
        self.hnw, self.hframe = hnw, hfr
        # frame rows r_0
        self.r0 = {i: zmul(self.twist, self.rpt[i], Z) for i in (1, 2)}
        # eta (C4'), eta^-1
        d = inv(A * B % P * self.e1 % P * self.e2 % P)
        self.eta = [[E * inv(a1 * a2 % P) % P * d % P, d], [d, (-a0 * d) % P]]
        det = (self.eta[0][0] * self.eta[1][1] - self.eta[0][1] * self.eta[1][0]) % P
        di = inv(det)
        self.etainv = [[self.eta[1][1] * di % P, (-self.eta[0][1] * di) % P], [(-self.eta[1][0] * di) % P, self.eta[0][0] * di % P]]
        # Delta_i = 1/eta(e_i, e_i), e_i = (-L_j, 1)/(L_i - L_j) in the basis (1, H)
        self.Delta = {}
        for i, j in ((1, 2), (2, 1)):
            Li, Lj = self.L[i], self.L[j]
            num = self.eta[0][0] * Lj * Lj - (2 * self.eta[0][1]) * Lj + self.eta[1][1]
            den = (Li - Lj) * (Li - Lj)
            self.Delta[i] = den / num
        # rho
        self.rho = s_sqrt1(self.Delta[1] * self.t * self.t * inv((-self.e1 * self.e2 % P * pow(E, 3, P)) % P)).inverse()
        # H rows
        self.r1 = {}
        for i in (1, 2):
            dd = self.Delta[i].der() / self.Delta[i] * inv(2)
            Li = self.L[i]
            row = {}
            for m in range(Z + 1):
                x = Li * zcoef(self.r0[i], m, prec)
                if m >= 1:
                    f = zcoef(self.r0[i], m - 1, prec)
                    x = x + f.der() - dd * f
                row[m] = x
            self.r1[i] = row
        self.v = {i: (self.r0[i], self.r1[i]) for i in (1, 2)}

    # ---- Wick expansion (D1): Lnorm = Sum_n (2n-1)!! (-z/W_2)^n [v^(2n)] Exp[Sum_{k>=3} W_k v^k/(k! z)]
    def wick_row(self, Wk):
        Z = self.Z
        N = self.N
        nv = 6 * Z + 1  # v-degrees 0..6Z
        # Ex[n] = dict j -> S : coefficient of v^n y^j in Exp[Sum_k (W_k/k!) v^k y]
        Ex = [dict() for _ in range(nv)]
        Ex[0][0] = const(1, N)
        fk = {k: Wk[k].scale(inv(math.factorial(k - 1))) for k in range(3, self.kmax + 1)}  # k W_k/k! = W_k/(k-1)!
        for n in range(1, nv):
            acc = {}
            for k in range(3, min(n, self.kmax) + 1):
                for j, s in Ex[n - k].items():
                    p = fk[k] * s
                    acc[j + 1] = acc[j + 1] + p if (j + 1) in acc else p
            invn = inv(n)
            Ex[n] = {j: s.scale(invn) for j, s in acc.items()}
        W2inv = Wk[2].inverse()
        W2pow = [const(1, N)]
        for n in range(1, 3 * Z + 1):
            W2pow.append(W2pow[-1] * W2inv)
        row = {}
        for m in range(Z + 1):
            acc = None
            for n in range(m, 3 * m + 1):
                j = n - m
                if j not in Ex[2 * n]:
                    continue
                df = 1
                for q in range(1, 2 * n, 2):
                    df = df * q % P
                term = Ex[2 * n][j].scale(df if n % 2 == 0 else (-df) % P)
                term = term * W2pow[n] if n > 0 else term
                acc = term if acc is None else acc + term
            row[m] = acc if acc is not None else const(0, N)
        return row

    # ---- pairing of covectors a = (a_0, a_1) (S-valued), <a, b> = Sum etainv^{ab} a_a b_b
    def pair(self, a, b):
        M = self.etainv
        return M[0][0] * a[0] * b[0] + M[0][1] * (a[0] * b[1] + a[1] * b[0]) + M[1][1] * a[1] * b[1]

    def vcomp(self, i, p):
        prec = self.rowprec
        return (zcoef(self.r0[i], p, prec), zcoef(self.r1[i], p, prec))

    def Pij(self, i, j, p, q):
        """P^ij_pq = -(-1)^(p+q) <v_i[p], v_j[q]>"""
        x = self.pair(self.vcomp(i, p), self.vcomp(j, q))
        return -x if (p + q) % 2 == 0 else x

    def edge(self, i, j, m, n):
        """E_ij(m, n) = Sum_{k=0}^{n} (-1)^k P^ij_{m+1+k, n-k}; needs m + n + 1 <= Z"""
        assert m + n + 1 <= self.Z
        acc = None
        for k in range(n + 1):
            term = self.Pij(i, j, m + 1 + k, n - k)
            if k % 2:
                term = -term
            acc = term if acc is None else acc + term
        return acc

    def dil(self, i, d):
        x = zcoef(self.r0[i], d - 1, self.rowprec)
        return -x if d % 2 else x

    def symp(self, i, j, m):
        """[w^m] <v_i(w), v_j(-w)>"""
        acc = None
        for p in range(m + 1):
            q = m - p
            term = self.pair(self.vcomp(i, p), self.vcomp(j, q))
            if q % 2:
                term = -term
            acc = term if acc is None else acc + term
        return acc

    # ---- closed-form ingredients
    def dphi(self):
        """dphi = phi_tau(s*) - phi_0(A/E), s* = (1 - U_1)/tau"""
        A, a1, a2, E = self.A, self.a1, self.a2, self.E
        sstar = (1 - self.U[1]) / self.tau
        self.sstar = sstar
        l1 = s_log1(sstar.scale(E * inv(A) % P))
        l2 = s_log1((1 - sstar).scale(E * inv(a1) % P))
        l3 = s_log1(self.U[1])
        return A * l1 + a1 * l2 - a2 * l3

    def key_rhs(self, zcap, Apar=None, a2par=None, Epar=None):
        """(E/(E+z)) (tau/t) 2F1(1 + a2/z, 1 + A/z; 2 + E/z; tau) Exp[-dphi/z], as a z-Laurent dict truncated at zcap.
        Optional replacements of the hypergeometric parameters (controls)."""
        A = self.A if Apar is None else Apar
        a2 = self.a2 if a2par is None else a2par
        E = self.E if Epar is None else Epar
        N = self.N
        prec = self.rowprec
        dphi = self.dphi()
        zbig = zcap + 2 * N + 2
        # Ephi = Sum_l (-dphi)^l/(l! z^l)
        Ephi = {}
        pw = const(1, N)
        for l in range(N):
            Ephi[-l] = pw.scale(inv(math.factorial(l)))
            pw = pw * (-dphi)
        # F = Sum_n F_n tau^n, F_n = prod_{j<=n}(a2 + j z)(A + j z) / (n! z^n prod_{j=2}^{n+1}(E + j z))
        Fz = {0: const(1, N)}
        taun = const(1, N)
        for n in range(1, N):
            taun = taun * self.tau
            num = {0: const(1, N)}
            for j in range(1, n + 1):
                num = zmul(num, {0: const(a2, N), 1: const(j, N)}, zbig)
                num = zmul(num, {0: const(A, N), 1: const(j, N)}, zbig)
            den = {0: const(inv(math.factorial(n)), N)}
            for j in range(2, n + 2):
                den = zmul(den, geom(E, j, zbig, N), zbig)
            Fn = zmul(num, den, zbig)
            Fn = {k - n: x for k, x in Fn.items()}
            Fz = zadd(Fz, zscale(Fn, taun))
        pref = zscale(geom(E, 1, zbig, N), (self.tau / self.t).scale(E))
        rhs = zmul(zmul(pref, Fz, zbig), Ephi, zbig)
        return {k: x for k, x in rhs.items() if k <= zcap}

    def A_j_closed(self, j, dphi):
        """A_j(t) = (-1)^j C_j t^-1 zeta^j 2F1(j B/E, j a1/E; j; tau), zeta = tau e^(dphi/E), C_j = P_{j-1}/(j!(j-1)! E^(2j-4))"""
        a1, a2, A, B, E = self.a1, self.a2, self.A, self.B, self.E
        N = self.N
        Pk = 1
        k = j - 1
        for i in range(1, k + 1):
            Pk = Pk * (((k + 1) * a2 - i * E) % P) % P * (((k + 1) * A - i * E) % P) % P
        Cj = Pk * inv(math.factorial(j) * math.factorial(j - 1) % P) % P * pow(inv(E), 2 * j - 4, P) % P if 2 * j - 4 >= 0 else \
            Pk * inv(math.factorial(j) * math.factorial(j - 1) % P) % P * pow(E, 4 - 2 * j, P) % P
        zeta = self.tau * s_exp(dphi.scale(inv(E)))
        al, be = j * B % P * inv(E) % P, j * a1 % P * inv(E) % P
        hyp = const(0, N)
        term = const(1, N)
        taun = const(1, N)
        coef = 1
        for n in range(N):
            hyp = hyp + taun.scale(coef)
            coef = coef * ((al + n) % P) % P * ((be + n) % P) % P * inv(((j + n) % P) * (n + 1) % P) % P
            taun = taun * self.tau
        zj = const(1, N)
        for _ in range(j):
            zj = zj * zeta
        out = (zj * hyp).scale(Cj if j % 2 == 0 else (-Cj) % P)
        return S(out.v - 1, out.c)  # divide by t

# ---------------------------------------------------------------- leaf decomposition tools (mod p linear algebra)
def solve(Mat, rhs):
    """solve square system mod p (Gaussian elimination)"""
    n = len(Mat)
    M = [row[:] + [rhs[i]] for i, row in enumerate(Mat)]
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c] % P:
                piv = r
                break
        if piv is None:
            raise ZeroDivisionError("singular")
        M[c], M[piv] = M[piv], M[c]
        ic = inv(M[c][c])
        M[c] = [x * ic % P for x in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [(x - f * y) % P for x, y in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]

def leaf_fit(coeffs, k, E, sign=1):
    """coeffs[a] = [psi^a] f for a = 0..M.  Fit f = c_0 lam_0 + Sum_{j=1}^{k+1} c_j lam_j with lam_j = 1/(E - sign j psi),
    lam_0 = 1/E (constant), from a = 0..k+1; return (c list, first a > k+1 with mismatch or None, M)."""
    n = k + 2
    nodes = [0] + list(range(1, k + 2))
    iE = inv(E)
    def lamc(j, a):
        if j == 0:
            return iE if a == 0 else 0
        return pow((sign * j) % P, a, P) * pow(iE, a + 1, P) % P
    Mat = [[lamc(j, a) for j in nodes] for a in range(n)]
    c = solve(Mat, [coeffs[a] % P for a in range(n)])
    bad = None
    for a in range(n, len(coeffs)):
        pred = sum(c[idx] * lamc(j, a) for idx, j in enumerate(nodes)) % P
        if pred != coeffs[a] % P:
            bad = a
            break
    return c, bad, len(coeffs) - 1

def leaf_fit2(C, k, E, sign=1):
    """C[a][b] = [psi^a psi'^b] f, a, b = 0..M.  Fit f = Sum_{i,l=0}^{k+1} c_il lam_i(psi) lam_l(psi') from the (k+2)x(k+2) block;
    return (c matrix, first (a,b) mismatch or None)."""
    n = k + 2
    nodes = list(range(0, k + 2))
    iE = inv(E)
    def lamc(j, a):
        if j == 0:
            return iE if a == 0 else 0
        return pow((sign * j) % P, a, P) * pow(iE, a + 1, P) % P
    V = [[lamc(j, a) for j in nodes] for a in range(n)]
    # Vinv
    Vinv = [solve(V, [1 if r == i else 0 for r in range(n)]) for i in range(n)]  # columns
    Vinv = [[Vinv[j][i] for j in range(n)] for i in range(n)]  # Vinv[i][j]
    blk = [[C[a][b] % P for b in range(n)] for a in range(n)]
    # c = Vinv blk Vinv^T
    tmp = [[sum(Vinv[i][a] * blk[a][b] for a in range(n)) % P for b in range(n)] for i in range(n)]
    c = [[sum(tmp[i][b] * Vinv[l][b] for b in range(n)) % P for l in range(n)] for i in range(n)]
    M = len(C) - 1
    bad = None
    for a in range(M + 1):
        for b in range(M + 1):
            if a < n and b < n:
                continue
            pred = sum(c[i][l] * lamc(i, a) % P * lamc(l, b) for i in range(n) for l in range(n)) % P
            if pred != C[a][b] % P:
                bad = (a, b)
                break
        if bad:
            break
    return c, bad

def polymul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] = (out[i + j] + x * y) % P
    return out
