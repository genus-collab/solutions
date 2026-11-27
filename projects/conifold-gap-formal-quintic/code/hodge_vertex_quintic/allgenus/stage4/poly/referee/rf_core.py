# -*- coding: utf-8 -*-
"""rf_core.py -- REFEREE of Stage 4 / poly: independent arithmetic and conifold-side data.  No project module and no
module of the worker (pp_*) is imported.  Everything mod a prime p (2^61-1, 2^89-1; both are 1 mod 5, so zeta_5 lies in F_p),
plus a few exact helpers over Q (Bernoulli numbers by Akiyama-Tanigawa, Q(zeta_5) traces).

Contents
  md(x, p)                       Fraction -> F_p
  Ser(p, K)                      truncated power series mod p (Kronecker multiplication), inv, power, exp, log, deriv, integ
  bernoulli(N)                   exact B_0..B_N (B_1 = +1/2), only even ones are used
  zeta5(p)                       a primitive 5th root of unity in F_p
  hcoef_mod(p, N), hcoef_exact(N)  [z^n] h_{w_0}(z), w_0 = {1 - zeta^j (j = 1..4)} U {-5}
  periods(S)                     f0, f1, f2 of u D^4 - (1-u)(2D^3 + 7/5 D^2 + 2/5 D + 24/625), D = -(1-u) d/du; gamma = 7/10 checked
  frame_data(S, al, be)          C0 = f0 + al g + be f2, t_c = g/C0, tau, v = u tau^(-5/4), the cyclic systems b (quintic), bt (model)
  model_rows(p, N)               colour-0 rows F_b(t) of the LG model -(x^4/4 + v/x) (own derivation of the cycle, see docstring there)
"""
from fractions import Fraction as Fr

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1
CEXP = (0, -1, 3, 2, 1)


def md(x, p):
    x = Fr(x)
    return x.numerator % p * pow(x.denominator % p, -1, p) % p


class Ser(object):
    """power series mod p truncated at u^K (lists of length K)"""
    def __init__(self, p, K):
        self.p, self.K = p, K
        self.nb = (2 * p.bit_length() + K.bit_length() + 2) // 8 + 1

    def z(self): return [0] * self.K
    def one(self): return [1] + [0] * (self.K - 1)
    def var(self): return [0, 1] + [0] * (self.K - 2)
    def cst(self, c): return [md(c, self.p)] + [0] * (self.K - 1)

    def _pack(self, a):
        nb = self.nb
        return int.from_bytes(b''.join(x.to_bytes(nb, 'little') for x in a), 'little')

    def mul(self, a, b):
        p, K, nb = self.p, self.K, self.nb
        la = len(a) if len(a) < K else K
        lb = len(b) if len(b) < K else K
        if la == 0 or lb == 0: return [0] * K
        c = (self._pack(a[:la]) * self._pack(b[:lb])).to_bytes(nb * (la + lb), 'little')
        n = min(K, la + lb - 1)
        out = [int.from_bytes(c[i * nb:(i + 1) * nb], 'little') % p for i in range(n)]
        return out + [0] * (K - n)

    def add(self, a, b, c=1):
        p = self.p
        return [(x + c * y) % p for x, y in zip(a, b)]

    def scale(self, a, c):
        p = self.p; c = md(c, p) if isinstance(c, Fr) else c % p
        return [x * c % p for x in a]

    def shift(self, a, m):
        """u^m a (m may be negative: then the dropped coefficients must vanish)"""
        if m >= 0: return ([0] * m + list(a))[:self.K]
        assert all(x == 0 for x in a[:-m]), "not divisible by u^%d" % (-m)
        return list(a[-m:]) + [0] * (-m)

    def inv(self, a):
        p, K = self.p, self.K
        i0 = pow(a[0], -1, p); out = [i0] + [0] * (K - 1)
        for n in range(1, K):
            s = 0
            for k in range(1, n + 1):
                if a[k]: s += a[k] * out[n - k]
            out[n] = -s % p * i0 % p
        return out

    def deriv(self, a):
        p = self.p
        return [(i + 1) * a[i + 1] % p for i in range(len(a) - 1)] + [0]

    def integ(self, a):
        p = self.p
        return [0] + [a[i] * pow(i + 1, -1, p) % p for i in range(self.K - 1)]

    def log1(self, a):
        """log a for a[0] = 1"""
        assert a[0] == 1
        return self.integ(self.mul(self.deriv(a), self.inv(a)))

    def exp0(self, a):
        """exp a for a[0] = 0 (recurrence n E_n = Sum k a_k E_(n-k))"""
        p, K = self.p, self.K
        assert a[0] == 0
        E = [1] + [0] * (K - 1)
        for n in range(1, K):
            s = 0
            for k in range(1, n + 1):
                if a[k]: s += k * a[k] * E[n - k]
            E[n] = s % p * pow(n, -1, p) % p
        return E

    def powr(self, a, r):
        """a^r, a[0] = 1, r rational (via exp(r log a))"""
        return self.exp0(self.scale(self.log1(a), md(r, self.p)))


def bernoulli(N):
    """exact B_0..B_N by the Akiyama-Tanigawa algorithm (B_1 = +1/2)"""
    A = [Fr(0)] * (N + 1); B = []
    for m in range(N + 1):
        A[m] = Fr(1, m + 1)
        for j in range(m, 0, -1): A[j - 1] = j * (A[j - 1] - A[j])
        B.append(A[0])
    return B


def zeta5(p):
    assert (p - 1) % 5 == 0
    for g in range(2, 200):
        z = pow(g, (p - 1) // 5, p)
        if z != 1:
            assert pow(z, 5, p) == 1
            return z
    raise ValueError


def hcoef_mod(p, N):
    """[z^n] h_{w_0}(z), n = 0..N, mod p, with zeta_5 in F_p"""
    zt = zeta5(p); B = bernoulli(N + 2)
    w = [(1 - pow(zt, j, p)) % p for j in range(1, 5)] + [(-5) % p]
    S = Ser(p, N + 1); lg = [0] * (N + 1)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        ps = sum(pow(pow(a, -1, p), m, p) for a in w) % p
        lg[m] = -md(B[2 * k] / (2 * k * (2 * k - 1)), p) * ps % p
    return S.exp0(lg)


def _cyc_mul(a, b):
    """Q(zeta_5) in the basis 1, z, z^2, z^3 (z^4 = -1 - z - z^2 - z^3)"""
    c = [Fr(0)] * 7
    for i in range(4):
        for j in range(4): c[i + j] += a[i] * b[j]
    for k in (6, 5, 4):                                   # z^k = z^(k-4) z^4
        x = c[k]; c[k] = Fr(0)
        for t in range(4): c[k - 4 + t] -= x
    return c[:4]


def hcoef_exact(N):
    """[z^n] h_{w_0}(z) over Q: Sum_(j=1..4) (1 - zeta^j)^(-m) = Tr(y^m), y = 1/(1-zeta) = (4 + 3z + 2z^2 + z^3)/5"""
    B = bernoulli(N + 2); y = [Fr(4, 5), Fr(3, 5), Fr(2, 5), Fr(1, 5)]
    assert _cyc_mul(y, [Fr(1), Fr(-1), Fr(0), Fr(0)]) == [1, 0, 0, 0]
    yp = [Fr(1), Fr(0), Fr(0), Fr(0)]; tr = [Fr(4)]
    for m in range(1, N + 1):
        yp = _cyc_mul(yp, y); tr.append(4 * yp[0] - yp[1] - yp[2] - yp[3])
    lg = [Fr(0)] * (N + 1)
    for k in range(1, N // 2 + 2):
        m = 2 * k - 1
        if m > N: break
        lg[m] = -B[2 * k] / (2 * k * (2 * k - 1)) * (tr[m] + Fr(-5) ** (-m))
    E = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1): E[n] = sum(k * lg[k] * E[n - k] for k in range(1, n + 1)) / n
    return E


def _D(S, f):
    """D f = -(1-u) f'"""
    d = S.deriv(f)
    return [(-d[i] + (d[i - 1] if i else 0)) % S.p for i in range(S.K)]


def pf4(S, f):
    d1 = _D(S, f); d2 = _D(S, d1); d3 = _D(S, d2); d4 = _D(S, d3); p = S.p
    c2, c3, c4 = md(Fr(7, 5), p), md(Fr(2, 5), p), md(Fr(24, 625), p)
    inner = [(2 * d3[i] + c2 * d2[i] + c3 * d1[i] + c4 * f[i]) % p for i in range(S.K)]
    return [((d4[i - 1] if i else 0) - inner[i] + (inner[i - 1] if i else 0)) % p for i in range(S.K)]


def periods(S):
    """f0, f1, f2 = u^k + O(u^3) (analytic solutions), mod u^K; the indicial pivot n(n-1)^2(n-2) is checked numerically"""
    p, K = S.p, S.K; big = Ser(p, K + 4); out = []
    for k in range(3):
        a = [0] * (K + 4); a[k] = 1
        for n in range(3, K + 4):
            a[n] = 0; r = pf4(big, a)[n - 3]
            a[n] = 1; piv = (pf4(big, a)[n - 3] - r) % p
            assert piv == n * (n - 1) ** 2 * (n - 2) % p
            a[n] = -r * pow(piv, -1, p) % p
        assert not any(pf4(big, a)[:K])
        out.append(a[:K])
    return out


def frame_data(S0, al=Fr(0), be=Fr(0), gamma=Fr(7, 10)):
    """C0, t_c, tau, v, b, bt (b[1], bt[1] returned as u b_1, u bt_1); computed at precision K+2, returned exact mod u^K"""
    K0 = S0.K; S = Ser(S0.p, K0 + 2)
    d = _frame_data(S, al, be, gamma)
    cut = lambda x: x[:K0] if isinstance(x[0], int) else [y[:K0] for y in x]
    return {k: cut(x) for k, x in d.items()}


def _frame_data(S, al, be, gamma):
    p = S.p; f0, f1, f2 = periods(S)
    g = S.add(f1, S.scale(f2, md(gamma, p)))
    C0 = S.add(S.add(f0, S.scale(g, md(al, p))), S.scale(f2, md(be, p)))
    tc = S.mul(g, S.inv(C0)); tcp = S.deriv(tc)
    tau = S.inv([md(Fr(4, 5 * n + 4), p) for n in range(S.K)])
    v = S.shift(S.powr(tau, Fr(-5, 4)), 1); vp = S.deriv(v)
    om = S.inv([1, p - 1] + [0] * (S.K - 2))                 # 1/(1-u)
    om3 = S.mul(om, S.mul(om, om))
    ub1 = S.scale(S.inv(S.mul(S.mul(C0, C0), S.mul(tcp, tcp))), p - 1)
    ub1 = S.mul(ub1, om3)
    b = [tcp, ub1, tcp, S.scale(S.mul(C0, om), p - 1), S.scale(S.mul(C0, om), p - 1)]
    ubt1 = S.scale(S.mul(vp, S.inv(S.shift(v, -1))), p - 1)   # u bt_1 = -u v'/v
    bt = [vp, ubt1, vp, S.scale(vp, p - 1), S.scale(vp, p - 1)]
    return dict(C0=C0, tc=tc, g=g, tau=tau, v=v, b=b, bt=bt, f0=f0, f1=f1, f2=f2)


def model_rows(p, N):
    """F_b(t), b = 0..4, colour 0, lists of length N+1 mod p.  Own derivation: Y_b = e^(U/z) sig^(c_b-1) F_b(z/sig^4),
    U = -(5/4) sig^4, sig^5 = v, in z dY_b/dv = bmod_b Y_(b+1), bmod = (1,-1/v,1,-1,-1), gives
    bmod_b sig^(c_(b+1)-1) F_(b+1) = sig^(c_b-2) T_(c_b - 1) F_b,  T_k = -1 + (t/5)(k - 4 t d/dt), i.e.
    F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4 (closure).  F_0(0) = 1."""
    M = N + 2; i5 = pow(5, -1, p)
    def T(k, F):
        out = [0] * M
        for n in range(M):
            out[n] = -F[n] % p
            if n >= 1: out[n] = (out[n] + (k - 4 * (n - 1)) * i5 * F[n - 1]) % p
        return out
    def neg(F): return [(-x) % p for x in F]
    def cyc(F0):
        F1 = T(-1, F0); F2 = neg(T(-2, F1)); F3 = T(2, F2); F4 = neg(T(1, F3))
        return [F0, F1, F2, F3, F4], neg(T(0, F4))
    F0 = [0] * M; F0[0] = 1
    for n in range(1, N + 1):
        F0[n] = 0; r = (cyc(F0)[1][n + 1] - F0[n + 1]) % p
        F0[n] = 1; piv = (cyc(F0)[1][n + 1] - F0[n + 1] - r) % p
        assert piv == 4 * n % p, "pivot"
        F0[n] = -r * pow(piv, -1, p) % p
    rows, back = cyc(F0)
    assert back[:N + 1] == F0[:N + 1]
    return [R[:N + 1] for R in rows]
