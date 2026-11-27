# -*- coding: utf-8 -*-
"""rf_core.py -- REFEREE of the Stage-4 task 'evid': independent core.  Imports NOTHING from the project or from the worker (stdlib only).
Fields: QF() exact (Fraction), PF(p) mod p.  u-series / polynomials are lists (index = power of u).
 periods(F, K, al, be) -> (C0, C1, tc, gamma): quintic PF  D^4 - (1-u)(D+1/5)(D+2/5)(D+3/5)(D+4/5),  D = q d/dq = (u-1) d/du, solved at
     u = 0 by Frobenius on terms u^m (log u)^l; the vanishing period g is FOUND as the analytic solution (initial data a_0, a_1, a_2) whose
     log partner g log u + h exists (3x3 obstruction matrix at u^-3, u^-2, u^-1 of rank 2); C0 = f0 + al g + be f2, tc = g/C0, C1 = D tc.
 itau(F, K) = 1/tau from the canonical one-form (5/4) u^(4/5)/tau = Int_0^u t^(-1/5) dt/(1-t);  spow(F, f, a, K) = f^a (f(0) = 1).
 qrr(F, N, zi, twist): [z^n] exp(-Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_(w) w^(1-2k)), weights zi*{1 - zeta^j} U {-5 zi}; B_n by
     Akiyama-Tanigawa; exact power sums by the Galois trace of (1 - zeta)^-m in Q(zeta5) (basis 1, z, z^2, z^3).
 model(F, N, zi): rows F_b of the LG model x^4/4 + v/x for colour zi, from F_1 = T_-1 F_0, F_2 = -T_-2 F_1, F_3 = T_2 F_2,
     F_4 = -T_1 F_3, F_0 = -T_0 F_4 (closure solved order by order), T_k = -zi + (t/5)(k - 4 t d/dt).
 Wkb(F, zi): psi = Sum_n z^n psi_n(s) DIRECTLY from the PF operator P = M^5 - 1 - q prod_(k=1..5)(5M + k z), M = zD + zi L, L = 1/s,
     u = s^5 = 1 - 3125 q, acting on Laurent monomials (M s^e = zi s^(e-1) + (e/5) z (s^e - s^(e-5))); psi_0 = 1/s; psi_n is a dict
     {e: c}; the free constant (s^-1, 5 | n) is set by step(free); ph_n(u) = Sum_k c_(-1-4n+5k) u^k.
 A_at(F, e, zi, Xn, Xp, n, K): [zhat^n] A_e X, used for the hatted rows Qh_0 = zi Phi/C0, Qh_(a+1) = kh_a A_(e_a) Qh_a.
"""
from fractions import Fraction as Fr

class QF(object):
    name = "Q"
    def of(self, x): return Fr(x)
    def red(self, x): return x
    def inv(self, x): return 1 / Fr(x)

class PF(object):
    def __init__(self, p): self.p = p; self.name = ("2^%d-1" % p.bit_length()) if (p + 1) & p == 0 else str(p)
    def of(self, x): x = Fr(x); return x.numerator % self.p * pow(x.denominator, -1, self.p) % self.p
    def red(self, x): return x % self.p
    def inv(self, x): return pow(x % self.p, -1, self.p)
    def zeta5(self):
        for g in range(2, 500):
            z = pow(g, (self.p - 1) // 5, self.p)
            if z != 1: return z

def mul(F, a, b, K):
    out = [0] * K
    for i, x in enumerate(a[:K]):
        if x:
            for j, y in enumerate(b[:K - i]):
                if y: out[i + j] += x * y
    return [F.red(c) for c in out]

def inv(F, a, K):
    i0 = F.inv(a[0]); out = [F.red(i0)]
    for k in range(1, K):
        s = 0
        for j in range(1, min(k, len(a) - 1) + 1): s += a[j] * out[k - j]
        out.append(F.red(-s * i0))
    return out

def spow(F, f, al, K):
    al = Fr(al); f = list(f[:K]) + [0] * max(0, K - len(f)); y = [F.of(1)]
    for n in range(1, K):
        s = 0
        for k in range(1, n + 1):
            if f[k]: s += F.of((al + 1) * k - n) * f[k] * y[n - k]
        y.append(F.red(s * F.inv(n)))
    return y

def deriv(F, a): return [F.red(k * a[k]) for k in range(1, len(a))]
def one_minus(F, a, K): return [F.red((a[k] if k < len(a) else 0) - (a[k - 1] if 0 < k <= len(a) else 0)) for k in range(K)]

# ------------------------------------------------------------------ periods (Frobenius with log terms)
def _D(F, t):
    out = {}
    for (m, l), c in t.items():
        for key, v in (((m, l), m * c), ((m - 1, l), -m * c)) + ((((m, 0), c), ((m - 1, 0), -c)) if l else ()):
            out[key] = out.get(key, 0) + v
    return {k: F.red(v) for k, v in out.items()}

def _L(F, t):
    d = [t]
    for _ in range(4): d.append(_D(F, d[-1]))
    out = {}
    for (m, l), c in d[4].items(): out[(m + 1, l)] = out.get((m + 1, l), 0) + c
    for j, cf in ((3, Fr(2)), (2, Fr(7, 5)), (1, Fr(2, 5)), (0, Fr(24, 625))):
        c0 = F.of(cf)
        for (m, l), c in d[j].items():
            out[(m, l)] = out.get((m, l), 0) - c0 * c; out[(m + 1, l)] = out.get((m + 1, l), 0) + c0 * c
    return {k: F.red(v) for k, v in out.items()}

def periods(F, K, al=0, be=0):
    N = K + 4; ell = [_L(F, {(m, 0): F.of(1)}) for m in range(N + 1)]
    assert all(k[0] >= 0 for m in range(3) for k, v in ell[m].items() if v), "L(power series) has no negative powers"
    def sol(a0):
        a = list(a0)
        for n in range(3, N):
            piv = ell[n].get((n - 3, 0), 0); assert F.red(piv - n * (n - 1) ** 2 * (n - 2)) == 0
            s = sum(a[m] * ell[m].get((n - 3, 0), 0) for m in range(max(0, n - 4), n))
            a.append(F.red(-s * F.inv(piv)))
        return a
    f = [sol([F.of(int(i == k)) for i in range(3)]) for k in range(3)]
    ob = [[_L(F, {(k, 1): F.of(1)}).get((N_, 0), 0) for k in range(3)] for N_ in (-3, -2, -1)]
    assert all(_L(F, {(k, 1): F.of(1)}).get((N_, 1), 0) == 0 for k in range(3) for N_ in (-3, -2, -1))
    # kernel with a_1 = 1: rows of ob must have rank 2; solve ob . (a0, 1, a2) = 0
    pairs = [(r, s) for i, r in enumerate(ob) for s in ob[i + 1:] if F.red(r[0] * s[2] - r[2] * s[0]) != 0]
    assert pairs, "no non-singular 2x2 minor in (a_0, a_2)"
    r, s = pairs[0]; det = F.red(r[0] * s[2] - r[2] * s[0])
    a0 = F.red((-r[1] * s[2] + r[2] * s[1]) * F.inv(det)); a2 = F.red((-r[0] * s[1] + r[1] * s[0]) * F.inv(det))
    assert all(F.red(r[0] * a0 + r[1] + r[2] * a2) == 0 for r in ob), "obstruction matrix rank 3: no log solution"
    g = [F.red(a0 * x + y + a2 * w) for x, y, w in zip(*f)]
    C0 = [F.red(x + F.of(al) * y + F.of(be) * w) for x, y, w in zip(f[0], g, f[2])][:K + 1]
    tc = mul(F, g[:K + 1], inv(F, C0, K + 1), K + 1)
    C1 = [F.red(-x) for x in one_minus(F, deriv(F, tc), K)]          # D tc = (u - 1) tc'
    return C0[:K], C1, tc[:K], (a0, a2)

def itau(F, K): return [F.of(Fr(4, 5 * n + 4)) for n in range(K)]

# ------------------------------------------------------------------ QRR constants
def bernoulli(F, n):
    A = [0] * (n + 1); B = []
    for m in range(n + 1):
        A[m] = F.of(Fr(1, m + 1))
        for j in range(m, 0, -1): A[j - 1] = F.red(j * (A[j - 1] - A[j]))
        B.append(A[0])
    return B

def _zmul(x, y):
    c = [0] * 7
    for i in range(4):
        for j in range(4): c[i + j] += x[i] * y[j]
    c[0] += c[5]; c[1] += c[6]
    return [c[k] - c[4] for k in range(4)]

def powersums_Q(N):
    """p_m = Sum_(j=1..4) (1 - zeta^j)^-m, m = 0..N, exact (Galois trace in Q(zeta5))"""
    y = [Fr(4, 5), Fr(3, 5), Fr(2, 5), Fr(1, 5)]; x = [Fr(1), 0, 0, 0]; out = []
    for m in range(N + 1):
        out.append(4 * x[0] - x[1] - x[2] - x[3]); x = _zmul(x, y)
    return out

def qrr(F, N, zi=1, twist=True, B=None, zeta=None):
    """h_n of colour zi (zi = zeta^i in F_p; exact path only for zi = 1); B: optional list of Bernoulli numbers (for controls)"""
    B = bernoulli(F, N + 1) if B is None else B
    if isinstance(F, QF): S = [F.of(x) + (F.of(Fr(-1, 5) ** m) if twist else 0) for m, x in enumerate(powersums_Q(N))]
    else:
        zeta = F.zeta5() if zeta is None else zeta
        w = [F.red(zi - zi * pow(zeta, j, F.p)) for j in range(1, 5)] + ([F.red(-5 * zi)] if twist else [])
        wi = [F.inv(x) for x in w]; S = [F.red(sum(pow(x, m, F.p) for x in wi)) for m in range(N + 1)]
    lg = [0] * (N + 1)
    for k in range(1, N // 2 + 2):
        if 2 * k - 1 <= N: lg[2 * k - 1] = F.red(-B[2 * k] * F.inv(2 * k * (2 * k - 1)) * S[2 * k - 1])
    h = [F.of(1)]
    for n in range(1, N + 1): h.append(F.red(sum(k * lg[k] * h[n - k] for k in range(1, n + 1)) * F.inv(n)))
    return h

# ------------------------------------------------------------------ the model rows
OPS = ((-1, 1), (-2, -1), (2, 1), (1, -1), (0, -1))          # F_(b+1) = sign T_k F_b, b = 0..4 (F_5 = F_0)

def model(F, N, zi=1):
    def T(k, t):
        out = {}
        for m, c in t.items():
            out[m] = out.get(m, 0) - zi * c; out[m + 1] = out.get(m + 1, 0) + F.of(Fr(k - 4 * m, 5)) * c
        return out
    o = []
    for m in range(N + 1):
        t = {m: F.of(1)}
        for k, sg in OPS: t = {e: F.red(sg * c) for e, c in T(k, t).items()}
        assert F.red(t[m] - 1) == 0; o.append(t)
    F0 = [F.red(zi)]                                              # normalisation F_0^(i)(0) = zeta^i = the quintic's Qh_0^(i)(0, 0)
    for n in range(2, N + 2):
        piv = o[n - 1][n]; assert F.red(piv - 4 * (n - 1) * zi ** 4) == 0, "model pivot"
        s = sum(o[m].get(n, 0) * F0[m] for m in range(max(0, n - 5), n - 1))
        F0.append(F.red(-s * F.inv(piv)))
    rows = [F0]
    for k, sg in OPS:
        t = T(k, dict(enumerate(rows[-1]))); rows.append([F.red(sg * t.get(m, 0)) for m in range(N + 1)])
    assert rows[5] == rows[0], "model closure"
    return rows[:5]

# ------------------------------------------------------------------ WKB directly from the PF operator
class Wkb(object):
    CF = (1, 3, Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625))   # (1/3125) prod (5M + kz) = Sum_j CF_j z^j M^(5-j)

    def __init__(self, F, zi=1):
        self.F = F; self.zi = zi; self.zinv = F.inv(zi); self.psi = [{-1: F.of(1)}]; self.cache = {}; self.res = []
        self.cf = [F.of(c) for c in self.CF]

    def T(self, e):
        if e in self.cache: return self.cache[e]
        F = self.F; pw = [{(0, e): F.of(1)}]
        for _ in range(5):
            nx = {}
            for (j, E), c in pw[-1].items():
                q = F.of(Fr(E, 5)) * c
                for key, v in (((j, E - 1), self.zi * c), ((j + 1, E), q), ((j + 1, E - 5), -q)): nx[key] = nx.get(key, 0) + v
            pw.append({k: F.red(v) for k, v in nx.items()})
        P = dict(pw[5]); P[(0, e)] = P.get((0, e), 0) - 1
        for j in range(6):
            for (jj, E), c in pw[5 - j].items():
                v = self.cf[j] * c; P[(jj + j, E)] = P.get((jj + j, E), 0) - v; P[(jj + j, E + 5)] = P.get((jj + j, E + 5), 0) + v
        out = [dict() for _ in range(7)]
        for (j, E), c in P.items():
            c = F.red(c)
            if c: out[j][E] = c
        assert not out[0] and not out[6], "characteristic part"
        w = F.red(self.zinv * (e + 1))
        assert out[1] == ({} if w == 0 else {e + 1: w, e - 4: F.red(-w)}), "transport operator"
        self.cache[e] = out
        return out

    def step(self, free=None):
        F = self.F; n = len(self.psi); R = {}
        for j in range(2, 6):
            if n + 1 - j < 0: continue
            for e, c in self.psi[n + 1 - j].items():
                for E, t in self.T(e)[j].items(): R[E] = R.get(E, 0) + c * t
        Rp = {E: F.red(-self.zi * v) for E, v in R.items()}; Rp = {E: v for E, v in Rp.items() if v}
        W = {}
        if Rp:
            lo, hi = min(Rp), max(Rp); assert all((E - lo) % 5 == 0 for E in Rp)
            for E in range(lo, hi - 4, 5): W[E] = F.red(W.get(E - 5, 0) - Rp.get(E, 0))
            assert F.red(Rp[hi] - W.get(hi - 5, 0)) == 0, "not divisible by s^5 - 1 at n = %d" % n
        c = {}
        for E, w in W.items():
            if E + 4 == -1:
                self.res.append((n, w == 0)); assert w == 0, "resonance obstruction at n = %d" % n
            elif w: c[E + 4] = F.red(w * F.inv(E + 5))
        assert all((e + 1 + 4 * n) % 5 == 0 for e in c)
        if n % 5 == 0:
            v = free(n, c) if callable(free) else (free if free is not None else 0)
            if F.red(v): c[-1] = F.red(v)
        self.psi.append(c)
        return c

    def ph(self, n):
        c = self.psi[n]; d = max([(e + 1 + 4 * n) // 5 for e in c] or [-1])
        return [c.get(-1 - 4 * n + 5 * k, 0) for k in range(d + 1)]

def qrr_rule(F, h): return lambda n, c: F.red(h[n] - sum(c.values()))

def A_at(F, e, zi, Xn, Xp, n, K):
    """[zhat^n] A_e X = zi X_n + (1-u)((4(n-1) - e)/5 - theta) X_(n-1)   (A_e: M (s^e X) = s^(e-1) A_e X in hatted variables)"""
    t = [F.red(zi * x) for x in Xn[:K]] + [0] * max(0, K - len(Xn))
    if n:
        c = F.of(Fr(4 * (n - 1) - e, 5)); y = [F.red(c * x - k * x) for k, x in enumerate(Xp)]
        t = [F.red(a + b) for a, b in zip(t, one_minus(F, y, K))]
    return t
