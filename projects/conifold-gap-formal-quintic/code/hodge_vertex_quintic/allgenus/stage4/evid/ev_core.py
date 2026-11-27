# -*- coding: utf-8 -*-
"""ev_core.py -- Stage 4 evidence: exact COLOUR-0 machinery for the reduced criterion of the connection problem (EVID.md s.2).

Fields QF() = Q (Fractions), PF(p) = F_p (ints); F.red reduces, F.of converts a rational, F.inv inverts.  u-series and polynomials are
lists (index = power of u).  Writes no file.  Independent of the fq_* / ref_* modules (used only as cross-checks in ev_validate.py).
  periods(F, K, al, be) -> (C0, C1, tc, g): C0 = varpi_0 = f0 + al g + be f2, t_c = g/varpi_0, C1 = D t_c = -(1-u) t_c' (K terms);
        f_k = u^k + O(u^3) and g = u + gamma u^2 + .. solve u D^4 - (1-u)(2D^3 + (7/5)D^2 + (2/5)D + 24/625), D = (u-1) d/du
        (the log solution g log u + h is found with dual numbers u^(n+eps), exactly as in natural_form/referee/ref_rows.py)
  tau_series(F, K) = 1/2F1(1,4/5;9/5;u) = 1/Sum 4u^n/(5n+4)
  qrr_h(F, N): h_n = [z^n] h_{w_0}(z), w_0 = {1 - zeta^j (j = 1..4)} U {-5}; y_j = 1/(1 - zeta^j) are the roots of
        5y^4 - 10y^3 + 10y^2 - 5y + 1, so Sum_{a in w_0} a^-m = p_m(y) + (-1/5)^m with p_m from Newton's identities (exact over Q)
  model_F(F, N): F_0..F_4 (t^0..t^N), the homogeneous WKB rows of x^4/4 + v/x, colour 0 (UNIFORM_REDUCTION s.2 Prop. S)
  Wkb(F): ph_n(u), colour 0, from (HPF) of ref_rows.py, INCREMENTAL; free coefficient [u^(4n/5)] ph_n (5 | n) supplied or set later
  Rows(F, C0, C1, KU): hatted rows Qh_a (a = 0..4) of colour 0 from the ph_n, truncated at u^KU
  pair(F, Qa, Fb, taupow, N, KU) = Sum_{n<=N} Qa[n] (-tau)^(N-n) F_{b,N-n}: the [zhat^N] coefficient of the colour-reduced pairing
"""
from fractions import Fraction as Fr
import math


class QF(object):
    name = "Q"
    def of(self, x): return Fr(x)
    def red(self, x): return x
    def inv(self, x): return 1 / Fr(x)


class PF(object):
    def __init__(self, p): self.p = p; self.name = ("2^%d-1" % ((p + 1).bit_length() - 1)) if (p + 1) & p == 0 else "F_%d" % p
    def of(self, x): x = Fr(x); return x.numerator % self.p * pow(x.denominator % self.p, self.p - 2, self.p) % self.p
    def red(self, x): return x % self.p
    def inv(self, x): return pow(x % self.p, self.p - 2, self.p)


# ------------------------------------------------------------------ series / polynomials
def smul(F, a, b, K):
    la = min(len(a), K); lb = min(len(b), K); out = [0] * K
    for i in range(la):
        x = a[i]
        if not x: continue
        for j in range(min(lb, K - i)):
            y = b[j]
            if y: out[i + j] += x * y
    return [F.red(c) for c in out]


def sinv(F, a, K):
    i0 = F.inv(a[0]); out = [F.red(i0)]
    for k in range(1, K):
        s = 0
        for j in range(1, min(k, len(a) - 1) + 1): s += a[j] * out[k - j]
        out.append(F.red(-s * i0))
    return out


def padd(F, a, b, c=1):
    n = max(len(a), len(b)); out = [0] * n
    for k, x in enumerate(a): out[k] = x
    for k, x in enumerate(b): out[k] = out[k] + c * x
    out = [F.red(x) for x in out]
    while out and not out[-1]: out.pop()
    return out


def one_minus(F, a, K=None):                                      # (1 - u) a; polynomial (K None) or truncated series
    n = len(a) + 1 if K is None else K
    return [F.red((a[k] if k < len(a) else 0) - (a[k - 1] if 0 < k <= len(a) else 0)) for k in range(n)]


def Bop(F, c, a, K=None):
    """(1 - u)(c a - theta a), theta = u d/du"""
    return one_minus(F, [F.red(c * x - k * x) for k, x in enumerate(a)], K)


# ------------------------------------------------------------------ periods, tau
def _dm(F, a, b): return (F.red(a[0] * b[0]), F.red(a[0] * b[1] + a[1] * b[0]))
def _da(F, a, b): return (F.red(a[0] + b[0]), F.red(a[1] + b[1]))

def _pf_mono(F, n):
    """the PF operator applied to u^(n+eps): {k: dual coefficient of u^(k+eps)}"""
    def D(c):
        out = {}
        for k, x in c.items():
            for kk, sg in ((k, 1), (k - 1, -1)):
                out[kk] = _da(F, out.get(kk, (0, 0)), _dm(F, (F.red(sg * k), F.red(sg)), x))
        return out
    jets = [{n: (F.of(1), 0)}]
    for _ in range(4): jets.append(D(jets[-1]))
    out = {}
    def acc(c, cf, sh):
        for k, x in c.items(): out[k + sh] = _da(F, out.get(k + sh, (0, 0)), (F.red(x[0] * cf), F.red(x[1] * cf)))
    acc(jets[4], F.of(1), 1)
    for j, cf in ((3, Fr(2)), (2, Fr(7, 5)), (1, Fr(2, 5)), (0, Fr(24, 625))):
        c = F.of(cf); acc(jets[j], F.red(-c), 0); acc(jets[j], c, 1)
    return out


def periods(F, K, al=0, be=0):
    N = K + 3; mono = {n: _pf_mono(F, n) for n in range(N + 1)}
    def solve(b):
        for k in range(-3, 0):
            r = (0, 0)
            for n in range(3): r = _da(F, r, _dm(F, b[n], mono[n].get(k, (0, 0))))
            assert r == (0, 0), ("obstruction", k)
        for n in range(3, N):
            k = n - 3; r = (0, 0)
            for m in range(max(0, k - 1), n): r = _da(F, r, _dm(F, b[m], mono[m].get(k, (0, 0))))
            ind = mono[n][k]; assert F.red(ind[0] - n * (n - 1) ** 2 * (n - 2)) == 0
            i0 = F.inv(ind[0]); b[n] = _dm(F, (F.red(-r[0]), F.red(-r[1])), (i0, F.red(-ind[1] * i0 * i0)))
        return b
    fs = []
    for k in range(3):
        b = solve({n: (0, F.of(1) if n == k else 0) for n in range(3)}); fs.append([b[n][1] for n in range(N)])
    x = mono[1][-1]; y = mono[2][-1]; assert x[0] == 0 and y == (0, F.of(2))
    gam = F.red(-x[1] * F.inv(2))
    b = solve({0: (0, 0), 1: (F.of(1), 0), 2: (gam, 0)}); g = [b[n][0] for n in range(N)]
    a_, b_ = F.of(al), F.of(be)
    v0 = [F.red(fs[0][k] + a_ * g[k] + b_ * fs[2][k]) for k in range(K + 2)]
    tc = smul(F, g[:K + 2], sinv(F, v0, K + 2), K + 2)
    dtc = [F.red((k + 1) * tc[k + 1]) for k in range(K + 1)]
    C1 = [F.red(-x) for x in one_minus(F, dtc, K + 1)][:K]
    return v0[:K], C1, tc[:K], g[:K]


def tau_series(F, K):
    return sinv(F, [F.of(Fr(4, 5 * n + 4)) for n in range(K)], K)


# ------------------------------------------------------------------ QRR constants of colour 0
_B = [Fr(1)]
def bern(n):
    while len(_B) <= n:
        m = len(_B); _B.append(-sum(math.comb(m + 1, j) * _B[j] for j in range(m)) / Fr(m + 1))
    return _B[n]


def qrr_h(F, N, twist=True):
    p = [F.of(4), F.of(2), F.of(0), F.of(-1), F.of(Fr(-4, 5))]
    for m in range(5, N + 2): p.append(F.red(2 * p[m - 1] - 2 * p[m - 2] + p[m - 3] - F.of(Fr(1, 5)) * p[m - 4]))
    lg = [0] * (N + 1)
    for k in range(1, N // 2 + 2):
        if 2 * k - 1 > N: break
        S = F.red(p[2 * k - 1] + (F.of(Fr(-1, 5) ** (2 * k - 1)) if twist else 0))   # twist: the weight -5 of the line bundle
        lg[2 * k - 1] = F.red(-F.of(bern(2 * k) / (2 * k * (2 * k - 1))) * S)
    out = [F.of(1)]
    for n in range(1, N + 1):
        out.append(F.red(sum(k * lg[k] * out[n - k] for k in range(1, n + 1)) * F.inv(n)))
    return out


# ------------------------------------------------------------------ model rows (colour 0)
def model_F(F, N):
    """F_b, b = 0..4, t^0..t^N; G_j = T_{k_j} G_(j-1), k = (-1,-2,2,1,0), T_k G = -G + (t/5)(k - 4 t d/dt) G; closure F_0 + G_5 = 0"""
    ks = (-1, -2, 2, 1, 0); G = [[F.of(1)]] + [[F.of((-1) ** j)] for j in range(1, 6)]
    def nxt(j, m, Gp_m, Gp_m1): return F.red(-Gp_m + F.of(Fr(ks[j - 1] - 4 * (m - 1), 5)) * Gp_m1)
    for n in range(1, N + 1):
        vals = []
        for X in (0, 1):                              # affine in X = F_0[n]
            g = [X]; h = [0]                          # g = G_j[n], h = G_j[n+1] with F_0[n+1] = 0
            for j in range(1, 6):
                g.append(nxt(j, n, g[j - 1], G[j - 1][n - 1])); h.append(nxt(j, n + 1, h[j - 1], g[j - 1]))
            vals.append(h[5])
        piv = F.red(vals[1] - vals[0]); assert F.red(piv + 4 * n) == 0, "model pivot"
        X = F.red(-vals[0] * F.inv(piv)); G[0].append(X); g = [X]
        for j in range(1, 6): g.append(nxt(j, n, g[j - 1], G[j - 1][n - 1])); G[j].append(g[j])
    assert all(F.red(G[0][n] + G[5][n]) == 0 for n in range(N + 1)), "model closure"
    return [G[0], G[1], [F.red(-x) for x in G[2]], [F.red(-x) for x in G[3]], G[4]]


# ------------------------------------------------------------------ the colour-0 WKB polynomials
class Wkb(object):
    """ph[n] (polynomials in u) for colour 0.  X[j][N] = [zhat^N] Phi_j, Phi_j = A_{-j} Phi_{j-1}.  step(f): next n; for 5 | n the free
    coefficient [u^(4n/5)] is f (callable(n, ph_without_f) -> value, or a value, or None = 0 for now; then set_free)."""
    HP = ((2, Fr(17, 5), 1, 3), (3, Fr(9, 5), 2, 2), (4, Fr(274, 625), 3, 1), (5, Fr(24, 625), 4, 0))

    def __init__(self, F):
        self.F = F; self.ph = [[F.of(1)]]; self.X = [[[F.of(1)]] for _ in range(6)]; self.resonance = []; self.free = {}
        self.hp = [(m, F.of(cf), up, j) for (m, cf, up, j) in self.HP]

    def step(self, f=None):
        F = self.F; n = len(self.ph); X = self.X
        X0 = [[]]
        for j in range(1, 6): X0.append(padd(F, X0[j - 1], Bop(F, F.of(Fr(4 * (n - 1) + j, 5)), X[j - 1][n - 1])))
        X1 = [[]]
        for j in range(1, 6): X1.append(padd(F, X1[j - 1], Bop(F, F.of(Fr(4 * n + j, 5)), X0[j - 1])))
        br = [F.red(3 * x) for x in X0[4]]
        for (m, cf, up, j) in self.hp:
            if n + 1 - m >= 0: br = padd(F, br, [0] * up + [F.red(cf * x) for x in X[j][n + 1 - m]])
        R = padd(F, X1[5], one_minus(F, br), -1)
        assert F.red(sum(R)) == 0, "R not divisible by 1-u at n = %d" % n
        rho = []; s = 0
        for x in R[:-1]: s = F.red(s + x); rho.append(F.red(-s))          # rho = -R/(1-u)
        kf = 4 * n // 5 if n % 5 == 0 else None
        c = [0] * max(len(rho), (kf + 1) if kf is not None else 1)
        if kf is not None:                                                # resonance: rho_kf must vanish (record if beyond deg rho)
            ok = kf >= len(rho) or F.red(rho[kf]) == 0; self.resonance.append((n, ok, kf >= len(rho)))
            assert ok, "resonance obstruction at n = %d" % n
        for k, x in enumerate(rho):
            if k != kf: c[k] = F.red(x * F.inv(4 * n - 5 * k))
        self.ph.append(c)
        for j in range(6): X[j].append(padd(F, X0[j], c))
        if kf is not None:
            val = f(n, c) if callable(f) else f
            if val is not None: self.set_free(n, val)
        return self.ph[n]

    def set_free(self, n, val):
        F = self.F; kf = 4 * n // 5; cur = self.ph[n][kf] if kf < len(self.ph[n]) else 0
        d = F.red(val - cur); self.free[n] = F.red(val)
        if d:
            mono = [0] * kf + [d]; self.ph[n] = padd(F, self.ph[n], mono)
            for j in range(6): self.X[j][n] = padd(F, self.X[j][n], mono)


def qrr_free(F, h):
    """free-coefficient rule of the q = 0 (QRR) normalisation: ph_n(1) = h_n"""
    return lambda n, c: F.red(h[n] - sum(c) + c[4 * n // 5])


# ------------------------------------------------------------------ rows and the colour-reduced pairing
class Rows(object):
    """hatted rows of colour 0, incremental in n: Qh[0][n] = ph_n/C0, Qh[a+1] = kh_a A_{e_a} Qh[a] (ref_hat.py conventions)"""
    E = (-1, -2, 2, 1)

    def __init__(self, F, C0, C1, KU, amax=2):
        self.F = F; self.KU = KU; self.amax = amax; self.Q = [[] for _ in range(amax + 1)]
        C0i = sinv(F, C0, KU); C1i = sinv(F, C1, KU)
        self.kh = [C1i, smul(F, smul(F, C0, C0, KU), smul(F, C1, C1, KU), KU), C1i, C0i]; self.C0i = C0i

    def set(self, n, ph_n):
        """(re)compute all rows at zhat^n (needs zhat^(n-1) final); drops anything above n"""
        F = self.F; K = self.KU
        for a in range(self.amax + 1): del self.Q[a][n:]
        self.Q[0].append(smul(F, ph_n, self.C0i, K))
        for a in range(self.amax):
            cur = self.Q[a][n]
            if n >= 1:
                e = self.E[a]; cur = padd(F, cur, Bop(F, F.of(Fr(4 * (n - 1) - e, 5)), self.Q[a][n - 1], K))
            self.Q[a + 1].append(smul(F, cur + [0] * (K - len(cur)), self.kh[a], K))


def pair(F, Qa, Fb, taupow, N, KU):
    out = [0] * KU
    for n in range(N + 1):
        m = N - n; c = Fb[m] if m % 2 == 0 else F.red(-Fb[m])
        if not c: continue
        t = smul(F, Qa[n], taupow[m], KU)
        for k in range(KU): out[k] += c * t[k]
    return [F.red(x) for x in out]
