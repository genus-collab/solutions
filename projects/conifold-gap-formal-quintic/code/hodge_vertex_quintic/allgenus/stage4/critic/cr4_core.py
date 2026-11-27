# -*- coding: utf-8 -*-
"""cr4_core.py -- Stage-4 critic (allgenus/stage4/critic/): exact colour-0 machinery, written from scratch.

Imports nothing from the project and nothing from the Stage-4 workers or referees (stdlib only).
Field(None) = Q (Fraction), Field(p) = Z/p.  Notation of UNIFORM_REDUCTION (UR) s.1-s.5 and of the Stage-4 notes:
u = 1 - 5^5 q, w = zhat, z = u^(4/5) w, colour 0, hatted WKB function Phi = Sum_n w^n ph_n(u), Phi_j = A_(-j) Phi_(j-1),
A_e = 1 + w (1-u) ((4/5) th_w - e/5 - th_u).  The PF operator M^5 - 1 - q prod_k (5M + kz), M = zD + L, gives (HPF); its
w^(n+1) coefficient is (1-u) [ (4n - 5 th_u) ph_n + r_n ] (re-derived by the critic; the left factor (1-u) is manifest).
Contents
  bernoulli(N); hcoef(F, N, w5)       QRR constants [z^n] h_(w_0), w_0 = {1 - zeta^j, j=1..4} U {w5}; w5 = -5 is the project's
  tau_ser(F, K)                       tau = 1/2F1(1,4/5;9/5;u) mod u^K
  periods(F, K, al, be, dg)           gamma, C0 = f0 + al g + be f2, t_c = g/C0, C1 = D t_c (dg: deliberately wrong g, a probe)
  model_rows(F, N)                    F_b(t), b = 0..4 (UR Prop. S, colour 0), closure F_0 = -T_0 F_4 asserted
  Wkb(F, N, h, mode, fv, ka)          ph_n, n <= N; free coefficient f_m = [u^4m] ph_5m by mode qrr / zero / list / can
  KA, Gauge.KE, Gauge.KP              the three readings of [z^5q] c(z):  KA = evid referee A (period free),
                                      KE = evid Lemma N/P ([u^4q zhat^5q] Ghat_11 through row 1), KP = poly s.3 (first row)
"""
from fractions import Fraction as Fr
from math import comb


class Field:
    def __init__(self, p=None): self.p = p
    def c(self, x):
        x = Fr(x)
        return x if self.p is None else x.numerator % self.p * pow(x.denominator % self.p, -1, self.p) % self.p
    def r(self, x): return x if self.p is None else x % self.p
    def inv(self, x): return 1 / Fr(x) if self.p is None else pow(x % self.p, -1, self.p)
    def name(self): return 'Q' if self.p is None else '2^%d-1' % self.p.bit_length()


# ------------------------------------------------------------------ truncated series / polynomials in one variable
def smul(F, a, b, K):
    out = [0] * K
    for i, x in enumerate(a[:K]):
        if x:
            for j, y in enumerate(b[:K - i]):
                if y: out[i + j] += x * y
    return [F.r(v) for v in out]

def sadd(F, a, b, cb=1):
    out = list(a) + [0] * (len(b) - len(a))
    for i, y in enumerate(b): out[i] = F.r(out[i] + cb * y)
    return out

def sinv(F, a, K):
    i0 = F.inv(a[0]); out = [i0] + [0] * (K - 1)
    for n in range(1, K):
        out[n] = F.r(-sum(a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1)) * i0)
    return out

def spow(F, a, e, K):
    """a^e for a[0] = 1 (n P_n = Sum_k ((e+1)k - n) a_k P_(n-k))"""
    e = F.c(e); out = [F.c(1)] + [0] * (K - 1)
    for n in range(1, K):
        s = sum(((e + 1) * k - n) * a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1))
        out[n] = F.r(s * F.inv(F.c(n)))
    return out

def oneminus(F, a, K=None):
    """(1-u) a (polynomial, or truncated mod u^K)"""
    out = [F.r(x - (a[k - 1] if k else 0)) for k, x in enumerate(a)] + ([F.r(-a[-1])] if a else [])
    return out if K is None else (out + [0] * K)[:K]

def thop(F, cst, a):
    """(cst - th_u) a, cst a rational"""
    c = F.c(cst); return [F.r(c * x - k * x) for k, x in enumerate(a)]

def peval(F, a, x):
    s = 0
    for c in reversed(a): s = F.r(s * x + c)
    return s


# ------------------------------------------------------------------ Bernoulli numbers, QRR constants, tau
def bernoulli(N):
    B = [Fr(1)]
    for m in range(1, N + 1): B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / Fr(m + 1))
    return B

def cyclo_powersums(M):
    """p_m = Sum_(j=1..4) (1 - zeta^j)^(-m), m = 0..M: Newton sums of y^4 - 2y^3 + 2y^2 - y + 1/5 (y_j = 1/(1-zeta^j))"""
    e = [None, Fr(2), Fr(2), Fr(1), Fr(1, 5)]; p = [Fr(4)]
    for m in range(1, M + 1):
        s = sum((-1) ** (i - 1) * e[i] * p[m - i] for i in range(1, min(m - 1, 4) + 1))
        if m <= 4: s += (-1) ** (m - 1) * m * e[m]
        p.append(s)
    return p

def hcoef(F, N, w5=-5):
    """[z^n] h_(w_0)(z), h_w = exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k)), b_k = B_2k/(2k(2k-1)), n = 0..N"""
    B = bernoulli(N + 2); p = cyclo_powersums(N + 1); lg = [Fr(0)] * (N + 1)
    for k in range(1, (N + 1) // 2 + 1):
        m = 2 * k - 1
        lg[m] = -B[2 * k] / (2 * k * (2 * k - 1)) * (p[m] + Fr(w5) ** (-m))
    h = [Fr(1)]
    for n in range(1, N + 1): h.append(sum(m * lg[m] * h[n - m] for m in range(1, n + 1)) / n)
    return [F.c(x) for x in h]

def tau_ser(F, K):
    return sinv(F, [F.c(Fr(4, 5 * n + 4)) for n in range(K)], K)


# ------------------------------------------------------------------ conifold periods (PF4 = u D^4 - (1-u)(D+1/5)..(D+4/5) + ..)
def _pf4_mono(m):
    """P(u^(m+eps)), P = u D^4 - (1-u)(2D^3 + (7/5)D^2 + (2/5)D + 24/625), D = -(1-u) d/du: {k: (re, eps-part)} of u^(k+eps)"""
    def D(c):
        out = {}
        for k, (x0, x1) in c.items():
            y0, y1 = k * x0, x0 + k * x1                       # (k + eps)(x0 + eps x1)
            for kk, s in ((k, 1), (k - 1, -1)):
                a = out.get(kk, (0, 0)); out[kk] = (a[0] + s * y0, a[1] + s * y1)
        return out
    jets = [{m: (Fr(1), Fr(0))}]
    for _ in range(4): jets.append(D(jets[-1]))
    out = {}
    def acc(c, cf, sh):
        for k, (x0, x1) in c.items():
            a = out.get(k + sh, (0, 0)); out[k + sh] = (a[0] + cf * x0, a[1] + cf * x1)
    acc(jets[4], 1, 1)
    for j, cf in ((3, Fr(2)), (2, Fr(7, 5)), (1, Fr(2, 5)), (0, Fr(24, 625))):
        acc(jets[j], -cf, 0); acc(jets[j], cf, 1)
    return out

def periods(F, K, al=0, be=0, dg=0):
    """returns dict: gamma, f0, f2, g, C0, tc, C1 (series mod u^K) in the conifold frame (al, be); dg adds dg*f2 to g (probe)"""
    Kp = K + 2; mono = [_pf4_mono(m) for m in range(Kp + 1)]
    x, y = mono[1][-1], mono[2][-1]
    assert x[0] == 0 and y == (0, 2) and mono[1].get(-2, (0, 0)) == (0, 0), (x, y)
    gam = -Fr(x[1]) / 2                                       # log solution exists iff the eps-part at u^(-1+eps) vanishes
    def sol(init):
        a = [Fr(v) for v in init] + [Fr(0)] * (Kp - 3)
        for n in range(3, Kp):
            a[n] = -sum(a[m] * mono[m].get(n - 3, (0, 0))[0] for m in range(max(0, n - 4), n)) / mono[n][n - 3][0]
        return a
    f0, f1, f2 = sol((1, 0, 0)), sol((0, 1, 0)), sol((0, 0, 1))
    g = [f1[i] + (gam + Fr(dg)) * f2[i] for i in range(Kp)]
    C0 = [f0[i] + Fr(al) * g[i] + Fr(be) * f2[i] for i in range(Kp)]
    cv = lambda s: [F.c(v) for v in s]
    tc = smul(F, cv(g), sinv(F, cv(C0), Kp), Kp)
    d = [F.r((k + 1) * tc[k + 1]) for k in range(Kp - 1)]
    C1 = [F.r(-(d[k] - (d[k - 1] if k else 0))) for k in range(K)]
    return dict(gamma=gam, f0=cv(f0)[:K], f2=cv(f2)[:K], g=cv(g)[:K], C0=cv(C0)[:K], tc=tc[:K], C1=C1)


# ------------------------------------------------------------------ model rows F_b(t) (UR Prop. S, colour 0)
def model_rows(F, N):
    def T(k, f): return [F.r(-f[0])] + [F.r(-f[m] + F.c(Fr(k - 4 * (m - 1), 5)) * f[m - 1]) for m in range(1, len(f))]
    F0 = [F.c(1)]
    for n in range(1, N + 1):
        g = F0 + [0, 0]
        for k in (-1, -2, 2, 1, 0): g = T(k, g)
        F0.append(F.r(g[n + 1] * F.inv(F.c(4 * n))))           # pivot -4n of (1 + T0 T1 T2 T-2 T-1) at t^(n+1)
    F1 = T(-1, F0); F2 = [F.r(-v) for v in T(-2, F1)]; F3 = T(2, F2); F4 = [F.r(-v) for v in T(1, F3)]
    assert [F.r(-v) for v in T(0, F4)] == F0, "model closure F_0 = -T_0 F_4 fails"
    return [F0, F1, F2, F3, F4]


# ------------------------------------------------------------------ colour-0 WKB polynomials and the period-free functional
def taupows(F, tau, M, K, sign):
    base = [F.r(sign * x) for x in tau[:K]]; out = [[F.c(1)] + [0] * (K - 1)]
    for _ in range(M): out.append(smul(F, out[-1], base, K))
    return out

def KA(F, P1, q, tpm, F2):
    """referee-A reading: K_q = -[u^4q] S_q, S_q = Sum_(n<=5q) (A_(-1)Phi)_n (-tau)^(5q-n) F_(2,5q-n); REG: S_q = O(u^4q)"""
    K = 4 * q + 1; S = [0] * K
    for n in range(5 * q + 1):
        c = F2[5 * q - n]
        if c: S = sadd(F, S, smul(F, P1[n], tpm[5 * q - n], K), c)
    return F.r(-S[4 * q]), all(x == 0 for x in S[:4 * q])

class Wkb:
    """ph_n(u), n = 0..N (colour 0).  mode 'qrr': ph_5m(1) = h_5m; 'zero': f_m = 0; 'list': f_m = fv[m]; 'can': f_m fixed by
    KA = 0 (needs ka = (tpm, F2)).  Records: H_fail (n with 5 not | n and ph_n(1) != h_n), res_ok, deg_ok, deg_eq."""
    def __init__(self, F, N, h=None, mode='qrr', fv=None, ka=None):
        self.F, self.N, self.f = F, N, {}
        self.ph = [[F.c(1)]]; self.P = [[[F.c(1)]] for _ in range(5)]      # P[j][m] = Phi_j^(m) = [w^m] Phi_j
        self.H_fail, self.res_ok, self.deg_ok, self.deg_eq = [], True, True, []
        self.shift = fv if (mode == 'qrr' and fv) else {}                    # qrr mode: ph_5m(1) = h_5m + fv[m] (referee-A rule)
        for n in range(1, N + 1): self._step(n, h, mode, fv, ka)

    def _step(self, n, h, mode, fv, ka):
        F, P = self.F, self.P
        Pb = [[]]
        for j in range(1, 5): Pb.append(sadd(F, Pb[j - 1], oneminus(F, thop(F, Fr(4 * (n - 1) + j, 5), P[j - 1][n - 1]))))
        r = []
        for i in range(1, 6): r = sadd(F, r, thop(F, Fr(4 * n + i, 5), Pb[i - 1]))
        r = sadd(F, r, Pb[4], -3)
        for cf, sh, j, m in ((Fr(17, 5), 1, 3, n - 1), (Fr(9, 5), 2, 2, n - 2), (Fr(274, 625), 3, 1, n - 3), (Fr(24, 625), 4, 0, n - 4)):
            if m >= 0: r = sadd(F, r, [0] * sh + P[j][m], -F.c(cf))
        ph = [0] * max(len(r), 4 * n // 5 + 1)
        for k, x in enumerate(r):
            if 5 * k == 4 * n: self.res_ok &= (x == 0)
            elif x: ph[k] = F.r(x * F.inv(F.c(5 * k - 4 * n)))
        if n % 5 == 0:
            m, k0 = n // 5, 4 * n // 5
            if mode == 'qrr': ph[k0] = F.r(h[n] - sum(ph) + F.c(self.shift.get(m, 0)))
            elif mode == 'list': ph[k0] = F.c(fv[m])
        while len(ph) > 1 and ph[-1] == 0: ph.pop()
        self.deg_ok &= len(ph) - 1 <= 4 * n // 5
        self.deg_eq.append(len(ph) - 1 == 4 * n // 5)
        self.ph.append(ph)
        for j in range(5): P[j].append(sadd(F, Pb[j], ph))
        if n % 5 == 0 and mode == 'can':
            m, k0 = n // 5, 4 * n // 5
            K0, _ = KA(F, P[1], m, ka[0], ka[1])               # slope 1 in f_m: shift every Phi_j^(n) at u^(4m)
            for lst in [self.ph[n]] + [P[j][n] for j in range(5)]:
                lst.extend([0] * (k0 + 1 - len(lst))); lst[k0] = F.r(lst[k0] - K0)
        if n % 5 == 0: self.f[n // 5] = self.ph[n][4 * n // 5] if len(self.ph[n]) > 4 * n // 5 else 0
        if h is not None and n % 5 and F.r(sum(self.ph[n]) - h[n]) != 0: self.H_fail.append(n)


class Gauge:
    """the two readings with periods: KE (row 1, Lemma P pairing) and KP (first-row triangular expansion, poly s.3)"""
    def __init__(self, F, W, Q, per, tau, Fb):
        self.F, self.Q, K = F, Q, 4 * Q + 1; self.K = K
        iC0, iC1 = sinv(F, per['C0'], K), sinv(F, per['C1'], K)
        Psi = [smul(F, W.ph[n], iC0, K) for n in range(5 * Q + 1)]
        self.Qh1 = [smul(F, sadd(F, Psi[n], oneminus(F, thop(F, Fr(4 * n - 3, 5), Psi[n - 1]), K) if n else []), iC1, K)
                    for n in range(5 * Q + 1)]
        self.tpm, self.F2 = taupows(F, tau, 5 * Q, K, -1), Fb[2]
        tpp = taupows(F, tau, 5 * Q, K, 1); g = []; self.div_ok = True
        for N in range(5 * Q + 1):
            s = Psi[N]
            for k in range(N):
                c = Fb[(-k) % 5][N - k]
                if c: s = sadd(F, s, smul(F, g[k], tpp[N - k], K), -c)
            g.append([F.r(x * F.inv(Fb[(-N) % 5][0])) for x in s])
            mu = 4 * (N // 5) + (0, 1, 2, 3, 3)[N % 5]
            self.div_ok &= all(x == 0 for x in g[N][:min(mu, K)])
        self.g = g

    def KE(self, q):
        F, K = self.F, 4 * q + 1; S = [0] * K
        for n in range(5 * q + 1):
            c = self.F2[5 * q - n]
            if c: S = sadd(F, S, smul(F, self.Qh1[n][:K], self.tpm[5 * q - n], K), c)
        return F.r(S[4 * q]), all(x == 0 for x in S[:4 * q])

    def KP(self, q):
        return self.F.r(self.g[5 * q][4 * q] + self.g[5 * q - 1][4 * q])
