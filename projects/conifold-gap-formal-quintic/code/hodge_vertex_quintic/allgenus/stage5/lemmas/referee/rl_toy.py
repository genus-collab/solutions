# -*- coding: utf-8 -*-
"""rl_toy.py -- referee check of the lemma statements L1b, L2, L3, L1c as the chart argument uses them (FORMAL_SP_LEMMAS s.8),
on a TOY chart with rational data (new code, stdlib only; uses rl_core only).

Toy Mellin-Barnes integrand (two x-colours, twist degree k): 4-variable formal integrand
   e^(Phi/z) dsigma dy dx_0 dx_1/(x_0 x_1),  Phi = lam sigma - y + k sigma log y + sum_j (-x_j + (c_j - sigma) log x_j),
critical point (sigma, y, x_j) = (L, kL, c_j - L), lam fixed by the critical equation (it never enters the Taylor data).
 (T) total 4-variable SP^norm (full 4x4 covariance, no lemma used);
 (A) order (a), L1b + L3: one-variable SP in sigma of e^(S_0/z) A_St, A_St = e^(S_1) e^(sum_k z^(2k-1) S_2k) (full Stirling);
 (B) order (b), L2 + L1b: chart SP in (y, x_1), phase psi = f(X_0(w), w), amplitude 1/x_1 (the Jacobian X_0 cancels 1/x_0);
 (D) continuation, s.8 warning (2): c_n(L) as Laurent series in d = L - c_0: no pole at L = c_0 and c(z; c_0) = h_toy(z),
     h_toy = exp(sum_k b_k z^(2k-1) [(k c_0)^(1-2k) + (c_1 - c_0)^(1-2k)])  (the toy (Gamma));  (B) at L = c_0 directly.
Controls (must fail): (A) without the factor (Delta_in(sigma)/Delta_in(L))^(-1/2); (B) without the Jacobian; (B) with dy/y
against (T) with dy; (D) with b_k -> -b_k; (D) with the colour-0 Stirling terms dropped.
Usage: py -I -B rl_toy.py [K] [N] [PADD]   (T, A, B through z^K, default 4; D through z^N, default 10; Ser margin PADD)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rl_core import Ser, ser_geom, comb_gen, sp_norm, sp_lb, bk_list, ser_exp

K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
N = int(sys.argv[2]) if len(sys.argv) > 2 else 10
PADD = int(sys.argv[3]) if len(sys.argv) > 3 else 3 * N + 10
KT, C0, C1 = 2, Fr(1), Fr(-2)                         # twist degree k and the colours c_0, c_1
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
REL = 'logs/rl_toy_%d_%d_%d.log' % (K, N, PADD); LOG = open(os.path.join(HERE, REL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + ('   ' + info if info else ''))
    if not ok:
        FAILS.append(name)


def first_diff(a, b):
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return i
    return None


BK = bk_list(max(K, N) + 3)


def total_data(L, dyy=False):
    b, a = KT * L, [C0 - L, C1 - L]; D = 2 * K + 2; tay = {}
    for n in range(1, D + 1):
        if n >= 2:
            tay[(0, n, 0, 0)] = Fr((-1) ** (n + 1), n) / b ** (n - 1)
            for j in range(2):
                e = [0, 0, 0, 0]; e[2 + j] = n; tay[tuple(e)] = Fr((-1) ** (n + 1), n) / a[j] ** (n - 1)
        if n + 1 <= D:
            tay[(1, n, 0, 0)] = KT * Fr((-1) ** (n + 1), n) / b ** n
            for j in range(2):
                e = [1, 0, 0, 0]; e[2 + j] = n; tay[tuple(e)] = Fr((-1) ** n, n) / a[j] ** n
    amp = {}
    for p in range(2 * K + 1):
        for q in range(2 * K + 1 - p):
            for r in (range(2 * K + 1 - p - q) if dyy else [0]):
                amp[((0, r, p, q), 0)] = Fr((-1) ** (p + q + r)) / (a[0] ** p * a[1] ** q * b ** r)
    return tay, amp


def sigma_data(L, M, J, ring=Fr, dyy=False, flip=1, drop0=False, nodelta=False):
    """Taylor data of S_0 at L (orders 2..M+2) and amplitude A_St coefficients A[j][p] (j <= J, p <= M), over 'ring'.
    L, a_j and their inverse powers are supplied by the helpers below (rational, or Ser in d = L - c_0)."""
    Lp, ap = ring['Lp'], ring['ap']                      # Lp(n) = L^n, ap(j, n) = a_j^n  (n any integer)
    f = [0, 0] + [(KT * (-1) ** m * Lp(1 - m) + ap(0, 1 - m) + ap(1, 1 - m)) * Fr(1, m * (m - 1)) for m in range(2, M + 3)]
    e1 = [comb_gen(Fr(1, 2), p) * Lp(-p) for p in range(M + 1)]                 # (1 + t/L)^(1/2)
    if dyy:
        e1 = [sum(e1[i] * (-1) ** (p - i) * Lp(-(p - i)) for i in range(p + 1)) for p in range(M + 1)]   # times (1+t/L)^-1
    if nodelta:                                                                      # control: inner series alone
        e1 = [Fr(int(p == 0)) for p in range(M + 1)]
    for j in range(2):
        g = [comb_gen(Fr(-1) if nodelta else Fr(-1, 2), p) * (-1) ** p * ap(j, -p) for p in range(M + 1)]   # (1 - t/a_j)^(-1/2)
        e1 = [sum(e1[i] * g[p - i] for i in range(p + 1)) for p in range(M + 1)]
    B = [[0] * (M + 1) for _ in range(J + 1)]                                        # exponent sum_k z^(2k-1) S_2k(L+t)
    for k in range(1, (J + 1) // 2 + 1):
        for p in range(M + 1):
            v = (Fr(KT) ** (1 - 2 * k)) * Lp(1 - 2 * k) * comb_gen(1 - 2 * k, p) * Lp(-p)
            for j in ((1,) if drop0 else (0, 1)):
                v = v + ap(j, 1 - 2 * k) * comb_gen(1 - 2 * k, p) * (-1) ** p * ap(j, -p)
            B[2 * k - 1][p] = flip * BK[k] * v
    A = [e1]                                                                         # exp in z, t-series coefficients
    for j in range(1, J + 1):
        Mj = max(M - 2 * j, 0); acc = [0] * (Mj + 1)
        for i in range(1, j + 1):
            if any(x != 0 for x in B[i]):
                for p in range(Mj + 1):
                    acc[p] = acc[p] + i * sum(B[i][q] * A[j - i][p - q] for q in range(p + 1))
        A.append([x * Fr(1, j) for x in acc])
    return f, A


def rat_ring(L):
    a = [C0 - L, C1 - L]
    return {'Lp': lambda n: Fr(L) ** n, 'ap': lambda j, n: a[j] ** n}


def chart_data(L, M, jac=True, dyy=False):
    b, a0, a1 = KT * L, C0 - L, C1 - L; tay = {}
    for n in range(M + 3):
        for p in range(M + 3 - n):
            if n + p < 2:
                continue
            v = Fr(0)
            if p == 0:
                v += C0 * KT * Fr((-1) ** (n + 1), n) / b ** n
            if n == 0:
                v += (C1 - C0) * Fr((-1) ** (p + 1), p) / a1 ** p
            v += -a0 * comb(KT, n) / b ** n * Fr((-1) ** p) / a1 ** p
            if v:
                tay[(n, p)] = v
    amp = {}
    for n in range(M + 1):
        for p in range(M + 1 - n):
            if jac:                          # 1/x_1 (times 1/y for dy/y)
                v = (Fr((-1) ** p) / a1 ** p) * ((Fr((-1) ** n) / b ** n) if dyy else (1 if n == 0 else 0))
            else:                            # control: 1/(X_0 x_1), X_0 = -q y^k/x_1:  (1 + eta/b)^(-k)
                v = (comb_gen(-KT, n) / b ** n) if p == 0 else 0
            if v:
                amp[((n, p), 0)] = v
    return tay, amp


def ser_ring():
    """L = c_0 + d, a_0 = c_0 - L = -d, a_1 = (c_1 - c_0) - d, as Ser in d"""
    cL, c1 = {}, {}

    def Lp(n):
        if n not in cL:
            cL[n] = ser_geom(C0, 1, n)
        return cL[n]

    def ap(j, n):
        if j == 0:
            return Ser({n: Fr((-1) ** (n % 2))})
        if n not in c1:
            c1[n] = ser_geom(C1 - C0, -1, n)
        return c1[n]
    return {'Lp': Lp, 'ap': ap}


def h_toy(N, flip=1):
    Lg = [Fr(0)] * (N + 1)
    for k in range(1, (N + 1) // 2 + 1):
        Lg[2 * k - 1] = flip * BK[k] * ((Fr(KT) * C0) ** (1 - 2 * k) + (C1 - C0) ** (1 - 2 * k))
    return ser_exp(Lg, N)


def sig_sp(f, A, K, how='wick'):
    if how == 'lb':
        return sp_lb(f, A, K)
    return sp_norm({(m,): f[m] for m in range(2, len(f))},
                   {((p,), j): A[j][p] for j in range(len(A)) for p in range(len(A[j]))}, 1, K)


def run_D(T, **kw):
    Ser.T = T
    f, A = sigma_data(None, 2 * N, N, ser_ring(), **kw)
    return sp_lb(f, A, N)


t0 = time.time()
out('rl_toy: toy chart, k = %d, (c_0, c_1) = (%s, %s); (T), (A), (B) through z^%d; (D) through z^%d (Ser margin %d); log %s'
    % (KT, C0, C1, K, N, PADD, REL))
H = h_toy(N)
out('     h_toy = 1 + %s z + %s z^2 + %s z^3 + ...' % (H[1], H[2], H[3]))
for L in (Fr(3), Fr(-3, 2)):   # generic points (at L = 2 the no-Jacobian control is blind, see REFEREE)
    t1 = time.time()
    NT = sp_norm(*total_data(L), 4, K)
    f, A = sigma_data(L, 2 * K, K, rat_ring(L))
    NA, NAlb = sig_sp(f, A, K), sig_sp(f, A, K, 'lb')
    NB = sp_norm(*chart_data(L, 2 * K), 2, K)
    check('L = %s: (T) total 4-variable SP = (A) sigma-SP of the full-Stirling form (L1b + L3), z^0..z^%d' % (L, K), NT == NA,
          'z^1, z^2 = %s, %s  (%.0f s)' % (NT[1], NT[2], time.time() - t1))
    check('L = %s: (A) by Wick (scaled exp, generating-function moments) = (A) by Lagrange-Buermann' % L, NA == NAlb)
    check('L = %s: (T) = (B) chart SP in (y, x_1), amplitude 1/x_1 (L2 + L1b), z^0..z^%d' % (L, K), NT == NB)
    NTd = sp_norm(*total_data(L, dyy=True), 4, K); fd, Ad = sigma_data(L, 2 * K, K, rat_ring(L), dyy=True)
    NBd = sp_norm(*chart_data(L, 2 * K, dyy=True), 2, K)
    check('L = %s: measure dy/y throughout: (T) = (A) = (B)' % L, NTd == sig_sp(fd, Ad, K) == NBd)
    fn, An = sigma_data(L, 2 * K, K, rat_ring(L), nodelta=True)
    c1 = first_diff(NT, sig_sp(fn, An, K)); c2 = first_diff(NT, sp_norm(*chart_data(L, 2 * K, jac=False), 2, K))
    c3 = first_diff(NT, NBd)
    check('control (must fail): (A) without (Delta_in(sigma)/Delta_in(L))^(-1/2)', c1 is not None, 'first differs at z^%s' % c1)
    check('control (must fail): (B) without the Jacobian X_0', c2 is not None, 'first differs at z^%s' % c2)
    check('control (must fail): (B) with dy/y against (T) with dy', c3 is not None, 'first differs at z^%s' % c3)
t1 = time.time(); NB1 = sp_norm(*chart_data(Fr(1), 2 * N), 2, N)
check('L = c_0 = 1 (q = 0): chart SP = h_toy (step (d): L1c(ii), (G4), L3), z^0..z^%d' % N, NB1 == H, '(%.0f s)' % (time.time() - t1))
t1 = time.time(); D1 = run_D(2 + PADD); D2 = run_D(8 + PADD)
stable = all({p: v for p, v in a.c.items() if p <= 2} == {p: v for p, v in b.c.items() if p <= 2} for a, b in zip(D1, D2))
check('precision guard: (D) coefficients d^(<=2) identical with Ser cut-offs %d and %d' % (2 + PADD, 8 + PADD), stable)
reg = all(c.val() is None or c.val() >= 0 for c in D1); val = [c.coef(0) for c in D1]
check('(D) toy (Gamma): c_n(L) has no pole at L = c_0, n <= %d' % N, reg,
      'valuations in d: %s' % [c.val() for c in D1])
check('(D) toy (Gamma): c_n(c_0) = [z^n] h_toy, n <= %d' % N, val == H, '(%.0f s)' % (time.time() - t1))
Dd = run_D(2 + PADD, dyy=True)
check('(D) with dy/y (chart predicts the same value, St_1 = St_0): regular, value h_toy',
      all(c.val() is None or c.val() >= 0 for c in Dd) and [c.coef(0) for c in Dd] == H)
for name, kw in [('b_k -> -b_k', {'flip': -1}), ('colour-0 Stirling terms dropped', {'drop0': True})]:
    Dc = run_D(2 + PADD, **kw)
    regc = all(c.val() is None or c.val() >= 0 for c in Dc)
    check('control (must fail): (D) with %s' % name, not (regc and [c.coef(0) for c in Dc] == H),
          'regular: %s; first pole/mismatch at n = %s' % (regc, next((n for n, c in enumerate(Dc) if (c.val() or 0) < 0
                                                                       or c.coef(0) != H[n]), None)))
out('RL_TOY %d %d %d: %s   (%.0f s)' % (K, N, PADD, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - t0))
LOG.close()
