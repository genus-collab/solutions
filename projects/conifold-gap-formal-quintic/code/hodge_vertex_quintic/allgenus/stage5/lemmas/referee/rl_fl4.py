# -*- coding: utf-8 -*-
"""rl_fl4.py -- referee checks of FL4 (FORMAL_SP_LEMMAS s.4.5) and of Lemma T(b) with p = +z (the shift sigma -> sigma + z of
Theorem P(a)) on Mellin-Barnes-type data: the toy sigma-phase of rl_toy (k = 2, colours 1, -2) with its full-Stirling amplitude
A_St, which is q-independent; the q-dependence is lam -> lam + eps (eps = log q - log q_0, so d/deps = D_q = q d/dq).
New code, stdlib only; uses rl_core only.  Exact over Q; eps-series through eps^T.
 (a) FL4 iterated: (z d/deps)^e S_0 = S_e for e = 1, 2, where
     S_e(eps) := e^((F_eps(c(eps)) - F_0(c(0)))/z) (Delta(eps)/Delta(0))^(-1/2) SP^norm_(c(eps))[F_eps; sigma^e A_St],
     F_eps = (lam + eps) sigma + (S_0 - lam sigma), c(eps) the MOVING critical point (formal IFT over Q[[eps]], own code).
 (b) Lemma T(b), p = +z: SP^norm_L[S_0; e^((S_0(L+t+z) - S_0(L+t))/z) A_St(L+t+z, z)] = SP^norm_L[S_0; A_St].
Controls (must fail): (a) without e^((F_eps(c) - F_0(c_0))/z); without the Hessian ratio; P(sigma) = sigma + z for e = 1;
(b) without the exponential factor.
Usage: py -I -B rl_fl4.py [K] [T]   (z^K, eps^T; defaults 6, 6)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from rl_core import Ser, comb_gen, sp_norm, sp_scaled, bk_list, inv_mat

K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
T = int(sys.argv[2]) if len(sys.argv) > 2 else 6
KT, C0, C1, L = 2, Fr(1), Fr(-2), Fr(3)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
REL = 'logs/rl_fl4_%d_%d.log' % (K, T); LOG = open(os.path.join(HERE, REL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + ('   ' + info if info else ''))
    if not ok:
        FAILS.append(name)


BK = bk_list(K + T + 4); MM = 2 * K + T + 4               # Taylor orders needed at the moving point
a = [C0 - L, C1 - L]
f = [Fr(0), Fr(0)] + [(KT * (-1) ** m * L ** (1 - m) + a[0] ** (1 - m) + a[1] ** (1 - m)) / (m * (m - 1)) for m in range(2, MM + 3)]
e1 = [comb_gen(Fr(1, 2), p) / L ** p for p in range(MM + 1)]
for j in range(2):
    g = [comb_gen(Fr(-1, 2), p) * (-1) ** p / a[j] ** p for p in range(MM + 1)]
    e1 = [sum(e1[i] * g[p - i] for i in range(p + 1)) for p in range(MM + 1)]
Bz = [[Fr(0)] * (MM + 1) for _ in range(K + 1)]
for k in range(1, (K + 1) // 2 + 1):
    for p in range(MM + 1):
        v = Fr(KT) ** (1 - 2 * k) * L ** (1 - 2 * k) * comb_gen(1 - 2 * k, p) / L ** p
        for j in range(2):
            v += a[j] ** (1 - 2 * k) * comb_gen(1 - 2 * k, p) * (-1) ** p / a[j] ** p
        Bz[2 * k - 1][p] = BK[k] * v
A = [e1]
for j in range(1, K + 1):
    acc = [Fr(0)] * (MM + 1)
    for i in range(1, j + 1):
        for p in range(MM + 1):
            acc[p] += i * sum(Bz[i][q] * A[j - i][p - q] for q in range(p + 1))
    A.append([x / j for x in acc])


def moving_sp(e_pow, with_exp=True, with_hess=True, zshift=False):
    """S_e(eps) as {z power: Ser in eps}"""
    Ser.T = T
    eps = Ser({1: 1}); p = Ser({})
    for _ in range(T + 2):                                   # p = -(eps + sum_{m>=3} m f_m p^(m-1)) / (2 f_2)
        s = eps
        for m in range(3, T + 3):
            s = s + m * f[m] * p ** (m - 1)
        p = s * (-1 / (2 * f[2]))
    pp = [p ** r for r in range(MM + 1)]
    tay = {}
    for m in range(2, 2 * K + 3):
        tay[(m,)] = sum((comb(m + r, r) * f[m + r] * pp[r] for r in range(T + 1)), Ser({}))
    sig = [L + p, Ser.const(1)]                              # sigma = (L + p) + t
    amp = {}
    for j in range(K + 1):
        Aj = [sum((comb(m + r, r) * A[j][m + r] * pp[r] for r in range(T + 1) if m + r < len(A[j])), Ser({}))
              for m in range(2 * K + 1)]
        for _ in range(e_pow):
            Aj = [sig[0] * Aj[m] + (Aj[m - 1] if m else 0) for m in range(2 * K + 1)]
        for m in range(2 * K + 1):
            amp[((m,), j)] = Aj[m]
    if zshift:                                               # control: P(sigma) = sigma + z (e = 1): add z A
        for j in range(K):
            for m in range(2 * K + 1):
                key = ((m,), j + 1); amp[key] = amp.get(key, Ser({})) + sum(
                    (comb(m + r, r) * A[j][m + r] * pp[r] for r in range(T + 1) if m + r < len(A[j])), Ser({}))
    Nk = sp_norm(tay, amp, 1, K)
    X = eps * (L + p)
    for m in range(2, T + 3):
        X = X + f[m] * pp[m]
    ratio = tay[(2,)] * (1 / f[2])
    hr = sum((comb_gen(Fr(-1, 2), d) * (ratio - 1) ** d for d in range(T + 1)), Ser({})) if with_hess else Ser.const(1)
    S = {}
    term = Ser.const(1)
    for n_ in range(T + 1 if with_exp else 1):               # e^(X/z) = sum_n X^n / (n! z^n)
        for k, Nv in enumerate(Nk):
            S[k - n_] = S.get(k - n_, Ser({})) + term * Nv * hr
        term = term * X * Fr(1, n_ + 1)
    return S


def zd(S):
    """z d/deps"""
    return {k + 1: v.deriv() for k, v in S.items()}


def agree(S1, S2, emax, kmax):
    """compare the coefficients z^k eps^d with d <= emax and k + d <= kmax (the complete region: N_k' is known for k' <= K)"""
    keys = set(S1) | set(S2)
    return all(S1.get(k, Ser({})).coef(d) == S2.get(k, Ser({})).coef(d) for k in keys for d in range(emax + 1) if k + d <= kmax)


t0 = time.time()
out('rl_fl4: FL4 and T(b) (p = +z) on the toy Mellin-Barnes sigma-phase at L = %s, z^%d, eps^%d; log %s' % (L, K, T, REL))
S0, S1, S2 = moving_sp(0), moving_sp(1), moving_sp(2)
check('FL4: z d/deps S_0 = S_1 (= SP of sigma A_St at the moving critical point), all z^k eps^d with k + d <= %d, d <= %d' % (K, T - 1),
      agree(zd(S0), S1, T - 1, K))
check('FL4 iterated: (z d/deps)^2 S_0 = S_2 (P(zD_q) SP = SP[P(sigma) ..], P = sigma^2), all z^k eps^d with k + d <= %d, d <= %d' % (K, T - 2),
      agree(zd(zd(S0)), S2, T - 2, K))
check('control (must fail): FL4 without e^((F_eps(c) - F_0(c_0))/z)', not agree(zd(moving_sp(0, with_exp=False)),
      moving_sp(1, with_exp=False), T - 1, K))
check('control (must fail): FL4 without (Delta(eps)/Delta(0))^(-1/2)', not agree(zd(moving_sp(0, with_hess=False)),
      moving_sp(1, with_hess=False), T - 1, K))
check('control (must fail): P(sigma) = sigma + z claimed for z d/deps', not agree(zd(S0), moving_sp(1, zshift=True), T - 1, K))
# (b) Lemma T(b) with p = +z ----------------------------------------------------------------------------------------------
M2 = 2 * K
Q = [dict() for _ in range(M2 + 1)]
for m in range(3, M2 + 3):                                   # F_3(s eta)/s^2
    Q[m - 2][(m,)] = Q[m - 2].get((m,), 0) + f[m]
QE = [dict() for _ in range(M2 + 1)]                         # E = (S_0(L+t+z) - S_0(L+t))/z, scaled
for m in range(2, M2 + 2):
    for i in range(1, m + 1):
        o = m - 2 + i
        if o <= M2:
            QE[o][(m - i,)] = QE[o].get((m - i,), 0) + f[m] * comb(m, i)
Ash = [dict() for _ in range(M2 + 1)]                        # A_St(L + t + z, z), scaled
for j in range(K + 1):
    for p_ in range(M2 + 1):
        for i in range(p_ + 1):
            o = 2 * j + p_ + i
            if o <= M2:
                Ash[o][(p_ - i,)] = Ash[o].get((p_ - i,), 0) + A[j][p_] * comb(p_, i)
Qtot = [{**{k: v for k, v in Q[o].items()}} for o in range(M2 + 1)]
for o in range(M2 + 1):
    for key, v in QE[o].items():
        Qtot[o][key] = Qtot[o].get(key, 0) + v
Cv = [[-1 / (2 * f[2])]]
lhs = sp_scaled(1, Qtot, Ash, Cv, K)
rhs = sp_norm({(m,): f[m] for m in range(2, M2 + 3)}, {((p_,), j): A[j][p_] for j in range(K + 1) for p_ in range(M2 + 1)}, 1, K)
check('T(b), p = +z (Theorem P(a) shift): SP^norm_L[S_0; e^((S_0(L+t+z) - S_0(L+t))/z) A_St(L+t+z)] = SP^norm_L[S_0; A_St], z^0..z^%d'
      % K, lhs == rhs, 'z^1 = %s' % rhs[1])
bad = sp_scaled(1, Q, Ash, Cv, K)
check('control (must fail): T(b) without the exponential factor', bad != rhs,
      'first differs at z^%s' % next((i for i in range(K + 1) if bad[i] != rhs[i]), None))
out('RL_FL4 %d %d: %s   (%.0f s)' % (K, T, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - t0))
LOG.close()
