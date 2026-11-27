# -*- coding: utf-8 -*-
"""rl_l3.py -- referee checks of Lemma L3 (FORMAL_SP_LEMMAS s.9) and of the steps of its formal proof (Proof A), exact over Q.
New code, stdlib only; uses rl_core only.
 (a) L3(i): SP^norm_a[a log x - x; (x/a)^(m-1)] = St_m(z/a), St_m by its DEFINITION (s.9.1), the SP by Lagrange-Buermann
     (Morse coordinate; no Wick contraction), m = -6..6, through z^N at a = 1; by Wick (rl_core) through z^NA at a = 3/7, -5/2.
 (b) L3(iv): rho(Y) = 0 in Q[[1/Y]] through Y^-(N+2), and the three intermediate identities of s.9.2; St_(m+1) = (1+mw) St_m.
 (c) Step A2 = L1c(iv) (T-z) in the L3 instance: S(eps; z) computed by Wick at the MOVING critical point 1 + eps over Q[[eps]],
     times e^(((1+eps)log(1+eps) - eps)/z) (1+eps)^(1/2), then eps -> z, equals SP^norm_1[F; 1] (no closed form used).
 (d) L1a with the non-polynomial map x = e^u at x = 1: SP_u[e^(m u)] = SP_x[x^(m-1)] (the Jacobian e^u included).
Controls (must fail): b_3 perturbed in rho; St_m used for m+1; (c) without (1+eps)^(1/2), and with eps -> 2z; (d) without
the Jacobian.
Usage: py -I -B rl_l3.py [N] [NA] [M]   (defaults 40, 16, 8)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial
from rl_core import Ser, comb_gen, sp_norm, sp_lb, bk_list, bernoulli, stirling_m, ser_exp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NA = int(sys.argv[2]) if len(sys.argv) > 2 else 16
M = int(sys.argv[3]) if len(sys.argv) > 3 else 8
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
REL = 'logs/rl_l3_%d_%d_%d.log' % (N, NA, M); LOG = open(os.path.join(HERE, REL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + ('   ' + info if info else ''))
    if not ok:
        FAILS.append(name)


def fd(a, b):
    return next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), None)


t0 = time.time(); BK = bk_list(N + 4)
out('rl_l3: Lemma L3 and Proof A; LB through z^%d, Wick at a = 3/7, -5/2 through z^%d, T-z through z^%d; log %s' % (N, NA, M, REL))
# (a) -------------------------------------------------------------------------------------------------------------------
ST = {m: stirling_m(m, N, BK) for m in range(-6, 8)}
f1 = [Fr(0), Fr(0)] + [Fr((-1) ** (l + 1), l) for l in range(2, 2 * N + 3)]          # log(1+t) - t at a = 1
for m in range(-6, 7):
    g = [[comb_gen(m - 1, p) for p in range(2 * N + 1)]]                               # (1 + t)^(m-1)
    lb = sp_lb(f1, g, N); ok = lb == ST[m][:N + 1]
    check('L3(i) a = 1, m = %d (dx/x with x^m; = dx with x^(m-1)): Lagrange-Buermann SP = St_%d by definition, z^0..z^%d'
          % (m, m, N), ok, '' if ok else 'first diff z^%s' % fd(lb, ST[m]))
for a in (Fr(3, 7), Fr(-5, 2)):
    tay = {(l,): Fr((-1) ** (l + 1), l) / a ** (l - 1) for l in range(2, 2 * NA + 3)}
    for m in (-4, 0, 3):
        amp = {((p,), 0): comb_gen(m - 1, p) / a ** p for p in range(2 * NA + 1)}
        w = sp_norm(tay, amp, 1, NA); ref = [ST[m][k] / a ** k for k in range(NA + 1)]
        check('L3(i) a = %s, m = %d: Wick SP = St_m(z/a) (homogeneity (H)), z^0..z^%d' % (a, m, NA), w == ref)
c = fd(sp_lb(f1, [[comb_gen(1, p) for p in range(2 * N + 1)]], N), ST[1])
check('control (must fail): St_1 claimed for the amplitude x^1 (m = 2)', c is not None, 'first differs at z^%s' % c)
# (b) -------------------------------------------------------------------------------------------------------------------
R = N + 2                                                     # series in u = 1/Y through u^R; f[i] = [u^i]


def mu_shift(bk, shift):
    """mu(Y + shift) = sum_k b_k (Y + shift)^(1-2k) as a series in u = 1/Y"""
    s = [Fr(0)] * (R + 1)
    for k in range(1, R // 2 + 2):
        for i in range(R + 1):
            e = 2 * k - 1 + i
            if e <= R:
                s[e] += bk[k] * comb_gen(1 - 2 * k, i) * Fr(shift) ** i
    return s


def rho(bk):
    lg = [Fr(0)] + [Fr((-1) ** (j + 1), j) for j in range(1, R + 2)]                   # log(1 + u)
    r = [lg[i + 1] + Fr(1, 2) * lg[i] for i in range(R + 1)]                            # (Y + 1/2) log(1 + 1/Y)
    r[0] -= 1
    m1, m0 = mu_shift(bk, 1), mu_shift(bk, 0)
    return [r[i] + m1[i] - m0[i] for i in range(R + 1)]


check('L3(iv): rho(Y) = (Y+1/2) log(1+1/Y) - 1 + mu(Y+1) - mu(Y) = 0 through Y^-%d' % R, all(x == 0 for x in rho(BK)))
bad = list(BK); bad[3] = bad[3] + Fr(1, 1000)
c = next((i for i, x in enumerate(rho(bad)) if x != 0), None)
check('control (must fail): b_3 + 1/1000 in rho', c is not None, 'first non-zero at Y^-%s' % c)
B = bernoulli(R + 2)
D2mu = [Fr(0)] * (R + 3)                                       # D^2 mu = sum_k B_2k Y^(-2k-1)
for k in range(1, R // 2 + 1):
    D2mu[2 * k + 1] = B[2 * k]
Tser = [Fr(0)] * (R + 3)                                       # T = sum_n (B_n/n!) D^n (1/Y) = sum_n B_n (-1)^n Y^(-n-1)
for n_ in range(R + 2):
    Tser[n_ + 1] = B[n_] * (-1) ** n_
ok1 = all(D2mu[i] == Tser[i] - (1 if i == 1 else 0) - (Fr(1, 2) if i == 2 else 0) for i in range(R + 2))
check('s.9.2: D^2 mu = T - 1/Y - 1/(2Y^2) through Y^-%d (B_1 = -1/2 convention)' % (R + 1), ok1)


def shift1(f):
    """f(Y+1) for f = sum f_i Y^-i (series in 1/Y): (Y+1)^-i = Y^-i (1+1/Y)^-i"""
    g = [Fr(0)] * len(f)
    for i, v in enumerate(f):
        if v:
            for j in range(len(f) - i):
                g[i + j] += v * comb_gen(-i, j)
    return g


lhs = [x - y for x, y in zip(shift1(D2mu[:R + 1]), D2mu[:R + 1])]
rhs = [Fr(0)] * (R + 1)
for j in range(R - 3):                                         # -1/(2 Y^2 (Y+1)^2) = -(1/2) Y^-4 (1 + 1/Y)^-2
    rhs[4 + j] = -Fr(1, 2) * comb_gen(-2, j)
check('s.9.2: (e^D - 1) D^2 mu = -1/(2 Y^2 (Y+1)^2) through Y^-%d' % R, lhs == rhs)
q = [Fr(0)] + [Fr((-1) ** (j + 1), j) for j in range(1, R + 3)]                          # l = log(1 + u)
Yl = [q[i + 1] + Fr(1, 2) * q[i] for i in range(R + 1)]                                  # (Y + 1/2) l in powers of 1/Y
D2 = [Fr(0)] * (R + 1)
for i in range(1, R - 1):                                                                # D^2 Y^-i = i(i+1) Y^-(i+2)
    D2[i + 2] = Yl[i] * i * (i + 1)
check('s.9.2: D^2[(Y + 1/2) log(1 + 1/Y)] = +1/(2 Y^2 (Y+1)^2) through Y^-%d' % R, D2 == [-x for x in rhs])
ok = all([ST[m][k] + (m * ST[m][k - 1] if k else 0) for k in range(N + 1)] == ST[m + 1][:N + 1] for m in range(-6, 7))
check('L3(iv): St_(m+1)(w) = (1 + m w) St_m(w) for m = -6..6 through w^%d (St_m by definition)' % N, ok)
# (c) -------------------------------------------------------------------------------------------------------------------
Ser.T = 2 * M                                                  # eps-series through eps^(2M): enough, since weight >= d/2
one_e = Ser({0: 1, 1: 1})                                      # 1 + eps
inv1e = one_e.inv()
tay = {}
for l in range(2, 2 * M + 3):
    tay[(l,)] = Fr((-1) ** (l + 1), l) * (inv1e ** (l - 1))     # [t^l] ((1+eps) log x - x) at x = 1 + eps
amp = {((p,), 0): inv1e * (inv1e ** p) * (-1) ** p for p in range(2 * M + 1)}   # (x/a)^(-1) = (1 + eps + t)^(-1), a = 1
Nk = sp_norm(tay, amp, 1, M)                                   # SP^norm at the moving critical point, coefficients in Q[[eps]]
E = {d: Fr((-1) ** d, d * (d - 1)) for d in range(2, 2 * M + 1)}   # (1+eps) log(1+eps) - eps
expo = {(0, 0): Fr(1)}; term = {(0, 0): Fr(1)}                 # e^(E/z) as {(eps power, z power): c}
for n_ in range(1, M + 1):
    new = {}
    for (d, k), v in term.items():
        for d2, v2 in E.items():
            if d + d2 <= 2 * M:
                new[(d + d2, k - 1)] = new.get((d + d2, k - 1), 0) + v * v2 / n_
    term = new
    for key, v in term.items():
        expo[key] = expo.get(key, 0) + v
half = {d: comb_gen(Fr(1, 2), d) for d in range(2 * M + 1)}    # (1+eps)^(1/2) = (Delta_eps/Delta)^(-1/2)


def assemble(with_half=True, factor=1):
    S = {}
    for (d1, k1), v1 in expo.items():
        for k2, ser in enumerate(Nk):
            for d2, v2 in ser.c.items():
                for d3 in (half if with_half else {0: 1}):
                    d = d1 + d2 + d3
                    if d <= 2 * M:
                        hv = half[d3] if with_half else 1
                        S[(d, k1 + k2)] = S.get((d, k1 + k2), 0) + v1 * v2 * hv
    sub = [Fr(0)] * (M + 1)                                    # eps -> factor * z
    for (d, k), v in S.items():
        if 0 <= d + k <= M:
            sub[d + k] += v * Fr(factor) ** d
    return sub


Ser.T = 0
direct = sp_norm({(l,): Fr((-1) ** (l + 1), l) for l in range(2, 2 * M + 3)}, {((0,), 0): Fr(1)}, 1, M)  # SP^norm_1[F; 1]
got = assemble()
check('L1c(iv)/T-z in Step A2: [e^(..) (1+eps)^(1/2) SP^norm_(1+eps)[F_eps; (x/a)^-1]]_(eps=z) = SP^norm_1[F; 1], z^0..z^%d' % M,
      got == direct, 'z^1..z^3 = %s' % [str(x) for x in got[1:4]])
check('   and SP^norm_1[F; 1] = St_1 = St_0 (L3 with m = 1)', direct == ST[1][:M + 1] == ST[0][:M + 1])
c1 = fd(assemble(with_half=False), direct); c2 = fd(assemble(factor=2), direct)
check('control (must fail): T-z without (Delta_eps/Delta)^(-1/2)', c1 is not None, 'first differs at z^%s' % c1)
check('control (must fail): T-z with eps -> 2z', c2 is not None, 'first differs at z^%s' % c2)
# (d) -------------------------------------------------------------------------------------------------------------------
tayu = {(l,): -Fr(1, factorial(l)) for l in range(2, 2 * NA + 3)}   # u - e^u + 1 = -sum_{l>=2} u^l/l!
for m in (-3, -1, 0, 2, 3):
    w = sp_norm(tayu, {((p,), 0): Fr(m) ** p / factorial(p) for p in range(2 * NA + 1)}, 1, NA)    # amplitude e^(m u)
    check('L1a, map x = e^u: SP_u[e^(m u)] = SP_x[x^(m-1)] = St_%d, m = %d, z^0..z^%d' % (m, m, NA), w == ST[m][:NA + 1])
wb = sp_norm(tayu, {((p,), 0): Fr(1) / factorial(p) for p in range(2 * NA + 1)}, 1, NA)   # m = 2 without the Jacobian e^u
c = fd(wb, ST[2])
check('control (must fail): L1a with x = e^u, m = 2, Jacobian e^u omitted', c is not None, 'first differs at z^%s' % c)
out('     St_0 = 1 + %s w + %s w^2 + %s w^3 + ...;  St_-6: w^1 = %s' % (ST[0][1], ST[0][2], ST[0][3], ST[-6][1]))
out('RL_L3 %d %d %d: %s   (%.0f s)' % (N, NA, M, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - t0))
LOG.close()
