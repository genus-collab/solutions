# g4_w11.py -- evidence for gap G4 items 2 and 4 (exact over Q and Q(zeta5)). Usage: py -I -B g4_w11.py
# Model curve x = W_v = p^4/4 + v/p, y = p, B = dp1 dp2/(p1-p2)^2, engine kernel (1/2) Int_{sigma q}^q B/((y(q)-y(sigma q)) dx(q)).
# (A) local formula: omega_11(p0) = -Sum_a (1/(4 alpha y1)) [B0 beta00 + B2/12 - B0 y3/(4 y1)], x - x_a = alpha zeta^2;
#     <-log p>_1 = Int_ray omega_11 = 5 x (contribution at c_0) at v = 1; controls.
# (B) the local formula (all five branch points) against the closed form -(5/8) p^4 (p^10 + 18 v p^5 + 6 v^2)/(p^5 - v)^4 dp.
# (C) item 2: Res_{c_b} x^k omega_11 = 0 (each b), controls Res_{c_b} p^k omega_11, Res_0 x^5 omega_11 (circle vs residues).
# (D) <-log p>_1 from the closed form by u = p^5 (exact antiderivative) and homogeneity v -> 32.
import os, sys
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4_field import Qz, ZERO, ONE, smul, sinv, ssqrt1, srevert


def qpow(a, k):
    r = ONE
    for _ in range(k): r = r * a
    return r


def local(c, v, N=8):
    """t(zeta) = p - c with x - x(c) = alpha zeta^2, alpha = x''(c)/2; returns alpha, t-series"""
    xk = [ZERO] * (N + 3)
    for k in range(2, N + 3):
        a = Qz.of(Fr(comb(4, k), 4)) * qpow(c, 4 - k) if k <= 4 else ZERO
        b = (Qz.of(v) / c) * Qz.of((-1) ** k) * qpow(ONE / c, k)
        xk[k] = a + b
    alpha = xk[2]
    u = [ONE] + [xk[m + 2] / alpha for m in range(1, N + 1)]
    zeta_of_t = [ZERO] + ssqrt1(u, N)[:N]
    return alpha, srevert(zeta_of_t, N)


def contrib_log(c, v, orient=-1, drop_beta=False, drop_y3=False):
    """contribution of branch point c to Int_{0}^{oo} omega_11 (log leg); orient=-1 means Int_ray B = -dq/q"""
    alpha, t = local(c, v)
    y1, y3 = t[1], t[3]
    beta = (t[3] * 6 / t[1] - Fr(3, 2) * (t[2] * 2 / t[1]) * (t[2] * 2 / t[1])) * Fr(1, 6)
    if drop_beta: beta = ZERO
    if drop_y3: y3 = ZERO
    dt = [t[k + 1] * (k + 1) for k in range(len(t) - 1)]
    q = [c + t[0]] + t[1:]
    L = smul(dt, sinv(q, 6), 6)                       # dq/q = t'(zeta) dzeta/(c + t)
    L = [x * orient for x in L]
    return -(ONE / (alpha * y1 * 4)) * (L[0] * beta + L[2] * Fr(1, 12) - L[0] * y3 / (y1 * 4))


def omega11_local(p0, v, sig):
    """omega_11(p0)/dp0 from the local formula, summed over the five branch points"""
    tot = ZERO
    for b in range(5):
        c = Qz.zeta(-b) * sig
        alpha, t = local(c, v)
        y1, y3 = t[1], t[3]
        beta = (t[3] * 6 / t[1] - Fr(3, 2) * (t[2] * 2 / t[1]) * (t[2] * 2 / t[1])) * Fr(1, 6)
        D = p0 - c
        dt = [t[k + 1] * (k + 1) for k in range(len(t) - 1)]
        den = [D - t[0]] + [-x for x in t[1:]]
        inv = sinv(den, 6); inv2 = smul(inv, inv, 6)
        Bs = smul(dt, inv2, 6)                        # B(p0, q(zeta))/(dp0 dzeta) = t'(zeta)/(p0 - c - t)^2
        tot = tot - (ONE / (alpha * y1 * 4)) * (Bs[0] * beta + Bs[2] * Fr(1, 12) - Bs[0] * y3 / (y1 * 4))
    return tot


def closed(p, v):
    return Qz.of(Fr(-5, 8)) * qpow(p, 4) * (qpow(p, 10) + qpow(p, 5) * (18 * v) + Qz.of(6 * v * v)) / qpow(qpow(p, 5) - v, 4)


def res_at(c, v, fpoly, N=8):
    """Res_{p=c} fpoly(p) * omega_11(p), omega_11 the closed form; fpoly: function p-series -> p-series (in t = p - c)"""
    # p^5 - v = t g(t), g = sum_{k>=1} C(5,k) c^(5-k) t^(k-1)
    g = [Qz.of(comb(5, k)) * qpow(c, 5 - k) for k in range(1, 6)] + [ZERO] * N
    g = g[:N + 1]
    ginv = sinv(g, N); g4 = smul(smul(ginv, ginv, N), smul(ginv, ginv, N), N)
    p = [c, ONE] + [ZERO] * (N - 1)
    pk = lambda k: [qpow(c, k - m) * comb(k, m) if m <= k else ZERO for m in range(N + 1)]
    num = [x * Fr(-5, 8) for x in pk(4)]
    poly = [a + b * (18 * v) + (Qz.of(6 * v * v) if m == 0 else ZERO) for m, (a, b) in enumerate(zip(pk(10), pk(5)))]
    f = smul(smul(num, poly, N), g4, N)
    f = smul(f, fpoly(p, N), N)
    return f[3]                                     # t^-4 * f(t): residue = [t^3] f


def xser(p, N):
    pinv = sinv(p, N)
    return [a * Fr(1, 4) + b for a, b in zip(smul(smul(p, p, N), smul(p, p, N), N), [x for x in pinv])]


def main():
    ok = True
    v = 1
    tot = ZERO; per = []
    for b in range(5):
        cb = contrib_log(Qz.zeta(-b), v); per.append(cb); tot = tot + cb
    print('(A) contributions per branch point at v = 1:', [str(x.rational()) for x in per], ' total <-log p>_1 =', tot.rational())
    ok &= tot == Qz.of(Fr(-1, 12)) and all(x == Qz.of(Fr(-1, 60)) for x in per)
    c1 = sum((contrib_log(Qz.zeta(-b), v, orient=+1) for b in range(5)), ZERO)
    c2 = sum((contrib_log(Qz.zeta(-b), v, drop_beta=True) for b in range(5)), ZERO)
    c3 = sum((contrib_log(Qz.zeta(-b), v, drop_y3=True) for b in range(5)), ZERO)
    print('    controls: opposite ray orientation', c1.rational(), '; beta00 dropped', c2.rational(), '; y3 dropped', c3.rational())
    ok &= c1 == Qz.of(Fr(1, 12)) and not (c2 == tot) and not (c3 == tot)
    bad = 0
    for v, sig in ((1, 1), (32, 2)):
        for p0 in (Fr(3), Fr(-2), Fr(1, 2), Fr(7, 3)):
            a = omega11_local(Qz.of(p0), v, sig); b = closed(Qz.of(p0), v)
            if not (a == b): bad += 1
    print('(B) local formula (sum over five branch points) = closed form at 8 points (v = 1, 32): mismatches', bad)
    ok &= bad == 0
    v = 1; zx = []; zp = []
    for b in range(5):
        c = Qz.zeta(-b)
        for k in range(0, 7):
            def fx(p, N, k=k):
                r = [ONE] + [ZERO] * N; xs = xser(p, N)
                for _ in range(k): r = smul(r, xs, N)
                return r
            if not res_at(c, v, fx).iszero(): zx.append((b, k))
        for k in range(1, 10):
            def fp(p, N, k=k):
                r = [ONE] + [ZERO] * N
                for _ in range(k): r = smul(r, p, N)
                return r
            if not res_at(c, v, fp).iszero(): zp.append((b, k))
    print('(C) Res_{c_b} x^k omega_11 != 0 for (b,k), k <= 6:', zx, '; control Res_{c_b} p^k omega_11 != 0 for', len(zp), 'of 45 (b,k)')
    ok &= zx == [] and len(zp) > 0
    # Res_{p=0} x^5 omega_11: x^5 = p^-5 (1 + p^5/4)^5, omega_11 = -(5/8) p^4 (6 v^2 + ...)/(-v)^4 (1 + O(p^5)) dp
    r0 = Fr(-5, 8) * 6  # [p^4] of omega_11/dp at v = 1 is -(5/8)*6; times [p^-5] of x^5 = 1
    print('    Res_{p=0} x^5 omega_11 (v = 1) =', r0, '(non-zero: the circle pairing of x^5 is not zero)')
    ok &= r0 != 0
    # (D) Int_0^oo omega_11 = -(1/8) Int (u^2 + 18 v u + 6 v^2)/(u - v)^4 du, w = u - v: w^-2 + 20 v w^-3 + 25 v^2 w^-4
    for vv in (Fr(1), Fr(32), Fr(-3)):
        Fw = lambda w: -1 / w - 10 * vv / w ** 2 - Fr(25, 3) * vv ** 2 / w ** 3
        val = Fr(-1, 8) * (0 - Fw(-vv))
        print('(D) closed form: <-log p>_1 at v = %s:' % vv, val, ' (-1/(12 v) =', Fr(-1, 12) / vv, ')')
        ok &= val == Fr(-1, 12) / vv
    print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')


main()
