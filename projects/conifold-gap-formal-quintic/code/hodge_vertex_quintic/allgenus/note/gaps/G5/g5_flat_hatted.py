# g5_flat_hatted.py -- cluster G5, item 1 (thm:con-Pa step 5): evidence only, exact over Q.
# The FLAT full-Stirling saddle series c(z; L) (normalised SP at sigma = L of e^(S/z) A_St dsigma, S'' = -5/(sigma(sigma^5 - 1)),
# A_St = (5 sigma/(sigma^5-1))^(1/2) exp(sum_k z^(2k-1) bk_k [(5 sigma)^(1-2k) - sum_j (sigma - zeta^j)^(1-2k)]), divided by its value at
# t = 0, z = 0) is computed by one-variable Wick contraction at rational points L0, and compared with L0^(4n) ph_n(L0^-5), ph_n the
# polynomial solution of HPF (POLY s.1) with the free coefficient fixed once from c_5.  Expected (thm:con-Pa (a2), prop:gam-hatflat):
#   c_n(L0) = L0^(4n) ph^SP_n(L0^-5) for every n, and c_5(L0) - L0^20 ph^(0)_5(L0^-5) = f_1 = B_6(1/5)/6 = 4537/3937500 at every L0.
# Controls (must fail): (X1) no half-density factor; (X2) twist term with the wrong sign; (X3) scaling L0^(4n-1) instead of L0^(4n).
# Usage: py -I -B g5_flat_hatted.py [NMAX]
import sys
from fractions import Fraction as Fr
from math import comb, factorial

def bernoulli(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B

def smul(a, b, T):
    r = [Fr(0)] * (T + 1)
    for i, x in enumerate(a[:T + 1]):
        if x:
            for j, y in enumerate(b[:T + 1 - i]):
                if y:
                    r[i + j] += x * y
    return r

def sinv(a, T):
    r = [Fr(0)] * (T + 1); r[0] = 1 / a[0]
    for n in range(1, T + 1):
        s = sum(a[k] * r[n - k] for k in range(1, min(n, len(a) - 1) + 1))
        r[n] = -s / a[0]
    return r

def slog1(a, T):                      # log of a series with constant term 1
    da = [k * a[k] for k in range(1, len(a))] + [Fr(0)]
    q = smul(da, sinv(a, T), T)
    return [Fr(0)] + [q[k - 1] / k for k in range(1, T + 1)]

def sexp0(a, T):                      # exp of a series with constant term 0
    r = [Fr(0)] * (T + 1); r[0] = Fr(1)
    da = [k * a[k] for k in range(len(a))]
    for n in range(1, T + 1):
        r[n] = sum(da[k] * r[n - k] for k in range(1, n + 1)) / n
    return r

def poly_shift(coeffs, L0, T):        # polynomial p(sigma) -> series in t at sigma = L0 + t
    r = [Fr(0)] * (T + 1)
    for d, c in enumerate(coeffs):
        for i in range(0, min(d, T) + 1):
            r[i] += c * comb(d, i) * L0 ** (d - i)
    return r

def flat_c(L0, NMAX, variant=''):
    T = 6 * NMAX + 2
    TT = T + 2 * NMAX + 2
    sig = poly_shift([0, 1], L0, TT)                                   # sigma = L0 + t
    s5m1 = poly_shift([-1, 0, 0, 0, 0, 1], L0, TT)                       # sigma^5 - 1
    h = [-5 * x for x in sinv(smul(sig, s5m1, TT), TT)]                 # S''(sigma)
    H = h[0]; C = -1 / H
    a = [Fr(0)] * (T + 1)
    for k in range(3, T + 1):
        a[k] = h[k - 2] / (k * (k - 1))                                  # F_3 coefficients
    # half-density (sigma/L0)(L0^5 - 1)/(sigma^5 - 1), square root
    r = smul([x / L0 for x in sig], [x * (L0 ** 5 - 1) for x in sinv(s5m1, TT)], T)
    hd = sexp0([x / 2 for x in slog1(r, T)], T) if variant != 'X1' else [Fr(1)] + [Fr(0)] * T
    # Stirling exponent: c_i(t) for z^i, i = 2k-1 <= NMAX
    B = bernoulli(2 * NMAX + 2)
    g = smul([5 * x for x in smul(smul(sig, sig, TT), smul(sig, sig, TT), TT)], sinv(s5m1, TT), TT)   # sum_j 1/(sigma - zeta^j)
    expo = {}
    for k in range(1, NMAX + 1):
        i = 2 * k - 1
        if i > NMAX:
            break
        n = 2 * k - 1
        bk = B[2 * k] / (2 * k * (2 * k - 1))
        sj = [Fr((-1) ** (n - 1), factorial(n - 1)) * factorial(m + n - 1) / factorial(m) * g[m + n - 1] for m in range(T + 1)]
        inv = sinv(poly_shift([0, 5], L0, TT), T)                         # 1/(5 sigma)
        pw = [Fr(1)] + [Fr(0)] * T
        for _ in range(n):
            pw = smul(pw, inv, T)
        sign = -1 if variant == 'X2' else 1
        expo[i] = [bk * (sign * pw[m] - sj[m]) for m in range(T + 1)]
    # E(t, z) = exp(sum_i z^i expo_i): E_r by the recursion r E_r = sum_i i c_i E_(r-i)
    E = [[Fr(1)] + [Fr(0)] * T]
    for rr in range(1, NMAX + 1):
        acc = [Fr(0)] * (T + 1)
        for i, ci in expo.items():
            if i <= rr:
                acc = [x + i * y for x, y in zip(acc, smul(ci, E[rr - i], T))]
        E.append([x / rr for x in acc])
    A = [smul(hd, E[rr], T) for rr in range(NMAX + 1)]
    F3p = [[Fr(1)] + [Fr(0)] * T]
    for m in range(1, 2 * NMAX + 1):
        F3p.append(smul(F3p[-1], a, T))
    c = []
    for n in range(NMAX + 1):
        tot = Fr(0)
        for rr in range(n + 1):
            for m in range(0, 2 * (n - rr) + 1):
                j = n - rr + m
                if 2 * j > T:
                    continue
                prod = smul(F3p[m], A[rr], 2 * j)
                df = 1
                for q in range(1, 2 * j, 2):
                    df *= q
                tot += df * C ** j * prod[2 * j] / factorial(m)
        c.append(tot)
    return c

def hpf_polys(N, free):               # POLY s.1 recursion (same code as g5_degree_bound.py)
    ph = [[Fr(1)]]
    def theta(p): return [k * c for k, c in enumerate(p)]
    def padd(a, b, c=1):
        n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + c * (b[i] if i < len(b) else 0) for i in range(n)]
    def pmul_1mu(p): return padd(p + [0], [0] + p, -1)
    def phis(phl, nmax):
        Phi = [[list(phl[m]) if m < len(phl) else [] for m in range(nmax + 1)]]
        for j in range(1, 6):
            prev = Phi[-1]; cur = []
            for m in range(nmax + 1):
                t = list(prev[m])
                if m >= 1:
                    q = padd([Fr(4 * (m - 1) + j, 5) * c for c in prev[m - 1]], theta(prev[m - 1]), -1)
                    t = padd(t, pmul_1mu(q))
                cur.append(t)
            Phi.append(cur)
        return Phi
    coef = [Fr(3), Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625)]
    for n in range(1, N + 1):
        Phi = phis(ph + [[]], n + 1); m = n + 1
        X = padd(Phi[5][m], Phi[0][m], -1); br = []
        for i, cc in enumerate(coef):
            mm = m - (i + 1)
            if mm >= 0:
                br = padd(br, [Fr(0)] * i + [cc * x for x in Phi[4 - i][mm]])
        X = padd(X, pmul_1mu(br), -1)
        r = []; rem = list(X)
        for k in range(len(rem) - 1):
            r.append(rem[k]); rem[k + 1] += rem[k]
        assert not rem or rem[-1] == 0
        new = []
        for k, cc in enumerate(r):
            if 5 * k == 4 * n:
                assert cc == 0
                new.append(free.get(n, Fr(0)))
            else:
                new.append(cc / (5 * k - 4 * n))
        if 4 * n % 5 == 0 and 4 * n // 5 >= len(new):
            new += [Fr(0)] * (4 * n // 5 - len(new)) + [free.get(n, Fr(0))]
        ph.append(new)
    return ph

def ev(p, u): return sum(c * u ** k for k, c in enumerate(p))

def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    pts = [Fr(2), Fr(3), Fr(-2), Fr(1, 2), Fr(5, 3)]
    print('g5_flat_hatted: NMAX = %d, points L0 =' % NMAX, [str(p) for p in pts])
    ph0 = hpf_polys(min(NMAX, 9), {})
    cs = {L0: flat_c(L0, NMAX) for L0 in pts}
    for L0 in pts:
        ok = all(cs[L0][n] == L0 ** (4 * n) * ev(ph0[n], L0 ** -5) for n in range(min(NMAX, 4) + 1))
        print('L0 = %s: c_n(L0) = L0^(4n) ph_n(L0^-5) for n <= 4:' % L0, ok)
    if NMAX >= 5:
        f1s = [cs[L0][5] - L0 ** 20 * ev(ph0[5], L0 ** -5) for L0 in pts]
        print('c_5(L0) - L0^20 ph^(0)_5(L0^-5) at the points:', [str(x) for x in f1s], '; equals 4537/3937500:', all(x == Fr(4537, 3937500) for x in f1s))
        ph = hpf_polys(min(NMAX, 9), {5: f1s[0]})
        for L0 in pts:
            ok = all(cs[L0][n] == L0 ** (4 * n) * ev(ph[n], L0 ** -5) for n in range(min(NMAX, 9) + 1))
            print('L0 = %s: c_n(L0) = L0^(4n) ph^SP_n(L0^-5) for n <= %d (free coefficient f_1 from c_5):' % (L0, min(NMAX, 9)), ok)
    L0 = pts[0]
    for var in ('X1', 'X2'):
        cv = flat_c(L0, 2, var)
        print('control %s at L0 = 2: c_1 matches L0^4 ph_1(L0^-5):' % var, cv[1] == L0 ** 4 * ev(ph0[1], L0 ** -5), '(must be False)')
    print('control X3 at L0 = 2: c_1 = L0^3 ph_1(L0^-5):', cs[L0][1] == L0 ** 3 * ev(ph0[1], L0 ** -5), '(must be False)')
    print('ph_1 =', [str(x) for x in ph0[1]], ' c_1(2) =', cs[L0][1])

if __name__ == '__main__':
    main()
