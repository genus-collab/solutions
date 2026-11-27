"""g5_valid.py -- validation of the eps-engine (g5_core) by an independent fixed-point stationary-phase code.

FP(sigma_c): normalised formal SP at a fixed rational critical point sigma_c (nothing singular there), computed
directly in t = sigma - sigma_c by Lagrange-Buermann, from the t-series of S_0'', A_0' and A_n (no eps, no scaling,
no groups).  Checks, exact over Q:
 V1 quintic colour 0: the engine with eps-truncation K = 5N is a polynomial in eps of degree <= 5n at z^n
    (c_n in Q[L], deg <= 4n: PROVED in Stage 4); the polynomial evaluated at L = 2, -2/3, 3/2 equals FP(L); and
    c_5(2) = 11495493197/31500 (the value printed in STAGE4_SYNTHESIS s.3, E3).
 V2 O(5)-twisted P^4, weights W1, colour 0 (c_n(1 + eps) rational, not polynomial): the truncated eps-series of the
    engine at eps0 = 1/1000 agrees with FP(1 + eps0) to better than 10^-20 relative (truncation error only).
Usage: py -I -B g5_valid.py N K2      Log: logs/g5_valid_<N>_<K2>.log
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
import g5_core as C, g5_data as D, g5_theory as T

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
K2 = int(sys.argv[2]) if len(sys.argv) > 2 else 24
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/g5_valid_%d_%d.log' % (N, K2)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []


def out(s):
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok:
        FAILS.append(name)


def ser_div(num, den, J):
    q = [Fr(0)] * (J + 1)
    for m in range(J + 1):
        s = (num[m] if m < len(num) else 0) - sum(den[i] * q[m - i] for i in range(1, min(m, len(den) - 1) + 1))
        q[m] = s / den[0]
    return q


def recip_pow(w, c, m, J):
    """t-series of (w + c t)^-m"""
    return [D.gb(-m, k) * Fr(c) ** k * Fr(w) ** (-m - k) for k in range(J + 1)]


def fp_sp(S2, A0p, An, N):
    """normalised SP at a fixed point: S2 = t-series of S_0'', A0p of A_0', An[n] of A_n (n = 1..N)"""
    J = 2 * N
    Pc = {k: S2[k - 2] / (k * (k - 1)) for k in range(2, J + 3)}
    a = -2 * Pc[2]
    Th = [Fr(1)] + [-2 * Pc[k + 2] / a for k in range(1, J + 1)]
    H0 = [Fr(0)] + [A0p[m] / (m + 1) for m in range(J)]
    g = [None] * (N + 1)                     # g[n] = t-series of [z^n] exp(H0 + Sum z^n A_n)
    g[0] = C.sexp(H0, J)
    for n in range(1, N + 1):
        acc = [Fr(0)] * (J + 1)
        for i in range(1, n + 1):
            pr = C.smul(An[i], g[n - i], J)
            acc = [x + i * y for x, y in zip(acc, pr)]
        g[n] = [x / n for x in acc]
    res = [Fr(0)] * (N + 1); dfm = 1
    for m in range(N + 1):
        if m:
            dfm *= 2 * m - 1
        Ps = C.spow_T([[x] for x in Th], Fr(-(2 * m + 1), 2), 2 * m, 0)
        Ps = [p[0] for p in Ps]
        for p in range(N - m + 1):
            res[m + p] += dfm * sum(g[p][j] * Ps[2 * m - j] for j in range(2 * m + 1)) / a ** m
    return res


def quintic_fp(L, N):
    J = 2 * N + 2
    F = ser_div([Fr(5) * D.gb(4, i) * Fr(L) ** (4 - i) for i in range(5)],
                [Fr(L) ** 5 - 1] + [D.gb(5, i) * Fr(L) ** (5 - i) for i in range(1, 6)], J + N + 2)  # Sum_j 1/(s-zeta^j)
    def pw(m):                               # t-series of Sum_j (sigma - zeta^j)^-m
        fac = Fr((-1) ** (m - 1))
        for i in range(1, m):
            fac /= i
        return [fac * F[k + m - 1] * Fr(1).__class__(1) * _ff(k + m - 1, m - 1) for k in range(J + 1)]
    S2 = [5 * x - y for x, y in zip(recip_pow(L, 1, 1, J), F)]
    A0p = [(x - y) / 2 for x, y in zip(recip_pow(L, 1, 1, J), F)]
    An = [None]
    for n in range(1, N + 1):
        if n % 2 == 0:
            An.append([Fr(0)] * (J + 1)); continue
        bk = D.beta(0, n); tw = recip_pow(5 * Fr(L), 5, n, J); bs = pw(n)
        An.append([bk * (x - y) for x, y in zip(tw, bs)])
    return fp_sp(S2, A0p, An, N)


def _ff(a, b):
    r = 1
    for i in range(b):
        r *= (a - i)
    return r


def rational_fp(fac, sc, N):
    J = 2 * N + 2
    S2 = [Fr(0)] * (J + 1); A0p = [Fr(0)] * (J + 1); An = [None] + [[Fr(0)] * (J + 1) for _ in range(N)]
    for f in fac:
        w = f['c'] * sc + f['d']
        S2 = [x + f['c'] ** 2 * y for x, y in zip(S2, recip_pow(w, f['c'], 1, J))]
        A0p = [x + (f['e'] - Fr(1, 2)) * f['c'] * y for x, y in zip(A0p, recip_pow(w, f['c'], 1, J))]
        for n in range(1, N + 1):
            An[n] = [x + D.beta(f['e'], n) * y for x, y in zip(An[n], recip_pow(w, f['c'], n, J))]
    return fp_sp(S2, A0p, An, N)


def main():
    out('g5_valid: engine vs independent fixed-point SP, z^0..z^%d, exact over Q' % N)
    K = 5 * N
    gr = D.quintic_groups(N + K + 1); sg = D.sing_base(N)
    ct = C.engine(D.build_loc(gr, sg, N, K), N, K)
    check('V1 quintic: engine series in eps vanishes beyond eps^(5n) at z^n (polynomial c_n), n <= %d' % N,
          all(ct[n][r] == 0 for n in range(N + 1) for r in range(5 * n + 1, K + 1)))
    for L in (Fr(2), Fr(-2, 3), Fr(3, 2)):
        eng = [sum(ct[n][r] * (L - 1) ** (r - n) for r in range(n, K + 1)) for n in range(N + 1)]
        fp = quintic_fp(L, N)
        check('V1 quintic: engine polynomial at L = %s equals the fixed-point SP, n <= %d' % (L, N), eng == fp)
        if L == 2 and N >= 5:
            check('V1 c_5(2) = 11495493197/31500 (STAGE4_SYNTHESIS s.3, E3)', fp[5] == Fr(11495493197, 31500), fp[5])
    lam = [Fr(1), Fr(3), Fr(-2), Fr(5, 2), Fr(-7, 3)]; fac = T.base(lam) + [T.convex(5)]
    gr2, _ = T.groups_for(fac, 0, N + K2 + 1)
    ct2 = C.engine(D.build_loc(gr2, D.sing_base(N), N, K2), N, K2)
    e0 = Fr(1, 1000); fp2 = rational_fp(fac, 1 + e0, N)
    worst = 0.0
    for n in range(N + 1):
        eng = sum(ct2[n][r] * e0 ** (r - n) for r in range(n, K2 + 1))
        rel = abs(float((eng - fp2[n]) / fp2[n])) if fp2[n] else abs(float(eng))
        worst = max(worst, rel)
        out('   V2 W1 colour 0, n = %d: |engine(trunc. eps^%d) - FP(1 + 1/1000)| / |FP| = %.3e' % (n, K2, rel))
    check('V2 agreement to better than 10^-20 relative (sanity check; truncation error only), all n <= %d' % N,
          worst < 1e-20, 'worst %.3e' % worst)
    out('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS)


if __name__ == '__main__':
    T0 = time.time(); main(); out('total %.1f s; log %s' % (time.time() - T0, LOGREL))
