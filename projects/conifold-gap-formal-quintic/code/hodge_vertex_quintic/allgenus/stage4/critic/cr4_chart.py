# -*- coding: utf-8 -*-
"""cr4_chart.py -- Stage-4 critic: falsifier for the chart step of the formal route to (Gamma) (STAGE4_SYNTHESIS s.1.4).

New code (stdlib + cr4_core).  The 7-variable formal integral exp(Phi/z) z^-1 dsigma dy prod_(j=0..4) dx_j/x_j,
Phi = sigma log(-q) - y + 5 sigma log y + Sum_j (-x_j + (zeta^j - sigma) log x_j), integrated first over (sigma, x_0) (phase linear in
sigma, amplitude free of sigma: restriction to x_0 = X_0 = -q y^5/prod_(j>=1) x_j, Jacobian |x_0|) leaves the CHART integral
    Omega_ch = exp( [5 log y - y + Sum_(j=1..4) (-x_j + (zeta^j - 1) log x_j) + q y^5/prod_j x_j] / z ) dy prod_(j=1..4) dx_j/x_j
(up to (-q)^(1/z) and constants).  Colour-0 critical point x_j = zeta^j - L, y = 5L (u L^5 = 1); there q y^5/prod x_j = L - 1.
The outline predicts: the normalised formal SP of Omega_ch equals R_n(L) = L^(4n) ph_n(L^-5) (= c_n(L), the full-Stirling series), and
at L = 1 (q = 0) it factorises into Stirling series: h_(w_0)(z).  Computed here by Wick's theorem in 5 variables, mod p (zeta in F_p).
Controls: drop the Jacobian |x_0| (amplitude / X_0): fails at z^1; dy/y instead of dy: fails first at z^2 (needs N >= 2).
Usage: py -I -B cr4_chart.py [N]   (default 3; log logs/cr4_chart_<N>.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from cr4_core import Field, hcoef, Wkb, peval

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
LOGREL = 'logs/cr4_chart_%d.log' % N
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()

def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()

def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)

def pmul(A, B, D, p):
    """A, B: {deg: {mono: c}} in 5 variables, product truncated at total degree D"""
    res = {}
    for da, Ad in A.items():
        for db, Bd in B.items():
            if da + db > D: continue
            od = res.setdefault(da + db, {})
            for e1, c1 in Ad.items():
                for e2, c2 in Bd.items():
                    e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2], e1[3] + e2[3], e1[4] + e2[4])
                    od[e] = (od.get(e, 0) + c1 * c2) % p
    return res

def uni(i, s):
    """univariate series s (list) in variable i as a 5-variable polynomial"""
    return {k: {tuple(k if a == i else 0 for a in range(5)): c} for k, c in enumerate(s) if c}

def padd(A, B, p, cb=1):
    res = {d: dict(m) for d, m in A.items()}
    for d, m in B.items():
        od = res.setdefault(d, {})
        for e, c in m.items(): od[e] = (od.get(e, 0) + cb * c) % p
    return res

def inv_mod_matrix(M, p):
    n = len(M); A = [list(r) + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        r = next(r for r in range(c, n) if A[r][c] % p); A[c], A[r] = A[r], A[c]
        iv = pow(A[c][c], -1, p); A[c] = [x * iv % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]

def chart_sp(p, zeta, L, N, amp='dy'):
    """normalised formal SP of Omega_ch at the colour-0 critical point, z^0..z^N, mod p"""
    D = 2 * N + 2; L = Fr(L) if not isinstance(L, int) else Fr(L)
    Lp = L.numerator * pow(L.denominator, -1, p) % p
    xc = [(pow(zeta, j, p) - Lp) % p for j in range(1, 5)]; yc = 5 * Lp % p
    aj = [(pow(zeta, j, p) - 1) % p for j in range(1, 5)]
    K0 = (Lp - 1) % p                                            # q y_c^5 / prod x_c at the critical point
    inv = lambda x: pow(x % p, -1, p)
    def logser(c, s):                                            # c log(1 + t/s), t^1..t^D
        return [0] + [c * (-1) ** (k + 1) * inv(k) * pow(inv(s), k, p) % p for k in range(1, D + 1)]
    Phi = {}
    for i in range(4):
        g = logser(aj[i], xc[i]); g[1] = (g[1] - 1) % p; Phi = padd(Phi, uni(i, g), p)
    g = logser(5, yc); g[1] = (g[1] - 1) % p; Phi = padd(Phi, uni(4, g), p)
    cpl = {0: {(0,) * 5: K0}}
    fy = [1, 5, 10, 10, 5, 1] + [0] * D
    cpl = pmul(cpl, uni(4, [fy[k] * pow(inv(yc), k, p) % p for k in range(D + 1)]), D, p)
    geo = lambda s: [(-1) ** k * pow(inv(s), k, p) % p for k in range(D + 1)]      # (1 + t/s)^-1
    for i in range(4): cpl = pmul(cpl, uni(i, geo(xc[i])), D, p)
    cpl[0] = {}; Phi = padd(Phi, cpl, p)
    assert all(c == 0 for c in Phi.get(1, {}).values()), "not a critical point"
    H = [[0] * 5 for _ in range(5)]
    for e, c in Phi.get(2, {}).items():
        idx = [a for a in range(5) for _ in range(e[a])]
        if idx[0] == idx[1]: H[idx[0]][idx[0]] = 2 * c % p
        else: H[idx[0]][idx[1]] = H[idx[1]][idx[0]] = c
    C = [[(-x) % p for x in r] for r in inv_mod_matrix(H, p)]  # covariance (z-free): -H^-1
    P = {d: m for d, m in Phi.items() if d >= 3}
    a = {0: {(0,) * 5: 1}}
    if amp == 'nojac':                                           # A / X_0, normalised: (1 + eta/y_c)^-5
        for _ in range(5): a = pmul(a, uni(4, geo(yc)), 2 * N, p)
    else:
        for i in range(4): a = pmul(a, uni(i, geo(xc[i])), 2 * N, p)   # prod (1 + xi_j/x_j)^-1: amplitude prod 1/x_j (dy)
        if amp == 'dyy': a = pmul(a, uni(4, geo(yc)), 2 * N, p)       # control: dy/y
    memo = {(0,) * 5: 1}
    def M(g):
        if g in memo: return memo[g]
        i = next(k for k in range(5) if g[k]); b = list(g); b[i] -= 1; s = 0
        for j in range(5):
            if b[j] and C[i][j]:
                bb = list(b); bb[j] -= 1; s += C[i][j] * b[j] * M(tuple(bb))
        memo[g] = s % p; return memo[g]
    sp = [0] * (N + 1); E = a
    for m in range(0, 2 * N + 1):
        for d, md in E.items():
            if d % 2 == 0 and d // 2 - m <= N and d // 2 - m >= 0:
                for e, c in md.items(): sp[d // 2 - m] = (sp[d // 2 - m] + c * M(e)) % p
        E = pmul(E, P, 2 * (N + m + 1), p); im = inv(m + 1)
        E = {d: {e: c * im % p for e, c in md.items()} for d, md in E.items()}
    return sp

out('cr4_chart: N = %d (z^0..z^%d); formal SP of the 5-variable chart integral mod p; log %s' % (N, N, LOGREL))
for p in (2 ** 61 - 1, 2 ** 89 - 1):
    F = Field(p); c = 2
    while pow(c, (p - 1) // 5, p) == 1: c += 1
    zeta = pow(c, (p - 1) // 5, p)
    h = hcoef(F, 5 * N + 5); W = Wkb(F, max(N, 5), h, 'qrr')
    t = time.time(); s1 = chart_sp(p, zeta, 1, N)
    check('[%s] L = 1 (q = 0): chart SP factorises into Stirling series = h_(w_0)(z) through z^%d' % (F.name(), N), s1 == h[:N + 1], s1[1] == F.c(Fr(-3, 20)))
    for L in (Fr(2), Fr(3, 2), Fr(-2, 3)):
        Lp = F.c(L); iu = F.c(1 / L ** 5)
        R = [F.r(pow(Lp, 4 * n, p) * peval(F, W.ph[n], iu)) for n in range(N + 1)]
        s = chart_sp(p, zeta, L, N)
        check('[%s] L = %s: chart SP = R_n(L) = L^4n ph_n(L^-5) (= c_n(L)) for n <= %d' % (F.name(), L, N), s == R)
        if N >= 5:                                               # z^5 is the first order that sees (N): stored normalisation
            Ws = Wkb(F, 5, h, 'list', fv={1: F.r(W.f[1] + F.c(Fr(-43397, 807975000)))})
            R5s = F.r(pow(Lp, 20, p) * peval(F, Ws.ph[5], iu))
            check('[%s] L = %s: z^5 sees (N): the stored normalisation (f_1 + m5) gives R_5 - chart_5 = m5 exactly' % (F.name(), L),
                  F.r(R5s - s[5]) == F.c(Fr(-43397, 807975000)))
    sj = chart_sp(p, zeta, Fr(2), N, amp='nojac'); sd = chart_sp(p, zeta, Fr(2), N, amp='dyy')
    R2 = [F.r(pow(F.c(2), 4 * n, p) * peval(F, W.ph[n], F.c(Fr(1, 32)))) for n in range(N + 1)]
    check('[%s] control: no Jacobian |x_0| (amplitude / X_0) fails at z^1' % F.name(), sj[0] == 1 and sj[1] != R2[1])
    fd = [n for n in range(N + 1) if sd[n] != R2[n]]
    check('[%s] control: measure dy/y instead of dy fails (first at z^2; z^1 is blind to it)' % F.name(), fd[:1] == [2], fd)
    out('     [%s] done in %.1f s' % (F.name(), time.time() - t))
out('CR4_CHART %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
