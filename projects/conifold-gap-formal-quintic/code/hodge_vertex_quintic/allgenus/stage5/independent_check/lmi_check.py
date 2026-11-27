# lmi_check.py N -- independent check (main session, 2026-11-14) of the two-variable chart proof of (Gamma)
# (stage5/chart/referee_adversary/REFEREE.md s.3, Theorem LMI with s* = 1, alpha0 = -1, M0 = 0).
# Own code (no project, worker, referee or critic module). Exact over Q: every Taylor coefficient is a symmetric function of the
# zeta^j, j = 1..4, computed from Lam1(x) = 1 + x + .. + x^4, so no number field is needed.
#  (1) c(z; L): normalised formal SP in sigma of the full-Stirling summand (all six factors), 1 variable, at rational L != 0, 1.
#  (2) the 2-variable chart: phase g(sigma) + (1 - sigma) mu - Q e^mu, amplitude F(sigma, z) (five regular Stirling factors),
#      at its critical point (L, mu_c), Q e^(mu_c) = 1 - L. For L != 1 it must equal (1) (Fubini + Stirling for Gamma((1-sigma)/z)).
#  (3) the chart at L = 1 (Q = 0) must equal h_(w_0)(z) (the delta lemma, C_ss = 0).
# Controls (must fail): amplitude without the singular factor's Stirling series in (1); twist weight +5 in h; chart with
# Q e^mu replaced by Q e^(2 mu) (wrong Euler integrand).
import sys, time
from fractions import Fraction as Fr
from math import comb
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
T0 = time.time()
B = [Fr(1), Fr(-1, 2)]                                             # Bernoulli numbers B_0, B_1, ..
for m in range(2, 2 * N + 4):
    B.append(-sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1))
bk = {k: B[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, N + 2)}
K = 2 * N + 4                                                      # s-degree needed


def pser_div(num, den, n):                                         # power series num/den, n terms
    r = []
    for i in range(n):
        s = (num[i] if i < len(num) else 0) - sum(den[j] * r[i - j] for j in range(1, min(i, len(den) - 1) + 1))
        r.append(Fr(s) / den[0])
    return r


def shift_poly(c, L):                                              # coefficients of p(L + t) from those of p(x)
    out = [Fr(0)] * len(c)
    for i, a in enumerate(c):
        for j in range(i + 1):
            out[j] += a * comb(i, j) * Fr(L) ** (i - j)
    return out


def Pp(L, n):                                                      # P'_r(L) = Sum_(j=1..4) (zeta^j - L)^-r, r = 1..n
    lam = shift_poly([1, 1, 1, 1, 1], L)                           # Lam1(L + t)
    dlam = [(i + 1) * lam[i + 1] for i in range(4)] + [0]
    g = pser_div([-x for x in dlam], lam, n)                       # Sum_r P'_r t^(r-1) = -Lam1'(L+t)/Lam1(L+t)
    return {r: g[r - 1] for r in range(1, n + 1)}


def data(L, j0, twist=-5):
    """Phase Taylor coefficients S_k (k >= 2) and log-amplitude as {(i, n): coeff of s^i z^n}; j0 includes Gamma((1-sigma)/z)."""
    L = Fr(L); P = Pp(L, K + 2 * N + 4)
    w5 = -twist                                                    # the factor Gamma(1 + w5 sigma/z): phase a log a - a, a = w5 sigma
    S = {}
    for k in range(2, K + 1):
        v = w5 * Fr((-1) ** k) / L ** (k - 1) + P[k - 1]
        if j0: v += 1 / (1 - L) ** (k - 1)
        S[k] = v / (k * (k - 1))
    la = {}
    def addc(i, n, c):
        if c: la[(i, n)] = la.get((i, n), 0) + c
    for i in range(1, K + 1):                                       # (1/2) log(5 sigma) - (1/2) Sum_(j>=1) log(zeta^j - sigma)
        addc(i, 0, Fr((-1) ** (i + 1), 2 * i) / L ** i + P[i] / (2 * i))
        if j0: addc(i, 0, Fr(1, 2 * i) / (1 - L) ** i)
    for k in range(1, N + 1):                                       # b_k z^(2k-1) [(w5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]
        r = 2 * k - 1
        for i in range(0, K + 1):
            c = Fr(w5) ** (1 - 2 * k) * Fr((-1) ** i * comb(r + i - 1, i)) / L ** (r + i) + comb(r + i - 1, i) * P[r + i]
            if j0: c += comb(r + i - 1, i) / (1 - L) ** (r + i)
            addc(i, r, bk[k] * c)
    return S, la


def bexp(la, maxdeg):                                              # exp of a bivariate series (no constant term), weight i + 2n <= maxdeg
    res = {(0, 0): Fr(1)}; term = {(0, 0): Fr(1)}
    for v in range(1, maxdeg + 1):
        new = {}
        for (i, n), c in term.items():
            for (i2, n2), c2 in la.items():
                if i + i2 + 2 * (n + n2) <= maxdeg:
                    new[(i + i2, n + n2)] = new.get((i + i2, n + n2), 0) + c * c2 / v
        term = new
        for kk, c in term.items(): res[kk] = res.get(kk, 0) + c
        if not term: break
    return res


def dfact(n): return 1 if n <= 0 else n * dfact(n - 2)


def moment(a, b, C):                                               # <eta1^a eta2^b> for covariance C = (c11, c12, c22)
    c11, c12, c22 = C; tot = Fr(0)
    for k in range(min(a, b) + 1):
        if (a - k) % 2 or (b - k) % 2: continue
        f = comb(a, k) * comb(b, k)
        for t in range(1, k + 1): f *= t
        tot += f * dfact(a - k - 1) * dfact(b - k - 1) * c12 ** k * c11 ** ((a - k) // 2) * c22 ** ((b - k) // 2)
    return tot


def sp(vert, amp, C, n):
    """Normalised formal SP through z^n. vert: list of (k, i, j, coeff) for a vertex c eta_s^i eta_m^j of degree k = i + j >= 3
    (weight eps^(k-2), eps = z^(1/2)); amp: {(i, nz): coeff} meaning coeff s^i z^nz (weight eps^(i + 2 nz))."""
    E = 2 * n
    expo = {(0, 0, 0): Fr(1)}; term = {(0, 0, 0): Fr(1)}
    for V in range(1, E + 1):
        new = {}
        for (e, a, b), c in term.items():
            for (k, i, j, cv) in vert:
                e2 = e + k - 2
                if e2 <= E:
                    key = (e2, a + i, b + j); new[key] = new.get(key, 0) + c * cv / V
        term = new
        for kk, c in term.items(): expo[kk] = expo.get(kk, 0) + c
        if not term: break
    out = [Fr(0)] * (n + 1)
    for (e, a, b), c in expo.items():
        for (i, nz), ca in amp.items():
            e2 = e + i + 2 * nz
            if e2 <= E and e2 % 2 == 0 and (a + i + b) % 2 == 0:
                out[e2 // 2] += c * ca * moment(a + i, b, C)
    return [x / out[0] for x in out]


def c_one(L, n, ctrl=None):                                       # (1): one variable, all six factors
    S, la = data(L, True)
    if ctrl == 'nostirl0':                                         # keep the z^0 part; Stirling series of the five regular factors only
        lax = data(L, False)[1]
        la = {kk: v for kk, v in la.items() if kk[1] == 0}
        la.update({kk: v for kk, v in lax.items() if kk[1] > 0})
    C = (-1 / (2 * S[2]), Fr(0), Fr(0))
    vert = [(k, k, 0, S[k]) for k in range(3, 2 * n + 3)]
    return sp(vert, bexp(la, 2 * n), C, n)


def c_chart(L, n, ctrl=None):                                     # (2)/(3): the two-variable chart, variables (s, m)
    L = Fr(L); S, la = data(L, False)
    q = 2 if ctrl == 'e2mu' else 1                                # -Q e^(q mu); critical: (1 - sigma) = q Q e^(q mu)
    cm = -(1 - L) / q                                             # coefficient: -Q e^(q mu_c) = -(1 - L)/q
    h_ss, h_sm, h_mm = 2 * S[2], Fr(-1), cm * q * q
    det = h_ss * h_mm - h_sm ** 2
    C = (-h_mm / det, h_sm / det, -h_ss / det)                     # C = -H^-1
    vert = [(k, k, 0, S[k]) for k in range(3, 2 * n + 3)] + [(k, 0, k, cm * Fr(q) ** k / prod(k)) for k in range(3, 2 * n + 3)]
    return sp(vert, bexp(la, 2 * n), C, n)


def prod(k):
    r = 1
    for t in range(2, k + 1): r *= t
    return r


def hser(n, twist=-5):
    # h_w(z) = exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k)), w = {1 - zeta^j : j = 1..4} U {twist}, from the DEFINITION:
    # the 1/(1 - zeta^j) are the roots of w^4 - 2 w^3 + 2 w^2 - w + 1/5 (Newton sums), independent of the P' used above
    e = [Fr(1), Fr(2), Fr(2), Fr(1), Fr(1, 5)]                     # elementary symmetric functions e_0..e_4 of the roots
    p = {}
    for r in range(1, 2 * n + 2):
        v = sum((-1) ** (i - 1) * e[i] * p[r - i] for i in range(1, min(r - 1, 4) + 1))
        if r <= 4: v += (-1) ** (r - 1) * r * e[r]
        p[r] = v
    lg = {}
    for k in range(1, n + 1):
        if 2 * k - 1 <= n:
            lg[(0, 2 * k - 1)] = -bk[k] * (Fr(twist) ** (1 - 2 * k) + p[2 * k - 1])
    e = bexp(lg, 2 * n)
    return [e.get((0, m), Fr(0)) for m in range(n + 1)]


fails = []
def report(ok, msg, control=False):
    good = ok != control
    if not good: fails.append(msg)
    print('%s %s%s   (%.1f s)' % ('PASS' if good else 'FAIL', msg, '  [control: must fail]' if control else '', time.time() - T0))


h = hser(N)
print('h_(w_0): h_1..h_5 =', [str(x) for x in h[1:6]])
for L in (Fr(2), Fr(-2, 3), Fr(1, 2), Fr(3)):
    a = c_one(L, N); b = c_chart(L, N)
    report(a == b, 'L = %s: two-variable chart SP = c(z; L) (one variable, full Stirling) through z^%d' % (L, N))
    if L == 2: print('     c_1..c_3(2) =', [str(x) for x in a[1:4]])
b1 = c_chart(1, N)
report(b1 == h, 'L = 1 (Q = 0): two-variable chart SP = h_(w_0)(z) through z^%d' % N)
report(c_one(Fr(2), N, 'nostirl0') == c_chart(Fr(2), N), 'L = 2: c(z; L) with the singular factor\'s Stirling series dropped = chart', control=True)
report(hser(N, twist=5) == b1, 'L = 1: chart = h with twist weight +5', control=True)
report(c_chart(Fr(2), N, 'e2mu') == c_one(Fr(2), N), 'L = 2: chart with Q e^(2 mu) = c(z; L)', control=True)
print('LMI_CHECK %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not fails else 'FAILED: %s' % fails, time.time() - T0))
