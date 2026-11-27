"""gamma_interp.py -- independent check of identity (Gamma) of Stage 4 (formal quintic, colour 0), exact over Q.

c(z; L) := normalised formal stationary-phase series, at the critical point sigma = L (u L^5 = 1), of the Mellin-Barnes summand
    q^s Gamma(1+5s) / prod_(j=0..4) Gamma(1 + s - zeta^j/z),   sigma = z s,
with every Gamma(1+y) replaced by its Stirling form sqrt(2 pi) y^(y+1/2) e^-y exp(Sum_k b_k y^(1-2k)), b_k = B_2k/(2k(2k-1)):
    phase  S(sigma) = sigma log q + 5 sigma log(5 sigma) - Sum_j (sigma - zeta^j) log(sigma - zeta^j)   (log z terms cancel),
    amplitude sqrt(5 sigma/(sigma^5 - 1)) exp(Sum_k b_k z^(2k-1) [(5 sigma)^(1-2k) - Sum_j (sigma - zeta^j)^(1-2k)]).
Claims tested, with no WKB recursion, no gauge and no project or agent code:
  (P) c_n(L) is a polynomial in L with exponents 4n - 5k, 0 <= k <= floor(4n/5) (fitted on some rational L != 0, 1, confirmed
      on EXTRA further values);
  (G) its value at L = 1 (q = 0) equals [z^n] h_(w_0)(z), h_w = exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k)),
      w_0 = {1 - zeta^j (j = 1..4)} U {-5}.
Controls (must FAIL): the Stirling correction of Gamma(1+5s) dropped; b_1 -> b_1 + 1/7; b_3 -> b_3 + 1/7 (z^5 only).
Run: py -I -B gamma_interp.py NMAX [EXTRA] [noctrl]   (stdout; working directory outside the project)
"""
import sys, time
from fractions import Fraction as Fr
from math import comb

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
EXTRA = int(sys.argv[2]) if len(sys.argv) > 2 else 3
LS = [Fr(2), Fr(3), Fr(-2), Fr(1, 2), Fr(3, 2), Fr(-1, 2), Fr(5, 2), Fr(4, 3), Fr(-3, 2), Fr(5, 3), Fr(7, 3), Fr(3, 4),
      Fr(-4, 3), Fr(5, 4), Fr(7, 4), Fr(-5, 2), Fr(2, 3), Fr(6, 5), Fr(-2, 3), Fr(9, 4)]
LS += [x for x in sorted({Fr(a, b) for b in range(1, 8) for a in range(-12, 13)}, key=lambda x: (abs(x.numerator) + x.denominator, x))
       if x not in (0, 1, -1) and x not in LS]                                  # enough distinct sample points for any NMAX used


def bernoulli(N):
    B = [Fr(0)] * (N + 1); B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BER = bernoulli(2 * NMAX + 2)
BK = {k: BER[2 * k] / (2 * k * (2 * k - 1)) for k in range(1, NMAX + 2)}


def smul(a, b, K):
    c = [Fr(0)] * (K + 1)
    for i, x in enumerate(a[:K + 1]):
        if x:
            for j, y in enumerate(b[:K + 1 - i]):
                if y:
                    c[i + j] += x * y
    return c


def sinv(a, K):
    c = [Fr(0)] * (K + 1); c[0] = 1 / a[0]
    for m in range(1, K + 1):
        c[m] = -sum(a[i] * c[m - i] for i in range(1, min(m, len(a) - 1) + 1)) / a[0]
    return c


def spow(a, alpha, K):                     # a[0] = 1; J.C.P. Miller recurrence for a^alpha
    c = [Fr(0)] * (K + 1); c[0] = Fr(1)
    for m in range(1, K + 1):
        c[m] = sum(((alpha + 1) * i - m) * a[i] * c[m - i] for i in range(1, min(m, len(a) - 1) + 1)) / m
    return c


def pmul(p, q):                            # polynomials in tau (lists)
    if not p or not q:
        return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                if y:
                    r[i + j] += x * y
    return r


def padd(p, q):
    r = [Fr(0)] * max(len(p), len(q))
    for i, x in enumerate(p): r[i] += x
    for i, x in enumerate(q): r[i] += x
    return r


def bmul(X, Y, A):                         # bivariate: list over eps-power a of tau-polynomials, truncated at eps^A
    Z = [[] for _ in range(A + 1)]
    for i in range(A + 1):
        if X[i]:
            for j in range(A + 1 - i):
                if Y[j]:
                    Z[i + j] = padd(Z[i + j], pmul(X[i], Y[j]))
    return Z


def bexp(X, A):                            # X[0] = 0;  a Y_a = Sum_i i X_i Y_(a-i)
    Y = [[] for _ in range(A + 1)]; Y[0] = [Fr(1)]
    for a in range(1, A + 1):
        acc = []
        for i in range(1, a + 1):
            if X[i] and Y[a - i]:
                acc = padd(acc, [i * c for c in pmul(X[i], Y[a - i])])
        Y[a] = [c / a for c in acc]
    return Y


def cseries(L, ctrl=None):
    """[c_0..c_NMAX] at the rational point L (exact)."""
    A = 2 * NMAX; K = A + 2; KF = K + 2 * NMAX + 2
    D = [L ** 5 - 1] + [comb(5, i) * L ** (5 - i) for i in range(1, 6)]          # (L+t)^5 - 1
    N4 = [comb(4, i) * L ** (4 - i) for i in range(5)]                          # (L+t)^4
    F = [5 * c for c in smul(N4, sinv(D, KF), KF)]                              # Sum_j 1/(L+t-zeta^j)
    invL = [(-1) ** m / L ** (m + 1) for m in range(KF + 1)]                    # 1/(L+t)
    T = [5 * invL[m] - F[m] for m in range(K + 1)]                              # S''(L+t)
    V = -1 / T[0]                                                               # <tau^2>
    P = [[] for _ in range(A + 1)]
    for k in range(3, A + 3):                                                   # S_k eps^(k-2) tau^k
        P[k - 2] = [Fr(0)] * k + [T[k - 2] / ((k - 1) * k)]
    R = smul([Fr(1), 1 / L], [c * (L ** 5 - 1) for c in sinv(D, K)], K)         # (1+t/L)(L^5-1)/((L+t)^5-1)
    am = spow(R, Fr(1, 2), K)
    AM = [[Fr(0)] * a + [am[a]] for a in range(A + 1)]
    E = [[] for _ in range(A + 1)]
    for k in range(1, NMAX + 2):
        r = 2 * k - 1; e0 = 4 * k - 2
        if e0 > A:
            break
        bk = BK[k] + (Fr(1, 7) if ((ctrl == 'b1' and k == 1) or (ctrl == 'b3' and k == 3)) else 0)
        for m in range(A - e0 + 1):
            # [t^m] (5(L+t))^(-r) = (-1)^m C(r+m-1, m) / (5^r L^(r+m))
            first = Fr(0) if ctrl == 'drop5' else Fr((-1) ** m * comb(r + m - 1, m), 5 ** r) / L ** (r + m)
            second = (-1) ** (r - 1) * F[m + r - 1] * comb(m + r - 1, r - 1)    # Sum_j (L+t-zeta^j)^(-r), [t^m]
            val = bk * (first - second)
            if val:
                E[e0 + m] = padd(E[e0 + m], [Fr(0)] * m + [val])
    Y = bmul(bmul(bexp(P, A), AM, A), bexp(E, A), A)
    out = []
    for n in range(NMAX + 1):
        s = Fr(0); dfac = Fr(1); vp = Fr(1)
        for b, c in enumerate(Y[2 * n]):
            if b % 2 == 0:
                if b > 0:
                    dfac *= (b - 1); vp *= V
                if c:
                    s += c * dfac * vp
        out.append(s)
        assert all(x == 0 for x in Y[2 * n + 1][::2]) if 2 * n + 1 <= A and Y[2 * n + 1] else True
    return out


def hseries():
    e = [Fr(2), Fr(2), Fr(1), Fr(1, 5)]                                          # w^4 - 2w^3 + 2w^2 - w + 1/5, w = 1/(1-zeta^j)
    p = [Fr(4)]
    for r in range(1, 2 * NMAX + 2):
        s = sum((-1) ** (i - 1) * e[i - 1] * p[r - i] for i in range(1, min(r - 1, 4) + 1))
        if r <= 4:
            s += (-1) ** (r - 1) * r * e[r - 1]
        p.append(s)
    H = [Fr(0)] * (NMAX + 1)
    for k in range(1, NMAX + 2):
        if 2 * k - 1 <= NMAX:
            H[2 * k - 1] = BK[k] * (Fr(1, 5 ** (2 * k - 1)) - p[2 * k - 1])     # -Sum_(a in w_0) a^(1-2k)
    h = [Fr(1)] + [Fr(0)] * NMAX                                                 # exp(H)
    for m in range(1, NMAX + 1):
        h[m] = sum(i * H[i] * h[m - i] for i in range(1, m + 1)) / m
    return h


def solve(M, y):                           # exact Gaussian elimination, square system
    n = len(y); assert all(len(row) == n for row in M), 'non-square system'
    M = [row[:] + [y[i]] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * z for x, z in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def run(ctrl=None):
    need = (4 * NMAX) // 5 + 1 + EXTRA
    pts = LS[:need]; t0 = time.time()
    assert len(pts) == need, 'not enough sample points: %d < %d' % (len(pts), need)
    vals = {L: cseries(L, ctrl) for L in pts}
    h = hseries(); bad_poly = []; bad_val = []
    for n in range(NMAX + 1):
        ks = list(range((4 * n) // 5 + 1)); m = len(ks)
        coef = solve([[L ** (4 * n - 5 * k) for k in ks] for L in pts[:m]], [vals[L][n] for L in pts[:m]])
        for L in pts[m:]:
            if sum(a * L ** (4 * n - 5 * k) for a, k in zip(coef, ks)) != vals[L][n]:
                bad_poly.append((n, str(L)))
        if sum(coef) != h[n]:
            bad_val.append(n)
        if ctrl is None and n <= 10:
            print('   n = %2d: %d free exponents, fit on %d points, checked on %d more; c_n(1) = %s; h_n = %s'
                  % (n, m, m, len(pts) - m, short(sum(coef)), short(h[n])))
    return bad_poly, bad_val, time.time() - t0, h


def short(x):
    s = str(x)
    return s if len(s) <= 40 else s[:18] + '..(%d chars)' % len(s)


if __name__ == '__main__':
    print('gamma_interp: NMAX = %d (z^0..z^%d), EXTRA = %d check points, exact over Q' % (NMAX, NMAX, EXTRA))
    bp, bv, dt, h = run()
    print('PASS' if not bp else 'FAIL', '(P) c_n(L) polynomial with exponents 4n-5k on all extra points, n <= %d:' % NMAX, bp)
    print('PASS' if not bv else 'FAIL', '(G) c_n(1) = h_n for all n <= %d (incl. (N) at n = 5, 10, ..):' % NMAX, bv, ' (%.1f s)' % dt)
    CTRLS = (('drop5', 'Stirling correction of Gamma(1+5s) dropped'), ('b1', 'b_1 -> b_1 + 1/7'),
             ('b3', 'b_3 -> b_3 + 1/7, first seen at z^5'))
    if len(sys.argv) > 3 and sys.argv[3] == 'noctrl':
        CTRLS = ()
        print('controls skipped (noctrl)')
    for ctrl, name in CTRLS:
        bp2, bv2, dt2, _ = run(ctrl)
        first = min([n for n, _ in bp2] + bv2) if (bp2 or bv2) else None
        print('PASS' if first is not None else 'FAIL', 'control (%s) fails, first at z^%s  (%.1f s)' % (name, first, dt2))
