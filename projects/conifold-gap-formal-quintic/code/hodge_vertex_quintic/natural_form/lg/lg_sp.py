# -*- coding: utf-8 -*-
"""lg_sp.py -- exact formal stationary-phase (Gaussian) expansions over Q (fractions.Fraction), no third-party modules.

ENGINE.  For a phase  Phi(y) = (a/2) y^2 + Sum_{m>=3} phi_m y^m  (a != 0) and an amplitude g(y) = Sum_j g_j y^j the formal integral
      Int exp(-Phi(y)/h) g(y) dy  =  sqrt(2 pi h / a) . Sum_{n>=0} A_n h^n
is computed by y = w t, w^2 = h/a:  exp(-t^2/2) . exp(-Sum_{m>=3} (phi_m/a) w^{m-2} t^m) . g(w t), expanded in w with polynomial
coefficients in t, and the Gaussian moments <t^{2j}> = (2j-1)!!, <t^odd> = 0.  Only even powers of w survive: A_n = [w^{2n}] / a^n.
     gauss_series(a, phi, g, N) -> [A_0, ..., A_N]     (phi = {m: phi_m}, g = list of Taylor coefficients, both exact)
No analytic statement is made: this is the formal (asymptotic) series attached to a non-degenerate critical point.

THE LG MODEL  lambda(x) = x^4 + u/x.  Critical points x_c^5 = u/4, critical value 5 x_c^4, lambda'' = 20 x_c^2.  With x = x_c (1 + y),
eps := z / x_c^4:   lambda/z = (1/eps) phi(y),   phi(y) = (1+y)^4 + 4/(1+y) = 5 + 10 y^2 + 0 y^3 + 5 y^4 - 4 y^5 + 4 y^6 - ...
      Int exp(lambda/z) x^k dx = x_c^{k+1} exp(5/eps) sqrt(2 pi (-eps)/20) . J_k(eps),      J_k(eps) = 1 + Sum_{n>=1} j_{k,n} eps^n,
i.e. h = -eps, Phi = phi - 5, a = 20, g = (1+y)^k.     J_series(k, N) -> [j_{k,0} = 1, ..., j_{k,N}]   (coefficients of eps^n).

VALIDATION (run as a script): Stirling series of Gamma(n+1) and the Airy series, both to order 8, against closed forms / tabulated values.
"""
from fractions import Fraction as Fr


def _pmul(A, B, dmax=None):
    """product of polynomials in t (lists of Fractions)"""
    if not A or not B: return []
    out = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x == 0: continue
        for j, y in enumerate(B):
            if y != 0: out[i + j] += x * y
    return out


def _padd(A, B, c=1):
    n = max(len(A), len(B)); out = [Fr(0)] * n
    for i, x in enumerate(A): out[i] += x
    for i, y in enumerate(B): out[i] += c * y
    return out


def dfact_odd(j):
    """(2j-1)!!"""
    r = 1
    for k in range(1, 2 * j, 2): r *= k
    return r


def gauss_series(a, phi, g, N):
    """[A_0..A_N] with Int exp(-Phi/h) g dy = sqrt(2 pi h/a) Sum A_n h^n; Phi = (a/2) y^2 + Sum phi[m] y^m (m >= 3)"""
    a = Fr(a); W = 2 * N                                   # need w-orders 0..2N
    # S(w, t) = -Sum_{m>=3} (phi_m / a) w^{m-2} t^m : S[j] = polynomial in t at w^j, j = m - 2 >= 1
    S = [[] for _ in range(W + 1)]
    for m, c in phi.items():
        j = m - 2
        if 1 <= j <= W and c != 0:
            S[j] = [Fr(0)] * m + [-Fr(c) / a]
    # E = exp(S):  n E_n = Sum_{j=1}^n j S_j E_{n-j}
    E = [[Fr(1)]] + [[] for _ in range(W)]
    for n in range(1, W + 1):
        acc = []
        for j in range(1, n + 1):
            if S[j] and E[n - j]: acc = _padd(acc, _pmul(S[j], E[n - j]), j)
        E[n] = [x / n for x in acc]
    # amplitude g(w t) = Sum g_j w^j t^j
    out = []
    for n in range(N + 1):
        tot = Fr(0)
        for j in range(min(2 * n, len(g) - 1) + 1):
            if g[j] == 0: continue
            P = E[2 * n - j]
            for d, c in enumerate(P):
                if c != 0 and (d + j) % 2 == 0: tot += g[j] * c * dfact_odd((d + j) // 2)
        out.append(tot / a ** n)
    return out


def binom_series(k, n):
    """Taylor coefficients of (1+y)^k through y^n, k any integer (or Fraction)"""
    out = [Fr(1)]
    for j in range(1, n + 1): out.append(out[-1] * (Fr(k) - (j - 1)) / j)
    return out


def lg_phase(M):
    """coefficients phi_m (m = 3..M) of phi(y) = (1+y)^4 + 4/(1+y)"""
    b4 = binom_series(4, M); bm1 = binom_series(-1, M)
    return {m: b4[m] + 4 * bm1[m] for m in range(3, M + 1)}


def J_series(k, N):
    """J_k(eps) = Sum_n j_n eps^n through eps^N for lambda = x^4 + u/x, amplitude x^k dx (see the module docstring)"""
    M = 2 * N + 2
    phi = lg_phase(M)
    b4 = binom_series(4, 2); bm1 = binom_series(-1, 2)
    assert b4[1] + 4 * bm1[1] == 0 and b4[2] + 4 * bm1[2] == 10
    A = gauss_series(20, phi, binom_series(k, 2 * N), N)      # expansion in h = -eps
    return [A[n] * (-1) ** n for n in range(N + 1)]


def J_general(k, N, amp_extra=None, power=4, pole=1):
    """same for lambda = x^power + u/x^pole (critical point: power x_c^{power+pole} = pole u), amplitude x^k: returns
    (lambda_c/x_c^power, lambda''/x_c^{power-2}, [coefficients in eps = z/x_c^power])"""
    M = 2 * N + 2; r = Fr(power, pole)
    bp = binom_series(power, M); bq = binom_series(-pole, M)
    ph = [bp[m] + r * bq[m] for m in range(M + 1)]
    assert ph[1] == 0
    a = 2 * ph[2]
    A = gauss_series(a, {m: ph[m] for m in range(3, M + 1)}, binom_series(k, 2 * N), N)
    return ph[0], a, [A[n] * (-1) ** n for n in range(N + 1)]


# ------------------------------------------------------------------------------------------------ validation
def stirling(N):
    """Gamma(n+1) = n^{n+1} e^{-n} Int exp(-n (y - log(1+y))) dy = sqrt(2 pi n) (n/e)^n Sum_k A_k n^-k, h = 1/n, Phi = y - log(1+y)"""
    M = 2 * N + 2
    phi = {m: Fr((-1) ** m, m) for m in range(3, M + 1)}     # y - log(1+y) = y^2/2 - y^3/3 + y^4/4 - ...
    return gauss_series(1, phi, [Fr(1)], N)


def airy(N):
    """Int exp(-(y^2 + y^3/3)/h) dy (the Airy saddle, Phi = y^2 + y^3/3, a = 2): A_k = (-1)^k... closed form checked below"""
    return gauss_series(2, {3: Fr(1, 3)}, [Fr(1)], N)


def airy_closed(k):
    """direct evaluation: Int exp(-y^2/h) Sum_j (-y^3/(3h))^j/j! dy: only j = 2k contributes at h^k:
       A_k = (1/(2k)!) (1/9)^k <y^{6k}>/h^{3k} ... = (6k-1)!! / ((2k)! 9^k 2^{3k})"""
    return Fr(dfact_odd(3 * k), _fact(2 * k) * 9 ** k * 2 ** (3 * k))


def _fact(n):
    r = 1
    for k in range(2, n + 1): r *= k
    return r


STIRLING_TABLE = [Fr(1), Fr(1, 12), Fr(1, 288), Fr(-139, 51840), Fr(-571, 2488320), Fr(163879, 209018880),
                  Fr(5246819, 75246796800), Fr(-534703531, 902961561600), Fr(-4483131259, 86684309913600)]

if __name__ == "__main__":
    S = stirling(8)
    print("Stirling series through 1/n^8:", ", ".join(str(x) for x in S))
    assert S == STIRLING_TABLE, "Stirling mismatch"
    A = airy(8)
    assert all(A[k] == airy_closed(k) for k in range(9)), "Airy mismatch"
    print("Airy saddle series through h^8 = (6k-1)!!/((2k)! 72^k):", ", ".join(str(x) for x in A))
    # second independent check of the engine: a Gaussian with a linear change of variable must be exact:
    # Int exp(-(y^2/2)/h) (1+y)^2 dy = sqrt(2 pi h)(1 + h)
    assert gauss_series(1, {}, [Fr(1), Fr(2), Fr(1)], 4) == [1, 1, 0, 0, 0]
    # reparametrisation invariance: Phi(y) = Y^2/2 with Y = y + y^2 (dy = dY/(1+2y)): Int exp(-Y(y)^2/(2h)) (1 + 2y) dy = sqrt(2 pi h)
    phi = {3: Fr(1), 4: Fr(1, 2)}
    assert gauss_series(1, phi, [Fr(1), Fr(2)], 6) == [1, 0, 0, 0, 0, 0, 0]
    print("exact-Gaussian and reparametrisation checks ok")
    for k in range(-2, 6):
        print("J_%d(eps) =" % k, ", ".join(str(x) for x in J_series(k, 7)))
    print("lg_sp ok")
