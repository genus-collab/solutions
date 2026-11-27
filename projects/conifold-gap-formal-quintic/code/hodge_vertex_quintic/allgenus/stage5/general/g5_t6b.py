"""g5_t6b.py -- tests F and G of g5_t6.py (called from there; see its docstring).

F  regular part NOT of Gamma type: B = exp(Phi/z + A_0 + Sum_(n>=1) z^n A_n) with random rational Taylor data (A_n for
   EVERY n >= 1, including even n, which never occur for Gamma factors), singular factor Gamma_St((s* - sigma)/z).
   Prediction (local Mellin inversion, the first pole is s* itself): no pole, value exp(Sum_n z^n A_n(s*)).
   Control F': the same with b_1 of the singular factor perturbed: poles.
F2 the same with a shifted singular factor Gamma(e* + c*(sigma - s*)/z), e* != 0: value = B at the first pole
   sigma = s* - e* z/c*, normalised: exp(Sum_m Phi_m delta^m/z + A_0(s*+delta) - A_0(s*) + Sum_k z^k A_k(s*+delta)).
G  coincident weights: M = 2 or 3 equal base weights, i.e. the singular factor Gamma_St((s* - sigma)/z)^M.
"""
import random
from fractions import Fraction as Fr
from math import factorial
import g5_core as C, g5_data as D


def build_raw(Phi2, A0p, An, sing, N, K):
    Phi = [Fr(0), Fr(0)] + [Fr(Phi2[m]) / ((m + 1) * (m + 2)) for m in range(K)]
    H0 = [C.cv(0)] + [C.cv(Fr(A0p[m]) / (m + 1)) for m in range(K)]
    Bk = [C.sexp(H0, K)]
    for n in range(1, N + 1):
        acc = [0] * (K + 1)
        for i in range(1, n + 1):
            pr = C.smul([C.cv(x) for x in An[i]], Bk[n - i], K)
            for t in range(K + 1):
                acc[t] += i * pr[t]
        Bk.append([C.red(x * C.finv(C.cv(n))) for x in acc])
    loc = dict(sing); loc['Phi'] = Phi; loc['Bk'] = Bk
    return loc


def main(N, out, check, FAILS):
    K = N
    rnd = random.Random(20261114)
    rr = lambda: Fr(rnd.randint(-9, 9), rnd.randint(1, 7))
    for trial in range(3):
        Phi2 = [rr() for _ in range(K + 1)]; A0p = [rr() for _ in range(K + 1)]
        An = [None] + [[rr() for _ in range(K + 1)] for _ in range(N)]
        sg = D.sing_base(N)
        ct = C.engine(build_raw(Phi2, A0p, An, sg, N, K), N, K)
        first, vals, poles = C.laurent_report(ct, N, K)
        lg = [Fr(0)] + [An[n][0] for n in range(1, N + 1)]
        E = [Fr(1)] + [Fr(0)] * N
        for m in range(1, N + 1):
            E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
        bad = [n for n in range(N + 1) if vals[n] != E[n]]
        check('F random non-Gamma WKB regular part, trial %d: no pole and value exp(Sum_n z^n A_n(s*)), n <= %d'
              % (trial, N), first is None and not bad, 'first pole %s, mismatches %s' % (first, bad[:4]))
        if trial == 0:
            sg2 = D.sing_base(N, lambda n: D.beta(0, n) + (Fr(1, 7) if n == 1 else 0))
            ct2 = C.engine(build_raw(Phi2, A0p, An, sg2, N, K), N, K)
            f2, _, _ = C.laurent_report(ct2, N, K)
            check('F\' control: same random data, singular factor b_1 + 1/7: pole (first at z^%s)' % f2, f2 is not None)
            sg3 = D.sing_base(N, lambda n: D.beta(0, n) + (Fr(1, 7) if n == 2 else 0))
            ct3 = C.engine(build_raw(Phi2, A0p, An, sg3, N, K), N, K)
            f3, _, _ = C.laurent_report(ct3, N, K)
            check('F\' control: singular factor with an extra even Stirling term (1/7) Y^-2: pole (first at z^%s)' % f3,
                  f3 is not None)
    # F2: random regular data, singular factor Gamma(e* + c*(sigma - s*)/z) with e* != 0: first pole at delta = -e* z/c*
    for cs, es in ((-1, Fr(1, 3)), (-1, Fr(1)), (2, Fr(1, 3)), (5, Fr(1))):
        Phi2 = [rr() for _ in range(K + 1)]; A0p = [rr() for _ in range(K + 1)]
        An = [None] + [[rr() for _ in range(K + 1)] for _ in range(N)]
        sg = dict(M=1, cs=Fr(cs), es=es, bstar=[0] + [D.beta(es, n) for n in range(1, N + 1)])
        ct = C.engine(build_raw(Phi2, A0p, An, sg, N, K), N, K)
        first, vals, poles = C.laurent_report(ct, N, K)
        dl = -es / cs                                     # delta = dl * z
        lg = [Fr(0)] * (N + 1)
        for m in range(2, N + 2):                         # Phi_m delta^m / z
            lg[m - 1] += Phi2[m - 2] / ((m - 1) * m) * dl ** m
        for m in range(1, N + 1):                         # A_0(s* + delta) - A_0(s*)
            lg[m] += A0p[m - 1] / m * dl ** m
        for k in range(1, N + 1):
            for m in range(0, N - k + 1):
                lg[k + m] += An[k][m] * dl ** m
        E = [Fr(1)] + [Fr(0)] * N
        for m in range(1, N + 1):
            E[m] = sum(i * lg[i] * E[m - i] for i in range(1, m + 1)) / m
        bad = [n for n in range(N + 1) if vals[n] != E[n]]
        check('F2 random non-Gamma regular part, singular Gamma(%s + %s(sigma - s*)/z): no pole, value = normalised '
              'B(s* + delta) at the first pole delta = %s z, n <= %d' % (es, cs, dl, N), first is None and not bad,
              'first pole %s, mismatches %s' % (first, bad[:4]))
    # G: coincident weights
    hank = [Fr(1)] + [Fr(_dfact(2 * n - 1) ** 2, factorial(n) * 16 ** n) for n in range(1, N + 1)]
    cases = [('P^1, lambda = (1, 1) [non-equivariant limit]', 2, []),
             ('P^2, lambda = (1, 1, -2)', 2, [D.group_rational(-1, 0, [-3], N + K + 1)]),
             ('P^2, lambda = (1, 1, 1) [non-equivariant]', 3, []),
             ('O(5)-twisted P^4, lambda = (1, 1, 3, -2, 5/2)', 2,
              [D.group_rational(-1, 0, [2, -3, Fr(3, 2)], N + K + 1), D.group_rational(5, 1, [5], N + K + 1)])]
    for name, M, gr in cases:
        sg = dict(M=M, cs=-1, es=0, bstar=[0] + [D.beta(0, n) for n in range(1, N + 1)])
        ct = C.engine(D.build_loc(gr, sg, N, K), N, K)
        first, vals, poles = C.laurent_report(ct, N, K)
        lead = [ct[n][0] for n in range(N + 1)]
        orders = [max([n - r for r in poles.get(n, [])] or [0]) for n in range(N + 1)]
        out('   G %s, M = %d: pole order of c_n = %s' % (name, M, orders[:12]))
        out('      eps^-n coefficients of c_n (n = 1..6): %s' % [str(x) for x in lead[1:7]])
        check('G %s: the removable singularity FAILS (poles from z^%s on), as the degenerate chart predicts; pole order '
              'exactly n at every z^n, 1 <= n <= %d: %s' % (name, first, N, all(orders[n] == n for n in range(1, N + 1))),
              first is not None)
        if M == 2:
            check('G   M = 2: eps^-n coefficient of c_n = Hankel coefficient ((2n-1)!!)^2/(n! 16^n) of K_0, n <= %d' % N,
                  lead == hank)


def _dfact(m):
    r = 1
    while m > 1:
        r *= m; m -= 2
    return r
