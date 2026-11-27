"""lf_dlmf K N: DLMF 5.11.8 as an identity of formal series, symbolically in h; the step LMI value -> Ginf.

D1  With x = 1/w: (1/x + h - 1/2) log(1 + h x) - h + Sum_(k>=1) b_k x^(2k-1) (1 + h x)^(1-2k)
      = Sum_(k>=2) (-1)^k B_k(h) x^(k-1)/(k(k-1)),
    i.e. DLMF 5.11.1 at the argument w + h, re-expanded in 1/w, equals DLMF 5.11.8. Checked as polynomials in h,
    coefficientwise through x^K (exact over Q).
D2  The difference equation behind it: B_k(h + 1) - B_k(h) = k h^(k-1), so the right side changes by log(1 + h x)
    under h -> h + 1, as Stirling's functional equation requires (k <= K).
D3  The quintic at s* = 0 (alpha0 = 5, M0 = 1, delta = -z/5): the LMI value exp(gt(delta)/z) F(delta, z)/F(0, 0), with
    g, F the five factors Gamma((zeta^j - sigma)/z), computed from their Taylor data at sigma = 0, equals Ginf, z^N.
Failing controls: the sign (-1)^(k+1); (h + 1/2) for (h - 1/2); delta = +z/5; h = 0 (no shift).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb
from lf_core import bern, bk, gbinom, power_sums_all
from lf_models import lmi_predict
from lf_log import Log, first_diff
from lf_quintic import ginf


def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pscale(p, c, shift=0):
    """c h^shift p(h)."""
    return [Fr(0)] * shift + [c * a for a in p]


def bern_poly_coeffs(n):
    """B_n(h) as a coefficient list in h."""
    out = [Fr(0)] * (n + 1)
    for k in range(n + 1):
        out[n - k] += comb(n, k) * bern(k)
    return out


def lhs(K, half=Fr(-1, 2)):
    """coefficients (polynomials in h) of x^0..x^K of the re-expanded Stirling series."""
    out = [[Fr(0)] for _ in range(K + 1)]
    for m in range(2, K + 2):            # (1/x) log(1 + h x) - h
        out[m - 1] = padd(out[m - 1], pscale([Fr(1)], Fr((-1) ** (m + 1), m), m))
    for m in range(1, K + 1):            # (h + half) log(1 + h x)
        c = Fr((-1) ** (m + 1), m)
        out[m] = padd(out[m], padd(pscale([Fr(1)], c, m + 1), pscale([Fr(1)], c * half, m)))
    k = 1
    while 2 * k - 1 <= K:                # b_k x^(2k-1) (1 + h x)^(1-2k)
        for r in range(0, K - (2 * k - 1) + 1):
            out[2 * k - 1 + r] = padd(out[2 * k - 1 + r], pscale([Fr(1)], bk(k) * gbinom(1 - 2 * k, r), r))
        k += 1
    return out


def rhs(K, sign=1):
    out = [[Fr(0)] for _ in range(K + 1)]
    for k in range(2, K + 2):
        out[k - 1] = padd(out[k - 1], pscale(bern_poly_coeffs(k), sign * Fr((-1) ** k, k * (k - 1))))
    return out


def norm(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def same(a, b):
    return all(norm(x) == norm(y) for x, y in zip(a, b))


def main():
    K, N = int(sys.argv[1]), int(sys.argv[2])
    log = Log("lf_dlmf_%d_%d" % (K, N))
    L, R = lhs(K), rhs(K)
    log.check("D1 DLMF 5.11.1 at w + h re-expanded = DLMF 5.11.8, as polynomials in h, through x^%d" % K, same(L, R))
    log.say("D1 x^1: %s ; x^2: %s  (coefficients of h^0, h^1, ...)" % (norm(R[1]), norm(R[2])))
    log.control("D1 sign (-1)^(k+1) in 5.11.8", same(L, rhs(K, sign=-1)))
    log.control("D1 (h + 1/2) for (h - 1/2)", same(lhs(K, half=Fr(1, 2)), R))
    ok = True
    for k in range(1, K + 2):
        p, q = bern_poly_coeffs(k), [Fr(0)] * (k + 1)
        for i, a in enumerate(p):        # B_k(h + 1)
            for r in range(i + 1):
                q[r] += a * comb(i, r)
        diff = norm([x - y for x, y in zip(q, p)])
        ok &= (diff == norm([Fr(0)] * (k - 1) + [Fr(k)]))
    log.check("D2 B_k(h + 1) - B_k(h) = k h^(k-1), k <= %d" % (K + 1), ok)
    # D3: g, F of the five regular factors at sigma = 0
    P = power_sums_all(Fr(0), 2 * N + 4)
    G = {k: P[k - 1] / (k * (k - 1)) for k in range(2, N + 2)}
    LF = {}
    for d in range(1, N + 1):
        LF[(d, 0)] = P[d] / (2 * d)
    k = 1
    while 2 * k - 1 <= N:
        for d in range(0, N - (2 * k - 1) + 1):
            LF[(d, 2 * k - 1)] = LF.get((d, 2 * k - 1), 0) + bk(k) * comb(2 * k - 2 + d, d) * P[2 * k - 1 + d]
        k += 1
    V = lmi_predict(G, LF, 5, 1, N)
    Gi = ginf(N)
    log.check("D3 LMI value at s* = 0 (delta = -z/5) = Ginf, through z^%d" % N, V == Gi)
    Vx = lmi_predict(G, LF, 5, 1, N, sign=-1)
    log.control("D3 delta = +z/5", Vx == Gi, "first diff z^%s" % first_diff(Vx, Gi))
    V0 = lmi_predict(G, LF, 5, 0, N)
    log.control("D3 no shift (M0 = 0)", V0 == Gi, "first diff z^%s" % first_diff(V0, Gi))
    log.say("D3 [z^5], [z^10], [z^15] of Ginf: %s, %s, %s" % (Gi[5], Gi[10] if N >= 10 else "-", Gi[15] if N >= 15 else "-"))
    log.close()


if __name__ == "__main__":
    main()
