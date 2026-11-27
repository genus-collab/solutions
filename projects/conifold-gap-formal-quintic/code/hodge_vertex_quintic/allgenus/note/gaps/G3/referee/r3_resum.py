# -*- coding: utf-8 -*-
"""r3_resum.py -- referee check (independent code) of the scalar identities behind G3-2e (vertex resummation) and G3-2f
(dilaton re-expansion), with failing controls.  Laurent series in u are pairs (o, list): u^o * Sum list[i] u^i, mod u^(o+L).

(E)  v = u tau^(-5/4), delta = u^2/30 + (127/3150) u^3 + u^5 (any element of u^2 Q[[u]]), F(v) = v^(2-2g):
     S_J := Sum_(j<=J) delta^j F^(j)(v)/j!  ->  (v + delta)^(2-2g) u-adically, with ord(S_J - target) >= 3 - 2g + J, and
     S_(3g-3) - (v + delta)^(2-2g) in Q[[u]] with u-order >= g (the prover's "shorter route", tail bound);
     CONTROL E1: v = u^2 (v'(0) = 0), delta = u^2: every term has u-order 4 - 4g, the partial sums do not converge;
     CONTROL E2: v as above, delta = u/2 (delta not in u^2 Q[[u]]): every term has u-order 2 - 2g, no convergence.
(F)  kappa = 1 + u + u^2/3, c = 1 - 1/kappa: Sum_(m<=K) binom(N+m-1, m) c^m = kappa^N mod u^(K+1), N = 1, 2, 5;
     CONTROL F1: kappa = 2 + u: the u^0 coefficients of the partial sums change at every m (no u-adic limit);
     (F2) the series grouped by the total number m of insertions, T_m = Sum_(m2<=2) c^(m-m2) binom(N0+m2+m-m2-1, m-m2) A_(m2)/m2!,
     with c(0) = 1/2 and A = (A_0, A_1, A_2) non-zero constants: [u^0] T_m != 0 for every m <= 80 (the m-grouped Taylor series
     has no limit either; the prover's (f2) proves only the m1-series).
Usage (working directory: the scratch folder):  py -I -B r3_resum.py L   -- log on stdout, nothing written, no path printed.
"""
import sys
from fractions import Fraction as Fr
from math import comb, factorial


def mul(a, b):
    L = min(len(a), len(b))
    out = [Fr(0)] * L
    for i in range(L):
        if a[i]:
            for j in range(L - i):
                out[i + j] += a[i] * b[j]
    return out


def powr(a, r):
    """a^r for a[0] = 1, r in Q: via a g' = r a' g (power series ODE)."""
    L = len(a)
    da = [(i + 1) * a[i + 1] for i in range(L - 1)] + [Fr(0)]
    g = [Fr(1)]
    for n in range(L - 1):
        rhs = r * sum(da[i] * g[n - i] for i in range(n + 1))
        rhs -= sum(a[i] * (n - i + 1) * g[n - i + 1] for i in range(1, n + 1))
        g.append(rhs / ((n + 1) * a[0]))
    return g


def lmul(A, B):
    return (A[0] + B[0], mul(A[1], B[1]))


def lsub(A, B, top):
    """A - B known mod u^top, returned as (o, list) with o = min offset."""
    o = min(A[0], B[0])
    out = [Fr(0)] * (top - o)
    for (oo, s), sg in ((A, 1), (B, -1)):
        for i, x in enumerate(s):
            if oo + i < top:
                out[oo + i - o] += sg * x
    return (o, out)


def lord(A):
    o, s = A
    i = next((i for i, x in enumerate(s) if x), None)
    return None if i is None else o + i


def model_unit(L):
    F = [Fr(4, 5 * n + 4) for n in range(L)]
    return powr(F, Fr(5, 4))                               # v = u * w(u)


def resum(vo, w, delta, g, Jmax, top):
    """partial sums S_J (J <= Jmax) and the target (v + delta)^(2-2g) as Laurent series known mod u^top; v = u^vo w."""
    L = len(w)
    r = 2 - 2 * g
    dl = (delta[0] - vo, mul(delta[1], powr(w, Fr(-1))))     # delta / v
    target = None
    S, out = None, []
    for j in range(Jmax + 1):
        coef = Fr(1)
        for i in range(j):
            coef *= Fr(r - i, i + 1)                          # binom(r, j)
        term = (r * vo, [coef * x for x in powr(w, Fr(r))])
        for _ in range(j):
            term = lmul(term, dl)
        S = term if S is None else lsub(S, (term[0], [-x for x in term[1]]), top)
        out.append(S)
    # target: v^r (1 + delta/v)^r, valid as a series only if ord(delta/v) >= 1
    if dl[0] >= 1:
        base = [Fr(0)] * L
        base[0] = Fr(1)
        for i, x in enumerate(dl[1]):
            if dl[0] + i < L:
                base[dl[0] + i] += x
        target = (r * vo, mul(powr(w, Fr(r)), powr(base, Fr(r))))
    return out, target


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    print("r3_resum.py  L = %d (series length); exact over Q" % L)
    ok = True
    w = model_unit(L)
    delta = (2, [Fr(1, 30), Fr(127, 3150), Fr(0), Fr(1)] + [Fr(0)] * (L - 4))
    for g in (2, 3, 4):
        top = 2 - 2 * g + L - 2
        S, target = resum(1, w, delta, g, 3 * g - 3 + L, top)
        conv = all((lambda d: d is None or d >= 3 - 2 * g + J)(lord(lsub(S[J], target, top))) for J in range(len(S)))
        tail = lord(lsub(S[3 * g - 3], target, top))
        print("[E] g = %d: ord(S_J - (v+delta)^(2-2g)) >= 3 - 2g + J for all J: %s;  ord(S_(3g-3) - target) = %s (claim >= %d)"
              % (g, conv, tail, g))
        ok &= conv and tail is not None and tail >= g
    g = 2
    one = [Fr(1)] + [Fr(0)] * (L - 1)
    S1, _ = resum(2, one, (2, one[:]), g, 12, 6)              # v = u^2, delta = u^2
    c1 = [s[1][0] if s[0] == 4 - 4 * g else None for s in S1]
    print("[E1] CONTROL v = u^2, delta = u^2: lowest offset %d; the u^(4-4g) coefficient of S_J for J = 0..12: %s (no limit)"
          % (S1[0][0], [str(x) for x in c1]))
    ok &= len(set(c1)) == len(c1)
    S2, t2 = resum(1, w, (1, [Fr(1, 2)] + [Fr(0)] * (L - 1)), g, 12, 6)
    c2 = [x for x in (s[1][0] if s[0] == 2 - 2 * g else None for s in S2)]
    print("[E2] CONTROL delta = u/2: target defined: %s; the u^(2-2g) coefficient of S_J, J = 0..12: %s (no limit)"
          % (t2 is not None, [str(x) for x in c2]))
    ok &= t2 is None and len(set(c2)) == len(c2)
    # (F)
    kap = [Fr(1), Fr(1), Fr(1, 3)] + [Fr(0)] * (L - 3)
    ik = powr(kap, Fr(-1))
    c = [-x for x in ik]
    c[0] += 1                                                   # c = 1 - 1/kappa, ord >= 1
    for N in (1, 2, 5):
        acc, cm = [Fr(0)] * L, [Fr(1)] + [Fr(0)] * (L - 1)
        for m in range(L):
            acc = [a + comb(N + m - 1, m) * x for a, x in zip(acc, cm)]
            cm = mul(cm, c)
        kN = powr(kap, Fr(N))
        print("[F] N = %d: Sum_(m<L) binom(N+m-1,m) c^m = kappa^N mod u^%d: %s" % (N, L, acc == kN))
        ok &= acc == kN
    c0 = Fr(1) - Fr(1, 2)                                       # kappa(0) = 2
    parts, ps = [], Fr(0)
    for m in range(12):
        ps += comb(2 + m - 1, m) * c0 ** m
        parts.append(ps)
    print("[F1] CONTROL kappa = 2 + u, N = 2: u^0 parts of the partial sums: %s (all distinct: %s)"
          % ([str(x) for x in parts[:6]], len(set(parts)) == len(parts)))
    ok &= len(set(parts)) == len(parts)
    N0, A = 3, (Fr(2), Fr(-5), Fr(7))
    zeros = [m for m in range(81) if sum(c0 ** (m - m2) * comb(N0 + m2 + m - m2 - 1, m - m2) * A[m2] / factorial(m2)
                                         for m2 in range(3) if m2 <= m) == 0]
    print("[F2] grouped by m, c(0) = 1/2, A = (2, -5, 7), N0 = 3: m <= 80 with [u^0] T_m = 0: %s (expected none)" % zeros)
    ok &= not zeros
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


if __name__ == "__main__":
    main()
