"""lf_fit N: the global shape of c_n(L) for the quintic, from exact values at rational points (n <= N).

For each n, c_n(L) is computed exactly (one-variable full-Stirling series, lf_models.q_onevar) at rational L
and fitted by a Laurent polynomial with exponents e = 4n - 5k, k = -2..floor(4n/5) + 2: two slots ABOVE 4n and two
slots BELOW 0 are allowed. The fit uses as many points as unknowns and is confirmed exactly on three further points.
Then: the extra slots vanish (c_n in Q[L], deg <= 4n, exponents = 4n mod 5); c_n(1) = h_n (h from Q(zeta_5));
c_(5m)(0) = [z^(5m)] Ginf. This is evidence (VERIFIED for n <= N), by a method that uses neither a Laurent
expansion at L = 1 or 0 nor any chart. Failing controls: the series with the singular (j = 0) Stirling series
dropped admits no such fit; h with the twist weight +5 differs from c_n(1).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from lf_wick import wick_1
from lf_models import q_onevar
from lf_log import Log
from lf_quintic import h_definition, ginf

Z, O = Fr(0), Fr(1)
POINTS = [Fr(2), Fr(3), Fr(1, 2), Fr(-1, 2), Fr(3, 2), Fr(-2), Fr(1, 3), Fr(-3), Fr(5, 2), Fr(2, 3), Fr(-2, 3),
          Fr(4), Fr(-4), Fr(5, 3), Fr(-5, 2), Fr(3, 4), Fr(-3, 4), Fr(7, 2), Fr(4, 3), Fr(-4, 3), Fr(5), Fr(1, 4),
          Fr(-1, 4), Fr(6, 5), Fr(-6, 5), Fr(7, 3), Fr(-7, 3), Fr(5, 4), Fr(-5, 4), Fr(8, 3), Fr(-3, 2), Fr(7, 4)]


def solve(A, b):
    n = len(A)
    M = [row[:] + [bb] for row, bb in zip(A, b)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        inv = 1 / M[c][c]
        M[c] = [x * inv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[r][n] for r in range(n)]


def fit(n, vals):
    """vals: list of (L, c_n(L)). Returns (coeffs by exponent, confirmed?)."""
    K = (4 * n) // 5
    exps = [4 * n - 5 * k for k in range(-2, K + 3)]
    U = len(exps)
    A = [[L ** e for e in exps] for L, _ in vals[:U]]
    sol = solve(A, [v for _, v in vals[:U]])
    conf = all(sum(a * L ** e for a, e in zip(sol, exps)) == v for L, v in vals[U:U + 3])
    return dict(zip(exps, sol)), conf


def main():
    N = int(sys.argv[1])
    log = Log("lf_fit_%d" % N)
    need = (4 * N) // 5 + 5 + 3
    assert need <= len(POINTS), "not enough points"
    pts = POINTS[:need]
    series, series_x = [], []
    for L in pts:
        ph3, logA, C = q_onevar(L, N)
        series.append(wick_1(ph3, logA, C, N, Z, O))
    for L in pts[:8]:
        ph3, logA, C = q_onevar(L, 1, sing_stirling=False)
        series_x.append(wick_1(ph3, logA, C, 1, Z, O))
    log.say("evaluated c(z; L) through z^%d at %d rational points" % (N, need))
    h = h_definition(N)
    hp = h_definition(N, twist=5)
    G = ginf(N)
    ok_shape, ok_val, ok_inf, conf_all = True, True, True, True
    for n in range(1, N + 1):
        co, conf = fit(n, [(L, s[n]) for L, s in zip(pts, series)])
        conf_all &= conf
        extra = [e for e, a in co.items() if a != 0 and (e < 0 or e > 4 * n)]
        ok_shape &= conf and not extra
        c1 = sum(co.values())
        ok_val &= (c1 == h[n])
        c0 = co.get(0, Fr(0))
        if n % 5 == 0:
            ok_inf &= (c0 == G[n])
            log.say("n = %2d: c_n(0) = %s = [z^n] Ginf: %s" % (n, c0, c0 == G[n]))
        else:
            ok_inf &= (c0 == 0 and G[n] == 0)
        if n <= 3:
            log.say("n = %d: c_n(L) = %s" % (n, " + ".join("(%s) L^%d" % (a, e) for e, a in sorted(co.items()) if a)))
    log.check("fit confirmed on 3 extra points for every n <= %d" % N, conf_all)
    log.check("c_n in Q[L], exponents 4n - 5k, 0 <= 4n - 5k <= 4n (no negative or too-large exponent), n <= %d" % N, ok_shape)
    log.check("(Gamma): c_n(1) = h_n (h from its definition), n <= %d" % N, ok_val)
    log.check("orbifold point: c_n(0) = [z^n] Ginf (0 unless 5 | n), n <= %d" % N, ok_inf)
    co, conf = fit(1, [(L, s[1]) for L, s in zip(pts[:8], series_x)])
    extra = [e for e, a in co.items() if a != 0 and (e < 0 or e > 4)]
    log.control("singular Stirling series dropped: c_1 admits the Laurent fit", conf and not extra)
    log.control("twist weight +5: c_n(1) = h_n for all n", all(
        sum(fit(n, [(L, s[n]) for L, s in zip(pts, series)])[0].values()) == hp[n] for n in range(1, min(N, 3) + 1)))
    log.close()


if __name__ == "__main__":
    main()
