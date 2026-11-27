# Fresh Hodge engine (round 6d, backward read). No project code imported.
# WK numbers by DVV; ch-classes by the Mumford recursion. Exact rationals.
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, factorial
from itertools import product
sys.setrecursionlimit(100000)


def dfact(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


def bernoulli_list(N):
    B = [Fr(0)] * (N + 1)
    B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BERN = bernoulli_list(80)


def stable(g, n):
    return 2 * g - 2 + n > 0


def submultisets(ms):
    """ms sorted tuple. Yield (sub, complement, multiplicity)."""
    vals = sorted(set(ms))
    cnt = [ms.count(v) for v in vals]
    for choice in product(*[range(c + 1) for c in cnt]):
        mult = 1
        a = []
        b = []
        for v, c, k in zip(vals, cnt, choice):
            mult *= comb(c, k)
            a += [v] * k
            b += [v] * (c - k)
        yield tuple(a), tuple(b), mult


@lru_cache(maxsize=None)
def wk(g, ds):
    """<tau_{d_1}...tau_{d_n}>_g, ds sorted tuple."""
    n = len(ds)
    if not stable(g, n) or g < 0:
        return Fr(0)
    if sum(ds) != 3 * g - 3 + n:
        return Fr(0)
    if g == 0:
        r = Fr(factorial(n - 3))
        for d in ds:
            r /= factorial(d)
        return r
    if (g, ds) == (1, (1,)):
        return Fr(1, 24)
    if ds and ds[0] == 0:  # string
        rest = list(ds[1:])
        tot = Fr(0)
        for i in range(len(rest)):
            if rest[i] > 0:
                new = rest[:]
                new[i] -= 1
                tot += wk(g, tuple(sorted(new)))
        return tot
    if 1 in ds:  # dilaton
        l = list(ds)
        l.remove(1)
        return (2 * g - 2 + len(l)) * wk(g, tuple(l))
    if (g, ds) == (1, (1,)):
        return Fr(1, 24)
    if not ds:
        return Fr(0)
    # DVV on the largest
    l = list(ds)
    kp1 = l.pop()
    k = kp1 - 1
    rest = tuple(l)
    tot = Fr(0)
    for j in range(len(rest)):
        dj = rest[j]
        new = list(rest)
        new[j] = dj + k
        tot += Fr(dfact(2 * k + 2 * dj + 1), dfact(2 * dj - 1)) * wk(g, tuple(sorted(new)))
    for r in range(k):
        s = k - 1 - r
        c = dfact(2 * r + 1) * dfact(2 * s + 1)
        t = wk(g - 1, tuple(sorted(rest + (r, s))))
        for I, J, m in submultisets(rest):
            for g1 in range(g + 1):
                g2 = g - g1
                if stable(g1, len(I) + 1) and stable(g2, len(J) + 1):
                    if sum(I) + r != 3 * g1 - 3 + len(I) + 1:
                        continue
                    t += m * wk(g1, tuple(sorted(I + (r,)))) * wk(g2, tuple(sorted(J + (s,))))
        tot += Fr(c, 2) * t
    return tot / dfact(2 * k + 3)


@lru_cache(maxsize=None)
def hodge(g, chs, ds):
    """Int_{Mbar_{g,n}} prod ch_{c} prod psi^{d}. chs sorted odd degrees. ds sorted."""
    n = len(ds)
    if not stable(g, n):
        return Fr(0)
    if sum(chs) + sum(ds) != 3 * g - 3 + n:
        return Fr(0)
    if not chs:
        return wk(g, ds)
    if g == 0:
        return Fr(0)
    # string and dilaton (ch pulled back)
    if ds and ds[0] == 0 and stable(g, n - 1):
        rest = list(ds[1:])
        tot = Fr(0)
        for i in range(len(rest)):
            if rest[i] > 0:
                new = rest[:]
                new[i] -= 1
                tot += hodge(g, chs, tuple(sorted(new)))
        return tot
    if 1 in ds and stable(g, n - 1):
        l = list(ds)
        l.remove(1)
        return (2 * g - 2 + len(l)) * hodge(g, chs, tuple(l))
    # Mumford on the largest ch
    c = chs[-1]
    S = chs[:-1]
    l2 = (c + 1)  # 2l
    pref = BERN[l2] / factorial(l2)
    tot = hodge(g, S, tuple(sorted(ds + (l2,))))
    for i in range(n):
        new = list(ds)
        new[i] += c
        tot -= hodge(g, S, tuple(sorted(new)))
    bdry = Fr(0)
    for a in range(c):
        b = c - 1 - a
        sg = (-1) ** b
        t = hodge(g - 1, S, tuple(sorted(ds + (a, b)))) if g >= 1 else Fr(0)
        for I, J, m in submultisets(ds):
            for S1, S2, ms in submultisets(S):
                for g1 in range(g + 1):
                    g2 = g - g1
                    if not (stable(g1, len(I) + 1) and stable(g2, len(J) + 1)):
                        continue
                    if sum(I) + a + sum(S1) != 3 * g1 - 3 + len(I) + 1:
                        continue
                    v1 = hodge(g1, S1, tuple(sorted(I + (a,))))
                    if v1 == 0:
                        continue
                    t += m * ms * v1 * hodge(g2, S2, tuple(sorted(J + (b,))))
        bdry += sg * t
    tot += bdry / 2
    return pref * tot


def ch_partitions(deg):
    """multisets of odd numbers summing to deg, as sorted tuples"""
    def rec(rem, mx):
        if rem == 0:
            yield ()
            return
        for c in range(min(mx, rem), 0, -1):
            if c % 2 == 1:
                for tail in rec(rem - c, c):
                    yield (c,) + tail
    return [tuple(sorted(x)) for x in rec(deg, deg)]


if __name__ == '__main__':
    assert wk(1, (1,)) == Fr(1, 24)
    assert wk(2, (4,)) == Fr(1, 1152)
    assert wk(2, (2, 3)) == Fr(29, 5760)
    assert wk(3, (7,)) == Fr(1, 82944)
    assert wk(2, (2, 2, 2)) == Fr(7, 240)
    # lambda_1 = ch_1
    assert hodge(1, (1,), (0,)) == Fr(1, 24)
    # lambda_1^3 on Mbar_2 = 1/2880
    assert hodge(2, (1, 1, 1), ()) == Fr(1, 2880)
    # 3! ch_3 = -1/5760 on Mbar_2
    assert 6 * hodge(2, (3,), ()) == Fr(-1, 5760)
    # lambda_1 lambda_2 lambda_3 on Mbar_3 = 1/1451520 via 5! ch_5 = lambda_3 lambda_2 (FB A3.3, g=3)
    # (2g-1)! ch_{2g-1} = (-1)^{g-1} lambda_g lambda_{g-1}, so Int ch_1 ch_5 *120 = Int lambda_1 lambda_2 lambda_3
    assert 120 * hodge(3, (1, 5), ()) == Fr(1, 1451520)
    # lambda_g formula g=2: Int lambda_2 psi^2 on Mbar_{2,1} = 7/5760; lambda_2 = ch_1^2/2
    assert hodge(2, (1, 1), (2,)) / 2 == Fr(7, 5760)
    print('engine controls pass')
