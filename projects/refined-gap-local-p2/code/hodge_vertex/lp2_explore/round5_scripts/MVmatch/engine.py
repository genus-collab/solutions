# engine.py -- Hodge side, written from scratch for the MV matching (round 5, key MVmatch).
# Nothing here uses the MV formula or any Schur function.
#   wk(g, a)         Witten-Kontsevich numbers, by the DVV recursion
#   hodge(g, a, ch)  Int_{Mbar_{g,n}} prod psi_i^{a_i} prod ch_k, Mumford's GRR formula (Faber's algorithm)
#   omega_monos      Lambda(t1)Lambda(t2)Lambda(t3) as a polynomial in the ch_k
#   Mg, Qg, Xi       the integrals of the Lemma A setting, computed from the definition
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import combinations
from math import comb, factorial
import sys
sys.setrecursionlimit(100000)


def dfact(n):
    r = 1
    while n > 1:
        r *= n
        n -= 2
    return r


_bern = [Fr(1)]


def bern(n):
    while len(_bern) <= n:
        m = len(_bern)
        _bern.append(-sum(comb(m + 1, k) * _bern[k] for k in range(m)) / (m + 1))
    return _bern[n]



def split_ms(a):
    """all ways to split the sorted tuple a into (I, J), I a sub-multiset, with multiplicity prod binom(m_v, c_v)"""
    vals = {}
    for v in a:
        vals[v] = vals.get(v, 0) + 1
    out = [((), (), 1)]
    for v, m in vals.items():
        new = []
        for I, J, w in out:
            for c in range(m + 1):
                new.append((I + (v,) * c, J + (v,) * (m - c), w * comb(m, c)))
        out = new
    return out


def _stable(g, n):
    return 2 * g - 2 + n > 0

@lru_cache(maxsize=None)
def _wk(g, a):
    n = len(a)
    if g < 0 or 2 * g - 2 + n <= 0:
        return Fr(0)
    if sum(a) != 3 * g - 3 + n or (n and a[0] < 0):
        return Fr(0)
    if g == 0 and a == (0, 0, 0):
        return Fr(1)
    if g == 1 and a == (1,):
        return Fr(1, 24)
    if a[0] == 0 and _stable(g, n - 1):      # string
        rest0 = a[1:]
        tot0 = Fr(0)
        vals = {}
        for v in rest0:
            vals[v] = vals.get(v, 0) + 1
        for v, m in vals.items():
            if v > 0:
                l0 = list(rest0); l0.remove(v); l0.append(v - 1)
                tot0 += m * wk(g, tuple(l0))
        return tot0
    if 1 in a and _stable(g, n - 1):         # dilaton
        l0 = list(a); l0.remove(1)
        return (2 * g - 2 + n - 1) * wk(g, tuple(l0))
    # DVV with k + 1 = largest exponent, k >= -1
    d = a[-1]
    rest = a[:-1]
    k = d - 1
    tot = Fr(0)
    # first sum
    for j in range(len(rest)):
        dj = rest[j]
        if k == -1 and dj == 0:
            continue
        new = rest[:j] + (dj + k,) + rest[j + 1:]
        if dj + k < 0:
            continue
        tot += Fr(dfact(2 * k + 2 * dj + 1), dfact(2 * dj - 1)) * wk(g, new)
    # second and third sums
    for r in range(0, k):
        s = k - 1 - r
        c = Fr(dfact(2 * r + 1) * dfact(2 * s + 1), 2)
        tot += c * wk(g - 1, rest + (r, s))
        for I, J, w in split_ms(rest):
            for g1 in range(0, g + 1):
                g2 = g - g1
                if sum(I) + r != 3 * g1 - 3 + len(I) + 1:
                    continue
                v1 = wk(g1, I + (r,))
                if v1 == 0:
                    continue
                tot += w * c * v1 * wk(g2, J + (s,))
    return tot / dfact(2 * k + 3)


def wk(g, a):
    return _wk(g, tuple(sorted(a)))


@lru_cache(maxsize=None)
def _hodge(g, a, ch):
    n = len(a)
    if 2 * g - 2 + n <= 0:
        return Fr(0)
    if sum(a) + sum(ch) != 3 * g - 3 + n:
        return Fr(0)
    if not ch:
        return wk(g, a)
    if g == 0:
        return Fr(0)      # rank 0 Hodge bundle
    if 0 in a and _stable(g, n - 1):         # string: the ch_k are pulled back
        l0 = list(a); l0.remove(0)
        tot0 = Fr(0)
        vals = {}
        for v in l0:
            vals[v] = vals.get(v, 0) + 1
        for v, m in vals.items():
            if v > 0:
                l1 = list(l0); l1.remove(v); l1.append(v - 1)
                tot0 += m * hodge(g, tuple(l1), ch)
        return tot0
    if 1 in a and _stable(g, n - 1):         # dilaton
        l0 = list(a); l0.remove(1)
        return (2 * g - 2 + n - 1) * hodge(g, tuple(l0), ch)
    k = ch[-1]
    rest = ch[:-1]
    assert k % 2 == 1
    l = (k + 1) // 2
    tot = hodge(g, a + (2 * l,), rest)                       # kappa_{2l-1}
    for i in range(n):
        tot -= hodge(g, a[:i] + (a[i] + 2 * l - 1,) + a[i + 1:], rest)    # psi_i^{2l-1}
    bd = Fr(0)
    m = len(rest)
    for x in range(0, 2 * l - 1):
        y = 2 * l - 2 - x
        sgn = -1 if x % 2 else 1
        v = hodge(g - 1, a + (x, y), rest)                   # irreducible divisor (degree 2 map)
        for I, J, wI in split_ms(a):
            for S, T, wS in split_ms(rest):
                dI = sum(I) + sum(S) + x
                # g1 is forced by the dimension
                num = dI + 3 - (len(I) + 1)
                if num % 3:
                    continue
                g1 = num // 3
                if g1 < 0 or g1 > g:
                    continue
                v1 = hodge(g1, I + (x,), S)
                if v1 == 0:
                    continue
                v += wI * wS * v1 * hodge(g - g1, J + (y,), T)
        bd += sgn * v
    tot += bd / 2
    return tot * bern(2 * l) / factorial(2 * l)


def hodge(g, a, ch=()):
    return _hodge(g, tuple(sorted(a)), tuple(sorted(ch)))


# ---- Omega as a polynomial in the ch_k
def partitions_odd(d, mx=None):
    """partitions of d into odd parts <= mx, as tuples (increasing)."""
    if mx is None:
        mx = d
    if d == 0:
        yield ()
        return
    for k in range(min(d, mx), 0, -1):
        if k % 2 == 0:
            continue
        for p in partitions_odd(d - k, k):
            yield p + (k,)


def omega_monos(g, d, ts):
    """degree d part of prod_t Lambda_g(t) = exp(sum_{k odd} (k-1)! p_k ch_k), p_k = sum t^k.
    returns a list of (ch tuple, coefficient). Uses ch_k = 0 for k >= 2g (Mumford)."""
    out = []
    if g == 0:
        return [((), Fr(1))] if d == 0 else []
    for p in partitions_odd(d, 2 * g - 1):
        c = Fr(1)
        mult = {}
        for k in p:
            mult[k] = mult.get(k, 0) + 1
        for k, m in mult.items():
            pk = sum(Fr(t) ** k for t in ts)
            c *= (factorial(k - 1) * pk) ** m / factorial(m)
        if c != 0:
            out.append((tuple(sorted(p)), c))
    return out


def ts_of(tau):
    """t-parameters of the normalised class Omega = Lambda(-1) Lambda(-1/tau) Lambda(1/(1+tau))"""
    tau = Fr(tau)
    return (Fr(-1), -1 / tau, 1 / (1 + tau))


def compositions(total, n):
    if n == 0:
        if total == 0:
            yield ()
        return
    if n == 1:
        yield (total,)
        return
    for i in range(total + 1):
        for rest in compositions(total - i, n - 1):
            yield (i,) + rest


def group_exps(legs, total):
    """Sum over exponents (a_i), sum a_i = total, of prod legs_i^{a_i}: grouped by equal legs.
    yields (exponent tuple, weight)."""
    groups = {}
    for v in legs:
        groups[Fr(v)] = groups.get(Fr(v), 0) + 1
    items = list(groups.items())
    def rec(i, rem):
        if i == len(items):
            if rem == 0:
                yield (), Fr(1)
            return
        v, c = items[i]
        for t in range(rem + 1):
            for e in parts_le(t, c):
                w0 = v ** t * _perm(e, c)
                if w0 == 0:
                    continue
                for tail, w1 in rec(i + 1, rem - t):
                    yield e + tail, w0 * w1
    yield from rec(0, total)


def parts_le(t, c):
    """nonincreasing c-tuples of nonnegative integers with sum t"""
    if c == 0:
        if t == 0:
            yield ()
        return
    def r(rem, k, mx):
        if k == 1:
            if rem <= mx:
                yield (rem,)
            return
        for v in range(min(rem, mx), -1, -1):
            for tail in r(rem - v, k - 1, v):
                yield (v,) + tail
    yield from r(t, c, t)


def _perm(e, c):
    m = {}
    for v in e:
        m[v] = m.get(v, 0) + 1
    r = factorial(c)
    for k in m.values():
        r //= factorial(k)
    return r

def leg_integral(g, legs, ts, extra_psi=()):
    """Int Omega_ts prod_i 1/(1 - legs_i psi_i) prod psi^{extra_psi}, over Mbar_{g, len(legs)+len(extra_psi)}.
    legs are numbers."""
    n = len(legs) + len(extra_psi)
    if 2 * g - 2 + n <= 0:
        raise ValueError("unstable")
    dim = 3 * g - 3 + n - sum(extra_psi)
    if dim < 0:
        return Fr(0)
    tot = Fr(0)
    for d in range(0, min(dim, 3 * g - 3 + n) + 1):
        monos = omega_monos(g, d, ts)
        if not monos:
            continue
        for A, w in group_exps(legs, dim - d):
            AA = tuple(A) + tuple(extra_psi)
            for ch, c in monos:
                tot += w * c * hodge(g, AA, ch)
    return tot


def Mg(g, mu, tau=None, ts=None):
    """M_g(mu) = Int Omega prod 1/(1-mu_i psi_i), with the LLZ/OP convention for the unstable cases."""
    mu = tuple(mu)
    l = len(mu)
    if g == 0 and l <= 2:
        if l == 0:
            raise ValueError("no (0,0)")
        return Fr(sum(mu)) ** (l - 3)
    if ts is None:
        ts = ts_of(tau)
    return leg_integral(g, mu, ts)


def _exps(total, cnt, mn):
    """nonincreasing lists of cnt integers >= mn with sum total, with multinomial multiplicity."""
    if cnt == 0:
        if total == 0:
            yield ()
        return
    def rec(rem, c, mx):
        if c == 0:
            if rem == 0:
                yield ()
            return
        for v in range(min(rem - mn * (c - 1), mx), mn - 1, -1):
            for t in rec(rem - v, c - 1, v):
                yield (v,) + t
    yield from rec(total, cnt, total)


def _mult(e):
    m = {}
    for v in e:
        m[v] = m.get(v, 0) + 1
    r = factorial(len(e))
    for c in m.values():
        r //= factorial(c)
    return r


def Qk(g, R, k, ts):
    """<Omega lam_R Tt^k>_{g, r+k}, Tt = psi^2/(1-psi), computed from the definition"""
    r = len(R)
    n = r + k
    if 2 * g - 2 + n <= 0:
        return Fr(0)
    dim = 3 * g - 3 + n
    tot = Fr(0)
    for E in range(2 * k, dim + 1):
        for e in _exps(E, k, 2):
            tot += _mult(e) * leg_integral(g, R, ts, extra_psi=e)
    return tot


def Qg_coeffs(g, R, ts, kmax=None):
    """coefficients [s^k] Q_g(s;R), k = 0..kmax, from the definition"""
    r = len(R)
    N = 3 * g - 3 + r
    if kmax is None:
        kmax = N
    return [Qk(g, R, k, ts) / factorial(k) if k <= N else Fr(0) for k in range(kmax + 1)]


def Xi_direct(g, R, ts):
    return sum(c * (-1) ** k for k, c in enumerate(Qg_coeffs(g, R, ts)))


if __name__ == "__main__":
    tests = [
        ("<tau_4>_2", wk(2, (4,)), Fr(1, 1152)),
        ("<tau_2 tau_3>_2", wk(2, (2, 3)), Fr(29, 5760)),
        ("<tau_2^3>_2", wk(2, (2, 2, 2)), Fr(7, 240)),
        ("<tau_7>_3", wk(3, (7,)), Fr(1, 82944)),
        ("<tau_1>_1", wk(1, (1,)), Fr(1, 24)),
        ("<tau_0 tau_2>_1", wk(1, (0, 2)), Fr(1, 24)),
        ("lam1 on M11", hodge(1, (0,), (1,)), Fr(1, 24)),
        ("ch3 on M_{1,3} psi^0 (must vanish)", hodge(1, (0, 0, 0), (3,)), Fr(0)),
        ("ch3 on M_{1,4} psi (must vanish)", hodge(1, (1, 0, 0, 0), (3,)), Fr(0)),
        ("lam1^3 on M2 = 1/2880", hodge(2, (), (1, 1, 1)), Fr(1, 2880)),
        ("ch3 on M2 = lam1^3/... : ch3 = (lam1^3 - 3 lam1 lam2 + 3 lam3)/6, lam2 = lam1^2/2",
         hodge(2, (), (3,)), (Fr(1, 2880) - Fr(3, 2) * Fr(1, 2880)) / 6),
        ("ch5 on M_{2,3} (must vanish, ch_k = 0 for k >= 2g)", hodge(2, (0, 0), (5,)), Fr(0)),
    ]
    # lambda_g psi^{2g-2}, g = 1,2,3: b_g = 1/24, 7/5760, 31/967680 (Faber-Pandharipande).
    # lambda_g = coefficient: Lambda(t) = exp(sum_{k odd}(k-1)! ch_k t^k), read off degree g.
    def lam(g, A):
        # Int lambda_g prod psi^A
        tot = Fr(0)
        for ch, c in omega_monos(g, g, (Fr(1),)):
            tot += c * hodge(g, A, ch)
        return tot
    tests.append(("b_1", lam(1, (0,)), Fr(1, 24)))
    tests.append(("b_2", lam(2, (2,)), Fr(7, 5760)))
    tests.append(("b_3", lam(3, (4,)), Fr(31, 967680)))
    # A_g = Int lam_{g-2} lam_{g-1} lam_g, closed form quoted in OP s.0.6 from Faber-Pandharipande
    def Ag_closed(g):
        return Fr(1, 2 * factorial(2 * g - 2)) * abs(bern(2 * g - 2)) / (2 * g - 2) * abs(bern(2 * g)) / (2 * g)
    # Int prod Lambda^vee(w_a)/(w1 w2 w3)^h = (-1)^h A_h/(w1w2w3)^{h-1}: test as <Omega>_h = -A_h sigma^{1-h}
    for g, tau in ((2, Fr(2, 7)), (2, Fr(-7, 3)), (3, Fr(2, 7)), (3, Fr(3, 2)), (4, Fr(2, 7))):
        sig = tau * (1 + tau)
        ts = ts_of(tau)
        d = 3 * g - 3
        val = sum(c * hodge(g, (), ch) for ch, c in omega_monos(g, d, ts))
        tests.append((f"<Omega>_{g} at tau={tau} = -A_g sigma^(1-g)", val, -Ag_closed(g) * sig ** (1 - g)))
    bad = 0
    for name, got, exp in tests:
        ok = got == exp
        bad += not ok
        print(("ok  " if ok else "FAIL"), name, got, "" if ok else f"expected {exp}")
    print("engine self-test:", "PASS" if not bad else f"{bad} FAILURES")
