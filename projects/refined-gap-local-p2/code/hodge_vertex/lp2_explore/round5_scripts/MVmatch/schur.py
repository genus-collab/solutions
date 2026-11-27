# schur.py -- Schur side, written from scratch for the MV matching (round 5, key MVmatch).
# Nothing here uses Hodge integrals.
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial


def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for p in partitions(n - k, k):
            yield (k,) + p


def conj(nu):
    return tuple(sum(1 for v in nu if v > j) for j in range(nu[0])) if nu else ()


def hooks(nu):
    c = conj(nu)
    return [nu[i] - j + c[j] - i - 1 for i in range(len(nu)) for j in range(nu[i])]


def kappa(nu):
    # LLZ: |nu| + sum (nu_i^2 - 2 i nu_i), i from 1
    return sum(nu) + sum(v * v - 2 * (i + 1) * v for i, v in enumerate(nu))


@lru_cache(maxsize=None)
def chi(nu, mu):
    """character by Murnaghan-Nakayama on beta numbers. mu a tuple (any order); order is irrelevant"""
    if not mu:
        return 1 if not nu else 0
    k = mu[0]
    rest = mu[1:]
    L = len(nu)
    beta = [nu[i] + L - 1 - i for i in range(L)]
    S = set(beta)
    tot = 0
    for i, b in enumerate(beta):
        t = b - k
        if t < 0 or t in S:
            continue
        ht = sum(1 for bb in beta if t < bb < b)
        nb = sorted([x if x != b else t for x in beta], reverse=True)
        nl = tuple(v for v in (nb[j] - (L - 1 - j) for j in range(L)) if v > 0)
        tot += (-1) ** ht * chi(nl, rest)
    return tot


def zee(mu):
    z = 1
    m = {}
    for v in mu:
        m[v] = m.get(v, 0) + 1
    for k, c in m.items():
        z *= k ** c * factorial(c)
    return z


def aut(mu):
    a = 1
    m = {}
    for v in mu:
        m[v] = m.get(v, 0) + 1
    for c in m.values():
        a *= factorial(c)
    return a


def Pi(m, tau):
    r = Fr(1)
    for a in range(1, m):
        r *= (m * Fr(tau) + a)
    return r / factorial(m - 1)


def PiMu(mu, tau):
    r = Fr(1)
    for m in mu:
        r *= Pi(m, tau)
    return r


# truncated power series as python lists of Fractions (length P+1)
def smul(a, b, P):
    c = [Fr(0)] * (P + 1)
    for i, x in enumerate(a):
        if x == 0 or i > P:
            continue
        for j, y in enumerate(b):
            if i + j > P:
                break
            if y != 0:
                c[i + j] += x * y
    return c


def sexp_linear(e, P):
    """exp(e x) truncated"""
    out = [Fr(1)]
    for n in range(1, P + 1):
        out.append(out[-1] * e / n)
    return out


_invsh = {}


def inv_sinh_series(h, P):
    """x h / (2 sinh(h x /2)) truncated at x^P (even series, constant term 1)"""
    key = (h, P)
    if key in _invsh:
        return _invsh[key]
    # 2 sinh(hx/2)/(hx) = sum_k (hx/2)^{2k}/(2k+1)!
    a = [Fr(0)] * (P + 1)
    for k in range(0, P // 2 + 1):
        a[2 * k] = (Fr(h, 2)) ** (2 * k) / factorial(2 * k + 1)
    b = [Fr(0)] * (P + 1)
    b[0] = Fr(1)
    for n in range(1, P + 1):
        b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1))
    _invsh[key] = b
    return b


_hookpart = {}


def hookpart(nu, P):
    """x^{|nu|} prod_box 1/(2 sinh(h x/2)) truncated at x^P"""
    key = (nu, P)
    if key in _hookpart:
        return _hookpart[key]
    c = conj(nu) if nu else ()
    if key[0] != c and (c, P) in _hookpart:
        _hookpart[key] = _hookpart[(c, P)]
        return _hookpart[key]
    s = [Fr(0)] * (P + 1)
    s[0] = Fr(1)
    H = 1
    for h in hooks(nu):
        H *= h
        s = smul(s, inv_sinh_series(h, P), P)
    s = [v / H for v in s]
    _hookpart[key] = s
    return s


def Wseries(nu, tau, P, frame_shift=Fr(1, 2)):
    """x^{|nu|} e^{(tau+shift) kappa x/2} prod_box 1/(2 sinh(h x/2)) as a truncated series"""
    e = (Fr(tau) + frame_shift) * kappa(nu) / 2
    return smul(hookpart(nu, P), sexp_linear(e, P), P)


def zmu_series(mu, tau, P, frame_shift=Fr(1, 2)):
    """[p_mu] Z = x^{-|mu|} * series, Z = sum_nu e^{(tau+shift) kappa x/2} prod 1/(2 sinh) s_nu(p),
       s_nu = sum_mu chi_nu(mu) p_mu / z_mu."""
    n = sum(mu)
    tot = [Fr(0)] * (P + 1)
    zz = zee(mu)
    mu = tuple(sorted(mu, reverse=True))
    for nu in partitions(n):
        ch = chi(nu, mu)
        if ch == 0:
            continue
        w = Wseries(nu, tau, P, frame_shift)
        for i in range(P + 1):
            if w[i] != 0:
                tot[i] += ch * w[i]
    return [t / zz for t in tot]


def sub_multisets(mu):
    """distinct sub-multisets of the tuple mu (sorted decreasing)"""
    m = {}
    for v in mu:
        m[v] = m.get(v, 0) + 1
    keys = sorted(m, reverse=True)
    out = [()]
    for k in keys:
        out = [o + (k,) * c for o in out for c in range(m[k] + 1)]
    return out


def msub(mu, rho):
    l = list(mu)
    for v in rho:
        l.remove(v)
    return tuple(l)


class ConnSeries:
    """[p_mu] log Z through the Euler recursion |mu| c_mu = |mu| z_mu - sum_{rho, rho != empty, mu} |rho| c_rho z_{mu-rho}.
       c_mu = x^{-|mu|} * series, all series truncated at P."""
    def __init__(self, tau, P, frame_shift=Fr(1, 2)):
        self.tau = Fr(tau)
        self.P = P
        self.fs = frame_shift
        self.z = {}
        self.c = {}

    def Z(self, mu):
        mu = tuple(sorted(mu, reverse=True))
        if mu not in self.z:
            self.z[mu] = zmu_series(mu, self.tau, self.P, self.fs)
        return self.z[mu]

    def C(self, mu):
        mu = tuple(sorted(mu, reverse=True))
        if mu in self.c:
            return self.c[mu]
        n = sum(mu)
        tot = list(self.Z(mu))
        for rho in sub_multisets(mu):
            if not rho or rho == mu:
                continue
            t = msub(mu, rho)
            f = Fr(sum(rho), n)
            pr = smul(self.C(rho), self.Z(t), self.P)
            for i in range(self.P + 1):
                tot[i] -= f * pr[i]
        self.c[mu] = tot
        return tot

    def coeff(self, mu, k):
        """[x^k][p_mu] log Z"""
        idx = k + sum(mu)
        if idx < 0:
            return Fr(0)
        if idx > self.P:
            raise ValueError("precision")
        return self.C(mu)[idx]


if __name__ == "__main__":
    # orthogonality and hook-length formula
    for n in range(1, 9):
        ps = list(partitions(n))
        for mu in ps:
            assert sum(chi(nu, mu) ** 2 for nu in ps) == zee(mu), ("orth", mu)
        for nu in ps:
            H = 1
            for h in hooks(nu):
                H *= h
            assert Fr(factorial(n), H) == chi(nu, (1,) * n)
    print("characters: column orthogonality and hook length formula OK, n <= 8")
    # LLZ (3) versus the hook product, as a pure multiset identity (LLZ Prop 2.1, Macdonald I.1 Ex.1)
    bad = 0
    tot = 0
    for n in range(1, 13):
        for nu in partitions(n):
            l = len(nu)
            num = []
            den = []
            for a in range(1, l + 1):
                for b in range(a + 1, l + 1):
                    num.append(nu[a - 1] - nu[b - 1] + b - a)
                    den.append(b - a)
            for i in range(1, l + 1):
                for v in range(1, nu[i - 1] + 1):
                    den.append(v - i + l)
            # num/den must equal 1/prod hooks as a multiset identity:
            lhs = sorted(num + hooks(nu))
            rhs = sorted(den)
            tot += 1
            bad += lhs != rhs
    print(f"LLZ (3) = prod 1/(2 sin(h lambda/2)) as multisets of arguments: {tot} partitions, n <= 12, failures {bad}")
    # q^rho evaluation: s_nu(q^rho) = q^{kappa/4} prod 1/[h], tested by numeric principal specialisation via
    # the identity sum_nu s_nu(q^rho) chi_nu(mu)/z_mu = prod_k (1/(q^{k/2}-q^{-k/2})) / z_mu   (p_k(q^rho) = 1/[k]).
    # Here: at tau = 0, frame e^{kappa x/4}, [p_mu] Z must be prod_i 1/(2 sinh(mu_i x/2)) / z_mu.
    P = 16
    bad = 0
    tot = 0
    for n in range(1, 9):
        for mu in partitions(n):
            s = zmu_series(mu, 0, P)
            r = [Fr(0)] * (P + 1)
            r[0] = Fr(1)
            for m in mu:
                r = smul(r, [v * Fr(1, m) for v in inv_sinh_series(m, P)], P)
            r = [v / zee(mu) for v in r]
            shift = sum(mu) - len(mu)   # s = x^{|mu|} [p_mu]Z, r = x^{l} prod 1/(2 sinh)
            r = ([Fr(0)] * shift + r)[:P + 1]
            tot += 1
            bad += r != s
    print(f"tau=0: [p_mu]Z = prod 1/(2 sinh(mu_i x/2))/z_mu (p_k(q^rho) = 1/[k]): {tot} partitions, failures {bad}")
