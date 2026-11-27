# MV Schur side -> M_g(mu) -> X_g(U), exact over Q at rational tau (uses MVmatch form (MV) only; no Hodge data).
from fractions import Fraction as Fr
from math import factorial, comb
from functools import lru_cache
from tl_lib import *

class MV:
    def __init__(self, tau, L=44):
        self.tau = Fr(tau)
        self.sigma = self.tau * (1 + self.tau)
        self.L = L             # relative x-order kept in x^{|mu|} Z_mu
        self.B = bernoulli_list(2 * L + 2)
        self.wcache = {}
        self.Zc = {}
        self.fc = {}
        self.Mc = {}
        self.Xc = {}

    def w(self, nu):
        """x^{|nu|} * e^{tau x f2} s_nu(q^rho), truncated at x^L. Uses s_nu(q^rho)=q^{f2/2}/prod [h]."""
        if nu in self.wcache: return self.wcache[nu]
        L = self.L
        a = [Fr(0)] * (L + 1)
        c = f2(nu)
        if L >= 1: a[1] = (self.tau + Fr(1, 2)) * c
        hs = hooks(nu) if nu else []
        for k in range(1, L // 2 + 1):
            # log S(y) = - sum B_2k y^2k / (2k (2k)!)
            co = -self.B[2 * k] / (2 * k * factorial(2 * k))
            a[2 * k] += co * sum(h**(2 * k) for h in hs)
        e = ser_exp(a, L)
        pre = Fr(1)
        for h in hs: pre /= h
        r = [pre * v for v in e]
        self.wcache[nu] = r
        return r

    def Z(self, mu):
        """x^{|mu|} [p_mu] Z_MV as truncated series."""
        if mu in self.Zc: return self.Zc[mu]
        n = sum(mu)
        r = [Fr(0)] * (self.L + 1)
        z = zmu(mu)
        for nu in partitions(n):
            c = chi(nu, mu)
            if c == 0: continue
            w = self.w(nu)
            fc = Fr(c, z)
            for i in range(self.L + 1):
                if w[i] != 0: r[i] += fc * w[i]
        self.Zc[mu] = r
        return r

    def f(self, mu):
        """x^{|mu|} [p_mu] log Z, by E Z = Z E F (recursion over sub-multisets)."""
        if mu in self.fc: return self.fc[mu]
        n = sum(mu)
        r = list(self.Z(mu))
        # n f_mu = n Z_mu - sum_{mu'+mu''=mu, mu' nonempty, mu'' nonempty} |mu''| f_mu'' Z_mu'
        for sub in sub_multisets(mu):
            if not sub or sub == mu: continue
            comp = multiset_minus(mu, sub)
            s2 = sum(sub)
            pr = ser_mul(self.f(sub), self.Z(comp), self.L)
            for i in range(self.L + 1):
                r[i] -= Fr(s2, n) * pr[i]
        self.fc[mu] = r
        return r

    def M(self, g, mu):
        """M_g(mu) from (MV): [x^{2g-2+l}][p_mu] log Z = sigma^{g+l-1} Pi_mu M_g(mu)/|Aut mu|."""
        key = (g, mu)
        if key in self.Mc: return self.Mc[key]
        l = len(mu); n = sum(mu)
        idx = 2 * g - 2 + l + n
        assert idx <= self.L, "increase L"
        val = self.f(mu)[idx]
        den = self.sigma**(g + l - 1)
        for k in mu: den *= Pi(k, self.tau)
        self.Mc[key] = val * aut(mu) / den
        return self.Mc[key]

    def Qcoef(self, g, U, top):
        """[s^m] Q_g(s;U) for m=0..top from the a-resummation, parts of U >= 1."""
        r = len(U); J = sum(U); e0 = 2 * g - 2 + r
        out = [Fr(0)] * (top + 1)
        for a in range(0, top + 1):
            mu = tuple(sorted(U + (1,) * a, reverse=True))
            if not mu: continue
            Ma = self.M(g, mu)
            # s^a (1+s)^{-(e0+a)} exp(-(J+a) s/(1+s)) / a!
            ls = [Fr(0)] * (top + 1)
            # log series: -(e0+a) log(1+s) - (J+a) s/(1+s)
            for k in range(1, top + 1):
                ls[k] += -(e0 + a) * Fr((-1)**(k - 1), k) - (J + a) * (-1)**(k - 1)
            ex = ser_exp(ls, top)
            for m in range(a, top + 1):
                out[m] += Ma / factorial(a) * ex[m - a]
        return out

    def X(self, g, U, check_top=False):
        """X_g(U): g=0 closed form; g>=1 Q_g(-1;U). U tuple, parts>=1, nonempty."""
        U = tuple(sorted(U, reverse=True))
        key = (g, U)
        if key in self.Xc and not check_top: return self.Xc[key]
        r = len(U); J = sum(U)
        if g == 0:
            v = x0(U)
        else:
            N = 3 * g - 3 + r
            if check_top:
                q = self.Qcoef(g, U, N + 1)
                assert q[N + 1] == 0, ("degree bound fails", g, U)
            else:
                q = self.Qcoef(g, U, N)
            v = sum((-1)**m * q[m] for m in range(N + 1))
        self.Xc[key] = v
        return v

def x0(U):
    r = len(U); J = sum(U)
    # Gamma(J)/Gamma(J-r+3)
    v = Fr(1)
    if J - 1 >= J - r + 2:
        for t in range(J - r + 3, J): v *= t
        return v
    for t in range(J, J - r + 3): v /= t
    return v

@lru_cache(maxsize=None)
def sub_multisets(mu):
    vals = sorted(set(mu), reverse=True)
    res = [()]
    for v in vals:
        c = mu.count(v)
        res = [s + (v,) * i for s in res for i in range(c + 1)]
    return [tuple(sorted(s, reverse=True)) for s in res]

def multiset_minus(mu, sub):
    l = list(mu)
    for s in sub: l.remove(s)
    return tuple(l)

def falling(a, k):
    r = 1
    for i in range(k): r *= (a - i)
    return r

def D_of(U): return sum(U) - len(U)
def delta(g, U): return D_of(U) - 2 * g + 2
