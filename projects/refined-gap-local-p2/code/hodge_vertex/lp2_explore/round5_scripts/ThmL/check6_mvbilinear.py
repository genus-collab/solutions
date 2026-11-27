# Check 6: P_k tau_0.tau_1 = 0 directly on Z_MV, tau_1 = Z(e^{tau x} o p), as Laurent series in x (no logarithm), exact.
# Controls: tau_1 replaced by Z(e^{-tau x} o p), swapped order.
import sys, itertools
from fractions import Fraction as Fr
from math import factorial
from tl_lib import *
from mvx import MV

tau = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 7)
dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 7   # total weight
L = int(sys.argv[3]) if len(sys.argv) > 3 else 24     # relative x-order
kmax = int(sys.argv[4]) if len(sys.argv) > 4 else 5
mv = MV(tau, L=L)

def shifted(mu, sgn):
    n = sum(mu)
    a = [Fr(0)] * (L + 1)
    if L >= 1: a[1] = sgn * n * tau
    return ser_mul(mv.Z(mu), ser_exp(a, L), L)

pars = [mu for n in range(dmax + 1) for mu in partitions(n)]
def Zc(mu): return mv.Z(mu) if mu else [Fr(1)] + [Fr(0)] * L
T0 = {mu: Zc(mu) for mu in pars}
T1 = {mu: (shifted(mu, 1) if mu else Zc(mu)) for mu in pars}
Tm = {mu: (shifted(mu, -1) if mu else Zc(mu)) for mu in pars}

def deriv(mu, idx):
    """d/dt_i ... on p_mu: returns (coefficient, nu) or None."""
    l = list(mu); co = 1
    for k in idx:
        e = l.count(k)
        if e == 0: return None
        co *= e * k
        l.remove(k)
    return co, tuple(l)

def hir(f, g, alpha, target):
    """coefficient of p_target in D^alpha f.g, x^{|target|+|alpha|} normalised."""
    n = len(alpha)
    res = [Fr(0)] * (L + 1)
    w = sum(target) + sum(alpha)
    for bits in itertools.product([0, 1], repeat=n):
        S = [alpha[i] for i in range(n) if bits[i]]
        C = [alpha[i] for i in range(n) if not bits[i]]
        sgn = (-1)**len(C)
        for sub in set(itertools.chain.from_iterable(itertools.combinations(target, r) for r in range(len(target) + 1))):
            lf = list(sub); lg = list(target)
            for s in lf: lg.remove(s)
            # coefficient of p_lf in d^S f: f has monomial lf+S
            mf = tuple(sorted(lf + S, reverse=True)); mg = tuple(sorted(lg + C, reverse=True))
            df = deriv(mf, S); dg = deriv(mg, C)
            if df is None or dg is None: continue
            pr = ser_mul(f[mf], g[mg], L)
            c = sgn * df[0] * dg[0]
            for i in range(L + 1): res[i] += c * pr[i]
    return res

tot = bad = c1 = c2 = 0
for k in range(2, kmax + 1):
    P = member_Pk(k)
    for n in range(0, dmax - k + 1):
        for target in partitions(n):
            # NB multiplicities: the coefficient of p_target of a product is the sum over distinct sub-multisets
            r = [Fr(0)] * (L + 1); r1 = [Fr(0)] * (L + 1); r2 = [Fr(0)] * (L + 1)
            for alpha, c in P.items():
                h = hir(T0, T1, list(alpha), target)
                h1 = hir(T0, Tm, list(alpha), target)
                h2 = hir(T1, T0, list(alpha), target)
                for i in range(L + 1):
                    r[i] += c * h[i]; r1[i] += c * h1[i]; r2[i] += c * h2[i]
            tot += 1
            if any(v != 0 for v in r):
                bad += 1
                print("FAIL", k, target)
            c1 += any(v != 0 for v in r1); c2 += any(v != 0 for v in r2)
print(f"tau={tau} total weight<={dmax}, relative x-order {L}: {tot} coefficients, {bad} failures; controls nonzero: tau_-1 {c1}, swapped {c2}")
