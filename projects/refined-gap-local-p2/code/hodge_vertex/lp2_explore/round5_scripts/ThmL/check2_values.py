# Check 2: values X_g(U) from the MV side.
# (a) Hodge control of (MV): genus 1, 2 M_g(mu) against direct Hodge integrals (base tables, string, dilaton, Mumford)
# (b) degree bound deg_s Q <= N
# (c) string: X_g(U+1) = J_U X_g(U); X_g(1^v) = 0, g>=2; X_1(1^v)=(v-1)! X_1(1), X_1(1) = -(1+1/sigma)/24
# (d) Lemma A range: X_h(T)=0 for non-empty T parts>=2, D<=2h-3; control: nonzero at D=2h-2 and beyond; empty profile nonzero
# (e) closed form (6.1) for one part (R2), an independent value control
import sys
from fractions import Fraction as Fr
from math import factorial, comb
from tl_lib import *
from mvx import *

tau = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 7)
gmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
maxdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 14
mv = MV(tau, L=2 * maxdeg + 2 * gmax + 2)
sigma = mv.sigma
print("tau", tau, "sigma", sigma, "gmax", gmax, "max |mu|", maxdeg)

# (a) direct Hodge integrals, genus 1 and 2
from functools import lru_cache
BASE = {
 (2, (), ()): None,
}
PSI2 = {(4,): Fr(1,1152), (3,2): Fr(29,5760), (2,2,2): Fr(7,240)}
LAM1 = {(3,): Fr(1,480), (2,2): Fr(5,576)}
LAM2 = {(2,): Fr(7,5760)}
L1L2 = {(): Fr(1,5760)}

@lru_cache(maxsize=None)
def corr(g, cls, ds):
    """<cls prod tau_d>_g, cls in 'psi','l1','l2','l1l2' (g=2), 'psi','l1' (g=1)."""
    ds = tuple(sorted(ds, reverse=True))
    n = len(ds)
    deg = {'psi':0,'l1':1,'l2':2,'l1l2':3}[cls]
    if sum(ds) + deg != 3*g - 3 + n: return Fr(0)
    if ds and ds[-1] == 0:
        rest = ds[:-1]
        if not rest: return Fr(1, 24) if (g == 1 and cls == 'l1') else Fr(0)
        tot = Fr(0)
        for i in range(len(rest)):
            if rest[i] >= 1:
                nd = list(rest); nd[i] -= 1
                tot += corr(g, cls, tuple(nd))
        return tot
    if ds and ds[-1] == 1:
        rest = ds[:-1]
        if g == 1 and not rest:
            return Fr(1, 24) if cls == 'psi' else Fr(0)
        return (2*g - 2 + n - 1) * corr(g, cls, rest)
    if g == 2:
        tab = {'psi':PSI2,'l1':LAM1,'l2':LAM2,'l1l2':L1L2}[cls]
        return tab.get(ds, Fr(0))
    if g == 1:
        if cls == 'l1' and ds == (): return Fr(1, 24)
        return Fr(0)
    return Fr(0)

def hodge_M(g, mu):
    l = len(mu)
    # Omega = Lambda(1)Lambda(tau)Lambda(-1-tau)/(-sigma)^g as dict on classes
    ws = [Fr(1), tau, -1 - tau]
    if g == 1:
        Om = {'psi': Fr(1), 'l1': Fr(0)}
        c0 = ws[0]*ws[1]*ws[2]
        c1 = -(ws[0]*ws[1] + ws[0]*ws[2] + ws[1]*ws[2])
        Om = {'psi': c0, 'l1': c1}
    else:
        # polynomial in l1,l2 with l1^2 = 2 l2, l2^2 = 0, l1^3 = 2 l1l2, l1^2 l2 = 0
        def mul(a, b):
            r = {}
            tbl = {('psi','psi'):('psi',1),('psi','l1'):('l1',1),('psi','l2'):('l2',1),('psi','l1l2'):('l1l2',1),
                   ('l1','l1'):('l2',2),('l1','l2'):('l1l2',1),('l1','l1l2'):('l2l2',0),('l2','l2'):('l2l2',0),('l2','l1l2'):('l2l2',0),('l1l2','l1l2'):('l2l2',0)}
            for x, cx in a.items():
                for y, cy in b.items():
                    k = (x, y) if (x, y) in tbl else (y, x)
                    z, m = tbl[k]
                    if m: r[z] = r.get(z, 0) + cx * cy * m
            return r
        Om = {'psi': Fr(1)}
        for w in ws:
            Om = mul(Om, {'psi': w*w, 'l1': -w, 'l2': Fr(1)})
    den = (-sigma)**g
    tot = Fr(0)
    dim = 3*g - 3 + l
    def rec(i, rem, cur, ds):
        nonlocal tot
        if i == l:
            for cls, co in Om.items():
                tot += cur * co * corr(g, cls, tuple(ds))
            return
        for d in range(rem + 1):
            rec(i + 1, rem - d, cur * mu[i]**d, ds + [d])
    rec(0, dim, Fr(1), [])
    return tot / den

cnt = bad = 0
for g in (1, 2):
    for mu in [(1,),(2,),(3,),(5,),(7,),(1,1),(2,1),(2,2),(3,2),(4,3),(5,1),(1,1,1),(2,1,1),(3,2,1),(2,2,2),(3,3,2),(1,1,1,1),(2,1,1,1)]:
        if sum(mu) + 2*g - 2 + len(mu) > mv.L: continue
        a = mv.M(g, mu); b = hodge_M(g, mu)
        cnt += 1
        if a != b:
            bad += 1
            print("  MISMATCH", g, mu, a, b)
print(f"(a) MV side against direct Hodge integrals, g=1,2: {cnt} comparisons, {bad} mismatches")
# genus 0 shape
cnt = bad = 0
for mu in [(1,),(2,),(3,1),(2,2),(3,2,1),(2,2,2),(4,1,1,1)]:
    a = mv.M(0, mu); b = Fr(sum(mu))**(len(mu)-3)
    cnt += 1; bad += (a != b)
print(f"    genus 0: M_0(mu) = |mu|^(l-3): {cnt} comparisons, {bad} mismatches")
# x-power shape control: no x^(2g-2+l-1) term
cnt = bad = 0
for mu in [(3,),(2,1),(3,2),(2,2,1)]:
    n = sum(mu); l = len(mu)
    for i in range(0, 2*gmax+l-2+n+1):
        e = i - n  # power of x
        v = mv.f(mu)[i]
        ok = (v == 0) if (e < l - 2 or (e - l) % 2 != 0) else True
        cnt += 1; bad += (not ok)
print(f"    shape of log Z (only x^(2g-2+l)): {cnt} coefficients, {bad} violations")

# (b)-(c) on a list
cnt = bad = 0
for g in range(1, gmax + 1):
    for U in [(1,),(2,),(3,),(2,2),(3,2),(4,),(2,1),(3,1),(1,1)]:
        r = len(U); N = 3*g - 3 + r
        if sum(U) + N + 1 > maxdeg: continue
        mv.X(g, U, check_top=True)   # asserts the coefficient of s^(N+1) vanishes
        cnt += 1
print(f"(b) deg_s Q_g(s;U) <= N: {cnt} cases, 0 failures (assert)")
cnt = bad = 0
for g in range(1, gmax + 1):
    for U in [(2,),(3,),(2,2),(3,2),(4,),(5,)]:
        r = len(U) + 1; N = 3*g - 3 + r
        if sum(U) + 1 + N > maxdeg: continue
        a = mv.X(g, U + (1,)); b = sum(U) * mv.X(g, U)
        cnt += 1; bad += (a != b)
    for v in (1, 2, 3):
        N = 3*g - 3 + v
        if v + N > maxdeg: continue
        a = mv.X(g, (1,)*v)
        b = (factorial(v-1) * (-(1 + 1/sigma)/24)) if g == 1 else Fr(0)
        cnt += 1; bad += (a != b)
print(f"(c) string on X and the values at U=1^v: {cnt} checks, {bad} failures")

# (d)
from itertools import combinations_with_replacement
def multisets_level(D, minpart=2):
    # all multisets of parts >=2 with sum (part-1) = D
    for lam in partitions(D):
        yield tuple(p + 1 for p in lam)
inr = nz = skipped = 0
for h in range(2, gmax + 1):
    for D in range(1, 2*h + 2):
        for T in multisets_level(D):
            r = len(T); N = 3*h - 3 + r
            if sum(T) + N > maxdeg: skipped += 1; continue
            v = mv.X(h, T)
            if D <= 2*h - 3:
                inr += 1
                if v != 0: print("  VIOLATION", h, T, v)
            else:
                nz += (v != 0)
                print(f"    outside range: h={h} T={T} D={D}: {'nonzero' if v != 0 else 'zero'}")
print(f"(d) Lemma A range values vanish: {inr} points, skipped (size) {skipped}; outside-range nonzero count {nz}")
# empty profile
for h in range(2, gmax + 1):
    N = 3*h - 3
    if N > maxdeg: continue
    # empty U: Q_h(-1;empty) with G from a>=0, M_h(empty) unavailable from log Z; use the formula Xi_h(0-part) route: omit
    pass

# (e) closed form (6.1): sigma^h X_h(j) = prod_{l=1}^{2h-2} (j-l) [lam^{2h}] S(e1 lam)^{1+j e2} S(e2 lam)^{1+j e1}, e1 = -tau, e2 = 1+tau
def closed(h, j):
    e1 = -tau; e2 = 1 + tau
    L = 2*h
    B = bernoulli_list(2*L + 2)
    a = [Fr(0)] * (L + 1)
    for k in range(1, L//2 + 1):
        co = -B[2*k] / (2*k*factorial(2*k))
        a[2*k] = co * ((1 + j*e2) * e1**(2*k) + (1 + j*e1) * e2**(2*k))
    ex = ser_exp(a, L)
    p = Fr(1)
    for l in range(1, 2*h - 1): p *= (j - l)
    return p * ex[2*h]
cnt = bad = 0
for h in range(1, gmax + 1):
    for j in range(1, 12):
        N = 3*h - 3 + 1
        if j + N > maxdeg: continue
        a = sigma**h * mv.X(h, (j,)); b = closed(h, j)
        cnt += 1
        if a != b:
            bad += 1; print("   closed-form mismatch", h, j, a, b)
print(f"(e) one-part closed form (6.1) against MV side: {cnt} comparisons, {bad} mismatches")
