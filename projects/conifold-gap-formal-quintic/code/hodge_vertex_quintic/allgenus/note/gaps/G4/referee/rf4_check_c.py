# rf4_check_c.py -- referee check C (item 3, route 2: prop:mod-period / lem:mod-vanper). Exact over Q. Usage: py -I -B rf4_check_c.py H
# C1 S_h := Sum_j (-1)^j (2h-2+j)!/j! [T^h] phit^j = -4^h B_2h/(2h), h <= 20; control phit -> Sum T^k/(2k+2).
# C2 From ABDKS (7.10) as transcribed in GAPS_G4_4 s.4 Step 1: Phi_h = Sum_r delta^(r-1)(a_r/x'), x = A p^4 + V/p, for h <= H,
#    (A, V) in {(1/4, 1), (1/4, 3), (1, 2)}: Laurent expansion at p = 0 has no negative power and constant term
#    V (2V)^(-2h) S_h = -(B_2h/2h) V^(1-2h); expansion at p = oo has top exponent <= 5 - 10h.
# C3 omega_(1,1) = -d Phi_1 equals the closed form -(5/8) p^4 (p^10 + 18 v p^5 + 6 v^2)/(p^5 - v)^4 dp (A = 1/4), to p^80 at 0.
import os, sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sys
from fractions import Fraction as F
from math import factorial, comb
from rf4_field import bernoulli

H = int(sys.argv[1]) if len(sys.argv) > 1 else 5
ok = True
B = bernoulli(42)


def S_h(h, ctrl=False):
    phit = [F(0)] + [F(1, (2 * k + 2) if ctrl else (2 * k + 1)) for k in range(1, h + 1)]
    tot = F(0); pw = [F(1)] + [F(0)] * h
    for j in range(1, h + 1):
        pw = [sum(pw[a] * phit[n - a] for a in range(n + 1)) for n in range(h + 1)]
        tot += F((-1) ** j * factorial(2 * h - 2 + j), factorial(j)) * pw[h]
    return tot


badS = sum(1 for h in range(1, 21) if S_h(h) != -F(4) ** h * B[2 * h] / (2 * h))
ctrlS = sum(1 for h in range(1, 21) if S_h(h, True) == -F(4) ** h * B[2 * h] / (2 * h))
print('C1 S_h = -4^h B_2h/(2h) for h = 1..20: %d failures; control agrees for %d of 20 (must be 0)' % (badS, ctrlS))
ok &= badS == 0 and ctrlS == 0


class L(dict):
    """Laurent series {exponent: coeff} truncated: at 0 keep e <= TOP, at oo keep e >= BOT"""
    pass


def lmul(a, b, keep):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            if keep(e):
                r[e] = r.get(e, 0) + c1 * c2
    return {e: c for e, c in r.items() if c != 0}


def ladd(a, b, s=1):
    r = dict(a)
    for e, c in b.items():
        r[e] = r.get(e, 0) + s * c
    return {e: c for e, c in r.items() if c != 0}


def lder(a):
    return {e - 1: e * c for e, c in a.items() if e != 0}


def eps(k, A, V):
    """eps_k = -x^(2k)/(4^k (2k+1)!) as an exact Laurent polynomial"""
    d = {}
    if 2 * k <= 4:
        d[4 - 2 * k] = A * F(factorial(4), factorial(4 - 2 * k))
    d[-2 * k - 1] = d.get(-2 * k - 1, 0) + V * factorial(2 * k)
    return {e: -c / (F(4) ** k * factorial(2 * k + 1)) for e, c in d.items()}


def compositions(h, j):
    if j == 1:
        yield (h,); return
    for k in range(1, h - j + 2):
        for rest in compositions(h - k, j - 1):
            yield (k,) + rest


def Phi_terms(h, A, V, at0, NT):
    if at0:
        keep = lambda e: e <= NT
        ixp = {2 + 5 * n: -1 / V * (4 * A / V) ** n for n in range((NT + 10 * h + 10) // 5 + 2)}   # p^2/(4A p^5 - V)
    else:
        keep = lambda e: e >= -NT
        ixp = {-3 - 5 * n: (1 / (4 * A)) * (V / (4 * A)) ** n for n in range((NT + 10 * h + 10) // 5 + 2)}
    terms = []
    for j in range(1, h + 1):
        r = 2 * h - 1 + j
        a = {}
        for ks in compositions(h, j):
            prod = {0: F(1)}
            for k in ks:
                prod = lmul(prod, eps(k, A, V), keep)
            a = ladd(a, prod)
        a = {e: -c / factorial(j) for e, c in a.items()}
        t = lmul(a, ixp, keep)
        for _ in range(r - 1):
            t = lmul(ixp, lder(t), keep)
        terms.append((r, t))
    return terms


for (A, V) in ((F(1, 4), F(1)), (F(1, 4), F(3)), (F(1), F(2))):
    for h in range(1, H + 1):
        t0 = Phi_terms(h, A, V, True, 20)
        neg = sum(1 for r, t in t0 for e in t if e < 0)
        const = sum(t.get(0, F(0)) for r, t in t0)
        target = -B[2 * h] / (2 * h) * V ** (1 - 2 * h)
        pure = V * (2 * V) ** (-2 * h) * S_h(h)
        tinf = Phi_terms(h, A, V, False, 20 * h + 20)
        top = max((e for r, t in tinf for e in t), default=None)
        okh = neg == 0 and const == target == pure and (top is None or top <= 5 - 10 * h)
        print('C2 A=%s V=%s h=%d: negative powers at 0: %d; Phi_h(0) = %s (target %s); top exponent at oo: %s (bound %d)%s'
              % (A, V, h, neg, const, target, top, 5 - 10 * h, '' if okh else '  <-- FAIL'))
        ok &= okh

# C3: omega_11 = -d Phi_1 vs the closed form (A = 1/4), expansions at p = 0 up to p^80, v = 1 and v = 3
for V in (F(1), F(3)):
    NT = 85
    (r, t), = Phi_terms(1, F(1, 4), V, True, NT)
    om = {e: -c for e, c in lder(t).items() if e <= 80}
    # closed form at 0: -(5/8) p^4 (p^10 + 18 V p^5 + 6 V^2) (p^5 - V)^-4 ; (p^5 - V)^-4 = V^-4 Sum C(n+3,3) (p^5/V)^n
    cf = {}
    for n in range(0, 18):
        c = F(comb(n + 3, 3)) / V ** (4 + n)
        for e0, c0 in ((14, F(1)), (9, 18 * V), (4, 6 * V * V)):
            e = e0 + 5 * n
            if e <= 80:
                cf[e] = cf.get(e, 0) - F(5, 8) * c0 * c
    diff = {e: c for e, c in ladd(om, cf, -1).items() if e <= 80 and c != 0}
    print('C3 V=%s: omega_11 = -d Phi_1 vs closed form, differences up to p^80: %d' % (V, len(diff)))
    ok &= len(diff) == 0
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
