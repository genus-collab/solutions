# g4_period.py -- evidence for gap G4 item 3 (route 2 of the Gap Theorem; vanishing period). Exact over Q.
# Usage: py -I -B g4_period.py GMAX
# Curve x = A z^4 + V/z, y = z (A = 1: the XY normalisation; A = 1/4, V = v: the note's x = W_v). One-point closed formula of arXiv:2212.00320 (7.10):
#   omega^A_{g,1} = [hbar^2g] Sum_{r>=1} D_x^r [w^r] W_1,  W_1 = -exp(E) dz/w,  E = -Sum_k w^(2k+1) hbar^(2k) x^(2k)(z)/(4^k (2k+1)!),
# so omega^A_{g,1} = d Phi_g, Phi_g = Sum_r delta^(r-1)(a_r/x'), a_r = -[hbar^2g][w^(r+1)] e^E, delta = (1/x') d/dz.
# Checks: (a) S_g := Sum_j (-1)^j (2g-2+j)!/j! [T^g] phit(T)^j = -4^g B_2g/(2g), phit = Sum_k T^k/(2k+1) (exact, g <= 15);
#   (b) every term delta^(r-1)(a_r/x') is regular at z = 0 with constant term = the pure-pole value; Phi_g(0) = -(B_2g/2g) V^(1-2g);
#   (c) every term vanishes at z = oo; (d) controls.
import os, sys
from fractions import Fraction as Fr
from math import factorial, comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g4_field import bernoulli

KMAX = 4 * (int(sys.argv[1]) if len(sys.argv) > 1 else 4) + 8


def lmul(a, b):
    r = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= KMAX: r[i + j] = r.get(i + j, 0) + x * y
    return {k: v for k, v in r.items() if v != 0}


def ladd(a, b, s=1):
    r = dict(a)
    for k, v in b.items(): r[k] = r.get(k, 0) + s * v
    return {k: v for k, v in r.items() if v != 0}


def lscale(a, s):
    return {k: v * s for k, v in a.items() if v * s != 0}


def S_sum(g, coef=lambda k: Fr(1, 2 * k + 1)):
    phit = [Fr(0)] + [coef(k) for k in range(1, g + 1)]
    tot = Fr(0); pw = [Fr(1)] + [Fr(0)] * g
    for j in range(1, g + 1):
        new = [Fr(0)] * (g + 1)
        for a, x in enumerate(pw):
            for b, y in enumerate(phit):
                if a + b <= g: new[a + b] += x * y
        pw = new
        tot += Fr((-1) ** j * factorial(2 * g - 2 + j), factorial(j)) * pw[g]
    return tot


def phi_terms(g, V, at, pure=False, A=Fr(1)):
    """list of the r-terms delta^(r-1)(a_r/x') as Laurent series; at = '0' (variable z) or 'oo' (variable u = 1/z)"""
    def xder(k2):                                       # x^(2k)(z) as a Laurent series in the chosen variable
        poly = {4 - k2: A * Fr(factorial(4), factorial(4 - k2))} if (k2 <= 4 and not pure) else {}
        pole = {-(k2 + 1): Fr(factorial(k2)) * V}
        s = ladd(poly, pole)
        return s if at == '0' else {-e: c for e, c in s.items()}
    if at == '0':                                       # 1/x' = z^2/(4 z^5 - V) = -(z^2/V) Sum (4 z^5/V)^m
        ixp = {2 + 5 * m: -(4 * A) ** m / V ** (m + 1) for m in range(KMAX)} if not pure else {2: -Fr(1) / V}
        dd = lambda f: lmul(ixp, {e - 1: c * e for e, c in f.items() if e != 0})
    else:                                               # u = 1/z: 1/x' = (u^3/4) Sum (V u^5/4)^m, d/dz = -u^2 d/du
        ixp = {3 + 5 * m: V ** m / (4 * A) ** (m + 1) for m in range(KMAX)}
        dd = lambda f: lmul(ixp, {e + 1: -c * e for e, c in f.items() if e != 0})
    # E = Sum_k hbar^(2k) w^(2k+1) eps_k,  eps_k = -x^(2k)/(4^k (2k+1)!);  exp(E) at hbar^(2g): dict w-power -> Laurent series
    eps = {k: lscale(xder(2 * k), Fr(-1, 4 ** k * factorial(2 * k + 1))) for k in range(1, g + 1)}
    layer = {0: {0: {0: Fr(1)}}}                       # layer[h][wpow] = series, for h = hbar-order/2
    expo = {0: {0: {0: Fr(1)}}}
    for j in range(1, g + 1):                           # E^j / j!
        new = {}
        for h, wd in layer.items():
            for k in range(1, g + 1 - h):
                for wp, ser in wd.items():
                    t = lmul(ser, eps[k]); key = (h + k, wp + 2 * k + 1)
                    new.setdefault(key[0], {}); new[key[0]][key[1]] = ladd(new[key[0]].get(key[1], {}), t)
        layer = new
        for h, wd in layer.items():
            for wp, ser in wd.items():
                expo.setdefault(h, {}); expo[h][wp] = ladd(expo[h].get(wp, {}), lscale(ser, Fr(1, factorial(j))))
    terms = []
    for wp, ser in sorted(expo.get(g, {}).items()):
        r = wp - 1
        if r < 1: continue
        f = lmul(lscale(ser, -1), ixp)
        for _ in range(r - 1): f = dd(f)
        terms.append((r, f))
    return terms


def main():
    G = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    B = bernoulli(2 * 15 + 2); ok = True
    bad = [g for g in range(1, 16) if S_sum(g) != -Fr(4) ** g * B[2 * g] / (2 * g)]
    ctrl = [g for g in range(1, 16) if S_sum(g, lambda k: Fr(1, 2 * k + 2)) == -Fr(4) ** g * B[2 * g] / (2 * g)]
    print('(a) S_g = -4^g B_2g/(2g) for g = 1..15: failures', bad, '| control (coefficients 1/(2k+2)) agrees at', ctrl)
    ok &= bad == [] and len(ctrl) < 3
    for A, V in ((Fr(1), Fr(1)), (Fr(1), Fr(3)), (Fr(1, 4), Fr(1)), (Fr(1, 4), Fr(3))):
        for g in range(1, G + 1):
            t0 = phi_terms(g, V, '0', A=A); tp = phi_terms(g, V, '0', pure=True, A=A); ti = phi_terms(g, V, 'oo', A=A)
            neg = [r for r, f in t0 if any(e < 0 for e in f)]
            const = sum((f.get(0, 0) for r, f in t0), Fr(0)); constp = sum((f.get(0, 0) for r, f in tp), Fr(0))
            pos = [r for r, f in ti if any(e <= 0 for e in f)]
            want = -B[2 * g] / (2 * g) * V ** (1 - 2 * g)
            sg = S_sum(g) * V * (2 * V) ** (-2 * g)
            print('(b,c) A = %s V = %s g = %d: terms r = %s; negative powers at 0 in r = %s; Phi_g(0) = %s (pure pole %s, '
                  'S-formula %s, -(B_2g/2g) V^(1-2g) = %s); terms not vanishing at oo: %s'
                  % (A, V, g, [r for r, f in t0], neg, const, constp, sg, want, pos))
            ok &= neg == [] and const == want == constp == sg and pos == []
    print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')


main()
