# Xi_h(R) directly at (E,-e1,-e2), Theorem VI, and t_cf in floating point. Evidence only.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import math, time
from fractions import Fraction as Fr
from itertools import product
from gs import *


def xi_val(h, R, e1, e2, V):
    E = Fr(e1) + Fr(e2)
    r = len(R)
    dim_base = 3 * h - 3 + r
    tot = 0
    Efp = fp(E)
    Ei = inv(Efp)
    # sum over slot degrees
    for a in product(*[range(dim_base + 1) for _ in R]):
        if sum(a) > dim_base:
            continue
        w0 = 1
        for ai, j in zip(a, R):
            w0 = w0 * pow(j, ai, P) % P * pow(Ei, ai + 1, P) % P
        budget = dim_base - sum(a)
        for D in partitions_ge2(budget):
            val = V.VI(h, tuple(a) + D)
            if val == 0:
                continue
            w = w0 * val % P
            for d in set(D):
                w = w * inv(factorial(D.count(d))) % P
            for d in D:
                w = w * (P - pow(Ei, d - 1, P)) % P  # T = -sum psi^d E^(1-d)
            tot = (tot + w) % P
    return tot


def parts_with_D(D, mx=None):
    # multisets of integers >= 2 with sum (j-1) = D
    def rec(rem, mx):
        if rem == 0:
            yield ()
            return
        for j1 in range(min(rem, mx), 0, -1):
            for t in rec(rem - j1, j1):
                yield (j1 + 1,) + t
    return list(rec(D, D))


for (e1, e2) in [(Fr(3), Fr(7)), (Fr(1), Fr(1)), (Fr(5, 2), Fr(-3, 7))]:
    V = Vert(e1, e2)
    for h in range(2, 6):
        t0 = time.time()
        nz = 0
        n = 0
        for D in range(1, 2 * h - 2):
            for R in parts_with_D(D):
                if len(R) > (3 if h >= 5 else 4):
                    continue
                n += 1
                if xi_val(h, R, e1, e2, V) != 0:
                    nz += 1
                    print('  NONZERO in range', h, R)
        # outside range
        nout = 0
        nout0 = 0
        for R in parts_with_D(2 * h - 2):
            if len(R) > 2:
                continue
            nout += 1
            if xi_val(h, R, e1, e2, V) == 0:
                nout0 += 1
        # Theorem VI
        Z = xi_val(h, (), e1, e2, V)
        E = e1 + e2
        vi = (fp(-e1 * e2 * E ** 3) ** (h - 1) if False else pow(fp(-e1 * e2 * E ** 3), h - 1, P)) * Z % P == fp(beta(2 * h - 2, e1, e2))
        print('eps', e1, e2, 'h', h, 'in-range R tested', n, 'nonzero', nz, '| D=2h-2 controls', nout, 'zero among them', nout0, '| Theorem VI', vi, '| time', round(time.time() - t0, 1), flush=True)

# t_cf in floats at y = -1/54, -1/81
pi = math.pi
for y in (-1 / 54, -1 / 81, -1 / 40):
    u = -27 * y
    xi = 1 - u
    # series of A, g, g2 in (-y)
    A = 0.0
    gg = 0.0
    g2 = 0.0
    a = 1.0
    Hd = 0.0
    H3 = 0.0
    py = 1.0
    for d in range(1, 4000):
        a *= (3 * d) * (3 * d - 1) * (3 * d - 2) / d ** 3
        py *= (-y)
        Hd += 1 / d
        H3 += 1 / (3 * d) + 1 / (3 * d - 1) + 1 / (3 * d - 2)
        c = a * py
        if c < 1e-30:
            break
        # A = 1 + theta g => g_d = a_d/d (-y)^d ; h = g + theta g2 => (g2)_d = (h_d - g_d)/d
        gd = c / d
        hd = 3 * c * (H3 - Hd)
        gg += gd
        g2 += (hd - gd) / d
    l = math.log(-y)
    tcf = math.sqrt(3) / (2 * pi) * (l * l / 2 + gg * l + g2 - pi * pi / 2)
    # t_c by series in xi
    w = 1.0
    ws = []
    for k in range(0, 3000):
        ws.append(w)
        w *= (1 / 3 + k) * (2 / 3 + k) / (k + 1) ** 2
    # dt/dxi = wc/(1-xi): coefficients partial sums
    s = 0.0
    tc = 0.0
    xp = xi
    for k in range(0, 3000):
        s += ws[k]
        tc += s * xp / (k + 1)
        xp *= xi
    print('y', y, 't_cf (BS, branch +i pi)', tcf, 't_c', tc, 'diff', tcf - tc)
