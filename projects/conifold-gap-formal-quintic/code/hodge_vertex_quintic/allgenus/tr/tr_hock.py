# -*- coding: utf-8 -*-
"""tr_hock.py -- the free energies F_g of x = P(z) + v/z, y = z from the x-y swap formula of Hock (arXiv:2502.18115, Cor. 1.2,
eq. (1.3), the case y = z without logarithmic singularities of x), evaluated exactly (TR_MODEL.md s.4):

    F_g = 1/(2-2g) Sum_{b in {0, inf}} Res_{z=b} I_g(z) dz,
    I_g = Sum_{j=1}^{g} (1/j!) D^(2g-2+j)(z) Sum_{k_1+..+k_j = g, k_i >= 1} prod_i x^(2k_i)(z) / (4^k_i (2k_i+1)!),   D = (1/x'(z)) d/dz,

([hbar^(2g-1)][u^m] of the exponential (S(u hbar d_z) - 1) x u, S(t) = sinh(t/2)/(t/2), m = 2g-1+j).  The points b are the singular
points of x, y with non-zero residue of y dx or x dy: z = 0 (simple pole of x) and z = infinity.
Closed form of the residue at the simple pole (PROVED in TR_MODEL.md s.4):  Res_{z=0} I_g dz = v^(2-2g) A_g,
A_g = Sum_j (-1)^j (2g-2+j)!/j! [t^g] phi(t)^j,  phi = Sum_{k>=1} t^k/(4^k (2k+1)),  A_g = -B_2g/(2g)  (Lagrange inversion).
Run: py -I -B tr_hock.py  (exact Fractions; checks for several P and g <= 8; writes tr_hock.log)
"""
import os, sys, random, time
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_compare import bernoulli


class L(object):
    """Laurent series sum_k c[k] q^(v+k), truncated: coefficients known up to absolute degree top"""
    def __init__(self, v, c, top): self.v = v; self.c = list(c); self.top = top
    def coef(self, n):
        k = n - self.v
        return self.c[k] if 0 <= k < len(self.c) and n <= self.top else Fr(0)
    def __mul__(self, o):
        v = self.v + o.v; top = min(self.top + o.v, o.top + self.v); n = top - v + 1
        out = [Fr(0)] * max(n, 0)
        for i, a in enumerate(self.c):
            if a == 0 or i >= n: continue
            for j, b in enumerate(o.c):
                if i + j >= n: break
                if b: out[i + j] += a * b
        return L(v, out, top)
    def deriv(self):
        return L(self.v - 1, [c * (self.v + k) for k, c in enumerate(self.c)], self.top - 1)
    def scale(self, s): return L(self.v, [c * s for c in self.c], self.top)
    def add(self, o):
        v = min(self.v, o.v); top = min(self.top, o.top)
        return L(v, [self.coef(n) + o.coef(n) for n in range(v, top + 1)], top)
    def inv(self):
        """1/f, f with non-zero leading coefficient c[0]"""
        n = self.top - self.v; a0 = self.c[0]; out = [Fr(1) / a0] + [Fr(0)] * n
        for k in range(1, n + 1):
            s = sum(self.c[i] * out[k - i] for i in range(1, min(k, len(self.c) - 1) + 1))
            out[k] = -s / a0
        return L(-self.v, out, -self.v + n)


def integrand_residues(Pc, v, g, prec=12):
    """(Res_{z=0} I_g dz, Res_{z=inf} I_g dz) for x = Sum_k Pc[k] z^k + v/z (Pc[0] irrelevant)"""
    r = max(k for k, c in enumerate(Pc) if c)
    out = []
    for where in ('0', 'inf'):
        T = prec + 4 * g + 6 * r + 10
        if where == '0':                                            # q = z
            x = L(-1, [v] + [Fr(0)] + list(Pc[1:]), T)
            dz = lambda f: f.deriv()
            zf = L(1, [Fr(1)], T)
        else:                                                       # q = w = 1/z, d/dz = -w^2 d/dw
            c = [Fr(0)] * (r + 2)
            for k, a in enumerate(Pc):
                if k >= 1: c[r - k] = a                                 # w^-k at index r - k (valuation -r)
            c[r + 1] = v                                            # v w
            x = L(-r, c, T)
            m2 = L(2, [Fr(-1)], T + 10)
            dz = lambda f: m2 * f.deriv()
            zf = L(-1, [Fr(1)], T)
        dx = dz(x); idx = dx.inv()
        D = [zf]
        for n in range(1, 3 * g):
            D.append(idx * dz(D[-1]))
        X = {0: x}
        for k in range(1, 2 * g + 1): X[k] = dz(X[k - 1])
        I = None
        for j in range(1, g + 1):
            for ks in _compositions(g, j):
                term = D[2 * g - 2 + j].scale(Fr(1, factorial(j)))
                for k in ks: term = term * X[2 * k].scale(Fr(1, 4 ** k * factorial(2 * k + 1)))
                I = term if I is None else I.add(term)
        if where == '0': out.append(I.coef(-1))                     # Res_{z=0} f dz = [z^-1] f
        else: out.append(-I.coef(1))                                # Res_{z=inf} f dz = Res_{w=0} f (-dw/w^2) = -[w^1] f
    return out


def _compositions(g, j):
    if j == 1:
        yield (g,); return
    for a in range(1, g - j + 2):
        for rest in _compositions(g - a, j - 1): yield (a,) + rest


def A_closed(g):
    base = [Fr(0)] + [Fr(1, 4 ** k * (2 * k + 1)) for k in range(1, g + 1)]
    tot = Fr(0); pw = [Fr(1)] + [Fr(0)] * g
    for j in range(1, g + 1):
        new = [Fr(0)] * (g + 1)
        for a, x in enumerate(pw):
            if x:
                for b, y in enumerate(base):
                    if a + b <= g: new[a + b] += x * y
        pw = new; tot += Fr((-1) ** j * factorial(2 * g - 2 + j), factorial(j)) * pw[g]
    return tot


if __name__ == "__main__":
    f = open(os.path.join(HERE, "tr_hock.log"), "w"); t0 = time.time(); fails = 0
    def log(s): print(s, flush=True); f.write(s + "\n")
    log("tr_hock.py (Hock arXiv:2502.18115 Cor. 1.2 for x = P(z) + v/z, y = z), exact over Q")
    for g in range(2, 13):
        a = A_closed(g); ok = a == -bernoulli(2 * g) / (2 * g); fails += not ok
        log("  [%s] A_%d = %s = -B_%d/(%d)  ->  F_%d = A_g/(2-2g) = %s" % ("PASS" if ok else "FAIL", g, a, 2 * g, 2 * g, g, a / (2 - 2 * g)))
    rng = random.Random(7)
    cases = [("W: P = z^4/4, v = 1", [0, 0, 0, 0, Fr(1, 4)], Fr(1)), ("lambda: P = z^4, v = 4", [0, 0, 0, 0, 1], Fr(4)),
             ("HZ: P = z, v = 1", [0, 1], Fr(1)), ("P = z^2/2, v = 3", [0, 0, Fr(1, 2)], Fr(3))]
    for _ in range(2):
        cases.append(("random quartic", [0] + [Fr(rng.randrange(-9, 10), rng.randrange(1, 7)) for _ in range(3)] + [Fr(1, 4)], Fr(rng.randrange(1, 9), rng.randrange(1, 5))))
    for name, Pc, v in cases:
        for g in range(2, 7):
            r0, rinf = integrand_residues(Pc, v, g)
            Fg = (r0 + rinf) / (2 - 2 * g); pred = bernoulli(2 * g) / (2 * g * (2 * g - 2)) * v ** (2 - 2 * g)
            ok = rinf == 0 and r0 == A_closed(g) * v ** (2 - 2 * g) and Fg == pred; fails += not ok
            log("  [%s] %s, P coeffs %s, g = %d: Res_0 = %s, Res_inf = %s, F_g = %s (B_2g/(2g(2g-2)) v^(2-2g) = %s)" % (
                "PASS" if ok else "FAIL", name, [str(c) for c in Pc], g, r0, rinf, Fg, pred))
    log("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails else "SOME CHECKS FAILED", time.time() - t0, fails))
    f.close()
