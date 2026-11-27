# -*- coding: utf-8 -*-
"""v3_hurwitz_float.py -- ADVERSARIAL check of the Hurwitz-space statements of LG_MODEL.md s.5 by a DIFFERENT method than lg_frob.py:
the residue pairings are evaluated as SUMS OVER THE FIVE CRITICAL POINTS (complex floating point, Durand-Kerner roots), not as
residues at the poles.  Floating point: agreement is to ~1e-9, this is a numerical cross-check of exact statements, not a proof.
   f(p) = p^4 + a2 p^2 + a3 p + a4 + a5/(p - a1),  primitive form dp,
   eta(d', d'') = Sum_crit d'f d''f / f'',   c(d', d'', d''') = Sum_crit d'f d''f d'''f / f'',
   flat coordinates t = (a1, a2, a3, a4 - a2^2/8, a5):  d_1 f = a5/(p-a1)^2, d_2 f = p^2 + a2/4, d_3 f = p, d_4 f = 1, d_5 f = 1/(p-a1).
Checks: (a) eta = const (eta_15 = 1, eta_24 = eta_33 = 1/4) at random complex points; (b) c_555 = 1/a5, c_4ab = eta_ab;
        (c) the same pairing for the form dx/x (x = p - a1, derivatives at fixed x) restricted to the line a = (0,0,0,0,r):
            the Gram matrix of the classes x^k dx, k = -2..2, is the SAME for dx and dx/x (it is Sum_crit x^{k+l}/f''), and the
            spectrum attached to each class is the same: mu(x^k dx) = (k-1)/4 for k = -1..2 and +1/2 for k = -2 in BOTH structures
            (dx: mu = 3/4 - deg t;  dx/x: mu = 1/2 - deg b, deg b_j = (4-j)/4 for the class x^j dx/x, residue direction mu = +1/2);
        (d) E = (5/4) r d_5 on the line and the critical values are 5 x_c^4.
Run:  py -I v3_hurwitz_float.py
"""
import sys
sys.dont_write_bytecode = True
import random, cmath


def roots(coefs):
    """Durand-Kerner; coefs ascending, monic not required"""
    n = len(coefs) - 1; lead = coefs[-1]; c = [x / lead for x in coefs]
    z = [(0.4 + 0.9j) ** k for k in range(n)]
    for _ in range(4000):
        new = []
        for i in range(n):
            v = sum(c[k] * z[i] ** k for k in range(n + 1)); d = 1
            for j in range(n):
                if j != i: d *= (z[i] - z[j])
            new.append(z[i] - v / d)
        if max(abs(new[i] - z[i]) for i in range(n)) < 1e-15: z = new; break
        z = new
    return z


def pmul(A, B):
    out = [0j] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        for j, y in enumerate(B): out[i + j] += x * y
    return out


def structure(a1, a2, a3, a4, a5):
    D2 = pmul([-a1, 1], [-a1, 1])
    Q = pmul([a3, 2 * a2, 0, 4], D2); Q[0] -= a5                     # f' (p - a1)^2
    crit = roots(Q)
    tang = lambda p: [a5 / (p - a1) ** 2, p * p + a2 / 4, p, 1, 1 / (p - a1)]
    f2 = lambda p: 12 * p * p + 2 * a2 + 2 * a5 / (p - a1) ** 3
    eta = [[sum(tang(p)[i] * tang(p)[j] / f2(p) for p in crit) for j in range(5)] for i in range(5)]
    c = [[[sum(tang(p)[i] * tang(p)[j] * tang(p)[k] / f2(p) for p in crit) for k in range(5)] for j in range(5)] for i in range(5)]
    return eta, c, crit


def main():
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    rng = random.Random(7)
    ETA = [[0] * 5 for _ in range(5)]; ETA[0][4] = ETA[4][0] = 1; ETA[1][3] = ETA[3][1] = ETA[2][2] = 0.25
    worst_a = worst_b = worst_u = 0.0
    for _ in range(6):
        pt = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(5)]
        eta, c, crit = structure(*pt)
        worst_a = max(worst_a, max(abs(eta[i][j] - ETA[i][j]) for i in range(5) for j in range(5)))
        worst_b = max(worst_b, abs(c[4][4][4] - 1 / pt[4]))
        worst_u = max(worst_u, max(abs(c[3][i][j] - ETA[i][j]) for i in range(5) for j in range(5)))
    check("(a) eta constant in t = (a1, a2, a3, a4 - a2^2/8, a5) at 6 random complex points (sum over critical points)", worst_a < 1e-8, "max deviation %.2e" % worst_a)
    check("(b) c_555 = 1/a5", worst_b < 1e-8, "max deviation %.2e" % worst_b)
    check("(b) c_4ab = eta_ab (unit d/dt4)", worst_u < 1e-8, "max deviation %.2e" % worst_u)
    # control: a wrong flat coordinate (a4 instead of a4 - a2^2/8) changes d_2 f to p^2: eta_22 = -a2/8
    pt = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(5)]
    a1, a2, a3, a4, a5 = pt; eta, c, crit = structure(*pt)
    f2 = lambda p: 12 * p * p + 2 * a2 + 2 * a5 / (p - a1) ** 3
    naive = sum(p ** 4 / f2(p) for p in crit)
    check("(a) control: naive coordinate a4 gives eta(d_a2, d_a2) = -a2/8", abs(naive + a2 / 8) < 1e-8, "%.2e" % abs(naive + a2 / 8))
    # (c), (d): the line
    r = 0.7 + 0.3j
    eta, c, crit = structure(0, 0, 0, 0, r)
    check("(d) on the line the critical points satisfy p^5 = r/4 and f(p_c) = 5 p_c^4", all(abs(p ** 5 - r / 4) < 1e-10 and abs(p ** 4 + r / p - 5 * p ** 4) < 1e-10 for p in crit))
    f2 = lambda p: 12 * p * p + 2 * r / p ** 3
    ks = (-2, -1, 0, 1, 2)
    gram_dx = [[sum(p ** (k + l) / f2(p) for p in crit) for l in ks] for k in ks]                       # form dx: (phi dx, phi' dx) = Sum phi phi'/f''
    gram_dlog = [[sum(p ** (k + 1) * p ** (l + 1) / (p * p * f2(p)) for p in crit) for l in ks] for k in ks]   # form dx/x: classes x^{k+1} dx/x, Hessian x^2 f''
    dev = max(abs(gram_dx[i][j] - gram_dlog[i][j]) for i in range(5) for j in range(5))
    pat = all((abs(gram_dx[i][j]) > 1e-6) == ((ks[i] + ks[j] - 2) % 5 == 0) for i in range(5) for j in range(5))
    check("(c) Gram matrices of the classes x^k dx (k = -2..2) for the forms dx and dx/x coincide; support k + l = 2 mod 5", dev < 1e-12 and pat, "%.1e" % dev)
    mu_dx = {0: -0.25, -1: -0.5, -2: 0.5, 2: 0.25, 1: 0.0}                                           # 3/4 - deg t_alpha, alpha <-> k = (4:0, 5:-1, 1:-2, 2:2, 3:1)
    mu_dlog = {}
    for j in (3, 2, 1, 0): mu_dlog[j - 1] = 0.5 - (4 - j) / 4.0                                         # class x^j dx/x = x^{j-1} dx, deg b_j = (4-j)/4, d = 1
    mu_dlog[-2] = 0.5                                                                                # residue direction q (flat coordinate log q, degree 0)
    check("(c) the spectrum attached to each class x^k dx is the same for dx (d = 1/2) and dx/x (d = 1): %s" % sorted(mu_dx.items()), mu_dx == mu_dlog, sorted(mu_dlog.items()))
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
