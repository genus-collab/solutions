# -*- coding: utf-8 -*-
"""rf_t3_swapformula.py -- referee checks of the proof of Theorem G (TR_MODEL.md s.4.1), own code.

(a) A_g = Sum_j (-1)^j (2g-2+j)!/j! [t^g] phi^j, phi = Sum_{k>=1} t^k/(4^k (2k+1)), against -B_2g/(2g) AND against the series form
    (2g-2)! 4^-g [s^1] artanh(s)^(1-2g) computed by series arithmetic (no closed form for B_2g used), g = 2..14, exact.
(b) the y = z formula of arXiv:2502.18115 Cor. 1.2 (1.3), transcribed independently as
       F_g = (1/(2-2g)) Sum_b Res_b dz Sum_j (1/j!) D^(2g-2+j)(z) Sum_{k_1+..+k_j = g} prod x^(2k_i)/(4^k_i (2k_i+1)!),  D = d/dx,
    full integrand (no leading-term shortcut), for x = P(z) + Sum_poles v/(z - b), b running over the poles and infinity.
    One pole, several P, g = 2..5: Res_inf = 0 and total = B_2g/(2g(2g-2)) v^(2-2g) (exact over Q).
(c) end-to-end: TWO simple poles, x = z + v1/z + v2/(z - b), y = z: formula (b) against the referee's TR, mod 2^61-1 and 2^89-1
    (branch points in F_p by construction), g = 2, 3; prediction of the locality mechanism: B_2g/(2g(2g-2)) (v1^(2-2g) + v2^(2-2g)).
Writes rf_t3_swapformula.log.  Run: py -I -B rf_t3_swapformula.py
"""
import os, sys, time, random, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from rf_ser import Ls, QField, FpField, P61, P89, bernoulli, ps_power
from rf_tr import TR, Curve, Mmax


def A_sum(g):
    phi = [Fr(0)] + [Fr(1, 4 ** k * (2 * k + 1)) for k in range(1, g + 1)]
    tot = Fr(0); pw = [Fr(1)] + [Fr(0)] * g
    for j in range(1, g + 1):
        pw = [sum(pw[a] * phi[b - a] for a in range(b + 1)) for b in range(g + 1)]
        tot += Fr((-1) ** j * math.factorial(2 * g - 2 + j), math.factorial(j)) * pw[g]
    return tot


def A_series(g):
    """(2g-2)! 4^-g [s^(2g)] (artanh(s)/s)^(1-2g)"""
    N = 2 * g; F = QField()
    a = [Fr(1, k + 1) if k % 2 == 0 else Fr(0) for k in range(N + 1)]
    b = ps_power(a, Fr(1 - 2 * g), N, F)
    return Fr(math.factorial(2 * g - 2), 4 ** g) * b[2 * g]


def comps(g, j):
    if j == 1: yield (g,); return
    for a in range(1, g - j + 2):
        for r in comps(g - a, j - 1): yield (a,) + r


def I_g(F, x, zf, dz, g):
    dx = dz(x).inv(F); D = [zf]
    for _ in range(1, 3 * g): D.append(dz(D[-1]).mul(dx, F))
    X = [x]
    for _ in range(2 * g): X.append(dz(X[-1]))
    tot = None
    for j in range(1, g + 1):
        for ks in comps(g, j):
            term = D[2 * g - 2 + j].scale(F.of(Fr(1, math.factorial(j))))
            for k in ks: term = term.mul(X[2 * k].scale(F.of(Fr(1, 4 ** k * math.factorial(2 * k + 1)))), F)
            tot = term if tot is None else tot.add(term, F)
    return tot


def _pw(F, a, e):
    r = F.one
    for _ in range(e): r = r * a
    return r


def formula_Fg(F, Pc, poles, g, T=60):
    """returns (F_g, list of residues at the poles, residue at infinity); Pc, poles = [(b, v)] in F"""
    r = max(k for k, c in enumerate(Pc) if not F.iszero(F.of(c))); res = []
    for (b0, v0) in poles:                                           # local variable q = z - b0
        reg = [F.zero] * (T + 1)
        for k, ck in enumerate(Pc):
            for i in range(k + 1): reg[i] = reg[i] + F.of(ck) * F.of(math.comb(k, i)) * _pw(F, b0, k - i)
        for (b1, v1) in poles:
            if b1 == b0: continue
            d = F.one / (b0 - b1); pw = d
            for i in range(T + 1): reg[i] = reg[i] + v1 * pw * F.of((-1) ** i); pw = pw * d
        x = Ls(-1, [v0] + reg[:T], T - 1); zf = Ls(0, [b0, F.one] + [F.zero] * (T - 1), T)
        res.append(I_g(F, x, zf, lambda f: f.deriv(F), g).coef(-1, F))
    c = [F.zero] * (T + r + 2)                                       # infinity, w = 1/z, d/dz = -w^2 d/dw
    for k, ck in enumerate(Pc): c[r - k] = c[r - k] + F.of(ck)
    for (b0, v0) in poles:
        pw = F.one
        for i in range(T + 1): c[r + 1 + i] = c[r + 1 + i] + v0 * pw; pw = pw * b0
    x = Ls(-r, c, T); w2 = Ls(2, [-F.one] + [F.zero] * (T + 10), T + 12)
    rinf = -(I_g(F, x, Ls(-1, [F.one] + [F.zero] * (T + 1), T), lambda f: w2.mul(f.deriv(F), F), g).coef(1, F))
    return (sum(res, F.zero) + rinf) / F.of(2 - 2 * g), res, rinf


def two_pole_curve(F, rng):
    """x = z + v1/z + v2/(z - b), y = z, 4 branch points in F_p (choose c1..c3, solve the quadratic constraint for c4; p = 3 mod 4)"""
    p = F.p
    while True:
        c = [F.of(rng.randrange(1, 10 ** 6)) for _ in range(3)]
        s1 = c[0] + c[1] + c[2]; s2 = c[0] * c[1] + c[0] * c[2] + c[1] * c[2]; s3 = c[0] * c[1] * c[2]
        A, B, C = s2, s1 * s2 - s3 * 3, s1 * s3                       # 4 e4 = e1 e3:  s2 c4^2 + (s1 s2 - 3 s3) c4 + s1 s3 = 0
        disc = B * B - A * C * 4; r = F.E(pow(disc.v, (p + 1) // 4, p))
        if not F.iszero(r * r - disc): continue
        c4 = (r - B) / (A * 2); cs = c + [c4]
        e1 = s1 + c4; e2 = s2 + c4 * s1; e3 = s3 + c4 * s2; e4 = c4 * s3
        b = e1 / 2; v1 = -(e3 / e1); v2 = b * b - v1 - e2
        assert F.iszero(e4 + v1 * b * b)
        if len(set(x.v for x in cs)) == 4 and not F.iszero(v1) and not F.iszero(v2) and not F.iszero(b): break
    def xt(a, N):
        out = [a, F.one] + [F.zero] * (N - 1)
        for (bb, vv) in ((F.zero, v1), (b, v2)):
            d = F.one / (a - bb); pw = d
            for k in range(N + 1): out[k] = out[k] + vv * pw * F.of((-1) ** k); pw = pw * d
        return out
    def yxt(a, N):
        X = xt(a, N + 1); dX = [X[k + 1] * F.of(k + 1) for k in range(N + 1)]
        return [a * dX[k] + (dX[k - 1] if k >= 1 else F.zero) for k in range(N + 1)]
    for ci in cs:
        X = xt(ci, 2); assert F.iszero(X[1]) and not F.iszero(X[2])
    return Curve(F, cs, xt, yxt, "two poles"), b, v1, v2


def main():
    log = open(os.path.join(HERE, "rf_t3_swapformula.log"), "w"); t0 = time.time(); fails = [0]
    def say(s): print(s, flush=True); log.write(s + "\n"); log.flush()
    def check(ok, s):
        say("  [%s] %s" % ("PASS" if ok else "FAIL", s)); fails[0] += (not ok)
    say("rf_t3_swapformula.py  %s" % time.strftime("%Y-%m-%d %H:%M"))
    for g in range(2, 15):
        a1, a2, a3 = A_sum(g), A_series(g), -bernoulli(2 * g) / (2 * g)
        check(a1 == a2 == a3, "(a) g = %d: A_g (finite sum) = %s, series form = %s, -B_2g/(2g) = %s" % (g, a1, a2, a3))
    F = QField(); rng = random.Random(5)
    plist = [("z^4/4", [0, 0, 0, 0, Fr(1, 4)]), ("z", [0, 1]), ("z^2/2 - 3z", [0, -3, Fr(1, 2)]),
             ("random quartic", [0] + [Fr(rng.randrange(-9, 10), rng.randrange(1, 5)) for _ in range(3)] + [Fr(1, 4)]),
             ("random quintic", [0] + [Fr(rng.randrange(-9, 10), rng.randrange(1, 5)) for _ in range(4)] + [Fr(1, 5)])]
    for name, Pc in plist:
        v = Fr(rng.randrange(1, 9), rng.randrange(1, 4))
        for g in range(2, 6):
            Fg, res, rinf = formula_Fg(F, Pc, [(F.zero, F.of(v))], g)
            pred = bernoulli(2 * g) / (2 * g * (2 * g - 2)) * v ** (2 - 2 * g)
            check(rinf == 0 and Fg == pred, "(b) P = %s, v = %s, g = %d: Res_0 = %s, Res_inf = %s, F_g = %s (B_2g/(2g(2g-2)) v^(2-2g) = %s)" % (name, v, g, res[0], rinf, Fg, pred))
    for p in (P61, P89):
        Fp = FpField(p); C, b, v1, v2 = two_pole_curve(Fp, random.Random(p % 1000))
        tr = TR(C, 2 * Mmax(2, 2) + 6)
        for g in (2, 3):
            t = time.time(); ftr = tr.Fg(g); ff, res, rinf = formula_Fg(Fp, [0, 1], [(Fp.zero, v1), (b, v2)], g)
            loc = Fp.of(bernoulli(2 * g) / (2 * g * (2 * g - 2))) * (Fp.one / _pw(Fp, v1, 2 * g - 2) + Fp.one / _pw(Fp, v2, 2 * g - 2))
            check(Fp.iszero(ftr - ff) and Fp.iszero(ftr - loc) and Fp.iszero(rinf),
                  "(c) mod 2^%d-1, two simple poles, g = %d: TR F_g == formula (1.3) == B_2g/(2g(2g-2)) (v1^(2-2g) + v2^(2-2g)); Res_inf = 0 (%.1f s)" % (
                      p.bit_length(), g, time.time() - t))
        d = tr.Fg(2) - Fp.of(bernoulli(4) / 8) / (v1 * v1)
        check(not Fp.iszero(d), "(c) control: the one-pole value B_4/8 v1^-2 alone does NOT reproduce F_2 of the two-pole curve")
    say("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not fails[0] else "SOME CHECKS FAILED", time.time() - t0, fails[0]))
    log.close()


if __name__ == "__main__":
    main()
