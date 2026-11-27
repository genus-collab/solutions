# -*- coding: utf-8 -*-
"""rf_res.py -- REFEREE: Theorem A(i) (the residue at x = 0 of the full-Stirling 1-form is the quintic's vanishing solution, K = 1)
and the closed form of the vanishing period, re-derived WITHOUT the project's vanishing solution, gauge or periods.

Periods: own Frobenius solver for L4 = u D^4 - (1-u)(2D^3 + (7/5)D^2 + (2/5)D + 24/625), D = -(1-u) d/du (the 5th-order PF is
z^5 (D o L4) - 1, see REFEREE.md s.3), on u^(rho) with dual exponents rho = n + eps (eps^2 = 0).  f_k = u^k + O(u^3) (k = 0,1,2);
g = u + gamma u^2 + .. is the analytic solution that has a log partner g log u + h (gamma is FORCED by the eps-part at u^-1).
Frame (alpha, beta): C0 = f0 + alpha g + beta f2, t_c = g/C0, C1 = D t_c.
Residue: psi = Res_{x=0} e^((l/x + G(x))/z) A(x; z) dx = Sum_k l^k z^-k/k! [x^(k-1)](e^(G/z) A),  l = log(1-u),
G = -Sum_l x^(5l+4)/((l+1)(5l+4)),  A as in rf_sp (power sums p_m by the Newton recurrence, now as x-series).
Criterion (GAUGE Lemma 1, uniqueness of regular solutions): rows Y0 = psi/C0, zD Y_a = c_a Y_(a+1), c = (C1, 1/(u C0^2 C1^2), C1, C0, C0);
if the closing relation zD Y4 = C0 Y0 holds, all Y_a are regular with no positive power of z, and Y(0) = -e_1, then Y = Y^C_Q.
Series: dict {power of z: u-series (list of Fractions)}.
"""
from fractions import Fraction as Fr
from rf_ring import FQ, mul, inv, add, scal, deriv, log1, exp0, powr, stirling_b, pad

F = FQ()


def _dmul(a, b): return (a[0] * b[0], a[0] * b[1] + a[1] * b[0])


def _dinv(a): return (1 / a[0], -a[1] / a[0] ** 2)


def _L4_mono(e):
    """L4 u^(e+eps) as {exponent: dual coefficient} (exponents relative: u^(k+eps))"""
    def D(c):
        out = {}
        for k, x in c.items():                                  # D u^(k+eps) = -(k+eps)(u^(k-1+eps) - u^(k+eps))
            t = _dmul((Fr(k), Fr(1)), x)
            out[k - 1] = add2(out.get(k - 1), (-t[0], -t[1])); out[k] = add2(out.get(k), t)
        return out
    jets = [{e: (Fr(1), Fr(0))}]
    for _ in range(4): jets.append(D(jets[-1]))
    out = {}
    for k, x in jets[4].items(): out[k + 1] = add2(out.get(k + 1), x)                      # u D^4
    for j, cf in ((3, Fr(2)), (2, Fr(7, 5)), (1, Fr(2, 5)), (0, Fr(24, 625))):
        for k, x in jets[j].items():                                                        # -(1-u) cf D^j
            out[k] = add2(out.get(k), (-cf * x[0], -cf * x[1])); out[k + 1] = add2(out.get(k + 1), (cf * x[0], cf * x[1]))
    return out


def add2(a, b): return b if a is None else (a[0] + b[0], a[1] + b[1])


def frobenius(start, NU):
    """dual-coefficient solution Sum_n b_n u^(n+eps) with b_n = start[n] for n <= 2, solved for n >= 3; residues at u^-3..u^-1
    returned (must vanish).  start values are duals."""
    mono = {n: _L4_mono(n) for n in range(NU + 4)}
    b = {n: start.get(n, (Fr(0), Fr(0))) for n in range(3)}
    res = {}
    for k in range(-3, 0):
        r = (Fr(0), Fr(0))
        for n in range(3): r = add2(r, _dmul(b[n], mono[n].get(k, (Fr(0), Fr(0)))))
        res[k] = r
    for n in range(3, NU):
        k = n - 3; r = (Fr(0), Fr(0))
        for m in range(max(0, n - 4), n): r = add2(r, _dmul(b[m], mono[m].get(k, (Fr(0), Fr(0)))))
        ind = mono[n][k]; assert ind[0] == n * (n - 1) ** 2 * (n - 2)
        b[n] = _dmul((-r[0], -r[1]), _dinv(ind))
    return b, res


def periods(NU):
    """f0, f1, f2, g (u-series mod u^NU), gamma; asserts the Frobenius structure"""
    m1 = _L4_mono(1)[-1]                                              # coefficient of u^(-1+eps) in L4 u^(1+eps): K_1(1+eps)
    assert m1[0] == 0                                                 # real part: the solution starting at u exists
    i2 = _L4_mono(2)[-1]; assert i2[0] == 0 and i2[1] == 2            # K_0(2+eps) = 2 eps
    gamma = -m1[1] / 2
    fs = []
    for k in range(3):
        b, res = frobenius({k: (Fr(1), Fr(0))}, NU)
        assert all(r[0] == 0 for r in res.values()), ("analytic solution obstructed", k, res)
        fs.append([b[n][0] for n in range(NU)])
    b, res = frobenius({1: (Fr(1), Fr(0)), 2: (gamma, Fr(0))}, NU)
    assert all(r == (0, 0) for r in res.values()), res               # log partner exists
    return fs[0], fs[1], fs[2], [b[n][0] for n in range(NU)], gamma


def Dop(f): return add(F, scal(F, deriv(F, f), -1), [Fr(0)] + deriv(F, f))     # D f = -(1-u) f'


def zD(Y):
    out = {}
    for e, f in Y.items(): out[e + 1] = Dop(f)
    return out


def ymul(Y, f, n): return {e: mul(F, g, f, n) for e, g in Y.items()}


def A_x(NX, stirling=True):
    """A(x; z) as {z^j: x-series}"""
    x = [Fr(0), Fr(1)] + [Fr(0)] * (NX - 2)
    x5 = [Fr(0)] * NX
    if NX > 5: x5[5] = Fr(1)
    kap = inv(F, add(F, [Fr(1)], x5, -1), NX)
    p = {0: [Fr(5)] + [Fr(0)] * (NX - 1)}
    for r in range(1, 5): p[-r] = [Fr(5 * (-1) ** r)] + [Fr(0)] * (NX - 1)
    for m in range(1, NX):
        cmb = [5 * a + 10 * b + 10 * c + 5 * d + e for a, b, c, d, e in zip(p[m - 1], p[m - 2], p[m - 3], p[m - 4], p[m - 5])]
        p[m] = scal(F, mul(F, kap, cmb, NX), -1)
    terms = {}; xp = [Fr(1)] + [Fr(0)] * (NX - 1)
    for j in range(1, NX):
        xp = mul(F, xp, x, NX)
        if j % 2 and stirling:
            terms[j] = scal(F, mul(F, xp, add(F, [Fr(1, 5 ** j)], p[j]), NX), stirling_b((j + 1) // 2))
    G = {0: [Fr(1)] + [Fr(0)] * (NX - 1)}
    for m in range(1, NX):
        acc = [Fr(0)] * NX
        for j, t in terms.items():
            if j <= m: acc = add(F, acc, scal(F, mul(F, t, G[m - j], NX), j))
        G[m] = scal(F, acc, Fr(1, m))
    pre = powr(F, add(F, [Fr(1)], x5, -1), Fr(-1, 2), NX)
    return {m: mul(F, pre, G[m], NX) for m in G if any(G[m])}, p


def psi_res(NU, stirling=True):
    """psi = Res Omega as {-n: u-series mod u^NU}"""
    NX = NU + 1; A, _ = A_x(NX, stirling)
    Gx = [Fr(0)] * NX
    for l in range(NX):
        if 5 * l + 4 < NX: Gx[5 * l + 4] = Fr(-1, (l + 1) * (5 * l + 4))
    eG = {0: [Fr(1)] + [Fr(0)] * (NX - 1)}; i = 0                     # z^-i: G^i/i!
    while True:
        i += 1; t = scal(F, mul(F, eG[i - 1], Gx, NX), Fr(1, i))
        if not any(t): break
        eG[i] = t
    prod = {}                                                          # e^(G/z) A as {z-power: x-series}
    for i, a in eG.items():
        for j, b in A.items():
            prod[j - i] = add(F, prod.get(j - i, [Fr(0)] * NX), mul(F, a, b, NX))
    lu = log1(F, [Fr(1), Fr(-1)] + [Fr(0)] * (NU - 2), NU)             # log(1-u)
    psi = {}; lk = [Fr(1)] + [Fr(0)] * (NU - 1); fact = 1
    for k in range(1, NU):
        lk = mul(F, lk, lu, NU); fact *= k
        for e, xs in prod.items():
            c = xs[k - 1] if k - 1 < len(xs) else 0
            if c: psi[e - k] = add(F, psi.get(e - k, [Fr(0)] * NU), scal(F, lk, c / fact))
    return {e: f for e, f in psi.items() if any(f)}


def rows(psi, C0, C1, NU):
    iC0 = inv(F, C0, NU); iC1 = inv(F, C1, NU)
    Y0 = ymul(psi, iC0, NU); Y1 = ymul(zD(Y0), iC1, NU)
    w = [Fr(0)] + mul(F, mul(F, C0, C0, NU), mul(F, C1, C1, NU), NU)[:NU - 1]   # u C0^2 C1^2
    Y2 = ymul(zD(Y1), w, NU); Y3 = ymul(zD(Y2), iC1, NU); Y4 = ymul(zD(Y3), iC0, NU)
    close = {e: add(F, f, (ymul(Y0, C0, NU)).get(e, []), -1) for e, f in zD(Y4).items()}
    for e, f in ymul(Y0, C0, NU).items():
        if e not in close: close[e] = scal(F, f, -1)
    return [Y0, Y1, Y2, Y3, Y4], close


def frame(f0, f2, g, al, be, NU):
    C0 = add(F, add(F, f0, scal(F, g, al)), scal(F, f2, be))
    tc = mul(F, g, inv(F, C0, NU), NU)
    return C0, Dop(tc)
