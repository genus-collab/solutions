# -*- coding: utf-8 -*-
"""cy_rows.py N NZ -- referee cyc_adversary: the Corollary of CYC_PROOF s.0 in the two directions the programme uses it.
(G) cyclic system => PF5 (P(b) step 5, POLY Lemma 2 for Phi^can): the vanishing solution of the quintic cyclic system built
    from the cyclic system ALONE (iterated integrals of GI Lemma Q: y^(0) = -e_1, y^(n+1)_a = Int_0^u b_a y^(n)_(a+1)),
    then Psi = varpi0 Y_0 = Sum_n z^-n varpi0 y^(n)_0 must satisfy PF5 = z^5 D PF4 - 1, i.e. D PF4(Psi_n) = Psi_(n-5).
(H) PF5 + the four forward relations => closing relation 4 -> 0 (UR s.1, s.4): the colour-0 WKB solution psi of PF5 (own
    recursion in L), rows X_0 = psi/varpi0, X_(a+1) = (zD + L) X_a / c_a (a = 0..3, the order of fq_rows' kinv), then
    (zD + L) X_4 = c_4 X_0 at every z-order, in the conifold frame; c = (Dt, L^5/(varpi0^2 Dt^2), Dt, varpi0, varpi0).
Exact over Q; u-series through u^(N-1); z-orders 0..NZ in (H).  Usage: py -I -B cy_rows.py 40 8"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cy_core import Ser, ONE_M, DL, op_apply, op_D_left
from cy_frob import pf4_poly, pf4_ser, frobenius
from cy_log import Log

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 8
LG = Log('cy_rows_%d_%d' % (N, NZ))
LG('cy_rows: exact over Q, u-series through u^%d, WKB z-orders 0..%d' % (N - 1, NZ))
A = pf4_poly(); As = pf4_ser(A); FR = frobenius(A, N)
f0, f1, f2 = [Ser(0, x, N) for x in FR['f']]; g = Ser(0, FR['g'], N); U = Ser.poly({1: 1})
def lin(*ps):
    acc = Ser(0, [])
    for x, s in ps: acc = acc.add(s.scale(x))
    return acc


DPF4 = op_D_left(As)


def cycl(y0, y1, swap=False, five=False):
    """c = (Dt, 1/(u y0^2 Dt^2), Dt, y0, y0) (L^5 = 1/u); controls: swap c_1 <-> c_2, or c_1 -> 5 c_1"""
    Dt = y1.mul(y0.inv()).D(); c1 = U.mul(y0).mul(y0).mul(Dt).mul(Dt).inv()
    c = [Dt, c1.scale(5) if five else c1, Dt, y0, y0]
    if swap: c[1], c[2] = c[2], c[1]
    return c


# ------------------------------------------------------------------ (G) cyclic => PF5
def vanishing(c, nmax):
    b = [x.mul(ONE_M.inv(cap=N)).neg() for x in c]                       # th Y_a = b_a Y_(a+1), b_a = -c_a/(1-u)
    y = [[Ser(0, []), Ser(0, [-1]), Ser(0, []), Ser(0, []), Ser(0, [])]]
    for n in range(nmax):
        y.append([b[a].mul(y[-1][(a + 1) % 5]).integ0() for a in range(5)])
    return y


def check_G(label, y0, y1, nmax, expect=True, **kw):
    y = vanishing(cycl(y0, y1, **kw), nmax)
    Psi = [y0.mul(yn[0]) for yn in y]
    ok = True; pr = N; first_bad = None
    for n in range(nmax + 1):
        lhs = op_apply(DPF4, Psi[n]); rhs = Psi[n - 5] if n >= 5 else Ser(0, [])
        d = lhs.sub(rhs); pr = min(pr, d.prec)
        if not d.is_zero():
            ok = False; first_bad = n if first_bad is None else first_bad
    info = 'n = 0..%d, precision >= u^%d%s' % (nmax, pr, '' if ok else ', first failure at n = %d' % first_bad)
    LG.check('(G) %s: Psi = varpi0 (Y^C_Q)_0 solves PF5 (D PF4 Psi_n = Psi_(n-5))' % label, ok, info=info, control=not expect)
    return y


cc = [0, -1, 3, 2, 1]
for lab, y0 in (('frame (0,0)', f0), ('frame (1/3,-2/7)', lin((1, f0), (Fr(1, 3), g), (Fr(-2, 7), f2))),
                ('frame (7/11,-3/5)', lin((1, f0), (Fr(7, 11), g), (Fr(-3, 5), f2)))):
    y = check_G(lab, y0, g, 16)
    grad = all(y[n][a].is_zero() for n in range(17) for a in range(5) if (n + a - 1) % 5)
    vals = all(y[n][a].val() == (4 * n + cc[a] + 1) // 5 for n in range(17) for a in range(5) if (n + a - 1) % 5 == 0)
    LG.check('(G) %s: GI Lemma Q(a): y^(n)_a = 0 unless n = 1 - a mod 5, ord_u = (4n + c_a + 1)/5, c = (0,-1,3,2,1), n <= 16' % lab,
             grad and vals)
    LG.check('(G) %s: y^(1)_0 = -t_c (GI Lemma Q(d))' % lab, y[1][0].add(g.mul(y0.inv())).is_zero())
check_G('control, wrong gamma: y1 = f1 + f2/2', f0, lin((1, f1), (Fr(1, 2), f2)), 16, expect=False)
check_G('control, non-isotropic: y1 = f1', f0, f1, 16, expect=False)
check_G('control, c_1 -> 5 c_1 (product 5 L^5)', f0, g, 16, expect=False, five=True)
check_G('control, c_1 <-> c_2 swapped', f0, g, 16, expect=False, swap=True)


# ------------------------------------------------------------------ (H) PF5 + forward relations => closing relation
def lp_add(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        x = out.get(k, 0) + s * v
        if x: out[k] = x
        else: out.pop(k, None)
    return out


def lp_sh(a, k): return {e + k: v for e, v in a.items()}
def lp_sc(a, s): return {e: s * v for e, v in a.items() if s * v}


def lp_D(a):                                    # D L^k = (k/5)(L^(k+5) - L^k)
    out = {}
    for k, v in a.items():
        if k: out = lp_add(out, {k + 5: Fr(k, 5) * v, k: Fr(-k, 5) * v})
    return out


def lp_mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items(): out = lp_add(out, {i + j: x * y})
    return out


QL = {0: Fr(1, 3125), -5: Fr(-1, 3125)}         # q = (1 - L^-5)/3125


def pcal(F):
    """colour-0 PF5 conjugated by e^(u_0/z): Nn^5 - 1 - q prod_k (5 Nn + k z), Nn = zD + L, on a z-series F"""
    Nn = lambda G: [lp_add(lp_D(G[m - 1]) if m else {}, lp_sh(G[m], 1)) for m in range(len(G))]
    G = F
    for _ in range(5): G = Nn(G)
    H = F
    for k in range(1, 6):
        NH = Nn(H); H = [lp_add(lp_sc(NH[m], 5), lp_sc(H[m - 1], k) if m else {}) for m in range(len(H))]
    return [lp_add(lp_add(G[m], F[m], -1), lp_mul(QL, H[m]), -1) for m in range(len(F))]


def solve_D(r):
    """phi with D phi = r, zero constant term; None if no Laurent-polynomial solution (a logarithm)"""
    if not r: return {}
    m = min(r); rho = {e - m: v for e, v in r.items()}; sig = {}
    for e in range(max(rho), 4, -1):
        x = rho.get(e, 0)
        if x: sig[e - 5] = sig.get(e - 5, 0) + x; rho = lp_add(rho, {e: -x, e - 5: x})
    if rho: return None
    dphi = {e + m - 1: 5 * v for e, v in sig.items()}
    if dphi.get(-1, 0): return None
    phi = {e + 1: v / (e + 1) for e, v in dphi.items() if v}
    return phi if lp_add(lp_D(phi), r, -1) == {} else None


psi = [{1: Fr(1)}]; okw = all(not x for x in pcal(psi + [{}]))      # Pcal_0 = 0 (symbol) and Pcal_1 L = 0
for n in range(1, NZ + 1):
    T = pcal(psi + [{}, {}])[n + 1]; phi = solve_D(lp_sc(T, Fr(-1, 5)))
    if phi is None: okw = False; break
    psi.append(lp_sh(phi, 1)); okw = okw and all(not x for x in pcal(psi + [{}])[:n + 2])
LG.check('(H) colour-0 WKB solution of PF5: psi_0 = L, psi_n = L phi_n, phi_n Laurent polynomials (no log), Pcal psi = 0 through z^%d' % (NZ + 1),
         okw, info='top L-degree of psi_n: %s' % [max(p) for p in psi])


def to_cls(p, n):
    """psi_n(L) -> (rho, U) with psi_n = s^rho U(u), s = 1/L, rho = (n-1) mod 5"""
    rho = (n - 1) % 5; d = {}
    for k, v in p.items():
        assert (-k - rho) % 5 == 0; d[(-k - rho) // 5] = v
    return rho, Ser.poly(d)


def Lm(x): return (x[0] - 1, x[1]) if x[0] else (4, x[1].shift(-1))
def Dc(x): return x[0], x[1].D().add(x[1].mul(DL).scale(Fr(x[0], 5)))
def addc(x, y): assert x[0] == y[0]; return x[0], x[1].add(y[1])


def check_H(label, y0, y1, expect=True, **kw):
    c = cycl(y0, y1, **kw); ic = [x.inv() for x in c]; iy0 = y0.inv()
    X = [[(r, Uu.mul(iy0)) for r, Uu in (to_cls(psi[n], n) for n in range(NZ + 1))]]
    for a in range(4):
        X.append([(lambda t: (t[0], t[1].mul(ic[a])))(addc(Dc(X[a][n - 1]), Lm(X[a][n])) if n else Lm(X[a][0])) for n in range(NZ + 1)])
    bad = []; pr = N
    for n in range(NZ + 1):
        Z = addc(Dc(X[4][n - 1]), Lm(X[4][n])) if n else Lm(X[4][0])
        d = addc(Z, (X[0][n][0], X[0][n][1].mul(c[4]).neg())); pr = min(pr, d[1].prec)
        if not d[1].is_zero(): bad.append(n)
    LG.check('(H) %s: closing relation (zD + L) X_4 = c_4 X_0 at z^0..z^%d' % (label, NZ), not bad,
             info='precision u^%d%s' % (pr, ', fails at z^%s' % bad if bad else ''), control=not expect)


for lab, y0 in (('frame (0,0)', f0), ('frame (1/3,-2/7)', lin((1, f0), (Fr(1, 3), g), (Fr(-2, 7), f2))),
                ('frame (7/11,-3/5)', lin((1, f0), (Fr(7, 11), g), (Fr(-3, 5), f2))),
                ('rescaled frame 2 varpi0(1/3,-2/7), 3 g', lin((2, f0), (Fr(2, 3), g), (Fr(-4, 7), f2)))):
    check_H(lab, y0, g if 'rescaled' not in lab else g.scale(3))
check_H('control, wrong gamma: y1 = f1 + f2/2', f0, lin((1, f1), (Fr(1, 2), f2)), expect=False)
check_H('control, non-isotropic: y1 = f1', f0, f1, expect=False)
check_H('control, c_1 -> 5 c_1', f0, g, expect=False, five=True)
check_H('control, c_1 <-> c_2 swapped', f0, g, expect=False, swap=True)
LG.close()
