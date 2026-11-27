# -*- coding: utf-8 -*-
"""cy_run.py N -- referee cyc_adversary: the ingredients of CYC_PROOF s.1-s.4 and the PF5 identity, own code, exact over Q,
series through u^(N-1).  (A) closed form and Frobenius basis; (B) self-adjointness (C1) with controls; (C) the concomitant
(C2, C4); (D) transversality (C3); (G) PF5 = z^5 D PF4 - 1 and the u-form behind HPF.  Usage: py -I -B cy_run.py 40"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cy_core import Ser, Rel, ONE_M, op_add, op_D_left, op_compose, op_adjoint, op_apply, op_diff_zero
from cy_frob import pf4_poly, pf4_ser, closed_form, frobenius, peval
from cy_log import Log

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LG = Log('cy_run_%d' % N)
LG('cy_run: exact over Q, series through u^%d' % (N - 1))

# ------------------------------------------------------------------ (A) closed form, Frobenius basis
A = pf4_poly(); As = pf4_ser(A); P = closed_form(A)
LG('(A) closed form PF4 u^x = Sum_j P_j(x) u^(x+j): j in %s' % sorted(P))
LG.check('(A1) P_-3(x) = x (x-1)^2 (x-2) (indicial polynomial, exponents 0, 1, 1, 2)',
         all(peval(P[-3], x) == x * (x - 1) ** 2 * (x - 2) for x in range(-3, 9)))
FR = frobenius(A, N); F = FR['F']; f = FR['f']
LG.check('(A2) PF4 maps Q[[u]] to Q[[u]]: P_j(n) = 0 whenever n + j < 0, n in {0, 1, 2}',
         all(peval(P[j], n) == 0 for j in P for n in range(3) if n + j < 0))
LG.check('(A3) f0 = 1 + (2/625)u^3, f1 = u - (19/60)u^3, f2 = u^2 + (37/30)u^3 (+ O(u^4)): the values of fq_conifold/UR',
         [f[0][3], f[1][3], f[2][3]] == [Fr(2, 625), Fr(-19, 60), Fr(37, 30)])
ob = FR['ob']
LG('    log obstructions (u^-3, u^-2, u^-1) of f0, f1, f2: %s' % [[str(x) for x in o] for o in ob])
LG.check('(A4) f1, f2 unobstructed at u^-3, u^-2; f0 obstructed at u^-3 (so g has no f0 part); gamma FOUND = 7/10',
         ob[1][:2] == [0, 0] and ob[2][:2] == [0, 0] and ob[0][0] != 0 and FR['gamma'] == Fr(7, 10), 'gamma = %s' % FR['gamma'])
LG.check('(A5) g = u + (7/10)u^2 + (41/75)u^3 + (1133/2500)u^4 (FORMAL_QUINTIC_genus2 s.1)',
         FR['g'][:5] == [0, 1, Fr(7, 10), Fr(41, 75), Fr(1133, 2500)] and FR['obg'] == [0, 0, 0])
basis = [Rel.of(Ser(0, x, N)) for x in f] + [Rel([Ser(0, FR['h'], N), Ser(0, FR['g'], N)])]
res = [op_apply(As, y) for y in basis]
LG.check('(A6) cross-method: the operator applied to f0, f1, f2, y3 = g log u + h vanishes', all(r.is_zero() for r in res),
         'precision u^%d' % min(c.prec for r in res for c in r.cs))

# ------------------------------------------------------------------ (B) self-adjointness
ok, pr = op_diff_zero(op_adjoint(As), As)
LG.check('(B1) C1: PF4* = Sum (-D)^k A_k = PF4 exactly', ok)
D = lambda s: s.D()
ok3 = As[3].sub(D(As[4]).scale(2)).is_zero(); ok1 = As[1].sub(D(As[2]).sub(D(D(D(As[4]))))).is_zero()
LG.check('(B2) own criterion: L = L* iff A_3 = 2 D A_4 and A_1 = D A_2 - D^3 A_4 (A_0 free); PF4 satisfies both', ok3 and ok1)
bad = op_add(As, [Ser(0, []), ONE_M.scale(Fr(1, 7))])
LG.check('(B3) control: PF4 + (1/7)(1-u) D is self-adjoint', op_diff_zero(op_adjoint(bad), bad)[0], control=True)


def M_op(y0, y1, a4):
    """(1/y0) D (1/Dt) D (a4 y0^2 Dt^2) D (1/Dt) D (1/y0), t = y1/y0"""
    t = y1.mul(y0.inv()); Dt = t.D(); iy0 = y0.inv(); iDt = Dt.inv()
    E = [iy0]
    for fac in (iDt, a4.mul(y0).mul(y0).mul(Dt).mul(Dt), iDt, iy0):
        E = op_D_left(E); E = [fac.mul(e) for e in E]
    return E


y0r = Ser(0, [1, Fr(3, 7), Fr(-2, 9), 5, Fr(1, 11)] + [Fr(k % 7 - 3, k + 2) for k in range(N - 5)], N)
y1r = Ser(1, [Fr(-2), 1, Fr(4, 3)] + [Fr((k * k) % 5 - 2, k + 3) for k in range(N - 4)], N)
Mr = M_op(y0r, y1r, Ser.poly({1: 1}))
LG.check('(B4) M(y0, t) is formally self-adjoint for ARBITRARY units y0, Dt (not solutions): no non-self-adjoint operator is an M',
         op_diff_zero(op_adjoint(Mr), Mr)[0], 'precision u^%d' % op_diff_zero(op_adjoint(Mr), Mr)[1])
LG.check('(B5) control: that random M equals PF4', op_diff_zero(Mr, As)[0], control=True)

# ------------------------------------------------------------------ (C) the concomitant
def concomitant(Aser, fr, gr):
    fj = [fr]
    for _ in range(3): fj.append(fj[-1].D())
    B = None
    for k in range(1, 5):
        t = gr.mulser(Aser[k])
        for j in range(k):
            term = t.mul(fj[k - 1 - j])
            if j % 2: term = term.neg()
            B = term if B is None else B.add(term)
            t = t.D()
    return B


def const_of(r):
    """(is a log-free constant, value, precision)"""
    ok = all(c.is_zero() for c in r.cs[1:]) and all(x == 0 for i, x in enumerate(r.cs[0].c) if r.cs[0].v + i != 0)
    return ok, (r.cs[0].coef(0) if r.cs[0].prec > 0 else None), min(c.prec for c in r.cs)


Bm = [[None] * 4 for _ in range(4)]; okc = True; prc = N
for i in range(4):
    for j in range(4):
        ok_ij, val, p = const_of(concomitant(As, basis[i], basis[j])); okc = okc and ok_ij; prc = min(prc, p); Bm[i][j] = val
LG('    Bmat (basis f0, f1, f2, y3): %s' % [[str(x) for x in row] for row in Bm])
LG.check('(C1) C2(ii): B(y_i, y_j) log-free and constant on the basis', okc, 'precision u^%d' % prc)
crit = [[0, Fr(-7, 5), 2, Fr(-9, 10)], [Fr(7, 5), 0, 0, -1], [-2, 0, 0, 0], [Fr(9, 10), 1, 0, 0]]
LG.check('(C2) Bmat equals CYC_PROOF s.3 and is antisymmetric; Pfaffian = 2, det = 4',
         Bm == crit and all(Bm[i][j] == -Bm[j][i] for i in range(4) for j in range(4))
         and Bm[0][1] * Bm[2][3] - Bm[0][2] * Bm[1][3] + Bm[0][3] * Bm[1][2] == 2)
gR = Rel.of(Ser(0, FR['g'], N))
bg = [const_of(concomitant(As, y, gR)) for y in basis]
LG.check('(C3) C4: B(f_i, g) = 0 (i = 0, 1, 2) and B(y3, g) = 1', all(b[0] for b in bg) and [b[1] for b in bg] == [0, 0, 0, 1])
LG.check('(C4) formal Picard-Lefschetz consistency: the u^-1 obstructions of f1, f2 equal B(f0, f1), B(f0, f2)',
         [ob[1][2], ob[2][2]] == [Bm[0][1], Bm[0][2]], '(-7/5, 2): gamma = 7/10 <=> B(f0, g) = 0')
def btop(fr, gr):
    """order-3 part Sum_j (-1)^j u D^j g D^(3-j) f of the concomitant"""
    acc = None; gj = gr
    for j in range(4):
        fk = fr
        for _ in range(3 - j): fk = fk.D()
        t = gj.mul(fk).mulser(As[4]).scale((-1) ** j); acc = t if acc is None else acc.add(t); gj = gj.D()
    return acc


LG.check('(C5) control: the order-3 part of B alone is constant on (f0, f1)', const_of(btop(basis[0], basis[1]))[0], control=True)

# ------------------------------------------------------------------ (D) transversality
def minv(M):
    n = len(M); a = [list(map(Fr, r)) + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if a[r][c] != 0); a[c], a[piv] = a[piv], a[c]
        x = a[c][c]; a[c] = [y / x for y in a[c]]
        for r in range(n):
            if r != c and a[r][c] != 0:
                m = a[r][c]; a[r] = [y - m * z for y, z in zip(a[r], a[c])]
    return [row[n:] for row in a]


Bi = minv(Bm)
jets = [[y] for y in basis]
for y in jets:
    for _ in range(3): y.append(y[-1].D())
Tk = []
for k in range(4):
    acc = None
    for i in range(4):
        for j in range(4):
            if Bi[i][j]:
                t = basis[i].mul(jets[j][k]).scale(Bi[i][j]); acc = t if acc is None else acc.add(t)
    Tk.append(acc)
okT = all(t.is_zero() for t in Tk[:3])
t3 = Tk[3]; ok3 = all(c.is_zero() for c in t3.cs[1:]) and t3.cs[0].sub(Ser.poly({-1: 1})).is_zero()
LG.check('(D1) C3: T_0 = T_1 = T_2 = 0 and T_3 = 1/u (log-free), basis-independent form with Bmat^-1', okT and ok3,
         'precision u^%d' % min(c.prec for t in Tk for c in t.cs))

# ------------------------------------------------------------------ (G) PF5 = z^5 D PF4 - 1
def zc(X, Y):
    out = {}
    for a, x in X.items():
        for b, y in Y.items():
            t = op_compose(x, y); out[a + b] = t if a + b not in out else op_add(out[a + b], t)
    return out


def zadd(X, Y, s=1):
    out = dict(X)
    for k, y in Y.items():
        yy = [e.scale(s) for e in y]; out[k] = yy if k not in out else op_add(out[k], yy)
    return out


c = lambda x: Ser.poly({0: Fr(x)})
Mz = {1: [c(0), c(1)]}; Id = {0: [c(1)]}; qop = {0: [Ser.poly({0: Fr(1, 3125), 1: Fr(-1, 3125)})]}
M5 = Id
for _ in range(5): M5 = zc(Mz, M5)
prod = Id
for k in range(1, 6): prod = zc({1: [c(k), c(5)]}, prod)
PF5 = zadd(zadd(M5, Id, -1), zc(qop, prod), -1)
rhs = zadd({5: op_D_left(As)}, Id, -1)
diffz = zadd(PF5, rhs, -1)
LG.check('(G1) M^5 - 1 - q prod_(k=1..5)(5M + kz) = z^5 D PF4 - 1 (M = zD; z-graded operator identity, exact)',
         all(all(e.is_zero() for e in op) for op in diffz.values()))
HPFu = zadd({5: [c(0)] * 5 + [Ser.poly({1: 1})]}, Id, -1)          # u M^5 - 1  (M^5 = z^5 D^5)
for m, x in ((1, 3), (2, Fr(17, 5)), (3, Fr(9, 5)), (4, Fr(274, 625)), (5, Fr(24, 625))):
    HPFu = zadd(HPFu, {5: [c(0)] * (5 - m) + [ONE_M.scale(x)]}, -1)  # (1-u) x z^m M^(5-m) = (1-u) x z^5 D^(5-m)
LG.check('(G2) PF5 = u M^5 - 1 - (1-u)(3zM^4 + (17/5)z^2M^3 + (9/5)z^3M^2 + (274/625)z^4M + (24/625)z^5) (POLY s.1)',
         all(all(e.is_zero() for e in op) for op in zadd(PF5, HPFu, -1).values()))
LG.close()
