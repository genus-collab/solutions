# -*- coding: utf-8 -*-
"""cy_frames.py N -- referee cyc_adversary: the operator identity of Theorem CYC,
    PF4 = M := (1/y0) D (1/Dt) D (A4 y0^2 Dt^2) D (1/Dt) D (1/y0),   t = y1/y0,
checked COEFFICIENTWISE as operators (direct composition; the critic checked M on solutions, pb on monomials), in every
frame the programme uses and in frames outside it, with controls; then deformations of PF4 (task item (d)).
Exact over Q, series through u^(N-1).  Usage: py -I -B cy_frames.py 40"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cy_core import Ser, Rel, op_D_left, op_adjoint, op_diff_zero
from cy_frob import pf4_poly, pf4_ser, frobenius
from cy_log import Log

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LG = Log('cy_frames_%d' % N)
LG('cy_frames: exact over Q, series through u^%d' % (N - 1))
U = Ser.poly({1: 1})


def M_op(y0, y1, a4):
    t = y1.mul(y0.inv()); Dt = t.D(); iy0 = y0.inv(); iDt = Dt.inv()
    E = [iy0]
    for fac in (iDt, a4.mul(y0).mul(y0).mul(Dt).mul(Dt), iDt, iy0):
        E = op_D_left(E); E = [fac.mul(e) for e in E]
    return E


def cvals(y0, y1):
    """valuations of c = (Dt, 1/(u y0^2 Dt^2), Dt, y0, y0); Theorem U needs (0, -1, 0, 0, 0)"""
    Dt = y1.mul(y0.inv()).D(); c1 = U.mul(y0).mul(y0).mul(Dt).mul(Dt).inv()
    return (Dt.val(), c1.val(), Dt.val(), y0.val(), y0.val())


def lin(*pairs):
    acc = Ser(0, [])
    for x, s in pairs: acc = acc.add(s.scale(x))
    return acc


def concomitant(Aser, fr, gr):
    fj = [fr]
    for _ in range(3): fj.append(fj[-1].D())
    B = None
    for k in range(1, 5):
        t = gr.mulser(Aser[k])
        for j in range(k):
            term = t.mul(fj[k - 1 - j]); term = term.neg() if j % 2 else term
            B = term if B is None else B.add(term); t = t.D()
    return B


def frame_obj(A, N):
    FR = frobenius(A, N)
    f0, f1, f2 = [Ser(0, x, N) for x in FR['f']]; g = Ser(0, FR['g'], N)
    return FR, f0, f1, f2, g


A = pf4_poly(); As = pf4_ser(A)
FR, f0, f1, f2, g = frame_obj(A, N)
prog = lambda al, be: lin((1, f0), (Fr(al), g), (Fr(be), f2))          # fq_conifold: varpi0 = f0 + alpha g + beta f2
crit = lambda a, b: lin((1, f0), (Fr(a), f1), (Fr(b), f2))             # cr5_cyc_run: varpi0 = f0 + a f1 + b f2
LG('(E0) the label "(1/3,-2/7)": programme varpi0 = f0 + g/3 - 2f2/7 has u^2 coefficient %s; the critic\'s f0 + f1/3 - 2f2/7 has %s'
   % (prog(Fr(1, 3), Fr(-2, 7)).coef(2), crit(Fr(1, 3), Fr(-2, 7)).coef(2)))
LG.check('(E0) the critic\'s frame "(1/3,-2/7)" is NOT the programme\'s frame (1/3,-2/7) (a label clash, harmless: CYC covers both)',
         prog(Fr(1, 3), Fr(-2, 7)).coef(2) != crit(Fr(1, 3), Fr(-2, 7)).coef(2))
cases = []
for al, be in ((0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5)), (1, 0), (0, 1), (Fr(2, 3), Fr(-1, 7)), (-2, 0)):
    cases.append(('programme frame (alpha, beta) = (%s, %s): y0 = f0 + alpha g + beta f2, y1 = g' % (al, be), prog(al, be), g, True, True))
for a, b in ((Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5))):
    cases.append(('critic-convention frame (%s, %s): y0 = f0 + a f1 + b f2, y1 = g' % (a, b), crit(a, b), g, True, True))
cases += [('rescaled: y0 = 2 (f0 + g/3 - 2f2/7), y1 = 3 g', prog(Fr(1, 3), Fr(-2, 7)).scale(2), g.scale(3), True, True),
          ('rescaled: y0 = -f0/2, y1 = -5 g', f0.scale(Fr(-1, 2)), g.scale(-5), True, True),
          ('same Lagrangian, y1 = g + f0/3 (t -> t + 1/3)', f0, g.add(f0.scale(Fr(1, 3))), True, True),
          ('varpi0(0) = 0: y0 = f1, y1 = g', f1, g, True, False),
          ('varpi0(0) = 0: y0 = f2, y1 = g', f2, g, True, False),
          ('varpi0(0) = 0: y0 = f1 - 2 f2, y1 = g', lin((1, f1), (-2, f2)), g, True, False),
          ('other Lagrangian pair: y0 = f1, y1 = f2', f1, f2, True, False),
          ('other Lagrangian pair: y0 = f2, y1 = f1', f2, f1, True, False),
          ('non-isotropic: y0 = f0, y1 = f1 (B = -7/5)', f0, f1, False, None),
          ('non-isotropic: y0 = f0, y1 = f2 (B = 2)', f0, f2, False, None),
          ('non-isotropic: y0 = f0 + f1, y1 = f2', lin((1, f0), (1, f1)), f2, False, None),
          ('wrong gamma 1/2: y0 = f0 + g/3 - 2f2/7, y1 = f1 + f2/2', prog(Fr(1, 3), Fr(-2, 7)), lin((1, f1), (Fr(1, 2), f2)), False, None),
          ('wrong gamma 7/10 + 1/1000: y0 = f0, y1 = f1 + gamma\' f2', f0, lin((1, f1), (Fr(701, 1000), f2)), False, None)]
for label, y0, y1, expect, unitshape in cases:
    ok, pr = op_diff_zero(M_op(y0, y1, U), As)
    cv = cvals(y0, y1)
    info = 'precision u^%d; val(c_a) = %s%s' % (pr, cv, '' if unitshape is None else
                                                  ('  [Theorem U shape]' if cv == (0, -1, 0, 0, 0) else '  [NOT the Theorem-U shape]'))
    LG.check('(E) M = PF4, %s' % label, ok, info=info, control=not expect)
    if expect and unitshape is not None:
        LG.check('(E) unit hypotheses of Theorem U (val c = (0,-1,0,0,0)) hold exactly when expected, %s' % label,
                 (cv == (0, -1, 0, 0, 0)) == unitshape)

# ------------------------------------------------------------------ deformations (task item (d))
def deform(label, pert, a4_poly=None, expect_sa=True, expect_cyc=True):
    Ap = pf4_poly(pert); Aps = pf4_ser(Ap)
    sa = op_diff_zero(op_adjoint(Aps), Aps)[0]
    LG.check('(F) %s: the perturbed operator is self-adjoint' % label, sa, control=not expect_sa)
    FRp, p0, p1, p2, gp = frame_obj(Ap, N)
    LG('    %s: obstructions %s, gamma = %s, residual obstruction of g: %s'
       % (label, [[str(x) for x in o] for o in FRp['ob']], FRp['gamma'], [str(x) for x in FRp['obg']]))
    a4 = Aps[4]
    ok, pr = op_diff_zero(M_op(p0, gp, a4), Aps)
    LG.check('(F) %s: CYC with its OWN Frobenius frame (y0 = f0, y1 = g of the perturbed operator)' % label, ok,
             info='precision u^%d' % pr, control=not expect_cyc)
    LG.check('(F) %s: CYC with the UNPERTURBED frame' % label, op_diff_zero(M_op(f0, g, U), Aps)[0], control=True)
    if a4_poly is not None:
        LG.check('(F) %s: CYC with the literal middle factor u y0^2 Dt^2 (instead of A4 y0^2 Dt^2)' % label,
                 op_diff_zero(M_op(p0, gp, U), Aps)[0], control=True)
    return FRp, Aps


deform('A0 -> A0 + (1/7)(1-u)', {0: {0: Fr(1, 7), 1: Fr(-1, 7)}})
deform('A2 -> A2 + u^2/3, A1 -> A1 + D(u^2)/3 (self-adjoint)', {2: {2: Fr(1, 3)}, 1: {1: Fr(-2, 3), 2: Fr(2, 3)}})
FR3, A3s = deform('A2 -> A2 + u^2/3 alone (not self-adjoint)', {2: {2: Fr(1, 3)}}, expect_sa=False, expect_cyc=False)
b3 = [Rel.of(Ser(0, x, N)) for x in FR3['f']] + [Rel([Ser(0, FR3['h'], N), Ser(0, FR3['g'], N)])]
Bnc = concomitant(A3s, b3[0], b3[1])
LG.check('(F) A2 + u^2/3 alone: the concomitant B(f0, f1) is a constant (C2(ii) needs C1)',
         all(c.is_zero() for c in Bnc.cs[1:]) and all(x == 0 for i, x in enumerate(Bnc.cs[0].c) if Bnc.cs[0].v + i), control=True)
a4n = Ser.poly({1: 1, 2: Fr(1, 4)})
d3 = a4n.D().scale(2).sub(As[3]); d1 = As[2].D().sub(a4n.D().D().D()).sub(As[1])
todict = lambda s: {s.v + i: x for i, x in enumerate(s.c) if x}
deform('A4 -> u + u^2/4 with A3 = 2 D A4, A1 = D A2 - D^3 A4', {4: {2: Fr(1, 4)}, 3: todict(d3), 1: todict(d1)}, a4_poly=True)
# exponents changed (not self-adjoint): A3 -> A3 + 1/2, indicial x(x-1)(x-2)(x-3/2): no log solution in K[log u]
Ap = pf4_poly({3: {0: Fr(1, 2)}}); FRx = frobenius(Ap, N)
obx = FRx['ob']
LG('    A3 -> A3 + 1/2: obstructions of f0, f1, f2 = %s' % [[str(x) for x in o] for o in obx])
def rank(M):
    a = [list(map(Fr, r)) for r in M]; rk = 0
    for c in range(len(a[0])):
        piv = next((r for r in range(rk, len(a)) if a[r][c] != 0), None)
        if piv is None: continue
        a[rk], a[piv] = a[piv], a[rk]
        for r in range(len(a)):
            if r != rk and a[r][c] != 0:
                m = a[r][c] / a[rk][c]; a[r] = [x - m * y for x, y in zip(a[r], a[rk])]
        rk += 1
    return rk


rank = rank(obx)
LG.check('(F) A3 -> A3 + 1/2 (exponents 0, 1, 3/2, 2): some g in span(f0, f1, f2) has a solution g log u + h (s.1 Frobenius structure)',
         rank < 3, info='obstruction map has rank %d on span(f0, f1, f2)' % rank, control=True)
LG.close()
