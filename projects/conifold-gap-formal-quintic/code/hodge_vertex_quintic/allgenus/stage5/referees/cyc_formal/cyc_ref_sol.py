# -*- coding: utf-8 -*-
"""cyc_ref_sol.py N -- referee cyc_formal: the solution-level steps of CYC_PROOF.md, exact over Q, series to u^N.

  S1  Frobenius structure: f_0, f_1, f_2 from the symbolic recursion, checked by applying PF4 as an operator; first terms.
  S2  log obstruction on the WHOLE 3-space span(f_0, f_1, f_2): the 3x3 matrix of [u^m] r(f_i), m = -3, -2, -1; rank 2, kernel
      (0, 1, 7/10) (so g is unique up to scale); [u^-3] r(f) = -2 f(0) (referee's hand formula); y_3 = g l + h solves PF4.
      Controls: gamma = 1/2 and g + f_0 are obstructed.
  S3  C2(ii),(iv): B on the basis (f_0, f_1, f_2, y_3) from the exact jet matrix: log-free, constant, alternating, the matrix of
      CYC s.3, det 4, Pfaffian 2.  Control: the concomitant without its k = 1 term is not constant.
  S4  Wronskian: det(D^a y_i) = C/u^2 (l-free), C^2 = 4 = det Bmat (det Bmat = (det Wr)^2 det P, det P = u^4).
  S5  C3: Sum_ij Binv_ij D^a y_i D^b y_j = (P^-1)_ab for ALL 16 entries (referee's hand P^-1).  Control: perturbed Binv.
  S6  C4: tau (l -> l+1) commutes with D on random elements of R; tau y_3 = y_3 + g; B(tau y_i, tau y_j) = B(y_i, y_j) computed
      directly; B(f, g) = 0 for random power-series solutions f; B(y_3, g) = 1.  Controls: l -> l + u; Lemma 0 control D l' = 1-u.
Usage: py -I -B cyc_ref_sol.py N     (prints the log; stored as logs/cyc_ref_sol_N.log)
"""
import os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U
from cyc_ref_ops import Op, REl
from cyc_ref_jet import concomitant_jet, det_inv, pf4, det_R, jets_of
from cyc_ref_frob import PF4, solve, logpart, basis, pairing, bmat, combo

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NF = [0]


def check(name, ok, control=False):
    if control:
        print(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        NF[0] += bool(ok)
    else:
        print(('PASS ' if ok else 'FAIL ') + name)
        NF[0] += (not ok)


t0 = time.time()
bas = basis(N)
f = bas['f']
check('S1 f_0, f_1, f_2 solve PF4 (operator applied) below u^%d' % (N - 3), all(PF4(REl.of(x)).zero_below(N - 3) for x in f))
check('S1 f_k = u^k + O(u^3): f_0 = 1 + 2/625 u^3, f_1 = u - 19/60 u^3, f_2 = u^2 + 37/30 u^3 (CYC s.8, fq_conifold, UR)',
      [[x.coef(e) for e in range(4)] for x in f] == [[1, 0, 0, Fr(2, 625)], [0, 1, 0, Fr(-19, 60)], [0, 0, 1, Fr(37, 30)]])
print('   obstruction matrix [u^m] r(f_i) (rows m = -3, -2, -1; columns i = 0, 1, 2): %s' % [[str(x) for x in r] for r in bas['Ob']])
check('S2 the l-part of PF4(f_i l) is PF4(f_i) = 0', all(logpart(x)[1].zero_below(N - 3) for x in f))
check('S2 obstruction rank 2 on span(f_0, f_1, f_2), kernel spanned by (0, 1, 7/10): g = f_1 + (7/10) f_2 unique up to scale',
      bas['rank'] == 2 and bas.get('kvec') == [0, 1, Fr(7, 10)])
check("S2 referee's formula [u^-3] r(f) = -2 f(0) on the basis", [bas['Ob'][0][i] for i in range(3)] == [-2, 0, 0])
check('S2 the h-equation has no obstruction for g; h_0 = h_1 = h_2 = 0', bas['h_ob'] == [0, 0, 0])
y3 = bas['y3']; g = bas['g']
check('S2 y_3 = g l + h solves PF4 below u^%d' % (N - 3), PF4(y3).zero_below(N - 3))
print('   g = u + %s u^2 + %s u^3 + ...;  h = %s u^3 + %s u^4 + ...' % (g.coef(2), g.coef(3), bas['h'].coef(3), bas['h'].coef(4)))
r_half = bas['r'][1] + bas['r'][2].sc(Fr(1, 2))
check('S2 gamma = 1/2 unobstructed', all(r_half.coef(m) == 0 for m in (-3, -2, -1)), control=True)
r_f0 = bas['r'][0] + bas['r'][1] + bas['r'][2].sc(Fr(7, 10))
check('S2 g + f_0 unobstructed', all(r_f0.coef(m) == 0 for m in (-3, -2, -1)), control=True)

# S3 the pairing B
Y = bas['Y']
B, lf, cst, prec = bmat(Y)
print('   Bmat (basis f_0, f_1, f_2, y_3) = %s   (constancy checked below u^%d)' % ([[str(x) for x in r] for r in B], prec))
check('S3 B(y_i, y_j) log-free', lf)
check('S3 B(y_i, y_j) constant (every u^e, e != 0, vanishes)', cst)
check('S3 B alternating', all(B[i][j] == -B[j][i] for i in range(4) for j in range(4)))
CYCB = [[0, Fr(-7, 5), 2, Fr(-9, 10)], [Fr(7, 5), 0, 0, -1], [-2, 0, 0, 0], [Fr(9, 10), 1, 0, 0]]
check('S3 Bmat equals the matrix of CYC s.3', B == CYCB)
det, Binv = det_inv(B)
check('S3 det Bmat = 4 and Pf(Bmat) = 2 (only B(f_0,f_2) = 2 and B(f_1,y_3) = -1 enter the Pfaffian)', det == 4 and pf4(B) == 2)
gl = [REl.of(f[0]), REl.of(g), REl.of(f[2])]
B3 = bmat(gl)[0]
check('S3 fq_conifold basis (f_0, g, f_2): column g vanishes ("pairing vacuous")', all(B3[a][1] == 0 for a in range(3)))
Pk1 = concomitant_jet(PF4.A, skip=(1,))
val, lf1, c1, _ = pairing(Y[0], Y[1], Pk1)
check('S3 concomitant without its k = 1 term constant on (f_0, f_1)', c1, control=True)

# S4 Wronskian
W = [[jets_of(y)[a] for y in Y] for a in range(4)]
dW = det_R(W)
s0 = dW.F[0]
lf_W = all(s.zero_below(s.N) for k, s in dW.F.items() if k)
Cw = s0.coef(-2)
check('S4 det Wr is l-free and equals C u^-2 (no other power below u^%d), C^2 = det Bmat = 4' % s0.N,
      lf_W and all(x == 0 for i, x in enumerate(s0.c) if s0.v + i != -2) and Cw * Cw == det)
print('   det Wr = %s u^-2' % Cw)

# S5 C3 transversality, all 16 entries
Z = Ser(0, [])
iu, iu2 = Ser(-1, [Fr(1)]), Ser(-2, [Fr(1), Fr(-1)])
Pinv = [[Z, Z, Z, iu], [Z, Z, iu.sc(-1), iu2.sc(-1)], [Z, iu, Z, iu2.sc(Fr(2, 5))], [iu.sc(-1), iu2, iu2.sc(Fr(-2, 5)), Z]]


def trans_ok(Bi, rows=range(4)):
    ok, pmin = True, 10 ** 9
    for a in rows:
        for b in range(4):
            S = None
            for i in range(4):
                for j in range(4):
                    if Bi[i][j]:
                        t = (W[a][i] * W[b][j]).sc(Bi[i][j]); S = t if S is None else S + t
            d = S - REl.of(Pinv[a][b])
            p = min(s.N for s in d.F.values()); pmin = min(pmin, p)
            ok = ok and d.zero_below(p)
    return ok, pmin


ok, pmin = trans_ok(Binv)
check('S5 C3: Wr Binv Wr^T = P^-1 entrywise (16 entries; T_k = 0 for k <= 2, T_3 = 1/u) below u^%d' % pmin, ok)
Bp = [r[:] for r in Binv]; Bp[0][1] += 1; Bp[1][0] -= 1
check('S5 the same with Binv perturbed by E_01 - E_10', trans_ok(Bp, rows=[0])[0], control=True)

# S6 C4 formal Picard-Lefschetz
random.seed(5)
def rnd_ser():
    return Ser.poly({e: Fr(random.randint(-9, 9), random.randint(1, 7)) for e in range(-3, 4)})
ok_c = True
for _ in range(6):
    x = REl({k: rnd_ser() for k in range(4)})
    d = x.D().shift_l(ONE) - x.shift_l(ONE).D()
    ok_c = ok_c and all(not s.c for s in d.F.values())
check('S6 tau D = D tau on 6 random elements of R (l-degree 3, exact)', ok_c)
lel = REl({1: ONE})
check('S6 l -> l + u commutes with D (on l)', all(not s.c for s in (lel.D().shift_l(U) - lel.shift_l(U).D()).F.values()), control=True)
ty3 = y3.shift_l(ONE)
dif = ty3 - y3 - REl.of(g)
check('S6 tau y_3 = y_3 + g, tau f_i = f_i', all(not s.c for s in dif.F.values()))
check('S6 l -> l + u maps y_3 to a solution', PF4(y3.shift_l(U)).zero_below(N - 5), control=True)
tY = [y.shift_l(ONE) for y in Y]
Bt = bmat(tY)[0]
check('S6 B(tau y_i, tau y_j) = B(y_i, y_j), computed directly on the transformed elements', Bt == B)
ok_pl = True
for _ in range(5):
    cf = [Fr(random.randint(-9, 9), random.randint(1, 9)) for _ in range(3)] + [0]
    ok_pl = ok_pl and pairing(combo(bas, cf), REl.of(g))[0] == 0
check('S6 B(f, g) = 0 for 5 random power-series solutions f', ok_pl)
check('S6 B(y_3, g) = 1 (CYC s.5)', pairing(y3, REl.of(g))[0] == 1)
dl2 = Ser(0, [Fr(1), Fr(-1)])                                  # control derivation D l' = 1 - u = D(-u)
cst2 = REl({1: ONE, 0: U}, dl=dl2).D()
check("S6 Lemma 0's conclusion for D l' = 1 - u = D(-u) instead of 1 - 1/u (D l' in D(K)): l' + u is not a constant",
      any(s.c for s in cst2.F.values()), control=True)
print('CYC_REF_SOL %d: %s   (%.1f s)' % (N, 'ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
