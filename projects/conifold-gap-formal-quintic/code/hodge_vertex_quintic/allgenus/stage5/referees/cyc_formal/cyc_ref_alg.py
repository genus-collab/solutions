# -*- coding: utf-8 -*-
"""cyc_ref_alg.py -- referee cyc_formal: the EXACT (no truncation) algebraic steps of CYC_PROOF.md.

Checks (each separately; every identity is an identity of operators or jets with polynomial coefficients, so exact = proof
of that identity for this PF4):
  A1  PF4 = D^4 - 5q prod_(k=1..4)(5D+k), q = (1-u)/3125, D = -(1-u) d/du, has exactly the A_k of CYC s.0 / fq_conifold;
      UR s.1's u-form of PF5 (coefficients 3, 17/5, 9/5, 274/625, 24/625).
  A2  Corollary: D^5 - q o prod_(k=1..5)(5D+k) = D o PF4, i.e. PF5 = z^5 D PF4 - 1.  Controls: q on the right; k = 0..4.
  A3  recursion polynomials c_s(n) of PF4 u^n = Sum_s c_s(n) u^(n+s), symbolic in n: c_-3(n) = n(n-1)^2(n-2); no negative powers
      from u^0, u^1, u^2.
  A4  C1: PF4* = PF4; the two self-adjointness conditions A_3 = 2 D A_4, A_1 = D A_2 - D^3 A_4; P(-D-1) = P(D); D o q = q o (D+1).
      Controls: PF4 + uD; A_3 -> -3(1-u); q on the right.
  A5  C2(i) telescoping identity g L f - f L* g = D P(f, g) as a jet identity (PF4 and a random operator); control: no signs.
      C2 alternation P(f,g) + P(g,f) = 0 identically for PF4; control PF4 + uD.
  A6  C2(iii) the matrix P_ab: zero for a + b >= 4, P_(3-j,j) = (-1)^j u, equal to the referee's hand computation; det P = u^4;
      the referee's hand formula for P^-1 satisfies P P^-1 = 1 (so (P^-1)_ab = 0 for a + b <= 2, (P^-1)_03 = 1/u).
Usage: py -I -B cyc_ref_alg.py     (prints the log; the referee stores it as logs/cyc_ref_alg.log)
"""
import os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U
from cyc_ref_ops import Op, REl
from cyc_ref_jet import (poly, pf4_from_large_radius, ops_equal_exact, jet_add, jet_neg, jet_D, jet_iszero, jet_swap,
                         concomitant_jet, lagrange_lhs, sm_apply, pn_mul, pn_add, pn_eval, det_R)

NF = [0]


def check(name, ok, control=False):
    if control:
        print(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        NF[0] += bool(ok)
    else:
        print(('PASS ' if ok else 'FAIL ') + name)
        NF[0] += (not ok)


t0 = time.time()
v = poly({0: 1, 1: -1})                                      # 1 - u
A_cyc = {4: poly({1: 1}), 3: v.sc(-2), 2: v.sc(Fr(-7, 5)), 1: v.sc(Fr(-2, 5)), 0: v.sc(Fr(-24, 625))}
PF4 = pf4_from_large_radius()
check('A1 PF4 = D^4 - 5q prod(5D+k) has the A_k of CYC s.0 / fq_conifold (u, -2(1-u), -7/5(1-u), -2/5(1-u), -24/625(1-u))',
      ops_equal_exact(PF4, Op(A_cyc)))
check('A1 D u = -(1-u) = q d/dq (1 - 3125 q) (D is q d/dq)', not (U.D() - v.sc(-1)).c)
e = [Fr(1)]                                                  # elementary symmetric polynomials of 1..5
for k in range(1, 6):
    e = [(e[i] if i < len(e) else 0) + (k * e[i - 1] if i >= 1 else 0) for i in range(len(e) + 1)]
URc = [Fr(5) ** (5 - j) * e[j] / 3125 for j in range(6)]   # q prod(5M + kz) = (1-u) Sum_j URc_j z^j M^(5-j)
check('A1 UR s.1 u-form of PF5: u M^5 - (1-u)[3 zM^4 + 17/5 z^2M^3 + 9/5 z^3M^2 + 274/625 z^4M + 24/625 z^5] - 1',
      URc == [1, 3, Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625)])

# A2 Corollary operator identity
q = poly({0: Fr(1, 3125), 1: Fr(-1, 3125)})
def prod5(ks):
    P = Op.fn(ONE)
    for k in ks:
        P = P * Op({1: poly({0: 5}), 0: poly({0: k})})
    return P
lhs = Op.Dk(5) - Op.fn(q) * prod5(range(1, 6))
check('A2 D^5 - q o prod_(k=1..5)(5D+k) = D o PF4 (so PF5 = z^5 D PF4 - 1, z a constant), exact', ops_equal_exact(lhs, Op.Dk(1) * PF4))
check('A2 same with q on the right', ops_equal_exact(Op.Dk(5) - prod5(range(1, 6)) * Op.fn(q), Op.Dk(1) * PF4), control=True)
check('A2 same with k = 0..4', ops_equal_exact(Op.Dk(5) - Op.fn(q) * prod5(range(0, 5)), Op.Dk(1) * PF4), control=True)

# A3 symbolic recursion
c = sm_apply(PF4)
for s in sorted(c):
    print('   c_%d(n) = %s   (coefficients of n^0, n^1, ...)' % (s, [str(x) for x in c[s]]))
target = pn_mul(pn_mul([Fr(0), Fr(1)], pn_mul([Fr(-1), Fr(1)], [Fr(-1), Fr(1)])), [Fr(-2), Fr(1)])
check('A3 [u^(n-3)] PF4 u^n = n(n-1)^2(n-2) as a polynomial in n', c.get(-3) == target)
check('A3 PF4 u^n lies in span(u^(n-3), .., u^(n+1)) for every n', set(c) <= set(range(-3, 2)))
check('A3 no negative powers from u^0, u^1, u^2 (c_s(n) = 0 when n + s < 0)',
      all(pn_eval(c.get(s, []), n) == 0 for n in range(3) for s in range(-3, 0) if n + s < 0))
check('A3 pivot n(n-1)^2(n-2) != 0 for 3 <= n <= 10^4', all(pn_eval(target, n) != 0 for n in range(3, 10001)))

# A4 C1
check('A4 C1: PF4* = PF4 exactly (formal adjoint Sum (-D)^k o A_k)', ops_equal_exact(PF4.adjoint(), PF4))
A = PF4.A
check('A4 A_3 = 2 D A_4 and A_1 = D A_2 - D^3 A_4 (the two conditions for an order-4 self-adjoint operator)',
      not (A[3] - A[4].D().sc(2)).c and not (A[1] - A[2].D() + A[4].D().D().D()).c)
Px = [Fr(1)]
for k in range(1, 5):
    Px = pn_mul(Px, [Fr(k, 5), Fr(1)])
Pm = [Fr(1)]
for k in range(1, 5):
    Pm = pn_mul(Pm, [Fr(k, 5) - 1, Fr(-1)])                  # P(-x-1)
check('A4 P(-D-1) = P(D), P(D) = prod_(k=1..4)(D + k/5)', Pm == Px)
check('A4 D o q = q o (D + 1)', ops_equal_exact(Op.Dk(1) * Op.fn(q), Op.fn(q) * Op({1: ONE, 0: ONE})))
check('A4 PF4 + u D self-adjoint', ops_equal_exact((PF4 + Op({1: U})).adjoint(), PF4 + Op({1: U})), control=True)
Lb = Op(dict(A)); Lb.A[3] = v.sc(-3)
check('A4 PF4 with A_3 -> -3(1-u) self-adjoint', ops_equal_exact(Lb.adjoint(), Lb), control=True)
Lr = pf4_from_large_radius('right')
check('A4 D^4 - 5 prod(5D+k) o q (q on the right) self-adjoint', ops_equal_exact(Lr.adjoint(), Lr), control=True)

# A5 C2(i) and alternation, as jet identities
def tele_ok(L, signs=True):
    J = jet_add(lagrange_lhs(L.A, L.adjoint().A), jet_neg(jet_D(concomitant_jet(L.A, signs))))
    return jet_iszero(J)
random.seed(20261114)
Lran = Op({k: poly({i: Fr(random.randint(-9, 9), random.randint(1, 9)) for i in range(3)}) for k in range(5)})
check('A5 C2(i) g L f - f L* g = D P(f,g), L = PF4 (jet identity, exact)', tele_ok(PF4))
check('A5 C2(i) the same identity for a random order-4 operator with quadratic coefficients', tele_ok(Lran))
check('A5 C2(i) without the signs (-1)^j', tele_ok(PF4, signs=False), control=True)
P = concomitant_jet(PF4.A)
check('A5 C2(ii) alternation: P(f,g) + P(g,f) = 0 identically for PF4', jet_iszero(jet_add(P, jet_swap(P))))
L1 = PF4 + Op({1: U}); P1 = concomitant_jet(L1.A)
check('A5 alternation for PF4 + uD', jet_iszero(jet_add(P1, jet_swap(P1))), control=True)
Q1 = jet_add(P1, jet_swap(P1))
check('A5 consistency: for PF4 + uD the sum Q = P(f,g) + P(g,f) is non-zero and so is D Q', not jet_iszero(Q1) and not jet_iszero(jet_D(Q1)))

# A6 C2(iii)
def Pent(a, b):
    s = P.get((a, b)); return s if s is not None else Ser(0, [])
print('   P_ab (a = order on f, b = order on g):')
for a in range(4):
    print('     ', ['%s' % {e: str(x) for e, x in zip(range(Pent(a, b).v, 99), Pent(a, b).c)} for b in range(4)])
check('A6 P_ab = 0 for a + b >= 4', all(not s.c for (a, b), s in P.items() if a + b >= 4))
check('A6 P_(3-j,j) = (-1)^j u', all(not (Pent(3 - j, j) - U.sc((-1) ** j)).c for j in range(4)))
Z = Ser(0, [])
hand = [[Z, v.sc(Fr(2, 5)), v, U.sc(-1)], [v.sc(Fr(-2, 5)), Z, U, Z], [v.sc(-1), U.sc(-1), Z, Z], [U, Z, Z, Z]]
check("A6 P equals the referee's hand computation [[0, 2/5(1-u), 1-u, -u], [-2/5(1-u), 0, u, 0], [-(1-u), -u, 0, 0], [u, 0, 0, 0]]",
      all(not (Pent(a, b) - hand[a][b]).c for a in range(4) for b in range(4)))
dP = det_R([[REl.of(Pent(a, b)) for b in range(4)] for a in range(4)])
check('A6 det P = u^4', set(dP.F) <= {0} and not (dP.F[0] - Ser(4, [Fr(1)])).c)
iu, iu2 = Ser(-1, [Fr(1)]), Ser(-2, [Fr(1), Fr(-1)])        # 1/u, (1-u)/u^2
Pinv = [[Z, Z, Z, iu], [Z, Z, iu.sc(-1), iu2.sc(-1)], [Z, iu, Z, iu2.sc(Fr(2, 5))], [iu.sc(-1), iu2, iu2.sc(Fr(-2, 5)), Z]]
prod_ok = True
for a in range(4):
    for b in range(4):
        s = Z
        for k in range(4):
            s = s + Pent(a, k) * Pinv[k][b]
        prod_ok = prod_ok and not (s - (ONE if a == b else Z)).c
check("A6 the referee's P^-1 (zero for a + b <= 2, (P^-1)_03 = 1/u) satisfies P P^-1 = 1 exactly", prod_ok)
print('CYC_REF_ALG: %s   (%.1f s)' % ('ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
