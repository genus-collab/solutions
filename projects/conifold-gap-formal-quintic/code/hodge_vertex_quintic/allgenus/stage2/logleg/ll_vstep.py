# -*- coding: utf-8 -*-
"""ll_vstep.py -- pairing-level checks of the induction step of (V) (LOGLEG_PROOF.md s.3) on x = W = p^4/4 + u/p (u = 1).
   py -I -B ll_vstep.py <out.log> <field: 89 | 61 | exact> <h,a,k> [<h,a,k> ..]
P_k(m) := <p^m_1 .. p^m_a, (-log p)^k>_h  (a >= 1 L-legs, k log legs, n = a + k),  e(m) = 2 - 2h - n + Sum m/5 (homogeneity exponent).
(S) step (Theorem V, eq. (V1)):  P_{k+1}(m) = e(m) P_k(m) + Sum_i m_i Sum_{q >= 0, m_i - 5 - 5q >= 1} P_k(m with m_i -> m_i - 5 - 5q)
    [d_v at u = 1 is multiplication by e;  (xi d Phi)_+ = (m p^m/(p^5 - u))_+ = m Sum_q u^q p^(m-5-5q) for Phi = p^m]
(D) dilaton cross-check (TR dilaton equation, EO Thm 4.7, independent of the variational formula; Phi_ydx = p^5/5 - u log p):
    P_{k+1}(m) = (2 - 2h - n) P_k(m) - (1/5) <p^m_1 .. p^m_a, p^5, (-log p)^k>_h   [sign fixed by the computation, see log]
(X) dictionary: <-log p, d xif_{a,d}> = delta_{a,1} delta_{d,0}  (tr_xi.XiBasis, flat basis H^a, descendant d): -log p is the ancestor tau_0(H).
Checked on ALL tuples with Sum m <= 5(2h - 2 + n + 1) + 5 (singular and regular range).  Controls (must FAIL): (S) with the xi-term of
the opposite sign, and without it.
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from ll_core import QzField, make_fp, P89, P61, setup, Legs, pair_table, contract_last, rat, compositions
from ll_falsifier import Log


def step_rhs(F, P, m, e, sgn=1):
    r = P.get(m, F.zero) * F.of(e)
    if sgn == 0: return r
    for i, mi in enumerate(m):
        q = 0
        while mi - 5 - 5 * q >= 1:
            mm = m[:i] + (mi - 5 - 5 * q,) + m[i + 1:]
            x = P.get(mm)
            if x is not None: r = r + x * mi * sgn
            q += 1
    return r


def main(out, fld, triples):
    log = Log(out); log('ll_vstep.py %s %s  (allgenus/stage2/logleg, %s)' % (fld, ' '.join('%d,%d,%d' % t for t in triples), time.strftime('%Y-%m-%d %H:%M')))
    F = {'89': make_fp(P89), '61': make_fp(P61), 'exact': QzField()}[fld]
    kmx = max(6 * h - 5 + 2 * (a + k + 2) for h, a, k in triples); cur, tr = setup(F, 2 * kmx + 8); c = cur.c; L = Legs(F, c)
    # (X)
    from tr_xi import XiBasis
    XB = XiBasis(cur, 8); bad = 0
    for a in range(5):
        for d in range(9):
            v = contract_last({((key),): x for key, x in XB.vec[('f', a, d)].items()}, [L.logleg]).get((), F.zero)
            want = F.one if (a, d) == (1, 0) else F.zero
            if not (v - want).iszero(): bad += 1
    log.check(bad == 0, '(X) <-log p, d xif_{a,d}> = delta_{(a,d),(1,0)} for a = 0..4, d = 0..8: the log leg is the ancestor insertion tau_0(H)')
    for h, a, k in triples:
        n = a + k; M = 5 * (2 * h - 2 + n + 1) + 5; t0 = time.time()
        Pk = pair_table(F, c, tr.omega(h, n), n, k, M, L)
        Pk1 = pair_table(F, c, tr.omega(h, n + 1), n + 1, k + 1, M, L)
        Pd = pair_table(F, c, tr.omega(h, n + 1), n + 1, k, M + 5, L)
        tot = ok = okD = c1 = c2 = 0; nz = 0
        keys = [m for s_ in range(5, M + 1, 5) for m in compositions(s_, a)]      # ALL weight tuples in range
        for m in keys:
            e = Fr(2 - 2 * h - n) + Fr(sum(m), 5)
            lhs = Pk1.get(m, F.zero); tot += 1; nz += (not lhs.iszero())
            if (lhs - step_rhs(F, Pk, m, e)).iszero(): ok += 1
            if not (lhs - step_rhs(F, Pk, m, e, -1)).iszero(): c1 += 1
            if not (lhs - step_rhs(F, Pk, m, e, 0)).iszero(): c2 += 1
            dil = Pk.get(m, F.zero) * (2 - 2 * h - n) - Pd.get(m + (5,), F.zero) * F.of(Fr(1, 5))
            if (lhs - dil).iszero(): okD += 1
        log.check(ok == tot and c1 > 0 and c2 > 0, '(S) h=%d a=%d k=%d -> k+1: step identity holds on %d/%d weight tuples (Sum m <= %d, %d non-zero LHS); '
                  'controls fail on %d (opposite xi-sign) and %d (no xi-term) tuples (%.1f s)' % (h, a, k, ok, tot, M, nz, c1, c2, time.time() - t0))
        log.check(okD == tot, '(D) h=%d a=%d k=%d: dilaton form P_{k+1} = (2-2h-n) P_k - (1/5)<.., p^5, ..> holds on %d/%d tuples' % (h, a, k, okD, tot))
    log.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], [tuple(int(x) for x in s.split(',')) for s in sys.argv[3:]])
