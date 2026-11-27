# -*- coding: utf-8 -*-
"""pb_glue.py -- Stage-5 'pb' audit: is the Phi^can of Theorem P(b) literally EVID's Phi^can?  mod 2^61-1 and 2^89-1.

EVID s.2.3-2.4 defines Phi^can by K_q = 0 (Lemma N); POLY s.6 defines it as the hatted first row of P_can Q^m with P_11(0,z) = 1.
pb_gauge builds the latter from the Theorem-U construction.  Compared here with the exact values printed in EVID s.4 and GI s.1:
  f_5 = 4537/3937500, f_10 = -177640111151/31007812500000, f_15 = -284328315069009259339/6226756347656250000000
  (top coefficients [u^4m] ph_5m), and the stored normalisation: [u^4 w^5](Phi^stored - Phi^can) = m5 = -43397/807975000.
Frames (0,0) and (1/3,-2/7).  Log: logs/pb_glue.log.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C, pb_gauge as G, pb_sp as S
from fractions import Fraction as Fr

EVID = {5: Fr(4537, 3937500), 10: Fr(-177640111151, 31007812500000),
        15: Fr(-284328315069009259339, 6226756347656250000000)}
M5 = Fr(-43397, 807975000)


def main():
    lines = ['pb_glue.py']; ok = True
    NMAX = 15; M = 4 * NMAX // 5 + 6; N = NMAX + 25
    for p in (S.P61, S.P89):
        C.set_field(p)
        K = C.Conifold(N + 6); tau = C.tau_inv(N + 4).inv(); F = G.model_rows(NMAX)
        for fr in ((0, 0), (Fr(1, 3), Fr(-2, 7))):
            K.set_frame(*fr)
            b, bt = G.systems(K, N)
            can = G.phi_can(K, G.gauge(b, bt, NMAX, 'can'), F, tau, NMAX, M)
            st = G.phi_can(K, G.gauge(b, bt, NMAX, 'stored'), F, tau, NMAX, M)
            fs = {n: can[n].coef(4 * n // 5) for n in EVID}
            good = all((fs[n] - C.red(EVID[n])) % p == 0 for n in EVID)
            top = all(can[n].ok_zero() is False or True for n in EVID)
            above = all(all(can[n].coef(k) % p == 0 for k in range(4 * n // 5 + 1, M)) for n in range(NMAX + 1))
            d = (st[5].coef(4) - can[5].coef(4)) % p
            gm5 = d == C.red(M5)
            lines.append('  p = %s, frame (%s,%s): f_5, f_10, f_15 of Phi^can = EVID values: %s; deg ph^can_n <= floor(4n/5), '
                         'n <= %d: %s; [u^4 w^5](Phi^stored - Phi^can) = m5: %s'
                         % ('2^61-1' if p == S.P61 else '2^89-1', fr[0], fr[1], 'PASS' if good else 'FAIL',
                            NMAX, 'PASS' if above else 'FAIL', 'PASS' if gm5 else 'FAIL'))
            ok = ok and good and above and gm5
    lines.append('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    print('\n'.join(lines))
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    with open(os.path.join(HERE, 'logs', 'pb_glue.log'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('log written to logs/pb_glue.log')


if __name__ == '__main__':
    main()
