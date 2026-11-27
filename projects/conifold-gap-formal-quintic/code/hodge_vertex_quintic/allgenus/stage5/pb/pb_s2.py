# -*- coding: utf-8 -*-
"""pb_s2.py N -- Stage-5 'pb' audit: (S2) of the poly referee, used in step 4 of Theorem P(b), mod 2^61-1 and 2^89-1.

(S2): the normalised formal SP at eta = 1 of eta^(c_b) exp(-(1/eta + eta^4/4)/t) d eta equals eps_b F_b(t), c = (0,-1,3,2,1),
eps = (+,-,-,+,+), F_b the model rows of UR s.2 Prop. S (pb_gauge.model_rows, a recursion that uses no integral).  PB_AUDIT s.2.4
proves it (the SPs solve the model system, are homogeneous, and the homogeneous solution with F_0(0) = 1 is unique); this checks
it through t^N.  Control: the all-plus signs must fail.  Log: logs/pb_s2_<N>.log.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C, pb_gauge as G, pb_sp as S
from fractions import Fraction as Fr


def sp_model(m, N, p):
    T = 2 * N + 1
    inv = lambda x: pow(x % p, p - 2, p)
    # F(1+x) - F(1) for F = -(1/eta + eta^4/4):  -[sum_(k>=2) (-1)^k x^k + (6x^2 + 4x^3 + x^4)/4]
    Fs = [0] * (T + 3)
    for k in range(2, T + 3):
        Fs[k] = -((-1) ** k)
    Fs[2] -= Fr(6, 4); Fs[3] -= 1; Fs[4] -= Fr(1, 4)
    Fs = [C.red(Fr(x)) if not isinstance(x, int) else x % p for x in Fs]
    a2 = Fs[2]                                     # F''(1)/2 = -5/2
    Q = [Fs[k + 2] * inv(a2) % p for k in range(T + 1)]
    Qm = S.ser_pow_unit(Q, Fr(-1, 2), T + 1, p); Qi = S.ser_pow_unit(Q, Fr(-1), T + 1, p)
    amp = S.binom_series(1, m, T + 1, p) if m >= 0 else S.binom_series(1, m, T + 1, p)
    out = []; Qp = Qm; A = (-2 * a2) % p
    for r in range(N + 1):
        if r:
            Qp = S.ser_mul(Qp, Qi, T + 1, p)
        df = 1
        for i in range(1, r + 1):
            df = df * (2 * i - 1) % p
        inner = sum(amp[i] * Qp[2 * r - i] for i in range(2 * r + 1)) % p
        out.append(df * pow(inv(A), r, p) % p * inner % p)
    return out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    lines = ['pb_s2.py %d' % N]; ok = True
    for p in (S.P61, S.P89):
        C.set_field(p)
        F = G.model_rows(N)
        good = True; ctl = False
        for b in range(5):
            sp = sp_model(G.CW[b], N, p)
            if any((sp[n] - G.EPS[b] * F[b][n]) % p for n in range(N + 1)):
                good = False
            if any((sp[n] - F[b][n]) % p for n in range(N + 1)):
                ctl = True
        lines.append('  p = %s: SP(eta^(c_b) e^(-(1/eta+eta^4/4)/t)) = eps_b F_b(t), all b, t^0..t^%d: %s; '
                     'control [all signs +]: %s' % ('2^61-1' if p == S.P61 else '2^89-1', N, 'PASS' if good else 'FAIL',
                     'PASS (fails as it must)' if ctl else 'CONTROL DID NOT BEHAVE'))
        ok = ok and good and ctl
    lines.append('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    print('\n'.join(lines))
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    with open(os.path.join(HERE, 'logs', 'pb_s2_%d.log' % N), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('log written to logs/pb_s2_%d.log' % N)


if __name__ == '__main__':
    main()
