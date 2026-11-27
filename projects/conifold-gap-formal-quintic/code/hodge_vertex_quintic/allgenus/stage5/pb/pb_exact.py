# -*- coding: utf-8 -*-
"""pb_exact.py NMAX -- Stage-5 'pb' audit: exact-over-Q record of the gauge side, and Phi^can (exact) vs Phi^SP (mod p61*p89).

Prints the exact Frobenius data (gamma found from the log obstruction; f0, g, t_c, v), the canonical gauge's p_(5,b)(0) (UR s.4 prints
43397/807975000, ., ., 283/2835000, -12659/201993750), ph^can_5 (POLY s.1 prints it), then compares the EXACT Phi^can through w^NMAX
with Phi^SP computed mod 2^61-1 and 2^89-1 and combined by CRT: every exact coefficient must agree modulo p61*p89.
Log: logs/pb_exact_<NMAX>.log.
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C, pb_gauge as G, pb_sp as S
from fractions import Fraction as Fr


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    M = 4 * NMAX // 5 + 6; N = NMAX + 25; t0 = time.time(); L = ['pb_exact.py %d' % NMAX]
    C.set_field(0)
    K = C.Conifold(N + 6)
    L.append('  gamma = %s (found), c0 = %s; f0 = %s; g = %s; t_c = %s' % (K.gamma, K.c0,
             [str(K.f[0].coef(i)) for i in range(5)], [str(K.g.coef(i)) for i in range(5)], [str(K.tc.coef(i)) for i in range(5)]))
    V, v = C.model_v(8)
    L.append('  v(u) = %s' % [str(v.coef(i)) for i in range(5)])
    tau = C.tau_inv(N + 4).inv(); F = G.model_rows(NMAX)
    b, bt = G.systems(K, N)
    P = G.gauge(b, bt, NMAX, 'can')
    L.append('  canonical gauge p_(5,b)(0), b = 0..4: %s' % [str(P[5][j].coef(0)) for j in range(5)])
    can = G.phi_can(K, P, F, tau, NMAX, M)
    L.append('  ph^can_5 = %s' % [str(can[5].coef(k)) for k in range(5)])
    L.append('  exact gauge and row: %.1f s' % (time.time() - t0))
    # Phi^SP mod the two primes
    bern = C.bernoulli(2 * NMAX + 8); sp = {}
    for p in (S.P61, S.P89):
        C.set_field(p)
        xs = []; ys = []
        for Lp in range(2, 2 + 4 * NMAX // 5 + 8):
            u, ph = S.hatted_values(Lp, NMAX, p, bern); xs.append(u); ys.append(ph)
        sp[p] = [S.interpolate(xs, [y[n] for y in ys], p) for n in range(NMAX + 1)]
    C.set_field(0)
    PP = S.P61 * S.P89
    bad = []; cnt = 0
    for n in range(NMAX + 1):
        for k in range(M):
            ex = can[n].coef(k)
            r61 = ex.numerator * pow(ex.denominator, S.P61 - 2, S.P61) % S.P61
            r89 = ex.numerator * pow(ex.denominator, S.P89 - 2, S.P89) % S.P89
            a61 = sp[S.P61][n][k] if k < len(sp[S.P61][n]) else 0
            a89 = sp[S.P89][n][k] if k < len(sp[S.P89][n]) else 0
            cnt += 1
            if r61 != a61 or r89 != a89:
                bad.append((n, k))
    L.append('  exact Phi^can = Phi^SP mod 2^61-1 and mod 2^89-1 (i.e. mod their product) in all %d coefficients w^0..w^%d, '
             'u^0..u^%d: %s' % (cnt, NMAX, M - 1, 'PASS' if not bad else 'FAIL first %s' % (bad[0],)))
    L.append('  total %.1f s' % (time.time() - t0))
    L.append('ALL CHECKS PASSED' if not bad else 'SOME CHECK FAILED')
    print('\n'.join(L))
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    with open(os.path.join(HERE, 'logs', 'pb_exact_%d.log' % NMAX), 'w') as f:
        f.write('\n'.join(L) + '\n')
    print('log written to logs/pb_exact_%d.log' % NMAX)


if __name__ == '__main__':
    main()
