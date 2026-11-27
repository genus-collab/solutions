# -*- coding: utf-8 -*-
"""pb_run_phi.py NMAX NPTS [quick] -- Stage-5 'pb' audit: Phi^SP = Phi^can through w^NMAX, mod 2^61-1 and 2^89-1.

Phi^can: the regular gauge of Theorem U built from scratch (pb_gauge), canonical normalisation P_11(0,z) = 1, hatted first row
times the model rows; u-series to O(u^M), M = floor(4 NMAX/5) + 12.  Phi^SP: the full-Stirling saddle series in the flat variable
at NPTS numerical saddles L = 2, 3, ... (pb_sp), turned into u-polynomials by interpolation in u = L^-5 (degree < NPTS).
Checks: (A) deg ph^SP_n <= floor(4n/5); (B) [u^k] ph^can_n = [u^k] ph^SP_n for all n <= NMAX, k < M (so ph^can_n is the same
polynomial); (C) frame (1/3,-2/7) gives the same Phi^can (one prime); (D) bonus: ph^SP_n(1) = h_n ((Gamma), all n <= NMAX).
Controls (each must FAIL its test): stored normalisation; free constant shifted at K = 10; SP without the Stirling term of
Gamma(1+5s); b_k -> -b_k; v = u (tau -> 1) must break the gauge construction.
Log: logs/pb_run_phi_<NMAX>_<NPTS>.log next to this script (relative paths only).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pb_core as C, pb_gauge as G, pb_sp as S
from fractions import Fraction as Fr

LOG = []


def log(s=''):
    print(s); LOG.append(s)


def first_mismatch(can, sp, nmax, M):
    for n in range(nmax + 1):
        for k in range(M):
            a = can[n].coef(k); b = sp[n][k] if k < len(sp[n]) else 0
            if (a - b) % C.MOD:
                return (n, k)
    return None


def run_prime(p, NMAX, NPTS, M, bern, frames, quick):
    C.set_field(p)
    N = NMAX + 25
    res = {}
    t0 = time.time()
    K = C.Conifold(N + 6)
    log('  p = %s: gamma = %s (found from the log obstruction), c0 = %s' % ('2^61-1' if p == S.P61 else '2^89-1',
        'ok 7/10' if (K.gamma * 10 - 7) % p == 0 else 'WRONG', K.c0))
    tau = C.tau_inv(N + 4).inv()
    F = G.model_rows(NMAX)
    log('  model rows F_b(0) = %s' % ([F[b][0] if F[b][0] < 10 else F[b][0] - p for b in range(5)]))
    # stationary phase at NPTS points
    xs = []; ys = []
    for L in range(2, 2 + NPTS):
        u, ph = S.hatted_values(L, NMAX, p, bern)
        xs.append(u); ys.append(ph)
    sp = [S.interpolate(xs, [y[n] for y in ys], p) for n in range(NMAX + 1)]
    degs = [max([i for i, c in enumerate(sp[n]) if c] or [0]) for n in range(NMAX + 1)]
    okA = all(degs[n] <= 4 * n // 5 for n in range(NMAX + 1))
    eq = sum(1 for n in range(NMAX + 1) if degs[n] == 4 * n // 5)
    log('  (A) deg ph^SP_n <= floor(4n/5) for all n <= %d: %s (equality for %d of %d n; n = 2 has degree %d)'
        % (NMAX, 'PASS' if okA else 'FAIL', eq, NMAX + 1, degs[2] if NMAX >= 2 else -1))
    log('      SP: %d saddles, %.1f s' % (NPTS, time.time() - t0))
    res['A'] = okA
    for (al, be) in frames:
        t1 = time.time()
        K.set_frame(al, be)
        b, bt = G.systems(K, N)
        P = G.gauge(b, bt, NMAX, 'can')
        can = G.phi_can(K, P, F, tau, NMAX, M)
        mm = first_mismatch(can, sp, NMAX, M)
        log('  (B) frame (%s,%s): Phi^can = Phi^SP in all [w^n u^k], n <= %d, k < %d (%d coefficients): %s  [gauge+row %.1f s]'
            % (al, be, NMAX, M, (NMAX + 1) * M, 'PASS' if mm is None else 'FAIL at %s' % (mm,), time.time() - t1))
        res['B', al, be] = mm is None
        if (al, be) == frames[0] and not quick:
            for name, norm, want in (('stored normalisation p_(5q,0)(0)=0', 'stored', (5, 4)),
                                     ('free constant shifted by 1/7 at K=10', {10: Fr(1, 7)}, (10, 8))):
                Pc = G.gauge(b, bt, NMAX, norm)
                cc = G.phi_can(K, Pc, F, tau, NMAX, M)
                m2 = first_mismatch(cc, sp, NMAX, M)
                log('  control [%s]: first mismatch %s (predicted %s): %s'
                    % (name, m2, want, 'PASS (fails as it must)' if m2 == want else 'CONTROL DID NOT BEHAVE'))
                res['ctl', name] = (m2 == want)
    # (D) (Gamma): ph^SP_n(1) = h_n
    h = S.h_qrr(NMAX, p, bern)
    bad = [n for n in range(NMAX + 1) if (sum(sp[n]) - h[n]) % p]
    log('  (D) (Gamma): ph^SP_n(1) = h_n for all n <= %d: %s' % (NMAX, 'PASS' if not bad else 'FAIL at %s' % bad[:5]))
    res['D'] = not bad
    # SP controls at low order
    for var, want in (('drop5', 1), ('flip', 1)):
        xs2 = []; ys2 = []
        for L in range(2, 2 + 8):
            u, ph = S.hatted_values(L, 6, p, bern, var)
            xs2.append(u); ys2.append(ph)
        sp2 = [S.interpolate(xs2, [y[n] for y in ys2], p) for n in range(7)]
        first = None
        for n in range(7):
            if any((sp2[n][k] - sp[n][k]) % p for k in range(min(len(sp2[n]), len(sp[n])))):
                first = n; break
        log('  control [SP variant %s]: first w-order differing from Phi^SP: %s (predicted %d): %s'
            % (var, first, want, 'PASS (fails as it must)' if first == want else 'CONTROL DID NOT BEHAVE'))
        res['ctl', var] = (first == want)
    # tau -> 1 control: the gauge construction must break
    K.set_frame(*frames[0])
    b, bt = G.systems(K, 30, tau_one=True)
    try:
        G.gauge(b, bt, 12, 'can'); broke = None
    except (AssertionError, ZeroDivisionError) as e:
        broke = str(e.args[0] if e.args else e)[:60]
    log('  control [v = u, tau -> 1]: gauge construction %s: %s'
        % ('breaks (%s)' % broke if broke else 'did NOT break', 'PASS (fails as it must)' if broke else 'CONTROL DID NOT BEHAVE'))
    res['ctl', 'tau1'] = broke is not None
    return res


def main():
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    NPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 4 * NMAX // 5 + 6
    quick = len(sys.argv) > 3 and sys.argv[3] == 'quick'
    M = 4 * NMAX // 5 + 12
    assert NPTS >= 4 * NMAX // 5 + 4, 'need at least 3 spare interpolation points'
    log('pb_run_phi.py %d %d%s   (M = %d u-coefficients per w-order)' % (NMAX, NPTS, ' quick' if quick else '', M))
    bern = C.bernoulli(2 * NMAX + 8)
    t0 = time.time(); allres = {}
    for p, frames in ((S.P61, [(0, 0), (Fr(1, 3), Fr(-2, 7))]), (S.P89, [(0, 0)])):
        allres[p] = run_prime(p, NMAX, NPTS, M, bern, frames, quick)
    ok = all(all(r.values()) for r in allres.values())
    log('total %.1f s' % (time.time() - t0))
    log('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    name = 'pb_run_phi_%d_%d%s.log' % (NMAX, NPTS, '_quick' if quick else '')
    with open(os.path.join(HERE, 'logs', name), 'w') as f:
        f.write('\n'.join(LOG) + '\n')
    print('log written to logs/' + name)


if __name__ == '__main__':
    main()
