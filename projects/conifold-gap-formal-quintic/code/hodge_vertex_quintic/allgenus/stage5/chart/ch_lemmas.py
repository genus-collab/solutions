# -*- coding: utf-8 -*-
"""ch_lemmas.py -- Stage 5 / chart: exact tests of the formal lemmas of CHART_PROOF s.2 on seeded random data over Q.

(S1)  N[Phi, z (d_i Phi/z) h + z d_i h] = 0 (formal integration by parts);
(L1a) N[Phi o chi, (A o chi) det D chi; 0] = det D chi(0) N[Phi, A; 0] for a random polynomial change of variables (n = 2);
(L1b) Fubini: N[Phi, A; 0] = N[phi, beta B; 0] with v_c(w), phi(w) = Phi(v_c(w), w), B = inner series (computed over the ring
      Q[[w]]/(w^D)), beta = (det H_vv(w)/det H_vv(0))^(-1/2); inner dimension 1 and 2, outer dimension 1;
(L2)  phase sigma G(x) + f(x), amplitude free of sigma: N = A(X_0, z) EXACTLY;
(L1c)(iii) one variable: the SP of exp((5 log y - y + q y^5)/z) dy at the q-dependent critical point, re-expanded in q, equals
      Sum_d (q/z)^d/d! N[5 log y - y, y^(5d); 5] coefficientwise (q^d z^n, n >= -d); also against (L3): 5^(5d) prod (1 + iz/5).
Controls (must FAIL): (L1a) without the Jacobian; (L1b) without beta; (L2) with a sigma-dependent amplitude and with a sigma^2 term.
Usage: py -I -B ch_lemmas.py [N] [SEEDS]   (default 4 3; log logs/ch_lemmas_<N>_<SEEDS>.log)
"""
import os, sys, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from ch_core import QF, RS
from ch_wick import gadd, gmul, gscale, gconst, guni, gtrunc, sp_series, sp1, matinv, det, hessian

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
LOGREL = 'logs/ch_lemmas_%d_%d.log' % (N, SEEDS)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8'); FAILS = []; T0 = time.time()
Q = QF()


def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()


def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)


def rf(rng): return Fr(rng.randint(-6, 6), rng.randint(1, 4))


def rpoly(rng, n, dmin, dmax):
    P = {}
    for d in range(dmin, dmax + 1):
        for _ in range(3 * n):
            e = [0] * n
            for _ in range(d): e[rng.randrange(n)] += 1
            P.setdefault(d, {})[tuple(e)] = rf(rng)
    return P


def rphase(rng, n, dmax):
    """random phase with zero linear part and a random invertible Hessian (all leading minors non-zero)"""
    while True:
        P = rpoly(rng, n, 3, dmax); H = [[rf(rng) for _ in range(n)] for _ in range(n)]
        H = [[(H[i][k] + H[k][i]) / 2 for k in range(n)] for i in range(n)]
        if all(det(Q, [r[:m] for r in H[:m]]) != 0 for m in range(1, n + 1)): break
    for i in range(n):
        for k in range(i, n):
            e = [0] * n; e[i] += 1; e[k] += 1
            P.setdefault(2, {})[tuple(e)] = H[i][i] / 2 if i == k else H[i][k]
    return P


def dpoly(P, i):
    R = {}
    for d, m in P.items():
        for e, c in m.items():
            if e[i]:
                f = list(e); f[i] -= 1; R.setdefault(d - 1, {})[tuple(f)] = R.get(d - 1, {}).get(tuple(f), 0) + c * e[i]
    return R


def compose(F, P, chis, D, n=2):
    """P(chi_1, .., chi_n) truncated at degree D, chi_i graded polynomials (n variables) without constant term"""
    pw = []
    for ch in chis:
        L = [gconst(F, n, F.one)]
        for _ in range(D): L.append(gmul(F, L[-1], ch, D))
        pw.append(L)
    R = {}
    for d, m in P.items():
        if d > D: continue
        for e, c in m.items():
            T = pw[0][e[0]]
            for i in range(1, len(e)): T = gmul(F, T, pw[i][e[i]], D)
            R = gadd(F, R, T, F.of(c))
    return R


def rpow(R, a, alpha):
    c = [Fr(0)] * R.D; c[0] = Fr(1)
    for m in range(1, R.D):
        c[m] = sum(((alpha + 1) * i - m) * a[i] * c[m - i] for i in range(1, m + 1)) / m
    return tuple(c)


def test_S1(seed):
    rng = random.Random(seed); n = 2; Phi = rphase(rng, n, 2 * N + 2); h0 = rpoly(rng, n, 0, 2 * N); h1 = rpoly(rng, n, 0, 2 * N)
    ok = True
    for i in range(n):
        dP = dpoly(Phi, i); A0 = gmul(Q, dP, h0, 2 * N); A1 = gadd(Q, dpoly(h0, i), gmul(Q, dP, h1, 2 * N)); A2 = dpoly(h1, i)
        s, H = sp_series(Q, Phi, [A0, A1, A2], N, n); ok = ok and all(x == 0 for x in s)
    return ok


def test_L1a(seed, nojac=False):
    rng = random.Random(seed); n = 2; D = 2 * N + 2
    Phi = rphase(rng, n, D); amps = [rpoly(rng, n, 0, 2 * (N - m)) for m in range(2)]
    amps[0].setdefault(0, {})[(0, 0)] = Fr(1)
    while True:
        M = [[rf(rng) for _ in range(n)] for _ in range(n)]
        if det(Q, M) != 0: break
    chis = [gadd(Q, {1: {(1, 0): M[i][0], (0, 1): M[i][1]}}, rpoly(rng, n, 2, 3)) for i in range(n)]
    J = gadd(Q, gmul(Q, dpoly(chis[0], 0), dpoly(chis[1], 1), D), gmul(Q, dpoly(chis[0], 1), dpoly(chis[1], 0), D), Fr(-1))
    if nojac: J = gconst(Q, n, det(Q, M))
    Pc = compose(Q, Phi, chis, D); ac = [gmul(Q, compose(Q, A, chis, 2 * (N - m)), J, 2 * (N - m)) for m, A in enumerate(amps)]
    s1, H1 = sp_series(Q, Pc, ac, N, n); s0, H0 = sp_series(Q, Phi, amps, N, n)
    return s1 == [det(Q, M) * x for x in s0]


def test_L1b(seed, n1, nobeta=False):
    rng = random.Random(seed); n = n1 + 1; D = 2 * N + 3; R = RS(D); wg = tuple(Fr(1) if i == 1 else Fr(0) for i in range(D))
    Phi = rphase(rng, n, 2 * N + 2); amps = [rpoly(rng, n, 0, 2 * (N - m)) for m in range(2)]
    amps[0].setdefault(0, {})[(0,) * n] = Fr(1)
    s_tot, Ht = sp_series(Q, Phi, amps, N, n)
    Hvv = [[Ht[i][k] for k in range(n1)] for i in range(n1)]; Hi = matinv(Q, Hvv); grads = [dpoly(Phi, i) for i in range(n1)]
    vc = [R.zero] * n1

    def ev(P, vals):
        tot = R.zero
        for d, m in P.items():
            for e, c in m.items():
                t = R.of(c)
                for i, k in enumerate(e):
                    for _ in range(k): t = R.mul(t, vals[i])
                tot = R.add(tot, t)
        return tot
    for _ in range(D):
        g = [ev(gr, vc + [wg]) for gr in grads]
        vc = [R.sub(vc[i], R.add(R.zero, tuple(sum(Hi[i][k] * g[k][a] for k in range(n1)) for a in range(D)))) for i in range(n1)]
    assert all(R.iszero(ev(gr, vc + [wg])) for gr in grads), 'v_c(w) not converged'
    sh = [{0: {(0,) * n1: vc[i]}, 1: {tuple(1 if a == i else 0 for a in range(n1)): R.one}} for i in range(n1)]

    def inner(P, Dx):
        Rg = {}
        for d, m in P.items():
            for e, c in m.items():
                t = gconst(R, n1, R.of(c))
                for i in range(n1):
                    for _ in range(e[i]): t = gmul(R, t, sh[i], Dx)
                wpow = R.one
                for _ in range(e[n1]): wpow = R.mul(wpow, wg)
                Rg = gadd(R, Rg, t, wpow)
        return Rg
    Pin = inner(Phi, 2 * N + 2); phiw = Pin.get(0, {}).get((0,) * n1, R.zero)
    Pin = {d: m for d, m in Pin.items() if d >= 1}
    B, Hin = sp_series(R, Pin, [inner(A, 2 * (N - m)) for m, A in enumerate(amps)], N, n1)
    beta = rpow(R, R.mul(det(R, Hin), R.inv(R.of(det(Q, Hvv)))), Fr(-1, 2))
    if nobeta: beta = R.one
    oamps = [list(R.mul(beta, b)) for b in B]
    s_it = sp1(Q, list(phiw[:2 * N + 3]), oamps, N)
    return s_it == s_tot


def test_L2(seed, ctrl=None):
    rng = random.Random(seed); g = [Fr(0)] + [rf(rng) or Fr(1) for _ in range(3)]; f = [Fr(0)] + [rf(rng) for _ in range(4)]
    sc = -f[1] / g[1]; n = 2                                      # variables (s, x), sigma = sc + s
    P = {}
    for k in range(1, 5):
        c = sc * (g[k] if k < 4 else 0) + f[k]
        if c: P.setdefault(k, {})[(0, k)] = c
        if k < 4 and g[k]: P.setdefault(k + 1, {})[(1, k)] = g[k]
    if ctrl == 'sigma2': P.setdefault(2, {})[(2, 0)] = Fr(1, 3)
    a0 = {0: {(0, 0): Fr(1)}, 1: {(0, 1): rf(rng)}, 2: {(0, 2): rf(rng)}, 3: {(0, 3): rf(rng)}}
    a1 = {0: {(0, 0): rf(rng)}, 1: {(0, 1): rf(rng)}}; a2 = {0: {(0, 0): rf(rng)}}
    if ctrl == 'sigamp': a0 = gmul(Q, a0, {0: {(0, 0): Fr(1)}, 1: {(1, 0): Fr(2)}}, 2 * N)
    s, H = sp_series(Q, P, [a0, a1, a2], N, n)
    return s == [Fr(1), a1[0][(0, 0)], a2[0][(0, 0)]] + [Fr(0)] * (N - 2)


def test_L1c(DQ, NZ, ctrl=None):
    D = DQ + 1; R = RS(D); qg = tuple(Fr(1) if i == 1 else Fr(0) for i in range(D)); yc = R.of(5)
    for _ in range(D + 2):                                         # Newton: F(y) = 5/y - 1 + 5 q y^4
        Fy = R.add(R.sub(R.mul(R.of(5), R.inv(yc)), R.one), R.mul(R.of(5), R.mul(qg, R.mul(R.mul(yc, yc), R.mul(yc, yc)))))
        dF = R.add(R.neg(R.mul(R.of(5), R.inv(R.mul(yc, yc)))), R.mul(R.of(20), R.mul(qg, R.mul(yc, R.mul(yc, yc)))))
        yc = R.sub(yc, R.mul(Fy, R.inv(dF)))
    iy = R.inv(yc); K = 2 * (NZ + DQ) + 2; phi = [R.zero] * (K + 1); p = R.one
    y5 = [R.one]
    for _ in range(5): y5.append(R.mul(y5[-1], yc))
    for k in range(1, K + 1):
        p = R.mul(p, iy); v = R.smul(Fr(5 * (-1) ** (k + 1), k), p)
        if k <= 5: v = R.add(v, R.mul(qg, R.smul(comb(5, k), y5[5 - k])))
        if k == 1: v = R.sub(v, R.one)                             # the -y term
        phi[k] = v
    assert R.iszero(phi[1]), 'critical point'
    NN = NZ + DQ; Nq = sp1(R, phi, [[R.one]], NN)
    beta = rpow(R, R.smul(-1, R.mul(R.smul(2, phi[2]), R.of(5))), Fr(-1, 2))
    u = R.sub(R.mul(yc, R.of(Fr(1, 5))), R.one); lg = R.zero; up = R.one
    for k in range(1, D + 1): up = R.mul(up, u); lg = R.add(lg, R.smul(Fr((-1) ** (k + 1), k), up))
    Dl = R.add(R.sub(R.smul(5, lg), R.sub(yc, R.of(5))), R.mul(qg, y5[5]))
    if ctrl == 'nodelta': Dl = R.zero                              # control: critical-value exponential dropped
    if ctrl == 'nobeta': beta = R.one                              # control: Hessian ratio dropped
    lhs = {}; Ej = R.one
    for j in range(DQ + 1):
        for m in range(NN + 1):
            t = R.mul(Ej, R.mul(beta, Nq[m]))
            for d in range(DQ + 1): lhs[(d, m - j)] = lhs.get((d, m - j), 0) + t[d]
        Ej = R.smul(Fr(1, j + 1), R.mul(Ej, Dl))
    ok = ok3 = True; fac = Fr(1)
    for d in range(DQ + 1):
        if d: fac /= d
        Rn = sp1(Q, [Fr(0), Fr(0)] + [Fr(5 * (-1) ** (k + 1), k) / 5 ** k for k in range(2, 2 * (NN + d) + 3)],
                 [[comb(5 * d, k) * Fr(5) ** (5 * d - k) for k in range(2 * (NN + d) + 1)]], NN + d)
        poch = [Fr(1)]
        for i in range(1, 5 * d + 1): poch = [(poch[k] if k < len(poch) else 0) + (Fr(i, 5) * poch[k - 1] if k else 0) for k in range(len(poch) + 1)]
        for n in range(-d, NZ + 1):
            ok = ok and lhs.get((d, n), 0) == fac * Rn[n + d]
            ok3 = ok3 and Rn[n + d] == Fr(5) ** (5 * d) * sum(poch[k] * _st(n + d - k) for k in range(min(len(poch), n + d + 1)))
    return ok, ok3


def _st(n, _c={}):
    if not _c:
        from ch_core import bk_list, sexp_q
        bk = bk_list(40); H = [Fr(0)] * 61
        for k in range(1, 31): H[2 * k - 1] = bk[k] / Fr(5) ** (2 * k - 1)
        _c['s'] = sexp_q(H, 60)
    return _c['s'][n]


out('ch_lemmas: exact over Q, through z^%d, %d seeds per test; log %s' % (N, SEEDS, LOGREL))
for sd in range(SEEDS):
    check('(S1) seed %d: N[Phi, d_i(e^(Phi/z) h)] = 0 for i = 1, 2, all orders <= %d' % (sd, N), test_S1(100 + sd))
    check('(L1a) seed %d: change of variables, N transforms by det D chi(0)' % sd, test_L1a(200 + sd))
    check('(L1b) seed %d: Fubini, inner dimension 1' % sd, test_L1b(300 + sd, 1))
    check('(L1b) seed %d: Fubini, inner dimension 2' % sd, test_L1b(400 + sd, 2))
    check('(L2) seed %d: delta lemma, N = A(X_0, z) exactly' % sd, test_L2(500 + sd))
check('control: (L1a) without the Jacobian fails', not test_L1a(200, nojac=True))
check('control: (L1b) without beta fails', not test_L1b(300, 1, nobeta=True))
check('control: (L2) with a sigma-dependent amplitude fails', not test_L2(500, 'sigamp'))
check('control: (L2) with a sigma^2 term in the phase fails', not test_L2(500, 'sigma2'))
ok, ok3 = test_L1c(4, N)
check('(L1c)(iii): q-dependent critical point re-expanded = Sum_d (q/z)^d/d! N[psi_0, y^(5d)], q^0..q^4, z^-d..z^%d' % N, ok)
check('(L1c)(iii) + (L3): N[psi_0, y^(5d); 5] = 5^(5d) prod_(i=1..5d) (1 + iz/5) exp(Sum b_k (z/5)^(2k-1))', ok3)
check('control: (L1c)(iii) without the critical-value exponential e^((psi_q(y_c) - psi_0(5))/z) fails', not test_L1c(4, N, 'nodelta')[0])
check('control: (L1c)(iii) without the Hessian ratio fails', not test_L1c(4, N, 'nobeta')[0])
out('CH_LEMMAS: %s   (%.0f s)' % ('ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
