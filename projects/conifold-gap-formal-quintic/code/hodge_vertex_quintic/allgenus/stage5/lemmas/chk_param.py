# -*- coding: utf-8 -*-
"""chk_param.py -- Lemmas L1c (perturbation, T-q and T-z), FL4 (d/dq) and FL3 for the SP (u-dependent change of
variables with moving critical points), exact over Q.  Parameters q, eps, u have weight 1 (doubled weight 2).

S(q)  := e^((F(c(q)) - f(0))/z) (Delta(q)/Delta(0))^(-1/2) SP^norm_{c(q)}[F(.; q); a]     (non-perturbative side)
P(q)  := SP^norm_0[f; e^((F - f)/z) a],  F = f + q g + q^2 g2, g(0) = g2(0) = 0             (perturbative side)
T-q : S = P.   T-z : [S(eps)]_(eps = z) = SP^norm_0[f; e^g a]  (phase f + z g, critical point moving with z).
FL4 : z dS/dq = e^(..)(..)^(-1/2) SP^norm_{c(q)}[F; (d_q F) a + z d_q a]  (FL4 times z).
FL3 (SP): Y = X (1 + u h(X)); SP_Y at Y_c(u) of e^(F(Y;u)/z) b(Y) dY = SP_X at X_c(u) of the pull-back, with the
Hessian transformation Delta_X = Delta_Y Y'(X_c)^2 (normalised at u = 0).
T(b): SP^norm_0[f; e^((f(x-z) - f(x))/z) a(x-z)] = SP^norm_0[f; a] (translation by -z; gamma's (S4)).
Controls (must fail): T-q without the determinant ratio; FL4 without the d_q F/z term; T-z and T(b) without the exponential
factor; FL3 without Y'(X).
Usage: py -I -B chk_param.py [K]     (doubled-weight cut 2K: all q^d z^k with d + k <= K; default 8)
"""
import os, sys, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from sp_core import Ctx, rand_poly, inv_q
from sp_iter import crit_point, expand_at, hessian_const
from sp_rand import quad_form, rand_sym
from chk_util import Log

K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
log = Log(HERE, 'chk_param_%d' % K)
log.out('chk_param: L1c (T-q, T-z), FL4, FL3-SP; all terms of weight <= %d (q^d z^k, d + k <= %d); log %s' % (K, K, log.rel))


def nonpert(c, cb, F, amp_fn, xs):
    """S = e^((G - G(0))/z) (det K/det K(0))^(-1/2) SP^norm at the moving critical point; amp_fn(subs) -> amplitude"""
    xi = crit_point(cb, F, xs); G, Km, R = expand_at(cb, F, xs, xi)
    subs = {i: cb.add(xi[i], cb.var(i)) for i in xs}
    Kc = [[c.trunc(x) for x in row] for row in Km]; C = [[c.scal(x, -1) for x in row] for row in c.inv(Kc)]
    b = c.sp_norm(R, c.trunc(amp_fn(subs)), xs, C)
    D = c.det(Kc); ratio = c.add(c.scal(D, 1 / c.cpart(D)), c.const(1), -1)
    Gs = c.zp(cb.add(G, cb.const(cb.cpart(G)), -1), -1)
    return c.mul(c.mul(c.exp(Gs), c.binom_pow(ratio, Fr(-1, 2))), b), b, D, xi


def pert(c, f, P, amp, xs):
    """SP^norm_0[f; e^(P/z) amp]: P = everything beyond the quadratic part of f (must have w2 >= 3)"""
    H = hessian_const(c, f, xs); C = [[c.const(-x) for x in row] for row in inv_q(H)]
    return c.sp_norm(P, amp, xs, C)


def c_eq(A, B):
    return all(A.get(k, 0) == B.get(k, 0) for k in set(A) | set(B))


for n in (1, 2):
    xs = list(range(n)); q = n; nv = n + 1
    c = Ctx(nv, [1] * n + [2], 2 * K); cb = Ctx(nv, [1] * n + [2], 2 * K + 2)
    rng = random.Random(32452843 + n); t0 = time.time()
    H = rand_sym(rng, n); f = cb.add(quad_form(H, xs, nv), rand_poly(cb, rng, xs, 3, 5, 3, 0.7))
    g = rand_poly(cb, rng, xs, 1, 3, 3, 0.8); g2 = rand_poly(cb, rng, xs, 1, 2, 3, 0.8)
    F = cb.add(cb.add(f, cb.mul(cb.var(q), g)), cb.mul(cb.var(q, 1), cb.mul(cb.var(q), g2)))
    a = c.add(c.const(1), rand_poly(c, rng, xs + [q], 1, 3, 2, 0.6))
    S, b, D, xi = nonpert(c, cb, F, lambda s: cb.compose(a, s), xs)
    P = pert(c, f, cb.add(cb.add(F, quad_form(H, xs, nv), -1), {}), c.trunc(a), xs)
    log.check('T-q n = %d: S(q) = P(q) (moving critical point versus perturbation) for weight <= %d' % (n, K),
              c_eq(S, P), '%d terms  (%.1f s)' % (len(S), time.time() - t0))
    S0 = c.mul(c.exp(c.zp(c.add(cb.compose(F, xi), {}), -1)), b)                      # control: no det ratio
    log.control('T-q n = %d: without (Delta(q)/Delta(0))^(-1/2)' % n, not c_eq(S0, P))
    # FL4, multiplied by z (so that the amplitude has weight >= 0): z dS/dq = S'(q) with amplitude
    # (d_q F) a + z d_q a; z dS/dq is complete for weight <= K (S is complete for weight <= K, d/dq lowers it by 1)
    dF = cb.deriv(F, q); da = c.deriv(a, q)
    amp1 = lambda s: cb.add(cb.mul(cb.compose(dF, s), cb.compose(a, s)), cb.zp(cb.compose(da, s), 1))
    S1, _, _, _ = nonpert(c, cb, F, amp1, xs)
    dS = c.zp(c.deriv(S, q), 1)
    log.check('FL4 n = %d: d/dq of the full SP = full SP of d/dq(integrand) (times z), weight <= %d' % (n, K),
              c_eq(dS, S1))
    amp1b = lambda s: cb.zp(cb.compose(da, s), 1)
    S1b, _, _, _ = nonpert(c, cb, F, amp1b, xs)
    log.control('FL4 n = %d: without the (d_q F)/z term' % n, not c_eq(dS, S1b))
    # T-z: substitute eps = q -> z in S and compare with SP^norm_0[f; e^g a] (a taken q-free here)
    a0 = {k: v for k, v in a.items() if k[0][q] == 0}
    F1 = cb.add(f, cb.mul(cb.var(q), g))
    S, _, _, _ = nonpert(c, cb, F1, lambda s: cb.compose(a0, s), xs)
    Sz = {}
    for (e, j), v in S.items():
        k = (tuple(x if i != q else 0 for i, x in enumerate(e)), j + e[q]); Sz[k] = Sz.get(k, 0) + v
    gz = c.trunc({k: v for k, v in g.items()})
    Pz = pert(c, f, {k: v for k, v in f.items() if sum(k[0]) >= 3}, c.mul(c.exp(gz), a0), xs)
    log.check('T-z n = %d: [S(eps)]_(eps = z) = SP^norm_0[f; e^g a] through z^%d' % (n, K), c_eq(c.trunc(Sz), Pz))
    Pz_bad = pert(c, f, {k: v for k, v in f.items() if sum(k[0]) >= 3}, a0, xs)
    log.control('T-z n = %d: e^g dropped from the amplitude' % n, not c_eq(c.trunc(Sz), Pz_bad))
    # Lemma T(b) (= gamma's (S4)): translation by p = -z, SP^norm_0[f; e^((f(x-z) - f(x))/z) a(x-z)] = SP^norm_0[f; a]
    f3 = {k: v for k, v in f.items() if sum(k[0]) >= 3}
    sh = {i: cb.add(cb.var(i), {(cb.zero_e(), 1): Fr(-1)}) for i in xs}
    E = c.trunc(cb.zp(cb.add(cb.compose(f, sh), f, -1), -1))
    aT = c.trunc(cb.compose(a0, sh))
    lhsT = pert(c, f, f3, c.mul(c.exp(E), aT), xs); rhsT = pert(c, f, f3, a0, xs)
    log.check('T(b) n = %d: translation by -z leaves SP^norm unchanged, through z^%d' % (n, K), c_eq(lhsT, rhsT))
    log.control('T(b) n = %d: translated amplitude without e^((f(x-z) - f(x))/z)' % n, not c_eq(pert(c, f, f3, aT, xs), rhsT))

# FL3 for the SP: variables t (local coordinate) and u (parameter)
c = Ctx(2, [1, 2], 2 * K); cb = Ctx(2, [1, 2], 2 * K + 2)
rng = random.Random(49979687); t, u = 0, 1
FY = cb.add(quad_form([[-3]], [t], 2), rand_poly(cb, rng, [t], 3, 5, 3, 0.8))
FY = cb.add(FY, cb.mul(cb.var(u), rand_poly(cb, rng, [t], 1, 3, 3, 0.8)))        # critical point moves with u
beta = c.add(c.const(1), rand_poly(c, rng, [t, u], 1, 3, 2, 0.7))
h = rand_poly(cb, rng, [t], 0, 2, 2, 0.8)                                       # Y = X (1 + u h(X)), X = 1 + t
X = cb.add(cb.const(1), cb.var(t)); Ymap = cb.add(X, cb.mul(cb.var(u), cb.mul(X, h)))
Ysub = {t: cb.add(Ymap, cb.const(1), -1)}                                      # t_Y = Y(X) - 1 as a series in t_X, u
FX = cb.compose(FY, Ysub); dY = cb.deriv(Ymap, t)
SY, bY, DY, _ = nonpert(c, cb, FY, lambda s: cb.compose(beta, s), [t])
SX, bX, DX, _ = nonpert(c, cb, FX, lambda s: cb.mul(cb.compose(cb.compose(beta, Ysub), s), cb.compose(dY, s)), [t])
log.check('FL3 (SP): SP_Y at Y_c(u) = SP_X at X_c(u) of the pull-back, weight <= %d (critical points move with u)' % K,
          c_eq(SY, SX), '%d terms' % len(SY))
SXb, _, _, _ = nonpert(c, cb, FX, lambda s: cb.compose(cb.compose(beta, Ysub), s), [t])
log.control('FL3 (SP): pull-back without the factor Y\'(X)', not c_eq(SY, SXb))
log.finish('CHK_PARAM %d' % K)
