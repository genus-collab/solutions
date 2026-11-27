# -*- coding: utf-8 -*-
"""ch_hess_sym.py -- Stage 5 / chart: symbolic part (polynomials in L over Q(zeta_5)) and controls of ch_hess.py.
Polynomials are lists of Q(zeta) elements (index = power of L); determinants by fraction-free Bareiss elimination with exact
polynomial division. Imported by ch_hess.py; not a stand-alone script.
"""
from fractions import Fraction as Fr
from ch_wick import hessian, linear_is_zero, det


def _trim(F, a):
    a = list(a)
    while a and F.iszero(a[-1]): a.pop()
    return a


def padd(F, a, b):
    n = max(len(a), len(b)); return _trim(F, [F.add(a[i] if i < len(a) else F.zero, b[i] if i < len(b) else F.zero) for i in range(n)])


def pneg(F, a): return [F.neg(x) for x in a]


def pmul(F, a, b):
    if not a or not b: return []
    c = [F.zero] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if F.iszero(x): continue
        for j, y in enumerate(b): c[i + j] = F.add(c[i + j], F.mul(x, y))
    return _trim(F, c)


def pdiv(F, a, b):
    """exact division a/b (asserts zero remainder)"""
    a = _trim(F, a); b = _trim(F, b); assert b, 'division by 0'
    q = [F.zero] * max(len(a) - len(b) + 1, 0); ib = F.inv(b[-1]); r = list(a)
    for k in range(len(a) - len(b), -1, -1):
        c = F.mul(r[k + len(b) - 1], ib); q[k] = c
        for i, y in enumerate(b): r[k + i] = F.sub(r[k + i], F.mul(c, y))
    assert not _trim(F, r), 'inexact division'
    return _trim(F, q)


def bareiss(F, M):
    n = len(M); A = [[_trim(F, x) for x in row] for row in M]; sign = 1; prev = [F.one]
    for k in range(n - 1):
        if not A[k][k]:
            r = next((r for r in range(k + 1, n) if A[r][k]), None)
            if r is None: return []
            A[k], A[r] = A[r], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = pdiv(F, padd(F, pmul(F, A[i][j], A[k][k]), pneg(F, pmul(F, A[i][k], A[k][j]))), prev)
        prev = A[k][k]
    return A[n - 1][n - 1] if sign == 1 else pneg(F, A[n - 1][n - 1])


def symbolic_part(F, check, out):
    one = F.one; Z = [F.zeta(j) for j in range(5)]; Lp = [F.zero, one]
    lin = [[Z[j], F.neg(one)] for j in range(5)]                                    # zeta^j - L
    Lam = [one]
    for j in range(5): Lam = pmul(F, Lam, lin[j])
    check('symbolic: prod_j (zeta^j - L) = 1 - L^5 in Q(zeta)[L]', Lam == _trim(F, [one, F.zero, F.zero, F.zero, F.zero, F.neg(one)]))
    Lj = []
    for j in range(5):
        p = [one]
        for i in range(5):
            if i != j: p = pmul(F, p, lin[i])
        Lj.append(p)
    S = []
    for p in Lj: S = padd(F, S, p)
    check('symbolic: Sum_j prod_(i!=j) (zeta^i - L) = 5 L^4, i.e. Sum_j 1/(zeta^j - L) = 5L^4/(1 - L^5)',
          S == _trim(F, [F.zero] * 4 + [F.of(5)]))
    c = lambda x: [F.of(x)]
    M = [[[] for _ in range(7)] for _ in range(7)]
    M[0][1] = Lam
    for j in range(5): M[0][2 + j] = pneg(F, pmul(F, Lp, Lj[j]))                   # row sigma x L prod(zeta^j - L)
    M[1][0] = c(5); M[1][1] = c(-1)                                                  # row y x 5L
    for j in range(5): M[2 + j][0] = c(-1); M[2 + j][2 + j] = c(-1)                  # row x_j x (zeta^j - L)
    d = bareiss(F, M)
    check('symbolic: det of H_7 with rows scaled by L prod(zeta^j - L), 5L, (zeta^j - L) is the constant 5, '
          'hence det H_7 = 5/(5 L^2 (1-L^5)^2) = 1/(L^2 (1-L^5)^2) for all L', d == [F.of(5)], [str(x) for x in (d[0] if d else [])])
    v = [5, -1, -1, -1, -1]; Km = [F.neg(one), one]                                  # K = L - 1
    Hl = [[pmul(F, Km, c(v[i] * v[k])) for k in range(5)] for i in range(5)]
    Hl[0][0] = padd(F, Hl[0][0], [F.zero, F.of(-5)])
    for j in range(1, 5): Hl[j][j] = padd(F, Hl[j][j], pneg(F, lin[j]))
    dl = bareiss(F, Hl)
    check('symbolic: chart Hessian in log coordinates, diag(-5L, -(zeta^j - L)) + (L-1) v v^T, has determinant the constant -25',
          dl == [F.of(-25)], [str(x) for x in (dl[0] if dl else [])])
    return Hl


def controls_part(F, check, out, phasech, q_of):
    one = F.one; Z = [F.zeta(j) for j in range(5)]
    for L in (Fr(2), Fr(-2, 3), Fr(7, 5)):                                           # log-coordinate Hessian, entrywise
        Pc, K = phasech(L); Hc = hessian(F, Pc, 5); Dg = [F.of(5 * L)] + [F.sub(Z[j], F.of(L)) for j in range(1, 5)]
        v = [5, -1, -1, -1, -1]; good = True
        for i in range(5):
            for k in range(5):
                e = F.of((L - 1) * v[i] * v[k])
                if i == k: e = F.sub(e, Dg[i])
                good &= F.eq(F.mul(F.mul(Dg[i], Dg[k]), Hc[i][k]), e)
        check('log coordinates: D H_ch D = diag(-5L, -(zeta^j - L)) + (L-1) v v^T entrywise at L = %s (D = diag(5L, zeta^j - L))' % L, good)
    P1, K1 = phasech(Fr(1)); H1 = hessian(F, P1, 5)
    diag = [F.of(Fr(-1, 5))] + [F.neg(F.inv(F.sub(Z[j], one))) for j in range(1, 5)]
    check('q = 0 (L = 1): K = 0, gradient 0, chart Hessian = diag(-1/5, -1/(zeta^j - 1)) exactly',
          F.iszero(K1) and linear_is_zero(F, P1) and all(F.eq(H1[i][k], diag[i] if i == k else F.zero) for i in range(5) for k in range(5)))
    check('q = 0: det of the chart Hessian = -1/25 (= lim of -(1-L)^2/(L^2(1-L^5)^2))', F.eq(det(F, H1), F.of(Fr(-1, 25))))
    Pm, Km = phasech(Fr(2), qsign=-1)
    check('control: q -> -q in the chart phase moves the critical point (gradient at p_ch(L) != 0, L = 2)',
          not linear_is_zero(F, Pm))
    L = Fr(2); c = [F.sub(Z[j], F.of(L)) for j in range(5)]; prodc = one
    for x in c: prodc = F.mul(prodc, x)
    check('control: the sigma-equation with +q instead of -q fails at L = 2', not F.eq(F.mul(q_of(L), F.of((5 * L) ** 5)), prodc))
    check('control: a wrong closed form det H_7 = 1/(L^2 (1-L^5)) is rejected at L = 2',
          not F.eq(F.of(1 / (L ** 2 * (1 - L ** 5) ** 2)), F.of(1 / (L ** 2 * (1 - L ** 5)))))
