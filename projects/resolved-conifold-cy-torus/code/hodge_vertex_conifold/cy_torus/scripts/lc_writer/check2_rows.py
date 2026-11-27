# Fresh exact check (Fractions): the code's Wick row (in y = log U, via the rescaled wr_k, Lemma 3.9) at generic CY weights
# against the KEY FORMULA; L1 membership from the Wick-derived G; A_j closed form; H-row, Row 2, K_1, E_11 at t^0; a0 = 0 rescaling.
# Evidence only (finite orders in tau/t and z).
import sys, itertools
from fractions import Fraction as Fr
from math import factorial
T = int(sys.argv[1]) if len(sys.argv) > 1 else 6      # tau/t order
Z = int(sys.argv[2]) if len(sys.argv) > 2 else 9      # z order

# ---------- truncated power series in tau (lists of length T+1) ----------
def ps(c): return [Fr(x) for x in c] + [Fr(0)] * (T + 1 - len(c))
ONE = ps([1]); TAU = ps([0, 1])
def padd(a, b): return [x + y for x, y in zip(a, b)]
def psub(a, b): return [x - y for x, y in zip(a, b)]
def psc(k, a): return [k * x for x in a]
def pmul(a, b):
    out = [Fr(0)] * (T + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(T + 1 - i): out[i + j] += x * b[j]
    return out
def pinv(a):
    assert a[0] != 0
    out = [Fr(0)] * (T + 1); out[0] = 1 / a[0]
    for n in range(1, T + 1):
        out[n] = -sum(a[j] * out[n - j] for j in range(1, n + 1)) / a[0]
    return out
def ppow(a, k):
    r = ONE
    for _ in range(k): r = pmul(r, a)
    return r
def psqrt1(a):          # square root with value 1 of a series with a[0] = 1
    assert a[0] == 1
    out = [Fr(0)] * (T + 1); out[0] = Fr(1)
    for n in range(1, T + 1):
        out[n] = (a[n] - sum(out[j] * out[n - j] for j in range(1, n))) / 2
    return out
def pexp(a):            # exp of a series with a[0] = 0
    assert a[0] == 0
    out = [Fr(0)] * (T + 1); out[0] = Fr(1)
    da = [(i + 1) * a[i + 1] for i in range(T)] + [Fr(0)]
    for n in range(1, T + 1):
        out[n] = sum(da[j] * out[n - 1 - j] for j in range(n)) / n
    return out
def pder(a): return [(i + 1) * a[i + 1] for i in range(T)] + [Fr(0)]
def pint(a):            # integral from 0
    return [Fr(0)] + [a[i] / (i + 1) for i in range(T)]
def pcompose(a, b):     # a(b), b[0] = 0
    assert b[0] == 0
    out = [Fr(0)] * (T + 1)
    pw = ONE
    for i in range(T + 1):
        out = padd(out, psc(a[i], pw)); pw = pmul(pw, b)
    return out

# ---------- Laurent series in t: (v, list) with list[i] = coeff of t^(v+i), length T+1 ----------
def lmul(a, b): return (a[0] + b[0], pmul(a[1], b[1]))
def ladd(a, b):
    v = min(a[0], b[0]); n = T + 1
    out = [Fr(0)] * n
    for i in range(n):
        if 0 <= a[0] - v + i < n: out[a[0] - v + i] += a[1][i]
        if 0 <= b[0] - v + i < n: out[b[0] - v + i] += b[1][i]
    # precision: lowest of the two absolute precisions
    return (v, out)
def lsc(k, a): return (a[0], psc(k, a[1]))
def lder(a):            # d/dt
    v, c = a
    out = [Fr(0)] * (T + 1)
    for i in range(T + 1):
        if v + i != 0 and i + 1 <= T: pass
    # t^(v+i) -> (v+i) t^(v+i-1): new valuation v-1, coefficients shift
    out = [(v + i) * c[i] for i in range(T + 1)]
    return (v - 1, out)      # note: last coefficient is unreliable at absolute order v+T; we only read low orders
def lcoef(a, k):
    v, c = a
    if k < v: return Fr(0)
    return c[k - v]
def linv(a):
    v, c = a
    return (-v, pinv(c))

# tau(t) and t(tau)
TAU_T = ps([0] + [Fr((-1) ** (k + 1), factorial(k)) for k in range(1, T + 1)])       # 1 - e^-t
TAU_OVER_T = ps([Fr((-1) ** k, factorial(k + 1)) for k in range(T + 1)])           # (1 - e^-t)/t
T_OVER_TAU = pinv(TAU_OVER_T)
def tau_to_t(f_tau, v=0):
    """f = tau^v * (regular series in tau); returns Laurent series in t"""
    g = pcompose(f_tau, TAU_T)
    if v >= 0:
        g = pmul(g, ppow(TAU_OVER_T, v))
    else:
        g = pmul(g, ppow(T_OVER_TAU, -v))
    return (v, g)

# ---------- theta-polynomials Pol_n(p) = (p^2 - p) Pol_{n-1}'(p), Pol_0 = p ----------
def thetapolys(nmax):
    P = [[Fr(0), Fr(1)]]
    for n in range(1, nmax + 1):
        prev = P[-1]
        dprev = [i * prev[i] for i in range(1, len(prev))]
        p2mp = [Fr(0), Fr(-1), Fr(1)]
        new = [Fr(0)] * (len(dprev) + 2)
        for i, x in enumerate(p2mp):
            for j, y in enumerate(dprev): new[i + j] += x * y
        P.append(new)
    return P

# ---------- Wick coefficients: wickZ_m as a function of W_2, ..., W_{2m+2} (series) ----------
def partitions_min3(total, parts):
    """partitions of total into exactly `parts` parts, each >= 3, nonincreasing"""
    def rec(rem, k, maxp):
        if k == 0:
            if rem == 0: yield ()
            return
        for p in range(min(maxp, rem - 3 * (k - 1)), 2, -1):
            for rest in rec(rem - p, k - 1, p): yield (p,) + rest
    yield from rec(total, parts, total)
def dfact(n):
    r = 1
    while n > 1: r *= n; n -= 2
    return r
def wick_terms(m):
    """list of (coefficient, n, [(k, multiplicity)]) : wickZ_m = sum coeff * W_2^-n prod W_k^mult"""
    out = []
    for n in range(m, 3 * m + 1):
        pref = Fr(factorial(2 * n) * (-1) ** n, factorial(n) * 2 ** n)
        for p in partitions_min3(2 * n, n - m):
            mult = {}
            for k in p: mult[k] = mult.get(k, 0) + 1
            c = pref
            for k, mu in mult.items(): c /= factorial(k) ** mu * factorial(mu)
            out.append((c, n, sorted(mult.items())))
    return out

class CYData:
    def __init__(self, a0, a1, a2, e1, e2):
        a0, a1, a2, e1, e2 = map(Fr, (a0, a1, a2, e1, e2))
        self.a0, self.a1, self.a2, self.e1, self.e2 = a0, a1, a2, e1, e2
        E = a0 + a1 + a2; A = a0 + a2; B = a0 + a1
        self.E, self.A, self.B = E, A, B
        self.u = (E - 2 * a1) * (E - 2 * a2)
        # s*(tau): Newton for B tau s^2 - (E + a0 tau) s + A = 0, s(0) = A/E
        s = ps([A / E])
        for _ in range(6):
            N = padd(psub(psc(B, pmul(TAU, pmul(s, s))), pmul(padd(ps([E]), psc(a0, TAU)), s)), ps([A]))
            dN = psub(psc(2 * B, pmul(TAU, s)), padd(ps([E]), psc(a0, TAU)))
            s = psub(s, pmul(N, pinv(dN)))
        self.s = s
        self.delta = psub(padd(ps([E]), psc(a0, TAU)), psc(2 * B, pmul(TAU, s)))
        self.U1 = psub(ONE, pmul(TAU, s))
        self.L2 = psub(psc(B, s), ps([a0]))
        # alg^2 = a1 A delta U1/(E^3 s (1-s)) ; alg with alg(0) = 1
        alg2 = psc(a1 * A / E ** 3, pmul(pmul(self.delta, self.U1), pinv(pmul(s, psub(ONE, s)))))
        assert alg2[0] == 1
        self.alg = psqrt1(alg2)
        # dphi: dphi_tau = a2 s/(1 - tau s)
        self.dphi = pint(psc(a2, pmul(s, pinv(self.U1))))
        # rescaled derivatives wr_k = tau^k W_k(U_1), k = 2..2Z+2
        Pol = thetapolys(2 * Z + 2)
        X1 = psc(-1, pmul(psub(ONE, TAU), pinv(psub(ONE, s))))   # tau p_Q(U_1) = -(1-tau)/(1-s*)
        X2 = pinv(s)                                                # tau p_1(U_1) = 1/s*
        self.wr = {}
        for k in range(2, 2 * Z + 3):
            pc = Pol[k - 1]
            tot = ps([0])
            for i in range(1, k + 1):
                if pc[i] == 0: continue
                term = padd(psc(a1, ppow(X1, i)), psc(A, ppow(X2, i)))
                term = pmul(term, ppow(TAU, k - i))
                tot = padd(tot, psc(pc[i], term))
            self.wr[k] = psc(-1, tot)
        # Wick row in y = log U (= in terms of wr_k by homogeneity): point-frame row without constants
        W2i = pinv(self.wr[2])
        self.rowpt = []
        for m in range(Z + 1):
            tot = ps([0])
            for c, n, mult in wick_terms(m):
                term = psc(c, ppow(W2i, n))
                for k, mu in mult: term = pmul(term, ppow(self.wr[k], mu))
                tot = padd(tot, term)
            self.rowpt.append(tot)

def bern(n):
    A = [Fr(0)] * (n + 1)
    for m in range(n + 1):
        A[m] = Fr(1, m + 1)
        for j in range(m, 0, -1): A[j - 1] = j * (A[j - 1] - A[j])
    return A[0] if n != 1 else Fr(-1, 2)
def h_class(ws):
    """h_v(z) = Exp[-Sum_k bk_k z^(2k-1) Sum_a v_a^(1-2k)] as a z-series (list length Z+1)"""
    ex = [Fr(0)] * (Z + 1)
    k = 1
    while 2 * k - 1 <= Z:
        bk = bern(2 * k) / (2 * k * (2 * k - 1))
        ex[2 * k - 1] = -bk * sum(Fr(w) ** (1 - 2 * k) for w in ws)
        k += 1
    out = [Fr(0)] * (Z + 1); out[0] = Fr(1)
    for n in range(1, Z + 1):
        out[n] = sum(j * ex[j] * out[n - j] for j in range(1, n + 1)) / n
    return out
# z-series (regular) helpers, length Z+1
def zmul(a, b):
    out = [Fr(0)] * (Z + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(Z + 1 - i): out[i + j] += x * b[j]
    return out
def zinv(a):
    out = [Fr(0)] * (Z + 1); out[0] = 1 / a[0]
    for n in range(1, Z + 1): out[n] = -sum(a[j] * out[n - j] for j in range(1, n + 1)) / a[0]
    return out
def zlin(c0, c1): return [Fr(c0), Fr(c1)] + [Fr(0)] * (Z - 1)

def key_formula_check(D, label):
    E, A, B, a0, a1, a2 = D.E, D.A, D.B, D.a0, D.a1, D.a2
    print('==', label, (D.a0, D.a1, D.a2, D.e1, D.e2))
    # LHS: frame-B row = h_(a1,A,-E)(z) * rowpt : double series [n][m], m >= 0
    hB = h_class([a1, A, -E])
    LHS = {}
    for m in range(Z + 1):
        for n in range(T + 1):
            LHS[(n, m)] = sum(hB[j] * D.rowpt[m - j][n] for j in range(m + 1))
    # RHS: (E/(E+z)) F Ephi alg, with Laurent z-powers m >= -T
    # F_n(z) as Laurent series: z^-n * prod(a2+jz)(A+jz)/(n! prod_{j=2}^{n+1}(E+jz))  -- regular part as z-series length Z+T+1
    ZZ = Z + T
    def zmul2(a, b):
        out = [Fr(0)] * (ZZ + 1)
        for i, x in enumerate(a):
            if x:
                for j in range(ZZ + 1 - i): out[i + j] += x * b[j]
        return out
    def zinv2(a):
        out = [Fr(0)] * (ZZ + 1); out[0] = 1 / a[0]
        for n in range(1, ZZ + 1): out[n] = -sum(a[j] * out[n - j] for j in range(1, n + 1)) / a[0]
        return out
    def lin2(c0, c1): return [Fr(c0), Fr(c1)] + [Fr(0)] * (ZZ - 1)
    Fn = []   # Fn[n] = dict zexp -> coeff
    for n in range(T + 1):
        num = lin2(1, 0)
        for j in range(1, n + 1): num = zmul2(zmul2(num, lin2(a2, j)), lin2(A, j))
        den = lin2(factorial(n), 0)
        for j in range(2, n + 2): den = zmul2(den, lin2(E, j))
        q = zmul2(num, zinv2(den))
        Fn.append({i - n: q[i] for i in range(ZZ + 1)})
    # Ephi = Exp[-dphi/z]: tau^n coefficient is a polynomial in 1/z : compute via exp of series with coefficients in Q[1/z]
    # represent as dict (n, m) ; use the recursion out_n = (1/n) sum_j j a_j out_{n-j}, a_j = -dphi_j / z
    X = {(0, 0): Fr(1)}
    for n in range(1, T + 1):
        acc = {}
        for j in range(1, n + 1):
            if D.dphi[j] == 0: continue
            for (n2, m2), c in X.items():
                if n2 == n - j:
                    key = (n, m2 - 1)
                    acc[key] = acc.get(key, 0) + j * (-D.dphi[j]) * c
        for key, c in acc.items(): X[key] = c / n
    # product F * Ephi * alg * E/(E+z), truncated m <= Z, n <= T
    def dmul(P, Q):
        out = {}
        for (n1, m1), c1 in P.items():
            for (n2, m2), c2 in Q.items():
                n = n1 + n2; m = m1 + m2
                if n <= T and m <= Z + T: out[(n, m)] = out.get((n, m), 0) + c1 * c2
        return out
    Fd = {(n, m): c for n in range(T + 1) for m, c in Fn[n].items() if m <= Z + T}
    algd = {(n, 0): D.alg[n] for n in range(T + 1)}
    geo = zinv(zlin(E, 1)); geod = {(0, m): E * geo[m] for m in range(Z + 1)}
    RHS = dmul(dmul(dmul(Fd, X), algd), geod)
    RHS = {k: c for k, c in RHS.items() if k[1] <= Z}
    neg = [((n, m), c) for (n, m), c in RHS.items() if m < 0 and c != 0]
    print('   negative z-powers of (E/(E+z)) F Ephi alg vanish (n<=%d):' % T, len(neg) == 0)
    diff = [(k, LHS[k] - RHS.get(k, 0)) for k in LHS if LHS[k] - RHS.get(k, 0) != 0]
    print('   KEY FORMULA (Thm 3.17_CY): h_(a1,A,-E) x Wick row == (E/(E+z)) F Ephi alg to tau^%d z^%d:' % (T, Z), len(diff) == 0)
    if diff: print('      first mismatches', diff[:3])
    # t = 0 value: h_(E,-a1,-A)^-1 * rowpt(z;0) should be E/(E+z)  (Prop 3.15 rescaled)
    hbeta = h_class([E, -a1, -A])
    row0 = [D.rowpt[m][0] for m in range(Z + 1)]
    print('   r^pt_{0,1}(z;0) == h_(E,-a1,-A) E/(E+z):', row0 == zmul(hbeta, zmul([Fr(E)] + [Fr(0)] * Z, geo)))
    # ---------- G = rho r_{0,1} = (E/(E+z)) (tau/t) F Ephi as t-series; L1 and A_j ----------
    rowB = {m: [LHS[(n, m)] for n in range(T + 1)] for m in range(Z + 1)}      # frame-B row, tau-series per z-power
    rho_tau = pinv(D.alg)   # rho = (tau/t)/alg ; handle tau/t after conversion
    Gt = {}
    for m in range(Z + 1):
        f = pmul(rowB[m], rho_tau)
        v, g = tau_to_t(f, 0)
        g = pmul(g, TAU_OVER_T)
        Gt[m] = g              # t-series of [z^m] G
    # R_k(z) = z G_k(z): coefficient of z^m is Gt[m-1][k]
    def Rk(k): return [Fr(0)] + [Gt[m - 1][k] for m in range(1, Z + 1)]
    # L1: R_k * prod_{j=1}^{k+1}(E + j z) has degree <= k (coefficients z^{k+1}..z^Z vanish), and R_k(oo)-> E delta_{k0}
    for k in range(0, min(T, Z - 1) + 1):
        pol = zlin(1, 0)
        for j in range(1, k + 2): pol = zmul(pol, zlin(E, j))
        prod = zmul(Rk(k), pol)
        top = max([m for m in range(Z + 1) if prod[m] != 0], default=-1)
        ok = (top <= k) if k >= 1 else (top <= 1)
        print('   L1 k=%d: [t^k](zG) prod_{j<=k+1}(E+jz) has z-degree %d (bound %d), tested to z^%d:' % (k, top, k if k else 1, Z), ok)
    # A_j closed form: A_j = (-1)^j C_j t^-1 zeta^j 2F1(jB/E, j a1/E; j; tau), C_j = P_{j-1}/(j!(j-1)! E^(2j-4))
    zeta = pmul(TAU, pexp(psc(1 / E, D.dphi)))
    def hyp(a, b, c):
        out = ps([1]); term = Fr(1)
        for n in range(1, T + 1):
            term *= (a + n - 1) * (b + n - 1) / ((c + n - 1) * n)
            out[n] = term
        return out
    Aj_t = {}
    for j in range(1, T + 2):
        Pk = Fr(1)
        for i in range(1, j): Pk *= (j * a2 - i * E) * (j * A - i * E)
        Cj = Pk / (factorial(j) * factorial(j - 1) * E ** (2 * j - 4))
        f_tau = psc((-1) ** j * Cj, pmul(ppow(zeta, j), hyp(j * B / E, j * a1 / E, j)))
        v, g = tau_to_t(f_tau, 0)
        Aj_t[j] = (-1, g)   # times t^-1
    okA = True
    for k in range(0, T):
        target = Rk(k)
        if k == 0: target = [x - (E if m == 0 else 0) for m, x in enumerate(target)]
        s = [Fr(0)] * (Z + 1)
        for j in range(1, k + 2):
            s = [x + y for x, y in zip(s, [lcoef(Aj_t[j], k) * c for c in zinv(zlin(E, j))])]
        okA &= (s == target)
    print('   A_j closed form reproduces [t^k](zG) - E delta_k0 for k <= %d:' % (T - 1), okA)
    # first correction and rho_1
    R1 = Rk(1); u = D.u
    tgt = [u / 2 * (x - y) for x, y in zip(zinv(zlin(E, 1)), zinv(zlin(E, 2)))]
    print('   [t^1](zG) == (u/2)(1/(E+z) - 1/(E+2z)):', R1 == tgt)
    rho_t = pmul(pcompose(rho_tau, TAU_T), TAU_OVER_T)
    print('   rho_1 == u/(2E^2):', rho_t[1] == u / (2 * E * E))
    # ---------- H-row from the definition, Row 2, K_1, E_11 ----------
    # frame-B rows as t-Laurent series: r0[m] = (0, series); L_1 = A/(tau s*) - a0 ; Delta_1 = a1 A e1e2 wr_2/tau^2
    r0 = {m: tau_to_t(rowB[m], 0) for m in range(Z + 1)}
    L1 = ladd(tau_to_t(psc(A, pinv(D.s)), -1), (0, ps([-a0])))
    Delta1 = tau_to_t(psc(a1 * A * D.e1 * D.e2, D.wr[2]), -2)
    dlog = lsc(Fr(1, 2), lmul(lder(Delta1), linv(Delta1)))
    r1 = {0: lmul(L1, r0[0])}
    for m in range(1, Z + 1):
        r1[m] = ladd(lmul(L1, r0[m]), ladd(lder(r0[m - 1]), lsc(-1, lmul(dlog, r0[m - 1]))))
    rho_l = (0, rho_t)
    # F2: rho r_{1,1} = E/t + a1 a2/(E+z) + O(t)
    geoz = zinv(zlin(E, 1))
    f2ok = True
    for m in range(Z + 1):
        v = lmul(rho_l, r1[m])
        f2ok &= (lcoef(v, -1) == (E if m == 0 else 0)) and (lcoef(v, 0) == a1 * a2 * geoz[m])
    print('   F2: rho r_{1,1} = E/t + a1a2/(E+z) + O(t) (all z^m <= %d):' % Z, f2ok)
    # Row 2: (L1 + z D1) r_{1,1} + (a1a2/tau) r_{0,1} - (E/tau - a0) r_{1,1} = 0 ; 1/tau as t-series
    inv_tau = tau_to_t(ONE, -1)
    r2 = {0: lmul(L1, r1[0])}
    for m in range(1, Z + 1):
        r2[m] = ladd(lmul(L1, r1[m]), ladd(lder(r1[m - 1]), lsc(-1, lmul(dlog, r1[m - 1]))))
    row2ok = True
    for m in range(Z + 1):
        expr = ladd(r2[m], ladd(lmul(lsc(a1 * a2, inv_tau), r0[m]), lsc(-1, lmul(ladd(lsc(E, inv_tau), (0, ps([-a0]))), r1[m]))))
        # check low orders only (precision loss from derivatives): t^-3 .. t^(T-4)
        for k in range(-3, T - 3):
            if lcoef(expr, k) != 0: row2ok = False
    print('   Row 2 of the QDE at U_1 (z^m <= %d, t^-3..t^%d):' % (Z, T - 4), row2ok)
    # K_1 = <v_1(z), v_1(-z)>/Delta_1 : pairing e1e2[a0a1a2 a0b0 + a1a2(a0b1+a1b0) - E a1b1]
    e12 = D.e1 * D.e2
    def pairz(p, q, sgn):   # <v[p], v[q]> with v[q] evaluated at -z (sign (-1)^q)
        a = (r0[p], r1[p]); b = (r0[q], r1[q])
        val = ladd(lsc(e12 * a0 * a1 * a2, lmul(a[0], b[0])), ladd(lsc(e12 * a1 * a2, ladd(lmul(a[0], b[1]), lmul(a[1], b[0]))), lsc(-e12 * E, lmul(a[1], b[1]))))
        return lsc(sgn, val)
    K1ok = True
    for mm in range(0, Z + 1):
        tot = None
        for p in range(mm + 1):
            q = mm - p
            term = pairz(p, q, (-1) ** q)
            tot = term if tot is None else ladd(tot, term)
        K = lmul(tot, linv(Delta1))
        for k in range(0, T - 3):
            if lcoef(K, k) != (1 if (mm == 0 and k == 0) else 0): K1ok = False
    print('   K_1 = <v_1(z),v_1(-z)>/Delta_1 == 1 (w^m <= %d, t^0..t^%d):' % (Z, T - 4), K1ok)
    # E_11(0,0) = <v[1], v[0]> ; E_11(1,0) = P_20 - P_11 with P_pq = -(-1)^(p+q) <v[p],v[q]>
    E00 = pairz(1, 0, 1)
    E10 = lsc(-1, pairz(2, 0, 1))
    tgt00 = -(e12 * a1 * a2 * A * B / 2) / E ** 2
    print('   E_11(0,0): no negative t-powers:', all(lcoef(E00, k) == 0 for k in range(-3, 0)), '; [t^0] == -e1e2 a1a2 A B/(2E^2):', lcoef(E00, 0) == tgt00)
    print('   E_11(1,0): no negative t-powers:', all(lcoef(E10, k) == 0 for k in range(-3, 0)), '; [t^0] == -e1e2 a1a2 A B/(2E^3):', lcoef(E10, 0) == tgt00 / E)
    # Theta first order via Delta_2/(L_2 - L_1): Theta_1 = -e1e2 a1a2 A B/E^2
    L2t = tau_to_t(D.L2, 0)
    Theta = lsc(e12, ladd((0, ps([a1 * a2])), lsc(-E, L2t)))
    print('   Theta(0) = 0, Theta_1 == -e1e2 a1a2 A B/E^2:', lcoef(Theta, 0) == 0 and lcoef(Theta, 1) == -e12 * a1 * a2 * A * B / E ** 2)
    # E_11 flow (Thm 7.2_CY): d/dt [t^2 rho^2 E_11(m,n)] == Theta [psi^m]Dg [psi'^n]Dg, Dg = d/dt (t g), [psi^m] g = (-1)^m G[m]
    def E11(m, n):
        tot = None
        for k in range(n + 1):
            p, q = m + 1 + k, n - k
            term = lsc(-(-1) ** (p + q) * (-1) ** k, pairz(p, q, 1))
            tot = term if tot is None else ladd(tot, term)
        return tot
    t2rho2 = lmul((2, ps([1])), lmul(rho_l, rho_l))
    Dg = {m: lder(lmul((1, ps([1])), (0, [(-1) ** m * c for c in Gt[m]]))) for m in range(Z + 1)}
    flowok = True
    for m in range(0, 4):
        for n in range(0, 4 - m):
            lhs = lder(lmul(t2rho2, E11(m, n)))
            rhs = lmul(Theta, lmul(Dg[m], Dg[n]))
            for k in range(0, T - 4):
                if lcoef(lhs, k) != lcoef(rhs, k): flowok = False
    print('   E_11 flow d/dt(t^2 rho^2 E_11) = Theta Dg Dg (m+n <= 3, t^0..t^%d):' % (T - 5), flowok)
    # L3 at low order: [t^k](rho^2 E_11)(psi,psi') prod_{i<=k+1}(E - i psi) prod_{l<=k+1}(E - l psi') has bidegree <= (k,k), tested on m,n <= 4 for k = 0,1,2
    rho2 = lmul(rho_l, rho_l)
    Em = {(m, n): lmul(rho2, E11(m, n)) for m in range(5) for n in range(5)}
    for k in range(0, 3):
        pol = [Fr(1)]
        for i in range(1, k + 2):
            new = [Fr(0)] * (len(pol) + 1)
            for d, c in enumerate(pol): new[d] += E * c; new[d + 1] -= i * c
            pol = new
        bad = 0
        for m in range(5):
            for n in range(5):
                if m <= k and n <= k: continue
                s = sum(pol[i] * pol[j] * lcoef(Em[(m - i, n - j)], k) for i in range(len(pol)) for j in range(len(pol)) if m - i >= 0 and n - j >= 0)
                if s != 0: bad += 1
        print('   L3 k=%d: [t^k](rho^2 E_11) x annihilators has bidegree <= (k,k) on psi and psi-prime orders <= 4: violations %d' % (k, bad))
    return rowB

samples = [(Fr(5, 3), Fr(-4, 7), Fr(13, 6), Fr(9, 11)), (Fr(-2, 5), Fr(7, 4), Fr(3), Fr(1, 6))]
for a0, a1, a2, e1 in samples:
    e2 = a0 + a1 + a2 - e1
    D = CYData(a0, a1, a2, e1, e2)
    key_formula_check(D, 'CY sample')
# a0 = 0 rescaling: frame-B row at (0, a1, a2, e1, e2) vs fibrewise row at x = a1/a2 with z -> z/a2
a1, a2, e1 = Fr(5, 2), Fr(-3, 7), Fr(4, 9); e2 = a1 + a2 - e1
D0 = CYData(0, a1, a2, e1, e2); x = a1 / a2
DF = CYData(0, x, 1, x, 1)
hB0 = h_class([a1, a2, -(a1 + a2)]); hBF = h_class([x, 1, -(1 + x)])
ok = True
for m in range(Z + 1):
    rowm0 = [sum(hB0[j] * D0.rowpt[m - j][n] for j in range(m + 1)) for n in range(T + 1)]
    rowmF = [sum(hBF[j] * DF.rowpt[m - j][n] for j in range(m + 1)) for n in range(T + 1)]
    ok &= (rowm0 == [c / a2 ** m for c in rowmF])
print('== a0 = 0: frame-B row at (0,a1,a2,e1,e2) == fibrewise frame-B row at x = a1/a2 with z -> z/a2 (tau^%d, z^%d):' % (T, Z), ok)
print('   also h_(a1,a2,-E)(z) == h_(c,-x,-1)(z/a2)^-1:', hB0 == [c / a2 ** m for m, c in enumerate(zinv(h_class([1 + x, -x, -1])))])
