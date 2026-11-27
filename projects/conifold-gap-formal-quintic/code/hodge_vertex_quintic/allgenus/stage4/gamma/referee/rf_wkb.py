# -*- coding: utf-8 -*-
"""rf_wkb.py -- REFEREE: the colour-0 WKB series of the formal quintic in the FLAT variable L, from the Picard-Fuchs operator
applied literally (no hatted variables, no (HPF), no worker or project code).

PF: [(zD)^5 - 1 - q prod_{k=1..5}(5 zD + k z)] psi = 0,  D = q d/dq,  3125 q = 1 - L^-5,  L = (1 - 3125 q)^(-1/5),  so on functions of L
D = delta := ((L^6 - L)/5) d/dL.  Ansatz psi = e^(u_0/z) L f,  D u_0 = L,  f = Sum_n R_n(L) z^n,  R_0 = 1.  Conjugation: zD -> M = L + z delta.
T_j(R) := [z^j] P(L R) (P = conjugated operator, z-degree <= 5).  Derived (and asserted numerically below): T_0 = 0 (symbol),
T_1(R) = (L^6 - L) dR/dL (transport).  Order z^(n+1):  (L^6 - L) R_n' = -Sum_{j=2..5} T_j(R_(n+1-j)).  R_n is fixed up to an additive
constant (the L^0 coefficient).  Convention: the constant is 0 unless 5 | n; for n = 5m it is the free constant x_m, fixed by the
q = 0 (QRR) rule R_n(1) = h_n, or prescribed (controls).  (H) <=> R_n(1) = h_n automatically for 5 not| n ("h_auto").
Laurent polynomials in L are dicts {exponent: coefficient}.
h(z) = exp(-Sum_k b_k z^(2k-1) Sum_{a in w_0} a^(1-2k)),  w_0 = {1 - zeta^j (j = 1..4)} U {-5},  computed EXACTLY in Q(zeta_5)
(trace of (1 - zeta)^-m) and then reduced to the field.
"""
from fractions import Fraction as Fr
from rf_ring import FQ, stirling_b, cz_inv, cz_mul, cz_trace


def lp_add(F, a, b, c=1):
    out = dict(a)
    for e, x in b.items():
        out[e] = F.red(out.get(e, F.zero) + c * x)
    return {e: x for e, x in out.items() if x != F.zero}


def lp_shift(a, k): return {e + k: x for e, x in a.items()}


def lp_scal(F, a, c): return {e: F.red(c * x) for e, x in a.items() if F.red(c * x) != F.zero}


def lp_delta(F, a):
    """((L^6 - L)/5) d/dL on L^e: (e/5)(L^(e+5) - L^e)"""
    i5 = F.inv(F.of(5)); out = {}
    for e, x in a.items():
        if e:
            c = F.red(e * x * i5)
            out[e + 5] = F.red(out.get(e + 5, F.zero) + c); out[e] = F.red(out.get(e, F.zero) - c)
    return {e: x for e, x in out.items() if x != F.zero}


def lp_q(F, a):
    """q a = (1 - L^-5) a / 3125"""
    return lp_scal(F, lp_add(F, a, lp_shift(a, -5), -1), F.inv(F.of(3125)))


def lp_eval(F, a, L):
    Lr = F.of(L)
    if F.p is None: return sum((x * Fr(L) ** e for e, x in a.items()), F.zero)
    return F.red(sum(x * (pow(Lr, e, F.p) if e >= 0 else pow(F.inv(Lr), -e, F.p)) for e, x in a.items()))


def _M(F, G):
    """M = L + z delta on a z-polynomial G (list of Laurent dicts, index = z-power)"""
    out = [lp_shift(g, 1) for g in G] + [{}]
    for j, g in enumerate(G): out[j + 1] = lp_add(F, out[j + 1], lp_delta(F, g))
    return out


def T_ops(F, R):
    """[T_0(R), .., T_5(R)], T_j = [z^j] of the conjugated PF operator applied to L R"""
    g = [lp_shift(R, 1)]; A = g
    for _ in range(5): A = _M(F, A)
    B = g
    for k in range(1, 6):
        MB = _M(F, B); B = [lp_add(F, lp_scal(F, MB[j], F.of(5)), B[j - 1] if 1 <= j <= len(B) else {}, k) if j > 0
                            else lp_scal(F, MB[0], F.of(5)) for j in range(len(MB))]
    out = []
    for j in range(6):
        t = lp_add(F, A[j] if j < len(A) else {}, lp_q(F, B[j] if j < len(B) else {}), -1)
        if j == 0: t = lp_add(F, t, g[0], -1)
        out.append(t)
    return out


def qrr_h(F, N):
    """[h_0..h_N] (colour 0), exact via Q(zeta_5), then reduced into F"""
    t = cz_inv((Fr(1), Fr(-1), Fr(0), Fr(0)))                           # 1/(1 - zeta)
    lg = [Fr(0)] * (N + 1); tp = (Fr(1), Fr(0), Fr(0), Fr(0))
    for m in range(1, N + 1):
        tp = cz_mul(tp, t)
        if m % 2 == 1:
            k = (m + 1) // 2
            lg[m] = -stirling_b(k) * (cz_trace(tp) + Fr(-5) ** (-m))       # Sum_{a in w_0} a^(-m), m = 2k - 1
    out = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        out[n] = sum(k * lg[k] * out[n - k] for k in range(1, n + 1)) / n
    return [F.of(x) for x in out]


def wkb(F, N, h=None, free=None, checks=True):
    """R_n(L), n <= N.  Returns (R, info); info: h_auto (5 not| n), free constants x_m, min exponents, transport check."""
    h = h if h is not None else qrr_h(F, N)
    R = [{0: F.one}]; T = [T_ops(F, R[0])]; info = dict(h_auto={}, free={}, minexp=[0], T01=True)
    for n in range(1, N + 1):
        S = {}
        for j in range(2, 6):
            if n + 1 - j >= 0: S = lp_add(F, S, T[n + 1 - j][j])
        S = lp_scal(F, S, F.of(-1))                                       # (L^5 - 1) L R_n' = S
        Q = {}; emax = max(S) if S else 0; emin = min(S) if S else 0
        for e in range(emax, emin + 4, -1):
            Q[e - 5] = F.red(S.get(e, F.zero) + Q.get(e, F.zero))
        for e in range(emin, emin + 5):
            assert F.red(S.get(e, F.zero) + Q.get(e, F.zero)) == F.zero, "not divisible by L^5 - 1 at n = %d" % n
        Rn = {}
        for e, x in Q.items():                                            # R_n' = Q/L: L^(e-1) -> L^e/e
            if x == F.zero: continue
            assert e != 0, "logarithm at n = %d" % n
            Rn[e] = F.red(x * F.inv(F.of(e)))
        if n % 5 == 0:
            m = n // 5
            x = free[m] if (free is not None and m in free) else F.red(h[n] - lp_eval(F, Rn, 1))
            Rn = lp_add(F, Rn, {0: x}); info['free'][m] = x
        else:
            info['h_auto'][n] = (lp_eval(F, Rn, 1) == h[n])
        R.append(Rn); info['minexp'].append(min(Rn) if Rn else 0)
        T.append(T_ops(F, Rn))
        if checks and n <= 12:                                            # T_0 = 0 and T_1 = (L^6 - L) d/dL, asserted
            dR = {e - 1: F.red(e * x) for e, x in Rn.items() if e}
            info['T01'] &= (not T[-1][0]) and T[-1][1] == lp_add(F, lp_shift(dR, 6), lp_shift(dR, 1), -1)
    return R, info


def R_at(F, R, L):
    return [lp_eval(F, r, L) for r in R]


def ph_poly(R, n):
    """hatted polynomial ph_n(u) = Sum_k [L^(4n-5k)] R_n u^k (exponents must be = 4n mod 5)"""
    assert all((4 * n - e) % 5 == 0 for e in R[n]), "exponent not = 4n mod 5"
    out = {}
    for e, x in R[n].items(): out[(4 * n - e) // 5] = x
    return [out.get(k, 0) for k in range(max(out) + 1)]
