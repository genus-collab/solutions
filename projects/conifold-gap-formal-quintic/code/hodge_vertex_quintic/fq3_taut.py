# -*- coding: utf-8 -*-
"""fq3_taut.py -- the tautological u-correction of the symbol: R~(zhat; u) = colour-colour form of the hatted rows at all orders in u,
R~_{ji}(zhat; u) = Sum_ab Qhat_{b,0}^{(j)} eta^{ba} Qhat_{a,n}^{(i)} zhat^n, M_0 = R~(zhat; 0); with M(zhat; u) := R~(zhat; u) mod zhat^7 the
primed rows R' = R^ M^{-1} are 1 + O(zhat^7) exactly in u, so every edge graph vanishes identically (the trivial factorisation).
Its first u-correction N_1^taut := [u^1] log(M_0^{-1} R~(zhat; u)) = M_0^{-1} [u^1] R~ (right placement) is a Z/5-covariant element of sp
and MUST solve the [u^1] systems of genus 2 and genus 3 (every edge graph is O(u^2) for it): a consistency check of J and g1 against
an independently known solution.  taut_N1(st) -> (N_1 as sp series, checks dict).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_symbol import mat_series_mul
from fq_sser import mat_mul, mat_eq
from fq3_m import in_sp, sp_params, NSP, class_columns, KMAX


def colour_series_u(H, k):
    """[u^k] of the colour matrices C_n(u), n = 0..6 (list of 7 matrices)"""
    R = H.R; p = R.p; out = []
    for n in range(KMAX + 1):
        C = [[0] * 5 for _ in range(5)]
        for j in range(5):
            for i in range(5):
                acc = 0
                for (a, b), e in H.eta_inv.items():
                    # [u^k] of Qhat_{b,0}^{(j)} Qhat_{a,n}^{(i)}
                    s = H.Q[b][j][0].mul(H.Q[a][i][n])
                    acc = (acc + e * (s.coef(k) if k < s.prec else 0)) % p
                C[j][i] = acc
        out.append(C)
    return out


def taut_N1(st):
    R = st['R']; p = R.p; H = st['H']; M0 = st['M0']
    C0 = colour_series_u(H, 0); C1 = colour_series_u(H, 1)
    ok0 = all(mat_eq(C0[n], M0.M[n], p) for n in range(KMAX + 1))
    N1 = mat_series_mul(M0.Minv, C1, p, KMAX)           # [u^1] of M_0^{-1} R~(u) = N_1 (right placement) since exp(u N_1 + ..) = M_0^{-1} R~
    return N1, dict(C0_is_M0=ok0, N1_in_sp=in_sp(N1, p), N1_0_zero=not any(x for r in N1[0] for x in r))


def circulant_coords(R, N):
    """coordinates of N in the circulant class (least squares-free: solve A t = n exactly); None if N is not circulant"""
    from fq_mfit import fp_solve
    p = R.p; cols, names = class_columns(R, 'circ'); n = sp_params(N, p)
    A = [[cols[j][v] for j in range(len(cols))] for v in range(NSP)]
    sol = fp_solve(A, n, p)
    return (sol[0], names) if sol[0] is not None else (None, names)


if __name__ == "__main__":
    import json
    from fractions import Fraction as Fr
    from fq_fp import P2, P1
    from fq3_par import setup
    from fq3_rec import Rec
    for prime in (P2, P1):
        st = setup((prime, "0", "0", 34, 6)); R = st['R']; p = R.p; rec = Rec(R)
        N1, chk = taut_N1(st); print("p = 2^%d-1:" % (89 if prime == P2 else 61), chk)
        t, names = circulant_coords(R, N1); print("   N_1^taut circulant coordinates:", None if t is None else dict(zip(names, [rec.fmt(x) for x in t])))
        tag = "89" if prime == P2 else "61"
        path = os.path.join(HERE, "g3_state_%s.json" % tag)
        d = json.load(open(path)) if os.path.exists(path) else {}
        if 'J' not in d or 'J2' not in d or 'A_M0' not in d: print("   (state file for this prime incomplete: skipping the equation check)"); continue
        J3 = {int(i): v for i, v in d['J'].items()}; A0 = {int(i): v for i, v in d['A_M0'].items()}
        J2 = {int(i): v for i, v in d['J2'].items()}
        n = sp_params(N1, p)
        for (name, row, rhs) in (("genus-2 graph 2", J2[2], (-R.of(Fr(1, 750))) % p), ("genus-2 graph 4", J2[4], (-R.of(Fr(-13, 8000))) % p),
                                 ("genus-3 graph 2", J3[2], (-A0[2][1]) % p), ("genus-3 graph 5", J3[5], (-A0[5][1]) % p)):
            lhs = sum(row[v] * n[v] for v in range(NSP)) % p
            print("   %-16s J . N_1^taut = %s   required -g1 = %s   %s" % (name, rec.fmt(lhs), rec.fmt(rhs), "OK" if lhs == rhs else "MISMATCH"))
