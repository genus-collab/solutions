# -*- coding: utf-8 -*-
"""ur_model_check.py -- UNIFORM_REDUCTION.md, checks (A) and (B).

(A) The flat-frame quantum differential system of the formal quintic is CYCLIC (indices mod 5):
        (z D + zeta^i L) Q_a^{(i)} = c_a Q_{a+1}^{(i)},   c = (C_1, L^5/(C_0^2 C_1^2), C_1, C_0, C_0),   D = q d/dq = -(1-u) d/du,
    including the closing relation a = 4 -> 0 (which fq_rows does not use).  Checked on the conifold-frame s-series rows of
    fq_sser.SRows (C_0 -> varpi_0^cf, C_1 -> D t_c), all colours, z-orders 1..6.
(B) The conifold symbol R^ is the homogeneous WKB solution of the MODEL system  z d/dv Y_b = b^mod_b Y_{b+1},
    b^mod = (1, -1/v, 1, -1, -1)  (the u -> 0, zhat fixed limit of (A); it is the Gauss-Manin system of lambda = -(x^4/4 + v/x),
    primitive form dx, in the basis dx, dx/x, v dx/x^2, x^2 dx, x dx up to constants; scalar equation v th^5 + 3 z th^4 + 1 = 0).
    Ansatz Y_b = exp(U/z) sigma^(c_b - 1) F_b(t), sigma = v^(1/5), t = z/sigma^4, U = -(5/4) zeta^i sigma^4:
        z d/dv [sigma^k e^(U/z) F] = sigma^(k-1) e^(U/z) T_k F,   T_k = -zeta^i + (t/5)(k - 4 t d/dt),
        F_1 = T_{-1} F_0, F_2 = -T_{-2} F_1, F_3 = T_2 F_2, F_4 = -T_1 F_3, F_0 = -T_0 F_4.
    The series F_b are computed over Q for colour 0 (F_0(0) = 1) and compared with sym.Rh[b][i][n] for all rows, colours (by the
    covariance R^_b^(i)(t) = zeta^(i(b+1)) R^_b^(0)(zeta^-i t)) and orders n <= 6, at both primes.
Run: py -I ur_model_check.py        (single process, about a minute per prime)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq_series import Ser
from fq3_par import setup
CEXP = (0, -1, 3, 2, 1)
OUT = open(os.path.join(HERE, "ur_model_check.log"), "w", encoding="utf-8")
FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


# ------------------------------------------------------------------ (B) the model series over Q (colour 0)
def T(kappa, F):
    """T_kappa F = -F + (t/5)(kappa - 4 t d/dt) F on a coefficient list (same length, truncated)"""
    out = [-x for x in F]
    for n in range(1, len(F)): out[n] += Fr(kappa - 4 * (n - 1), 5) * F[n - 1]
    return out


def closure(F):
    """F + T_0 T_1 T_2 T_{-2} T_{-1} F  (must vanish for the solution F_0)"""
    G = T(0, T(1, T(2, T(-2, T(-1, F)))))
    return [a + b for a, b in zip(F, G)]


def model_rows(N):
    """F_0..F_4 as coefficient lists t^0..t^N; F_0 solved from the closure equation order by order"""
    L = N + 3
    F = [Fr(0)] * L; F[0] = Fr(1)
    for n in range(1, N + 1):
        mono = [Fr(0)] * L; mono[n] = Fr(1)
        e = closure(mono)[n + 1]                    # coefficient of f_n in [t^(n+1)] of the closure
        assert e == -4 * n, (n, e)
        c = closure(F)[n + 1]
        F[n] = -c / e
    res = closure(F)
    assert all(res[k] == 0 for k in range(N + 2)), res[:N + 2]
    F0 = F
    F1 = T(-1, F0); F2 = [-x for x in T(-2, F1)]; F3 = T(2, F2); F4 = [-x for x in T(1, F3)]
    back = [-x for x in T(0, F4)]
    assert all(back[k] == F0[k] for k in range(N + 1))
    return [r[:N + 1] for r in (F0, F1, F2, F3, F4)]


def stationary_phase(N):
    """independent of the ODE: formal stationary phase of Int exp(-(x^4/4 + v/x)/z) dx at v = 1, critical point x = 1 (-> colour 0).
    lambda(1+y) - lambda(1) = -(5/2) y^2 - Sum_{k>=3} a_k y^k; Gaussian moments <y^(2m)> = (2m-1)!! (z/5)^m.  Returns the series in z."""
    K = 2 * N + 2
    # -(x^4/4 + 1/x) at x = 1 + y: -(1/4)(1+y)^4 - Sum_k (-1)^k y^k
    lam = [Fr(0)] * (K + 3)
    binom = [1, 4, 6, 4, 1]
    for k in range(5): lam[k] += Fr(-binom[k], 4)
    for k in range(K + 3): lam[k] += -Fr((-1) ** k)
    assert lam[1] == 0 and lam[2] == Fr(-5, 2), lam[:4]
    # exp( Sum_{k>=3} lam_k y^k / z ) as a polynomial in y and w = 1/z, graded so that y ~ z^(1/2): weight (k-2)/2 per vertex
    # E[j][m] = coefficient of y^j w^m; term y^k w contributes z-order (k-2)/2; keep total order <= N: j - 2m <= 2N
    E = {(0, 0): Fr(1)}
    cur = {(0, 0): Fr(1)}; fact = 1
    for r in range(1, 2 * N + 1):
        new = {}
        for (j, m), c in cur.items():
            for k in range(3, K + 3):
                if lam[k] == 0: continue
                jj = j + k; mm = m + 1
                if jj - 2 * mm > 2 * N: break
                new[(jj, mm)] = new.get((jj, mm), 0) + c * lam[k]
        fact *= r; cur = new
        if not cur: break
        for key, c in cur.items(): E[key] = E.get(key, 0) + c / fact
    out = [Fr(0)] * (N + 1)
    for (j, m), c in E.items():
        if j % 2: continue
        h = j // 2; n = h - m
        if n < 0 or n > N: continue
        dd = 1
        for q in range(1, j, 2): dd *= q                 # (j-1)!!
        out[n] += c * dd * Fr(1, 5) ** h
    return out


# ------------------------------------------------------------------ (A) the cyclic system on the s-series rows
def cyclic_check(st):
    R = st['R']; p = R.p; S = st['S']; C = st['C']; Nz = S.Nz
    C0 = C.to_s(C.varpi0); C1 = C.to_s(C.Dop(C.tc))                 # C_0 -> varpi_0^cf,  C_1 = D(I_1/I_0) -> D t_c
    L5 = Ser.mono(R, -5, 1, C0.prec - 5)
    c = [C1, L5.mul(C0.mul(C0).mul(C1).mul(C1).inv()), C1, C0, C0]
    Linv = lambda prec: Ser.mono(R, -1, 1, prec)
    bad = []; minrel = None
    for a in range(5):
        for i in range(5):
            for n in range(1, Nz + 1):
                Qn = S.Q[a][i][n]
                lhs = C.Ds(S.Q[a][i][n - 1]).add(Qn.mul(Linv(Qn.prec + 4)).scale(R.z(i)))
                rhs = c[a].mul(S.Q[(a + 1) % 5][i][n])
                d = lhs.sub(rhs)
                lead = min(lhs.v, rhs.v)
                rel = d.prec - lead                                 # number of s-orders on which the identity is tested
                minrel = rel if minrel is None else min(minrel, rel)
                if d.c: bad.append((a, i, n, d.v))
    return bad, minrel


def main():
    N = 10
    F = model_rows(N); sp = stationary_phase(N)
    log("model identity row F_0(t) = " + " + ".join("(%s) t^%d" % (F[0][n], n) for n in range(N + 1)))
    check("stationary phase of Int exp(-(x^4/4 + v/x)/z) dx (critical point x = 1, v = 1) = F_0 from the cyclic system, through t^%d (exact, Q)" % N, sp == F[0],
          [str(x) for x in sp[:4]])
    f_known = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000), Fr(52989974037, 25600000000)]
    check("F_0 = the identity row f of the conifold symbol as recorded in CONIFOLD_SYMBOL_genus2.md (through t^6, exact, Q)", F[0][:7] == f_known)
    log("leading coefficients F_b(0), b = 0..4: %s" % [str(F[b][0]) for b in range(5)])
    for b in range(5): log("   F_%d(t)/F_%d(0) = %s" % (b, b, " + ".join("(%s) t^%d" % (F[b][n] / F[b][0], n) for n in range(7))))
    for prime, tag in ((P2, "2^89-1"), (P1, "2^61-1")):
        t0 = time.time(); st = setup((prime, "0", "0", 34, 6)); R = st['R']; p = R.p; sym = st['sym']
        log("-" * 100); log("p = %s (setup %.0f s)" % (tag, time.time() - t0))
        bad, minrel = cyclic_check(st)
        check("(A) cyclic system (zD + zeta^i L) Q_a = c_a Q_{a+1}, a in Z/5 (incl. 4 -> 0), all colours, z^1..z^6, s-series (>= %d s-orders each)" % minrel, not bad, bad[:5])
        ok_all = True; badB = []
        for b in range(5):
            for i in range(5):
                for n in range(7):
                    pred = R.of(F[b][n]) * R.z(i * (b + 1 - n)) % p
                    if pred != sym.Rh[b][i][n] % p: ok_all = False; badB.append((b, i, n))
        check("(B) conifold symbol R^_b^(i)(zhat) = zeta^(i(b+1)) F_b(zeta^-i zhat) for ALL rows b, colours i, orders n <= 6 (model = Gauss-Manin system of -(x^4/4 + v/x))", ok_all, badB[:5])
        H = st['H']
        okH = all(H.Q[b][i][n].coef(0) == sym.Rh[b][i][n] for b in range(5) for i in range(5) for n in range(7))
        check("hatted rows are pure u-series (asserted by HatRowsU) with u^0 term = the symbol", okH)
    log("UR_MODEL_CHECK " + ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); OUT.close()


if __name__ == "__main__":
    main()
