# check3.py -- evidence only.
# T5  Prop MV2: Q_g(s;R) from the definition (explicit dilaton leaves) = (1+s)^-(2g-2+r) e^{-uJ} G_R(u e^{-u}), string/part-0 operators,
#     Laguerre formula; real and non-integer leaves.
# T7  Xi_h(R) from the Schur side alone (formula (STAR)) against the definition. Lemma A range: zeros; controls just outside.
# T8  The Mathematica tables of lemmaA/data (independent code), including c-normalisation.
# T9  E-scaling (7.2), string and lam_0 identities with E.
from fractions import Fraction as Fr
from math import factorial, comb
import re, sys, time, os
from engine import *
from schur import *


def sig(t):
    return t * (1 + t)


def lag(n, al, x):
    return sum(Fr((-1) ** i) * gbin(n + al, n - i) * x ** i / factorial(i) for i in range(n + 1))


def gbin(top, k):
    if k < 0:
        return Fr(0)
    r = Fr(1)
    for i in range(k):
        r *= Fr(top - i)
    return r / factorial(k)


def series_from_M(Ms, n0, J, m):
    """(1+s)^{-n0} e^{-uJ} G(u e^{-u}), truncated s^m; Ms[a] = M_a. Expanded as Sum_a s^a/a! (1+s)^{-(n0+a)} exp(-s(J+a)/(1+s)) M_a."""
    u = [Fr(0)] + [Fr((-1) ** (i - 1)) for i in range(1, m + 1)]
    out = [Fr(0)] * (m + 1)
    for a in range(0, m + 1):
        # (1+s)^{-(n0+a)}
        w = [Fr(1)]
        for i in range(1, m + 1):
            w.append(w[-1] * (-(n0 + a) - i + 1) / i)
        e = [Fr(1)]
        d = -(J + a)
        for n in range(1, m + 1):
            e.append(sum(i * d * u[i] * e[n - i] for i in range(1, n + 1)) / n)
        pr = smul(w, e, m)
        for i in range(m + 1 - a):
            out[i + a] += pr[i] * Ms[a] / factorial(a)
    return out


def T5():
    print("T5: Prop MV2")
    ok = bad = 0
    cases = [(1, (2,)), (1, (Fr(5, 2),)), (1, (3, 2)), (2, ()), (2, (2,)), (2, (3,)), (2, (Fr(1, 2),)), (2, (2, 2)), (2, (2, 3)),
             (2, (Fr(1, 2), Fr(5, 2))), (2, (1,)), (2, (1, 2)), (2, (0,)), (2, (0, 3)), (0, (2, 3, 4)), (0, (1, 2, 3, 2)), (3, (2,)), (3, ())]
    extra = {3: 1}
    for tau in (Fr(2, 7), Fr(-7, 3)):
        ts = ts_of(tau)
        for g, R in cases:
            r = len(R)
            N = 3 * g - 3 + r
            m = N + extra.get(g, 3)
            J = sum(R)
            Ms = []
            for a in range(m + 1):
                if g == 0 or 2 * g - 2 + r + a > 0:
                    Ms.append(leg_integral(g, tuple(R) + (Fr(1),) * a, ts))
            q = Qg_coeffs(g, tuple(R), ts, kmax=m)
            ser = series_from_M(Ms, 2 * g - 2 + r, J, m)
            b = sum(1 for x, y in zip(q, ser) if x != y)
            ok += m + 1 - b
            bad += b
            # Laguerre coefficient formula and value at s=-1
            n0 = 2 * g - 3 + r
            for mm in range(m + 1):
                val = sum(Ms[a] / factorial(a) * (-1) ** (mm - a) * lagr(mm - a, n0 + a, -(J + a)) for a in range(mm + 1))
                if val != q[mm]:
                    bad += 1
                else:
                    ok += 1
            xi = sum(c * (-1) ** k for k, c in enumerate(q))
            xi2 = sum(Fr((-1) ** a, factorial(a)) * lagr(N - a, 2 * g - 2 + r + a, -J - a) * Ms[a] for a in range(N + 1))
            if xi != xi2:
                bad += 1
                print("Xi Laguerre mismatch", g, R, tau, xi, xi2)
            else:
                ok += 1
    print(f"  dilaton-leaf elimination + Laguerre: {ok} checks pass, {bad} fail")
    # operators for parts 1 and 0 (polynomial identities in s)
    ok = bad = 0
    for tau in (Fr(2, 7),):
        ts = ts_of(tau)
        for g, R in ((2, ()), (2, (2,)), (2, (3,)), (2, (2, 3)), (1, (4,)), (3, (2,))):
            r = len(R)
            J = sum(R)
            n = 2 * g - 2 + r
            q = Qg_coeffs(g, R, ts, kmax=3 * g - 3 + r + 1 + 2)
            q1 = Qg_coeffs(g, R + (1,), ts, kmax=len(q) - 1)
            q0 = Qg_coeffs(g, R + (0,), ts, kmax=len(q) - 1)
            L = len(q)
            dq = [(i + 1) * q[i + 1] if i + 1 < L else Fr(0) for i in range(L)]
            for i in range(L - 2):
                # (1+s)^2 Q' + [J + (1+s) n] Q
                a = dq[i] + (2 * dq[i - 1] if i >= 1 else 0) + (dq[i - 2] if i >= 2 else 0) + J * q[i] + n * q[i] + (n * q[i - 1] if i >= 1 else 0)
                # (J + s n) Q + s(1+s) Q'
                b = J * q[i] + (n * q[i - 1] if i >= 1 else 0) + (i * q[i]) + (dq[i - 2] if i >= 2 else 0) * 1
                if i >= 1:
                    b += 0
                if a == q1[i]:
                    ok += 1
                else:
                    bad += 1
                    print("string op mismatch", g, R, i, a, q1[i])
                if b == q0[i]:
                    ok += 1
                else:
                    bad += 1
                    print("part0 op mismatch", g, R, i, b, q0[i])
    print(f"  Q(s;R+1) = (1+s)^2 Q' + [J+(1+s)n] Q and Q(s;R+0) = (J+sn) Q + s(1+s) Q': {ok} pass, {bad} fail")


def lagr(n, al, x):
    return lag(n, al, x)


class SchurXi:
    """Xi_h(R) from the Schur side only, by (STAR). R has parts >= 2, non-empty."""
    def __init__(self, tau):
        self.tau = Fr(tau)

    def M(self, h, R, a, C):
        mu = tuple(sorted(R + (1,) * a, reverse=True))
        l = len(mu)
        s = sig(self.tau)
        return C.coeff(mu, 2 * h - 2 + l) * aut(R) * factorial(a) / (s ** (h + l - 1) * PiMu(R, self.tau))

    def xi(self, h, R):
        r = len(R)
        J = sum(R)
        N = 3 * h - 3 + r
        P = 2 * h - 2 + r + N + J + N
        C = ConnSeries(self.tau, P)
        tot = Fr(0)
        for a in range(N + 1):
            tot += Fr((-1) ** a, factorial(a)) * lag(N - a, 2 * h - 2 + r + a, -J - a) * self.M(h, R, a, C)
        return tot


def multisets_level(D):
    """multisets of integers >= 2 with sum(j-1) = D"""
    out = []
    for p in partitions(D):
        out.append(tuple(sorted((v + 1 for v in p), reverse=True)))
    return out


def T7():
    print("T7: Xi_h(R) from Schur only (STAR) vs definition; range zeros and controls")
    for tau in (Fr(2, 7), Fr(-7, 3)):
        X = SchurXi(tau)
        ts = ts_of(tau)
        for h in (2, 3):
            for D in range(1, 2 * h + 1):
                for R in multisets_level(D):
                    if h == 3 and (len(R) > 2 and D > 3):
                        continue
                    if h == 3 and sum(R) > 6:
                        continue
                    t1 = time.time()
                    s = X.xi(h, R)
                    d = Xi_direct(h, R, ts) if (h == 2 or len(R) <= 3) else None
                    tag = "in range" if D <= 2 * h - 3 else "outside"
                    match = "(no direct)" if d is None else ("match" if d == s else f"MISMATCH direct {d}")
                    print(f"  tau={tau} h={h} R={R} D={D} [{tag}]: Xi = {s} {match}  t={time.time()-t1:.1f}s", flush=True)
    # genus 4, 5 from Schur only (one and two leaves, range and first point outside)
    tau = Fr(2, 7)
    X = SchurXi(tau)
    for h, Rs in ((4, [(2,), (3,), (4,), (5,), (6,), (2, 2), (3, 2), (7,)]), (5, [(2,), (3,), (9,)])):
        for R in Rs:
            t1 = time.time()
            D = sum(R) - len(R)
            print(f"  Schur only: tau={tau} h={h} R={R} D={D} [{'in range' if D <= 2*h-3 else 'outside'}]: Xi = {X.xi(h, R)}  t={time.time()-t1:.1f}s", flush=True)


# --------------- tables ---------------
# frozen copy (2026-11-27) of lemmaA/data/ of the resolved-conifold project, see conifold_lemmaA_data/README.md
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "conifold_lemmaA_data")


def load_expr(path):
    with open(path) as f:
        line = [l for l in f if l.strip()][0]
    g, r, expr = [t.strip() for t in line.split("|")]
    assert re.fullmatch(r"[0-9a-zA-Z_\s\+\-\*/\^\(\)]+", expr), "unexpected characters"
    expr = expr.replace("^", "**")
    expr = re.sub(r"(?<![a-zA-Z0-9_])(\d+)", r"Fr(\1)", expr)
    return int(g), int(r), compile(expr, path, "eval")


def T8():
    print("T8: Mathematica tables (independent code)")
    ok = bad = 0
    for g, r in ((2, 1), (2, 2), (3, 1), (3, 2)):
        gg, rr, code = load_expr(f"{DATA}/q_g{g}_r{r}.txt")
        assert (gg, rr) == (g, r)
        for xt in (Fr(-2, 9), Fr(-3, 5)):
            c = 1 + xt
            tau = -xt / (1 + xt)
            ts = ts_of(tau)
            legs = [(2,), (3,), (5,)] if r == 1 else [(2, 2), (2, 3), (4, 3)]
            if g == 3:
                legs = legs[:2]
            for R in legs:
                N = 3 * g - 3 + r
                q = Qg_coeffs(g, R, ts)
                for sval in (Fr(0), Fr(-1), Fr(2)):
                    env = {"Fr": Fr, "x": xt, "s": sval, "j1": Fr(R[0])}
                    if r == 2:
                        env["j2"] = Fr(R[1])
                    tab = eval(code, {"__builtins__": {}}, env)
                    mine = sum(cc * sval ** k for k, cc in enumerate(q)) * c ** (-(3 * g - 3 + 2 * r))
                    if tab == mine:
                        ok += 1
                    else:
                        bad += 1
                        print("  table mismatch", g, r, xt, R, sval, tab, mine)
    print(f"  q tables g2,g3 r1,r2 at c-normalisation: {ok} match, {bad} mismatch")
    # control: wrong c power
    gg, rr, code = load_expr(f"{DATA}/q_g2_r1.txt")
    xt = Fr(-2, 9); c = 1 + xt; tau = -xt / (1 + xt)
    q = Qg_coeffs(2, (3,), ts_of(tau))
    tab = eval(code, {"__builtins__": {}}, {"Fr": Fr, "x": xt, "s": Fr(2), "j1": Fr(3)})
    print("  control, exponent 3-2r instead of 3h-3+2r:", "fails as it must" if tab != sum(cc * Fr(2) ** k for k, cc in enumerate(q)) * c ** (-(3 - 2)) else "DID NOT FAIL")


def T9():
    print("T9: E-scaling and E-string identities, direct at E != 1")
    from engine import _exps, _mult
    def xiE(h, R, tau, E):
        # Omega_w = prod Lambda(-1/w_a), w = E(1,tau,-1-tau); lam_j = 1/(E - j psi); T = -psi^2/(E - psi)
        ts = (-1 / E, -1 / (E * tau), 1 / (E * (1 + tau)))
        r = len(R)
        tot = Fr(0)
        for k in range(0, 3 * h - 3 + r + 1):
            n = r + k
            dim = 3 * h - 3 + n
            for tE in range(2 * k, dim + 1):
                for e in _exps(tE, k, 2):
                    # T legs: -psi^n E^{-(n-1)}
                    wT = Fr((-1) ** k)
                    for v in e:
                        wT *= E ** (-(v - 1))
                    # lam legs: E^{-1} 1/(1-(j/E) psi)
                    val = leg_integral(h, tuple(Fr(j) / E for j in R), ts, extra_psi=e) * E ** (-r)
                    tot += _mult(e) * wT * val / factorial(k)
        return tot
    ok = bad = 0
    for E in (Fr(3), Fr(-5, 2)):
        for tau in (Fr(2, 7), Fr(-7, 3)):
            ts1 = ts_of(tau)
            for h, R in ((2, ()), (2, (2,)), (2, (3,)), (2, (2, 3)), (3, (3,))):
                r = len(R)
                a = xiE(h, R, tau, E)
                b = E ** (-(3 * h - 3 + 2 * r)) * Xi_direct(h, R, ts1)
                ok += a == b
                bad += a != b
                # (S_P2), (L0_P2)
                J = sum(R)
                D = J - r
                a1 = xiE(h, R + (1,), tau, E)
                a0 = xiE(h, R + (0,), tau, E)
                if a1 == J / E ** 2 * a: ok += 1
                else: bad += 1; print("S_P2 mismatch", E, tau, h, R)
                if a0 == (D - 2 * h + 2) / E ** 2 * a: ok += 1
                else: bad += 1; print("L0_P2 mismatch", E, tau, h, R, a0, (D - 2*h+2)/E**2*a)
    print(f"  Xi^(E) = E^-(3h-3+2r) Xi^(1), (S_P2), (L0_P2): {ok} pass, {bad} fail")


if __name__ == "__main__":
    t0 = time.time()
    T5()
    T8()
    T9()
    T7()
    print(f"total time {time.time()-t0:.1f}s")
