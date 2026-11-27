# check2.py -- evidence only. Direct test of the all-profile generating-function statement (Theorem MV3):
#   F(x,s;P) = log sum_nu e^{-u|nu|} e^{(tau+1/2) kappa x'/2} prod 1/(2 sinh(h x'/2)) s_nu(p0),
#   x' = x/(1+s), u = s/(1+s), p0_k = P_k/(sigma x Pi_k), p0_1 = s/(sigma x)
# computed from the Schur sum in K((x))[[s,P]] with NO use of MV, against the Hodge engine (dilaton leaves explicit).
from fractions import Fraction as Fr
from math import factorial, comb
import sys, time
from engine import *
from schur import *


class Ctx:
    def __init__(self, tau, m, Xbuf, P, opts):
        self.tau = Fr(tau)
        self.sig = self.tau * (1 + self.tau)
        self.m = m            # s-order
        self.Xbuf = Xbuf      # keep powers of x <= Xbuf
        self.P = P
        self.opts = opts
        self.fc = {}

    def f(self, k, d):
        """(1+s)^k e^{-u d} truncated at s^m. u = s/(1+s)."""
        key = (k, d)
        if key in self.fc:
            return self.fc[key]
        m = self.m
        # (1+s)^k
        w = [Fr(1)]
        for i in range(1, m + 1):
            w.append(w[-1] * (k - i + 1) / i)
        if self.opts.get("no_exp"):
            self.fc[key] = w
            return w
        # u = s - s^2 + ...
        u = [Fr(0)] + [Fr((-1) ** (i - 1)) for i in range(1, m + 1)]
        # exp(-d u)
        e = [Fr(1)]
        for n in range(1, m + 1):
            e.append(sum(i * (-d) * u[i] * e[n - i] for i in range(1, n + 1)) / n)
        r = smul(w, e, m)
        self.fc[key] = r
        return r


def add2(A, i, e, v):
    d = A[i]
    d[e] = d.get(e, Fr(0)) + v


def zero2(m):
    return [dict() for _ in range(m + 1)]


def mul2(A, B, ctx):
    m = ctx.m
    X = ctx.Xbuf
    C = zero2(m)
    for i in range(m + 1):
        if not A[i]:
            continue
        for j in range(m + 1 - i):
            if not B[j]:
                continue
            c = C[i + j]
            for e1, v1 in A[i].items():
                for e2, v2 in B[j].items():
                    e = e1 + e2
                    if e > X:
                        continue
                    c[e] = c.get(e, Fr(0)) + v1 * v2
    return C


def Zcoef(B, ctx):
    """coefficient of the monomial P_B in Z_s, as an s-array of Laurent dicts in x"""
    m = ctx.m
    tau = ctx.tau
    sig = ctx.sig
    A = zero2(m)
    JB = sum(B)
    rB = len(B)
    prefac = Fr(1)
    for k in B:
        prefac /= sig * (Pi(k, tau) if not ctx.opts.get("no_Pi") else 1)
    for a in range(0, m + 1):
        n = JB + a
        mu = tuple(sorted(B + (1,) * a, reverse=True))
        zz = zee(mu)
        sgn = -1 if ctx.opts.get("neg_p1") else 1
        cA = prefac * (Fr(sgn) / sig) ** a / zz
        for nu in partitions(n):
            ch = chi(nu, mu)
            if ch == 0:
                continue
            S = Wseries(nu, tau, ctx.P)
            for nn in range(0, ctx.P + 1):
                e = nn - n - a - rB
                if e > ctx.Xbuf:
                    break
                if S[nn] == 0:
                    continue
                k = n - nn
                if ctx.opts.get("x_plain"):
                    k = 0
                ff = ctx.f(k, n)
                base = cA * ch * S[nn]
                for i in range(a, m + 1):
                    v = ff[i - a]
                    if v != 0:
                        add2(A, i, e, base * v)
    return A


def inv2(Z0, ctx):
    """inverse of Z0, Z0[0] = {0:1}"""
    m = ctx.m
    assert Z0[0] == {0: Fr(1)} or all((v == 0 or e == 0) for e, v in Z0[0].items())
    G = zero2(m)
    G[0] = {0: Fr(1)}
    for i in range(1, m + 1):
        c = {}
        for j in range(1, i + 1):
            for e1, v1 in Z0[j].items():
                for e2, v2 in G[i - j].items():
                    e = e1 + e2
                    if e > ctx.Xbuf:
                        continue
                    c[e] = c.get(e, Fr(0)) - v1 * v2
        G[i] = c
    return G


def log2(Z0, ctx):
    """log of Z0 in K((x))[[s]]: i f_i = i z_i - sum_{j<i} j f_j z_{i-j}"""
    m = ctx.m
    F = zero2(m)
    for i in range(1, m + 1):
        c = {e: v * i for e, v in Z0[i].items()}
        for j in range(1, i):
            for e1, v1 in F[j].items():
                for e2, v2 in Z0[i - j].items():
                    e = e1 + e2
                    if e > ctx.Xbuf:
                        continue
                    c[e] = c.get(e, Fr(0)) - j * v1 * v2
        F[i] = {e: v / i for e, v in c.items()}
    return F


def logF(R, ctx):
    """s-array for [P_R] log Z_s"""
    R = tuple(sorted(R, reverse=True))
    Z0 = Zcoef((), ctx)
    # drop numerically zero entries at s^0 (only nu = empty contributes at s^0)
    Z0[0] = {e: v for e, v in Z0[0].items() if v != 0}
    if not R:
        return log2(Z0, ctx)
    Ginv = inv2(Z0, ctx)
    zh = {}
    for B in sub_multisets(R):
        if B:
            zh[B] = mul2(Zcoef(B, ctx), Ginv, ctx)
    c = {}
    for B in sorted(zh, key=lambda b: (sum(b), len(b))):
        t = [dict(d) for d in zh[B]]
        nB = sum(B)
        for rho in sub_multisets(B):
            if not rho or rho == B:
                continue
            pr = mul2(c[rho], zh[msub(B, rho)], ctx)
            f = Fr(sum(rho), nB)
            for i in range(ctx.m + 1):
                for e, v in pr[i].items():
                    t[i][e] = t[i].get(e, Fr(0)) - f * v
        c[B] = t
    return c[R]


def setup(tau, g, R, m, opts=None):
    opts = opts or {}
    J = sum(R)
    r = len(R)
    Xtarget = 2 * g - 2
    Xbuf = Xtarget + 2 * (m + J) + r + 2
    P = Xbuf + J + m + m + r + 2
    return Ctx(tau, m, Xbuf, P, opts), Xtarget


def extracted(tau, g, R, m, opts=None):
    ctx, X = setup(tau, g, R, m, opts)
    F = logF(R, ctx)
    sig = ctx.sig
    au = aut(R)
    series = [F[i].get(X, Fr(0)) * au / sig ** (g - 1) for i in range(m + 1)]
    # all powers of x below -2 and all odd powers must cancel, up to the target
    viol = 0
    for i in range(m + 1):
        for e, v in F[i].items():
            if v == 0 or e > X:
                continue
            if e < -2 or e % 2:
                viol += 1
    return series, viol


def ser_of_rational(num_den_fn, m):
    return [num_den_fn(i) for i in range(m + 1)]


def expected(g, R, tau, m):
    """the right side of Theorem MV3, from the Hodge engine or the closed unstable table"""
    tau = Fr(tau)
    sig = tau * (1 + tau)
    r = len(R)
    ts = ts_of(tau)
    if g == 0 and r == 0:
        # s + s^2/4 + s^3/(6(1+s))
        out = [Fr(0)] * (m + 1)
        if m >= 1: out[1] = Fr(1)
        if m >= 2: out[2] = Fr(1, 4)
        for i in range(3, m + 1):
            out[i] = Fr((-1) ** (i - 3), 6)
        return out
    if g == 0 and r == 1:
        j = R[0]
        out = [Fr(0)] * (m + 1)
        out[0] = Fr(1, j * j)
        if m >= 1:
            out[1] = Fr(1, j * j) - Fr(1, j * (j + 1))
        return out
    if g == 0 and r == 2:
        out = [Fr(0)] * (m + 1)
        out[0] = Fr(1, sum(R))
        return out
    if g == 1 and r == 0:
        # (1/24)[log(1+s) + (1+1/sigma)(1/(1+s) - 1)]
        out = [Fr(0)] * (m + 1)
        for i in range(1, m + 1):
            out[i] = (Fr((-1) ** (i - 1), i) + (1 + 1 / sig) * (-1) ** i) / 24
        return out
    N = 3 * g - 3 + r
    q = Qg_coeffs(g, R, ts, kmax=m)
    if r == 0:
        # subtract (1+s)^{-(2g-2)} <Omega>_g
        om = sum(c * hodge(g, (), ch) for ch, c in omega_monos(g, 3 * g - 3, ts))
        q = list(q)
        for i in range(m + 1):
            q[i] -= om * Fr((-1) ** i * comb(2 * g - 2 + i - 1, i))
    return q


def main():
    t0 = time.time()
    cases = [
        (0, (), 5), (0, (2,), 3), (0, (3,), 3), (0, (2, 3), 3), (0, (2, 2), 3), (0, (2, 2, 3), 3), (0, (2, 3, 4), 2),
        (1, (), 5), (1, (2,), 4), (1, (3,), 4), (1, (5,), 3), (1, (2, 2), 4), (1, (2, 3), 3), (1, (2, 2, 2), 3),
        (2, (), 5), (2, (2,), 5), (2, (3,), 5), (2, (4,), 5), (2, (2, 2), 6), (2, (2, 3), 5),
        (3, (), 6), (3, (2,), 8),
    ]
    tot_ok = tot_bad = viol_tot = 0
    for tau in (Fr(2, 7), Fr(3, 2)):
        for g, R, m in cases:
            if g == 3 and tau != Fr(2, 7):
                continue
            t1 = time.time()
            got, viol = extracted(tau, g, R, m)
            exp = expected(g, R, tau, m)
            bad = sum(1 for a, b in zip(got, exp) if a != b)
            tot_ok += len(got) - bad
            tot_bad += bad
            viol_tot += viol
            print(f"tau={tau} g={g} R={R} s^0..s^{m}: mismatches {bad}, x-power violations {viol}, t={time.time()-t1:.1f}s", flush=True)
            if bad:
                print("   got", got)
                print("   exp", exp)
    print(f"T6 Theorem MV3: {tot_ok} coefficients match, {tot_bad} mismatches, x-power violations {viol_tot}")
    print("Controls (each must fail at g=2, R=(2), tau=2/7, and g=1, R=(3)):")
    for name, opts in (("x' = x (no 1/(1+s))", {"x_plain": True}),
                       ("no e^{-u|nu|}", {"no_exp": True}),
                       ("p0_1 = -s/(sigma x)", {"neg_p1": True}),
                       ("p0_k without Pi_k", {"no_Pi": True})):
        for g, R, m in ((2, (2,), 4), (1, (3,), 3)):
            got, viol = extracted(Fr(2, 7), g, R, m, opts)
            exp = expected(g, R, Fr(2, 7), m)
            bad = sum(1 for a, b in zip(got, exp) if a != b)
            print(f"  {name}: g={g} R={R}: {bad} of {m+1} coefficients differ, x-power violations {viol}")
    print(f"total time {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
