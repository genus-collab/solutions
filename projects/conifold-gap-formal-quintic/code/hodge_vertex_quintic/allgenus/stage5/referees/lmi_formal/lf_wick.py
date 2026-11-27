"""lf_wick: the Wick functional N of CHART_PROOF s.1.2, implemented from its definition.

Scaled variables xi = t eta with t = z^(1/2). The normalised integrand exp(Phi_3(xi)/z) A(p + xi, z)/A_0(p)
equals exp(E), E = Sum_(k>=3) ph[k] t^(k-2) eta^k + Sum logA[(d, m)] t^(d + 2m) eta^d, and
N_n = < [t^(2n)] exp(E) >_C, with the Gaussian moments of Isserlis-Wick (covariance C = -H^-1).
Coefficients may be Fractions or lf_core.TS elements (pass zero and one of the ring).

* wick_1     one variable (ph3 = {k: coeff, k >= 3}, logA = {(d, m): coeff}, C = -1/Phi''(p)).
* wick_2sep  two variables (s, m) whose cubic-and-higher phase and log-amplitude are SEPARATED,
             Phi_3 = G(s) + K(m), A = A_s(s, z) A_m(m, z), with an arbitrary 2x2 covariance.
* wick_2gen  two variables, arbitrary (non-separated) data; used for cross-checks and controls.
* moments_2  <eta_s^a eta_m^b>: perfect matchings with k cross pairs, counted explicitly.
"""
from fractions import Fraction as Fr
from math import comb, factorial
from lf_core import dfact


def _powers(x, n, one):
    out = [one]
    for _ in range(n):
        out.append(out[-1] * x)
    return out


def exp_graded_1(E, R, zero, one):
    """E[r] (r = 1..R) polynomials in eta (lists); returns f[0..R], Sum_r f[r] t^r = exp(Sum_r E[r] t^r)."""
    f = [[one]]
    for r in range(1, R + 1):
        acc = []
        for k in range(1, r + 1):
            Ek, fk = E[k], f[r - k]
            for i, a in enumerate(Ek):
                if not a:
                    continue
                ak = a * k
                for j, b in enumerate(fk):
                    if b:
                        while len(acc) <= i + j:
                            acc.append(zero)
                        acc[i + j] = acc[i + j] + ak * b
        inv_r = Fr(1, r)
        f.append([c * inv_r for c in acc])
    return f


def _build_E1(ph3, logA, R, zero):
    E = [None] + [[] for _ in range(R)]

    def add(r, d, c):
        if 1 <= r <= R:
            row = E[r]
            while len(row) <= d:
                row.append(zero)
            row[d] = row[d] + c
    for k, c in ph3.items():
        assert k >= 3
        add(k - 2, k, c)
    for (d, m), c in logA.items():
        assert (d, m) != (0, 0) and d >= 0 and m >= 0
        add(d + 2 * m, d, c)
    return E


def wick_1(ph3, logA, C, nmax, zero, one):
    R = 2 * nmax
    f = exp_graded_1(_build_E1(ph3, logA, R, zero), R, zero, one)
    Cp = _powers(C, 3 * nmax + 2, one)
    N = []
    for n in range(nmax + 1):
        acc = zero
        for d, c in enumerate(f[2 * n]):
            if c and d % 2 == 0:
                acc = acc + c * (dfact(d - 1) * Cp[d // 2])
        N.append(acc)
    return N


def moments_2(Css, Csm, Cmm, amax, bmax, zero, one):
    """M[a][b] = <eta_s^a eta_m^b>_C: Sum_k C(a,k) C(b,k) k! (a-k-1)!! (b-k-1)!! Css^((a-k)/2) Cmm^((b-k)/2) Csm^k."""
    Pss = _powers(Css, amax // 2 + 1, one)
    Pmm = _powers(Cmm, bmax // 2 + 1, one)
    Psm = _powers(Csm, min(amax, bmax) + 1, one)
    M = [[zero] * (bmax + 1) for _ in range(amax + 1)]
    for a in range(amax + 1):
        for b in range(bmax + 1):
            if (a + b) % 2:
                continue
            acc = zero
            for k in range(a % 2, min(a, b) + 1, 2):
                coef = comb(a, k) * comb(b, k) * factorial(k) * dfact(a - k - 1) * dfact(b - k - 1)
                acc = acc + coef * (Pss[(a - k) // 2] * Pmm[(b - k) // 2] * Psm[k])
            M[a][b] = acc
    return M


def cov_2(Hss, Hsm, Hmm):
    """C = -H^-1 for H = [[Hss, Hsm], [Hsm, Hmm]]; returns (Css, Csm, Cmm, det H)."""
    det = Hss * Hmm - Hsm * Hsm
    idet = Fr(1) / det if not hasattr(det, "inv") else det.inv()
    return (-1 * Hmm) * idet, Hsm * idet, (-1 * Hss) * idet, det


def wick_2sep(phs3, logAs, phm3, logAm, cov, nmax, zero, one):
    """Separated two-variable data; cov = (Css, Csm, Cmm)."""
    R = 2 * nmax
    fs = exp_graded_1(_build_E1(phs3, logAs, R, zero), R, zero, one)
    fm = exp_graded_1(_build_E1(phm3, logAm, R, zero), R, zero, one)
    amax = max(len(x) for x in fs)
    bmax = max(len(x) for x in fm)
    M = moments_2(cov[0], cov[1], cov[2], amax, bmax, zero, one)
    N = []
    for n in range(nmax + 1):
        acc = zero
        for r1 in range(2 * n + 1):
            row_m = fm[2 * n - r1]
            for a, x in enumerate(fs[r1]):
                if not x:
                    continue
                Ma = M[a]
                for b, y in enumerate(row_m):
                    if y and (a + b) % 2 == 0:
                        acc = acc + x * y * Ma[b]
        N.append(acc)
    return N


def wick_2gen(Phi3, logA, cov, nmax, zero, one):
    """Generic two-variable data: Phi3 = {(i, j): c} (s^i m^j, i + j >= 3), logA = {(i, j, mz): c}."""
    R = 2 * nmax
    E = [None] + [dict() for _ in range(R)]
    for (i, j), c in Phi3.items():
        assert i + j >= 3
        r = i + j - 2
        if r <= R:
            E[r][(i, j)] = E[r].get((i, j), zero) + c
    for (i, j, mz), c in logA.items():
        assert (i, j, mz) != (0, 0, 0)
        r = i + j + 2 * mz
        if r <= R:
            E[r][(i, j)] = E[r].get((i, j), zero) + c
    f = [{(0, 0): one}]
    for r in range(1, R + 1):
        acc = {}
        for k in range(1, r + 1):
            for (i, j), a in E[k].items():
                ak = a * k
                for (p, q), b in f[r - k].items():
                    key = (i + p, j + q)
                    acc[key] = acc.get(key, zero) + ak * b
        inv_r = Fr(1, r)
        f.append({key: v * inv_r for key, v in acc.items()})
    amax = max([a for row in f for (a, b) in row] + [0])
    bmax = max([b for row in f for (a, b) in row] + [0])
    M = moments_2(cov[0], cov[1], cov[2], amax, bmax, zero, one)
    N = []
    for n in range(nmax + 1):
        acc = zero
        for (a, b), x in f[2 * n].items():
            if (a + b) % 2 == 0 and x:
                acc = acc + x * M[a][b]
        N.append(acc)
    return N
