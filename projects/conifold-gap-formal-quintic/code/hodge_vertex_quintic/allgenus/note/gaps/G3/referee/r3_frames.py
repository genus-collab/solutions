# -*- coding: utf-8 -*-
"""r3_frames.py -- referee check (independent code; nothing imported from the prover's scripts) of lem:tra-frames in the
graded setting (homogeneous elements over Q[v]).

Notation of NOTATION.md: H = Q^5, basis H^a (a in Z/5), nu = (0,-1,3,2,1), eps = (1,-1,-1,1,1), Y^C_b = eps_b I^C_(nu_b),
Pi_i(q) = -[z^i] <Y^C(-z), q(z)>.
Route of this script (different from both routes of g3_frames.py):
  (Y)  [z^-n] I^C_m is read off the RESIDUE DEFINITION Res_(x=0) x^m exp(-(x^4/4 + v/x)/z) dx by a direct double sum over the two
       exponential series, at two values v = 1 and v = 2 (so that homogeneity is tested, not assumed); then the model system
       z dY_b/dv = bmod_b Y_(b+1), bmod = (1, -1/v, 1, -1, -1), and Y^C(0) = -e_1 are checked coefficientwise;
  (F)  frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H is obtained by GENERIC Gaussian elimination (row pivoting) of the k x k linear
       system Pi_i(frame) = 0, 0 <= i < k, with Pi_i computed from the pairing itself (no shift rule, no Theta, no top-down rule);
       uniqueness = full rank; then Pi_i(frame) = 0 is re-checked for k <= i <= k + 12.
Checks: (T) trivial frames (all f_j = 0) are exactly the (k,a) with k < nu_a + 1 (ten), the source list of eight is incomplete;
(J) the top correction is f_(j0) = eps_a v^(nu_a+1)/(nu_a+1)!, j0 = k - nu_a - 1; (W) f_j(2)/f_j(1) = 2^e, e = (4(k-j)+nu_a+1)/5,
f_j = 0 unless e is an integer, e >= 1; (X) the outline's examples.
Controls: (C1) frames built from Y^C with the z^-n terms, n >= 5, dropped are NOT annihilated by the true tower (k >= 5);
(C2) the wrong exponent (4(k-j)+nu_a)/5 is refuted by (W); (C3) the mirrored rule "trivial iff k <= nu_a + 1" is refuted.
Usage (working directory: the scratch folder):  py -I -B r3_frames.py KMAX   -- log on stdout, no file written, no path printed.
"""
import sys
from fractions import Fraction as Fr
from math import factorial

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)


def ic_coeff(m, n, v):
    """[z^-n] Res_x x^m exp(-x^4/(4z)) exp(-v/(xz)) dx: terms x^(4i) z^-i (-1/4)^i/i! times x^-l z^-l (-v)^l/l!, m+4i-l = -1."""
    tot = Fr(0)
    for i in range(n + 1):
        l = m + 4 * i + 1
        if l < 0 or i + l != n:
            continue
        tot += Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-v) ** l, factorial(l))
    return tot


def yc(b, n, v):
    """[z^-n] Y^C_b at the value v (n >= 0); 0 for n < 0."""
    if n < 0:
        return Fr(0)
    return EPS[b] * ic_coeff(NU[b], n, v)


def yc_vpower(b, n):
    """the v-power of [z^-n] Y^C_b read off the double sum (l = n - i), or None if the coefficient is zero."""
    m = NU[b]
    for i in range(n + 1):
        l = m + 4 * i + 1
        if l >= 0 and i + l == n:
            return l
    return None


def pi(i, m, b, v, ycf=yc):
    """Pi_i(z^m H^b) = -[z^i] (Y^C_b(-z) z^m) = -(-1)^(m-i) [z^-(m-i)] Y^C_b, from the pairing."""
    n = m - i
    if n < 0:
        return Fr(0)
    return -((-1) ** n) * ycf(b, n, v)


def solve(M, rhs):
    """generic Gaussian elimination over Q with row pivoting; returns the unique solution or None if singular."""
    n = len(M)
    A = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return None
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fac = A[r][c] / A[c][c]
                A[r] = [x - fac * y for x, y in zip(A[r], A[c])]
    return [A[r][n] / A[r][r] for r in range(n)]


def frame(k, a, v, ycf=yc):
    """{j: f_j} at the value v, or None if the system is singular."""
    if k == 0:
        return {}
    M = [[pi(i, j, 1, v, ycf) for j in range(k)] for i in range(k)]
    rhs = [-pi(i, k, a, v, ycf) for i in range(k)]
    f = solve(M, rhs)
    if f is None:
        return None
    return {j: c for j, c in enumerate(f) if c != 0}


def annihilated(k, a, f, v, imax):
    for i in range(imax + 1):
        tot = pi(i, k, a, v) + sum((c * pi(i, j, 1, v) for j, c in f.items()), Fr(0))
        if tot != 0:
            return False
    return True


def yc_trunc(b, n, v):
    return yc(b, n, v) if n < 5 else Fr(0)


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    print("r3_frames.py  KMAX = %d  (exact over Q; Y^C from the residue definition at v = 1 and v = 2)" % K)
    ok = True
    # (Y) the vanishing solution: model system, initial value, closed form
    bad_ode = bad_init = bad_closed = 0
    for v in (Fr(1), Fr(2)):
        for b in range(5):
            if yc(b, 0, v) != (Fr(-1) if b == 1 else Fr(0)):
                bad_init += 1
            for n in range(1, K + 20):
                # z d/dv of c v^e z^-n contributes e c v^(e-1) at z^-(n-1); compare with bmod_b [z^-(n-1)] Y_(b+1)
                e = yc_vpower(b, n)
                lhs = Fr(0) if e is None or e == 0 else yc(b, n, v) * e / v
                rhs = yc((b + 1) % 5, n - 1, v) * (Fr(-1) / v if b == 1 else Fr((1, 1, 1, -1, -1)[b]))
                if lhs != rhs:
                    bad_ode += 1
                q = n // 5
                cl = Fr(EPS[b] * (-1) ** n, 4 ** q * factorial(q) * factorial(n - q)) * v ** (n - q) \
                    if (n - NU[b] - 1) % 5 == 0 else Fr(0)
                if cl != yc(b, n, v):
                    bad_closed += 1
    print("[Y] model system z dY_b/dv = bmod_b Y_(b+1): mismatches %d; Y^C(0) = -e_1: mismatches %d; closed form of "
          "lem:gau-vanshape: mismatches %d  (n < %d, v = 1, 2)" % (bad_ode, bad_init, bad_closed, K + 20))
    ok &= bad_ode == bad_init == bad_closed == 0
    # (F), (T), (J), (W)
    triv, nsing, nann, nexp, nj0, nfr = [], 0, 0, 0, 0, 0
    wrong_exp_refuted = 0
    for k in range(K + 1):
        for a in (0, 2, 3, 4):
            f1, f2 = frame(k, a, Fr(1)), frame(k, a, Fr(2))
            nfr += 1
            if f1 is None or f2 is None:
                nsing += 1
                continue
            if not annihilated(k, a, f1, Fr(1), k + 12) or not annihilated(k, a, f2, Fr(2), k + 12):
                nann += 1
            if not f1:
                triv.append((k, a))
            for j in set(f1) | set(f2):
                e = Fr(4 * (k - j) + NU[a] + 1, 5)
                c1, c2 = f1.get(j, Fr(0)), f2.get(j, Fr(0))
                if e.denominator != 1 or e < 1 or c1 == 0 or c2 != c1 * Fr(2) ** int(e):
                    nexp += 1
                ew = Fr(4 * (k - j) + NU[a], 5)
                if not (ew.denominator == 1 and c1 != 0 and c2 == c1 * Fr(2) ** int(ew)):
                    wrong_exp_refuted += 1
            j0 = k - NU[a] - 1
            if j0 >= 0 and f1.get(j0, Fr(0)) != Fr(EPS[a], factorial(NU[a] + 1)):
                nj0 += 1
    pred = sorted((k, a) for a in (0, 2, 3, 4) for k in range(NU[a] + 1))
    src8 = sorted([(0, 0), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3)])
    mirrored = sorted((k, a) for a in (0, 2, 3, 4) for k in range(NU[a] + 2))
    print("[F] %d frames (k <= %d, a != 1) by Gaussian elimination: singular systems %d; Pi_i(frame) != 0 for some i <= k+12: %d"
          % (nfr, K, nsing, nann))
    print("[T] trivial frames found: %s" % triv)
    print("[T] = {(k,a) : k < nu_a + 1} (%d elements): %s;  missing from the source list of eight: %s"
          % (len(pred), sorted(triv) == pred, sorted(set(triv) - set(src8))))
    print("[J] top correction f_(k - nu_a - 1) = eps_a v^(nu_a+1)/(nu_a+1)!: failures %d" % nj0)
    print("[W] exponent (4(k-j)+nu_a+1)/5 integral, >= 1, and f_j(2) = 2^e f_j(1): failures %d" % nexp)
    ex = {(1, 0): {0: Fr(1)}, (2, 4): {0: Fr(1, 2)}, (3, 3): {0: Fr(1, 6)}, (4, 2): {0: Fr(-1, 24)},
          (6, 0): {5: Fr(1), 0: Fr(-1, 120)}}
    exbad = [ka for ka, fv in ex.items() if frame(ka[0], ka[1], Fr(1)) != fv]
    print("[X] outline examples frame(1,0), (2,4), (3,3), (4,2), (6,0) at v = 1: mismatches %s" % (exbad or "none"))
    ok &= nsing == 0 and nann == 0 and sorted(triv) == pred and nj0 == 0 and nexp == 0 and not exbad
    # controls
    nctl = tot = 0
    for k in range(5, K + 1):
        for a in (0, 2, 3, 4):
            fc = frame(k, a, Fr(1), ycf=yc_trunc)
            if fc is None:
                continue
            tot += 1
            if not annihilated(k, a, fc, Fr(1), k + 12):
                nctl += 1
    print("[C1] control, frames from Y^C with z^-n (n >= 5) dropped: %d of %d frames (5 <= k <= %d) NOT annihilated (expected > 0)"
          % (nctl, tot, K))
    print("[C2] control, wrong exponent (4(k-j)+nu_a)/5: refuted on %d coefficients (expected > 0)" % wrong_exp_refuted)
    print("[C3] control, rule 'trivial iff k <= nu_a + 1' (%d elements) agrees with the computation: %s (expected False)"
          % (len(mirrored), sorted(triv) == mirrored))
    ok &= nctl > 0 and wrong_exp_refuted > 0 and sorted(triv) != mirrored
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


if __name__ == "__main__":
    main()
