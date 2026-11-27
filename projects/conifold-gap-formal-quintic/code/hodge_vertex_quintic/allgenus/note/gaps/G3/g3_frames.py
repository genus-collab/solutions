# -*- coding: utf-8 -*-
"""g3_frames.py -- exact check (over Q, homogeneous elements over Q[v]) of the frame construction of lem:tra-frames.

Gap cluster G3, item 1 (lem:tra-frames step 4).  Notation of NOTATION.md: H = Q^5, basis H^a (a in Z/5), charges
nu = (0,-1,3,2,1), signs eps = (1,-1,-1,1,1), weight wt(v) = 5, wt(z^k H^a) = 4k + nu_a + 1.  A homogeneous element of
H (x) Q[v][z] of weight w is stored as {(k, a): c}, meaning Sum c v^((w - 4k - nu_a - 1)/5) z^k H^a.

Two independent routes to the vanishing solution Y^C (lowered components, Y_b = eps_b I^C_(nu_b)):
  (I)  the closed form of lem:gau-vanshape: [z^-n] Y^C_b = eps_b (-1)^n v^(n-q) / (4^q q! (n-q)!), q = floor(n/5),
       if n = nu_b + 1 mod 5, else 0;
  (II) the regular solution of the model system z dY_b/dv = bmod_b Y_(b+1), bmod = (1, -1/v, 1, -1, -1), Y(0) = -e_1,
       by iterated integrals of monomials (the construction of lem:gau-uniq / lem:gau-vanq applied to the model).
Frames are built by two routes:
  (A)  Theta^-1 = Sum_r (-N)^r applied to (Pi_i(z^k H^a))_i  (the proof proposed in GAPS_G3_2.md, Step 1-3);
  (B)  top-down elimination of the non-negative powers of <Y^C(-z), n(z)> (the source's route), using Y^C from (II).
Checks: (I) = (II); (A) = (B); Pi_i(frame) = 0 for 0 <= i <= k + 12 with Y^C from (II); exponents integral and >= 1;
support rule; the source's examples.  Control C1: frames built from Y^C(+z) instead of Y^C(-z) are NOT annihilated.
Usage: py -I -B g3_frames.py KMAX     (prints a log to stdout; paths relative)
"""
import sys
from fractions import Fraction as Fr
from math import factorial

NU = (0, -1, 3, 2, 1)
EPS = (1, -1, -1, 1, 1)
BMOD = (1, -1, 1, -1, -1)          # coefficient of bmod_b; b = 1 carries the extra factor 1/v
PI0_H = None


def yc_closed(b, n):
    """[z^-n] Y^C_b at v = 1 (the v-power is n - floor(n/5) by homogeneity), closed form (I)."""
    if (n - NU[b] - 1) % 5 != 0 or n < 0:
        return Fr(0)
    q = n // 5
    if n - q < 0:
        return Fr(0)
    return Fr(EPS[b] * (-1) ** n, 4 ** q * factorial(q) * factorial(n - q))


def yc_ode(nmax):
    """{(b, n): (coefficient, v-power)} of the regular solution of z dY_b/dv = bmod_b Y_(b+1), Y(v=0) = -e_1, route (II)."""
    y = {(b, 0): ((Fr(-1), 0) if b == 1 else (Fr(0), 0)) for b in range(5)}
    for n in range(nmax):
        for b in range(5):
            c, j = y[((b + 1) % 5, n)]
            if c == 0:
                y[(b, n + 1)] = (Fr(0), 0)
                continue
            jj = j - (1 if b == 1 else 0)          # bmod_1 = -1/v lowers the power
            assert jj >= 0, "pole in the integrand"
            y[(b, n + 1)] = (BMOD[b] * c / (jj + 1), jj + 1)   # Int_0^v
    return y


def pi0(m, b):
    """Pi_0(z^m H^b) at v = 1 = (-1)^(m+1) [z^-m] Y^C_b (closed form)."""
    if m < 0:
        return Fr(0)
    return (-1) ** (m + 1) * yc_closed(b, m)


def frame_theta(k, a):
    """frame(k, a) by route (A): returns {j: f_j} (f_j at v = 1), frame = z^k H^a + Sum_j f_j z^j H."""
    rhs = [pi0(k - i, a) for i in range(k)]            # (Pi_i(z^k H^a))_(i<k); Pi_k = Pi_0(H^a) = 0 for a != 1
    piH = {5 * q: pi0(5 * q, 1) for q in range(1, k // 5 + 2)}

    def Nop(s):
        return [sum((piH[5 * q] * s[i + 5 * q] for q in range(1, (k - 1 - i) // 5 + 1)), Fr(0)) for i in range(k)]
    sigma, term, r = [Fr(0)] * k, rhs[:], 0
    while any(term):
        sigma = [x + y for x, y in zip(sigma, term)]
        term = [-x for x in Nop(term)]
        r += 1
        assert r <= k + 2, "N not nilpotent on finitely supported vectors?"
    return {j: -sigma[j] for j in range(k) if sigma[j] != 0}


def frame_elim(k, a, yode, sign=-1):
    """frame(k, a) by route (B): eliminate the powers z^(k-1), .., z^0 of <Y^C(sign z), n(z)> from the top, with Y^C from (II)."""
    def ycoef(b, n):               # coefficient of z^-n of Y^C_b(sign*z) at v = 1
        c, _ = yode[(b, n)]
        return c * (sign ** n)
    f = {}
    for i in range(k - 1, -1, -1):  # coefficient of z^i of <Y^C(sign z), z^k H^a + Sum_j f_j z^j H>
        tot = ycoef(a, k - i) + sum((f[j] * ycoef(1, j - i) for j in f if j >= i), Fr(0))
        # the z^i H term contributes f_i * ycoef(1, 0) = -f_i; choose f_i to kill tot
        f[i] = tot / (-ycoef(1, 0)) if i not in f else f[i]
    return {j: c for j, c in f.items() if c != 0}


def pairing_nonneg(k, a, f, yode, imax):
    """coefficients of z^i, 0 <= i <= imax, of <Y^C(-z), frame(z)> with Y^C from route (II)."""
    vec = {(k, a): Fr(1)}
    for j, c in f.items():
        vec[(j, 1)] = vec.get((j, 1), Fr(0)) + c
    out = []
    for i in range(imax + 1):
        tot = Fr(0)
        for (kk, b), c in vec.items():
            n = kk - i
            if n >= 0:
                tot += c * yode[(b, n)][0] * (-1) ** n
        out.append(tot)
    return out


def main():
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    print("g3_frames.py  KMAX = %d  (exact over Q; homogeneous, at v = 1 with v-powers from the weights)" % KMAX)
    yode = yc_ode(KMAX + 16)
    bad = 0
    for b in range(5):
        for n in range(KMAX + 16):
            c, j = yode[(b, n)]
            if c != yc_closed(b, n) or (c != 0 and j != n - n // 5):
                bad += 1
    print("[I=II] closed form of Y^C vs iterated-integral solution, b < 5, n < %d: mismatches %d" % (KMAX + 16, bad))
    nfr = nexp = nsupp = nab = nann = 0
    for k in range(KMAX + 1):
        for a in (0, 2, 3, 4):
            fa, fb = frame_theta(k, a), frame_elim(k, a, yode)
            nfr += 1
            if fa != fb:
                nab += 1
            for j, c in fa.items():
                e = Fr(4 * (k - j) + NU[a] + 1, 5)
                if e.denominator != 1 or e < 1:
                    nexp += 1
                if (j - (k - NU[a] - 1)) % 5 != 0:
                    nsupp += 1
            if any(pairing_nonneg(k, a, fa, yode, k + 12)):
                nann += 1
    print("[A=B] frames by Theta^-1 vs top-down elimination, k <= %d, a != 1: %d frames, mismatches %d" % (KMAX, nfr, nab))
    print("[ANN] Pi_i(frame) = 0 for 0 <= i <= k+12 (Y^C from the ODE): failures %d" % nann)
    print("[EXP] coefficients with non-integral or < 1 v-exponent: %d;  support rule j = k - nu_a - 1 mod 5 violated: %d"
          % (nexp, nsupp))
    ex = {(1, 0): {0: Fr(1)}, (2, 4): {0: Fr(1, 2)}, (3, 3): {0: Fr(1, 6)}, (4, 2): {0: Fr(-1, 24)},
          (6, 0): {5: Fr(1), 0: Fr(-1, 120)}}
    for (k, a) in [(0, 0), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3)]:
        ex[(k, a)] = {}
    exbad = [ka for ka, f in ex.items() if frame_theta(*ka) != f]
    print("[SRC] the source's examples (frame(1,0) = z1 + vH, .., frame(6,0) = z^6 1 + v z^5 H - (v^5/120) H, trivial frames): "
          "mismatches %s" % (exbad or "none"))
    # control C1: frames built from Y^C(+z) are not annihilated by the true tower
    nctl = tot = 0
    for k in range(1, KMAX + 1):
        for a in (0, 2, 3, 4):
            fc = frame_elim(k, a, yode, sign=+1)
            if fa_nontrivial(k, a):
                tot += 1
                if any(pairing_nonneg(k, a, fc, yode, k + 12)):
                    nctl += 1
    print("[C1] control, frames from Y^C(+z): %d of %d non-trivial frames FAIL to be annihilated (expected: most)" % (nctl, tot))
    ok = bad == 0 and nab == 0 and nann == 0 and nexp == 0 and nsupp == 0 and not exbad and nctl > 0
    print("RESULT:", "ALL CHECKS PASSED" if ok else "FAILURE")


def fa_nontrivial(k, a):
    return bool(frame_theta(k, a))


if __name__ == "__main__":
    main()
