# -*- coding: utf-8 -*-
"""ref_g1_model.py -- referee evidence (not proof) for gap cluster G1, item 5a (symplecticity of M_0).

Exact arithmetic in Q(zeta_5); stdlib only.  Run:  py -I -B ref_g1_model.py [M]   (M = t-order, default 24).
Writes ref_g1_model.log next to itself (relative name printed).

Model rows (OUTLINE_3 s02 item 10): z dY_b/dv = bm_b Y_(b+1), bm = (1, -1/v, 1, -1, -1), nu = (0,-1,3,2,1), eps = (1,-1,-1,1,1).
With Y_b = e^(U/z) varsigma^(nu_b - 1) Fs_b(t), t = z/varsigma^4, U = -(5/4) varsigma^4 (colour 0) the system becomes
T_(nu_b - 1) Fs_b = bm'_b Fs_(b+1), T_k = -1 + (t/5)(k - 4 t d/dt), bm' = (1,-1,1,-1,-1) (the referee's own derivation; it agrees with
prop:set-modrows (b)).  Fs_0 = 1 + ... is solved order by order (pivot 4n), Rhat^(i)_b(t) = zeta^(i(b+1)) Fs_b(zeta^-i t).
Claim (GAPS_G1_5 s.5.1): S_ij(t) := Sum_(a,b) Rhat^(i)_a(-t) eta^(ab) Rhat^(j)_b(t) = delta_ij, eta^(ab) = (1/5) delta_(a+b = 3 mod 5).
Controls that MUST fail: (m1) the diagonally rescaled system bm' = (2,-1,1,-1,-1/2) (same exponent, not self-adjoint for eta) with
the same eta; (m2) the pairing without t -> -t; (m3) the wrong pairing a+b = 2 mod 5.
Also printed: Fs_0 = 1 - (3/20) t + (9/800) t^2 + ... (the illustration of prop:set-modrows) and Fs_b(0) = eps_b.
"""
import os, sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
LOGNAME = "ref_g1_model.log"
M = int(sys.argv[1]) if len(sys.argv) > 1 else 24
NU = (0, -1, 3, 2, 1)


def zmul(a, b):                      # Q(zeta) elements as 5-tuples mod x^5 - 1
    r = [Fr(0)] * 5
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[(i + j) % 5] += x * y
    return r
def zadd(a, b):
    return [x + y for x, y in zip(a, b)]
def zred(a):                         # canonical form mod Phi_5
    return tuple(x - a[4] for x in a[:4])
def zc(x):
    return [Fr(x), Fr(0), Fr(0), Fr(0), Fr(0)]
def zp(k):
    r = [Fr(0)] * 5
    r[k % 5] = Fr(1)
    return r


def Top(k, f):
    """T_k f for a series f (list of Fractions): T_k t^n = -t^n + ((k - 4n)/5) t^(n+1)."""
    r = [Fr(0)] * (len(f) + 1)
    for n, c in enumerate(f):
        r[n] -= c
        r[n + 1] += c * Fr(k - 4 * n, 5)
    return r


def solve(bm):
    """colour-0 rows Fs_0..Fs_4 to order t^M for constants bm' with prod = -1."""
    f = [Fr(1)] + [Fr(0)] * M
    def chain(f0):
        Fs = [f0]
        for b in range(4):
            Fs.append([c / bm[b] for c in Top(NU[b] - 1, Fs[b])])
        back = [c / bm[4] for c in Top(NU[4] - 1, Fs[4])]
        return Fs, back
    for m in range(1, M + 1):
        f[m] = Fr(0)
        Fs, back = chain(f)
        res = back[m + 1] - (f[m + 1] if m + 1 <= M else 0)
        f[m] = -res / (4 * m)
    Fs, back = chain(f)
    assert all(back[n] == f[n] for n in range(M + 1)), "closing not solved"
    return [F[:M + 1] for F in Fs]


def Rhat(Fs, i):
    """Rhat^(i)_b(t) coefficients in Q(zeta): zeta^(i(b+1)) Fs_b(zeta^-i t)."""
    return [[zmul(zp(i * (b + 1) - i * n), zc(c)) for n, c in enumerate(Fs[b])] for b in range(5)]


def pairing(Ri, Rj, sign=-1, shift=3):
    """coefficients of Sum_a (1/5) Ri_a(sign t) Rj_(shift - a)(t), to order t^M."""
    out = []
    for n in range(M + 1):
        s = zc(0)
        for a in range(5):
            b = (shift - a) % 5
            for m in range(n + 1):
                x = Ri[a][m]
                if sign == -1 and m % 2:
                    x = [-y for y in x]
                s = zadd(s, zmul(x, Rj[b][n - m]))
        out.append(zred([y / 5 for y in s]))
    return out


def is_delta(P, i, j):
    return all(P[n] == ((Fr(1) if (n == 0 and i == j) else Fr(0)),) + (Fr(0),) * 3 for n in range(M + 1))


def main():
    lines = []
    def log(s):
        print(s); lines.append(s)
    bad = [0]
    def check(name, ok, must_fail=False):
        good = (not ok) if must_fail else ok
        if not good:
            bad[0] += 1
        tag = ("PASS (control fails as required)" if good else "FAIL (control did not fail)") if must_fail else ("PASS" if ok else "FAIL")
        log("%s  %s" % (tag, name))
    log("ref_g1_model: symplecticity of the model R-matrix, t-order M = %d, exact in Q(zeta)" % M)
    Fs = solve((Fr(1), Fr(-1), Fr(1), Fr(-1), Fr(-1)))
    log("   Fs_0 = %s + ..." % " ".join("%s t^%d" % (c, n) for n, c in enumerate(Fs[0][:4])))
    check("Fs_0 = 1 - (3/20) t + (9/800) t^2 + ..", Fs[0][:3] == [Fr(1), Fr(-3, 20), Fr(9, 800)])
    check("Fs_b(0) = eps_b = (1,-1,-1,1,1)", [F[0] for F in Fs] == [1, -1, -1, 1, 1])
    R = [Rhat(Fs, i) for i in range(5)]
    check("S_ij(t) = delta_ij for all i, j (to t^%d)" % M, all(is_delta(pairing(R[i], R[j]), i, j) for i in range(5) for j in range(5)))
    Fk = solve((Fr(2), Fr(-1), Fr(1), Fr(-1), Fr(-1, 2)))
    Rk = [Rhat(Fk, i) for i in range(5)]
    check("m1 rescaled system bm' = (2,-1,1,-1,-1/2), same eta", all(is_delta(pairing(Rk[i], Rk[j]), i, j) for i in range(5) for j in range(5)), must_fail=True)
    check("m2 pairing without t -> -t", all(is_delta(pairing(R[i], R[j], sign=1), i, j) for i in range(5) for j in range(5)), must_fail=True)
    check("m3 pairing a+b = 2 mod 5", all(is_delta(pairing(R[i], R[j], shift=2), i, j) for i in range(5) for j in range(5)), must_fail=True)
    log("TOTAL: %d unexpected result(s)" % bad[0])
    with open(os.path.join(HERE, LOGNAME), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("log written to " + LOGNAME)


if __name__ == "__main__":
    main()
