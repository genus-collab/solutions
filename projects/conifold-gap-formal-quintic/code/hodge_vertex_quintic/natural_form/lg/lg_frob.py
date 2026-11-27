# -*- coding: utf-8 -*-
"""lg_frob.py -- the Frobenius structure behind the conifold symbol: the Hurwitz space H_{0;1,4} of rational functions
        f(p) = p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)            (rank 5; primitive differential dp)
(notation of Carlet-van de Leur-Posthuma-Shadrin, arXiv:2211.12259, s.1, N = 5; Dubrovin's Hurwitz-space construction).  Exact rational checks:
 (a) FLAT COORDINATES / METRIC: with t = (t1..t5) = (a_1, a_2, a_3, a_4 - a_2^2/8, a_5) the residue metric
         (d', d'') = -(Res_{p=inf} + Res_{p=a_1}) d'f d''f / f' dp      [= sum over the critical points]
     is CONSTANT: eta_{15} = eta_{51} = 1, eta_{24} = eta_{42} = eta_{33} = 1/4, all else 0   (random rational points).
 (b) 3-point functions c(d', d'', d''') = -(Res_inf + Res_{a_1}) d'f d''f d'''f / f' dp:  c_555 = 1/a_5 at every point (the logarithm
     (1/2) t5^2 log t5 of the prepotential sits in the residue coordinate), unit e = d/dt4 (c_{4 alpha beta} = eta_{alpha beta}).
 (c) on the line a = (0,0,0,0,r):  d5 o: d4 -> d5 -> (1/r) d1,  d1 -> 4 d2,  d2 -> d3,  d3 -> d4  (cyclic, (d5 o)^5 = 4/r), E = (5/4) r d5,
     so U = E o has U^5 = (5/4)^5 4 r^4 = 5^5 (r/4)^4, eigenvalues 5 x_c^4 with x_c^5 = r/4 (the critical values).
 (d) charge and spectrum: deg t_alpha = alpha/4, d = 1/2, mu_alpha = (2-d)/2 - deg t_alpha = (3 - alpha)/4 (quoted formula (1.11) of the paper);
     dictionary flat index a of the quintic <-> alpha = a - 1 mod 5 (a = 0 <-> t4 = unit, a = 1 <-> t5 = residue), mu_a = (c_a - 1)/4.
 (e) the R-matrix equation in the canonical frame, for the symbol: (5/4) [C_{n+1}, Z] = (V - n) C_n, Z = diag(zeta^i),
     V_{ji} = (1/5) Sum_a zeta^{(i-j)(a+1)} mu_a (antisymmetric), checked mod p at both primes for the dumped symbol (n <= 5) and for the
     recursion output (n <= 30); together with C_0 = 1 this determines every C_n (the diagonal of C_{n+1} from the next order).
Run:  py -I lg_frob.py
"""
import os, sys, json, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from lg_gm import symbol_rows, CEXP


# ------------------------------------------------------------ polynomials in p (ascending lists of Fractions)
def pmul(A, B):
    if not A or not B: return []
    out = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x == 0: continue
        for j, y in enumerate(B): out[i + j] += x * y
    return out
def padd(A, B):
    n = max(len(A), len(B)); out = [Fr(0)] * n
    for i, x in enumerate(A): out[i] += x
    for i, x in enumerate(B): out[i] += x
    return out
def ptrim(A):
    A = list(A)
    while A and A[-1] == 0: A.pop()
    return A
def pshift(A, a):
    """A(v + a) as a polynomial in v (Horner)"""
    out = []
    for c in reversed(A):
        out = padd(pmul(out, [Fr(a), Fr(1)]), [c])
    return ptrim(out)
def res_inf(N, M):
    """Res_{p=inf} N/M dp = -[p^-1] of the expansion at infinity = -(coefficient of p^{deg M - 1} in N mod M)/lead(M)"""
    N = ptrim(N); M = ptrim(M); dM = len(M) - 1; R = list(N)
    while len(R) - 1 >= dM and R:
        c = R[-1] / M[-1]; k = len(R) - 1 - dM
        for i, m in enumerate(M): R[i + k] -= c * m
        R = ptrim(R)
    return -(R[dM - 1] if len(R) > dM - 1 else Fr(0)) / M[-1]
def res_at(N, Q, a, m):
    """Res_{p=a} N/((p-a)^m Q) dp = [v^{m-1}] N(v+a)/Q(v+a)   (Q(a) != 0)"""
    Nv = pshift(N, a); Qv = pshift(Q, a); out = []
    for k in range(m):
        s = (Nv[k] if k < len(Nv) else Fr(0)) - sum(out[j] * (Qv[k - j] if k - j < len(Qv) else 0) for j in range(k))
        out.append(s / Qv[0])
    return out[m - 1]


def structure(a1, a2, a3, a4, a5):
    """(eta 5x5, c 5x5x5) in the coordinates t = (a1, a2, a3, a4 - a2^2/8, a5) at the point"""
    D = [Fr(-a1), Fr(1)]; D2 = pmul(D, D)
    Q = padd(pmul([Fr(a3), Fr(2 * a2), Fr(0), Fr(4)], D2), [Fr(-a5)])                 # f' = Q / D^2
    tau = [[Fr(a5)], pmul([Fr(a2) / 4, Fr(0), Fr(1)], D2), pmul([Fr(0), Fr(1)], D2), D2, D]    # d_alpha f = tau_alpha / D^2
    def pair(N, m):                                                                  # -(Res_inf + Res_{a1}) N/(D^m Q)
        M = Q
        for _ in range(m): M = pmul(M, D)
        return -(res_inf(N, M) + res_at(N, Q, a1, m))
    eta = [[pair(pmul(tau[i], tau[j]), 2) for j in range(5)] for i in range(5)]
    c = [[[pair(pmul(pmul(tau[i], tau[j]), tau[k]), 4) for k in range(5)] for j in range(5)] for i in range(5)]
    return eta, c


def main():
    FAIL = []
    def check(name, ok, detail=""):
        print(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""), flush=True)
        if not ok: FAIL.append(name)
    ETA = [[Fr(0)] * 5 for _ in range(5)]
    ETA[0][4] = ETA[4][0] = Fr(1); ETA[1][3] = ETA[3][1] = ETA[2][2] = Fr(1, 4)
    rng = random.Random(20261110); pts = []
    for _ in range(4): pts.append(tuple(Fr(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(4)) + (Fr(rng.randint(1, 9), rng.randint(1, 7)),))
    oka = okb = oku = True
    for pt in pts:
        eta, c = structure(*pt)
        oka = oka and eta == ETA
        okb = okb and c[4][4][4] == 1 / pt[4]
        oku = oku and all(c[3][i][j] == eta[i][j] for i in range(5) for j in range(5))
    check("(a) metric constant in t = (a1, a2, a3, a4 - a2^2/8, a5): eta_15 = 1, eta_24 = eta_33 = 1/4 at 4 random rational points", oka)
    # control: without the correction a4 - a2^2/8 the metric is NOT constant (d/da2 at fixed a4 = d2 - (a2/4) d4)
    eta, c = structure(*pts[0]); a2 = pts[0][1]
    naive22 = eta[1][1] - 2 * (a2 / 4) * eta[1][3] + (a2 / 4) ** 2 * eta[3][3]
    check("(a) control: in the naive coordinate a4 the entry (d/da2, d/da2) = -a2/8 != 0", naive22 == -a2 / 8 and a2 != 0, str(naive22))
    check("(b) c_555 = 1/a5 at the 4 points (prepotential contains (1/2) t5^2 log t5)", okb)
    check("(b) e = d/dt4 is the unit: c_{4 alpha beta} = eta_{alpha beta}", oku)
    # (c) the line
    r = Fr(3, 7); eta, c = structure(0, 0, 0, 0, r)
    inv = [[Fr(0)] * 5 for _ in range(5)]
    inv[0][4] = inv[4][0] = Fr(1); inv[1][3] = inv[3][1] = inv[2][2] = Fr(4)
    m5 = [[sum(c[4][b][dd] * inv[dd][g] for dd in range(5)) for g in range(5)] for b in range(5)]      # d5 o d_b = Sum_g m5[b][g] d_g
    exp = [[Fr(0)] * 5 for _ in range(5)]
    exp[3][4] = Fr(1); exp[4][0] = 1 / r; exp[0][1] = Fr(4); exp[1][2] = Fr(1); exp[2][3] = Fr(1)
    check("(c) on the line: d5 o d4 = d5, d5 o d5 = (1/r) d1, d5 o d1 = 4 d2, d5 o d2 = d3, d5 o d3 = d4", m5 == exp)
    P = [[Fr(int(i == j)) for j in range(5)] for i in range(5)]
    for _ in range(5): P = [[sum(P[i][k] * m5[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
    check("(c) (d5 o)^5 = (4/r) id, hence U = (5/4) r d5 o has U^5 = 5^5 (r/4)^4: eigenvalues 5 x_c^4, x_c^5 = r/4",
          P == [[(4 / r) * int(i == j) for j in range(5)] for i in range(5)])
    # (d) degrees / spectrum bookkeeping
    deg = [Fr(al, 4) for al in range(1, 6)]; d = Fr(1, 2)
    okd = all(deg[i] + deg[j] == 2 - d for i in range(5) for j in range(5) if ETA[i][j] != 0)
    mu_al = [Fr(2 - d) / 2 - deg[i] for i in range(5)]
    mu_a = [Fr(CEXP[a] - 1, 4) for a in range(5)]
    okm = all(mu_a[a] == mu_al[(a - 1) % 5 - 1 if (a - 1) % 5 else 4] for a in range(5))
    check("(d) deg t_alpha = alpha/4 pair to 2 - d with d = 1/2 on the support of eta; mu_alpha = (3-alpha)/4; mu_a = (c_a-1)/4 under alpha = a-1 mod 5",
          okd and okm, "mu_a = %s" % [str(x) for x in mu_a])
    # (e) canonical-frame R-matrix equation mod p
    N = 30; rows = symbol_rows(N + 1)
    for tag in ("89", "61"):
        dmp = json.load(open(os.path.join(HERE, "lg_symbol_%s.json" % tag))); p = dmp['p']; z = dmp['zeta']; Rh = dmp['Rh']
        f = lambda x: Fr(x).numerator % p * pow(Fr(x).denominator % p, p - 2, p) % p
        inv5 = pow(5, p - 2, p); zp = [pow(z, k, p) for k in range(5)]
        V = [[sum(zp[((i - j) * (a + 1)) % 5] * f(mu_a[a]) for a in range(5)) * inv5 % p for i in range(5)] for j in range(5)]
        okV = all((V[j][i] + V[i][j]) % p == 0 for i in range(5) for j in range(5))
        def Cmat(n, src):
            if src == 'rec': return [[sum(zp[(j * (4 - a) + i * (a + 1 - n)) % 5] * f(rows[a][n]) for a in range(5)) * inv5 % p for i in range(5)] for j in range(5)]
            return [[sum(Rh[b][j][0] * e * Rh[a][i][n] for (a, b, e) in dmp['eta_inv']) % p for i in range(5)] for j in range(5)]
        for src, nmax in (('dump', 5), ('rec', N)):
            bad = 0
            for n in range(nmax + 1):
                C0 = Cmat(n, src); C1 = Cmat(n + 1, src)
                for j in range(5):
                    for i in range(5):
                        lhs = f(Fr(5, 4)) * C1[j][i] * (zp[i] - zp[j]) % p
                        rhs = (sum(V[j][k] * C0[k][i] for k in range(5)) - n * C0[j][i]) % p
                        if lhs != rhs: bad += 1
            check("(e) p = 2^%s-1, %s symbol: (5/4)[C_{n+1}, Z] = (V - n) C_n for n <= %d; V antisymmetric %s" % (tag, src, nmax, okV), bad == 0 and okV, "%d bad entries" % bad)
    print("ALL PASS" if not FAIL else "FAILURES: %s" % FAIL)


if __name__ == "__main__":
    main()
