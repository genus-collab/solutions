# ref5_item3.py -- referee of gaps/G5, item 3 (cor:mai-natfact) on actual data.  Evidence only.
# Hatted rows of the formal quintic built from the HPF polynomials (ref5_item4 recursion, free coefficients [z^5m] Ginf, i.e. D = 1),
# the conifold periods of PF4 in the frame (alpha, beta) (my own Frobenius recursion; the 7/10 of varpi_1 re-derived from the
# existence of the logarithmic solution), the forward relations Qh_(a+1) = kh_a A_(e_a) Qh_a (ref_rows docstring) and colour
# covariance.  Then Psi(u) := Qhat(u, 0), Rt := Psi^T eta^-1 Qhat, M = M_0(tau zhat) (model T-system), G := Rt M^-1.  Checked:
# Lemma 3.A, orthonormality, (b) Qhat(0, zhat) = Rhat^mod(zhat), (c) [u^n] G polynomial of zhat-degree <= floor((5n+4)/4),
# G(u,0) = G(0,zhat) = 1, Rt and G symplectic.  Arithmetic mod p = 2^61 - 1 (zeta in F_p), truncation u^U, zhat^K.
# Controls: M_0 on the left; M_0(zhat) without tau; gamma = 7/10 + 1/3; the OU3 relation; f_1 shifted (D = 1 + d z^5: invisible to
# polynomiality, seen by symplecticity at u^8 zhat^10 only).   Usage: py -I -B ref5_item3.py [U] [K]
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref5_item4 import build_P, flat_recursion
from ref5_item1 import ginf

p = (1 << 61) - 1
def inv(a): return pow(a % p, p - 2, p)
def fr(x): return x.numerator % p * inv(x.denominator) % p
ZETA = next(pow(h, (p - 1) // 5, p) for h in range(2, 100) if pow(h, (p - 1) // 5, p) != 1)
NU = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)

def qd(f):                                              # D = -(1-u) d/du on exact power series (lists)
    d = [k * f[k] for k in range(1, len(f))] + [Fr(0)]
    return [-d[k] + (d[k - 1] if k else 0) for k in range(len(f))]

def pf4(f):
    D1 = qd(f); D2 = qd(D1); D3 = qd(D2); D4 = qd(D3)
    Qf = [2 * a + Fr(7, 5) * b + Fr(2, 5) * c + Fr(24, 625) * d for a, b, c, d in zip(D3, D2, D1, f)]
    uD4 = [Fr(0)] + D4[:-1]
    return [x - y + (Qf[k - 1] if k else 0) for k, (x, y) in enumerate(zip(uD4, Qf))]

def frob(U):                                            # f_k = u^k + O(u^3), k = 0, 1, 2
    out = []
    for k0 in range(3):
        f = [Fr(0)] * (U + 4); f[k0] = Fr(1)
        for n in range(0, U + 1):
            r0 = pf4(f)[n]; f[n + 3] = Fr(1); r1 = pf4(f)[n]; f[n + 3] = -r0 / (r1 - r0)
        out.append(f[:U + 1])
    return out

def logcheck(f1, f2, gamma, U):                         # negative powers of PF4(phi ell) - ell PF4(phi), D ell = 1 - 1/u
    phi = {k: a + gamma * b for k, (a, b) in enumerate(zip(f1, f2)) if a + gamma * b}
    def D(el):                                          # el = (b, a) = b + a ell, b, a Laurent dicts
        b, a = el
        def d(f):
            r = {}
            for k, c in f.items():
                r[k - 1] = r.get(k - 1, 0) - k * c; r[k] = r.get(k, 0) + k * c
            return r
        nb = d(b)
        for k, c in a.items():
            nb[k] = nb.get(k, 0) + c; nb[k - 1] = nb.get(k - 1, 0) - c
        return (nb, d(a))
    e0 = ({}, phi); e1 = D(e0); e2 = D(e1); e3 = D(e2); e4 = D(e3)
    tot = {}
    for (b, _), c in ((e4, 'u'), (e3, -2), (e2, Fr(-7, 5)), (e1, Fr(-2, 5)), (e0, Fr(-24, 625))):
        for k, v in b.items():
            if c == 'u':
                tot[k + 1] = tot.get(k + 1, 0) + v
            else:
                tot[k] = tot.get(k, 0) + c * v; tot[k + 1] = tot.get(k + 1, 0) - c * v
    return {k: v for k, v in tot.items() if k < 0 and v != 0}

def umul(a, b, U):
    r = [0] * U
    for i, x in enumerate(a[:U]):
        if x:
            for j, y in enumerate(b[:U - i]):
                r[i + j] = (r[i + j] + x * y) % p
    return r
def uinv(a, U):
    r = [inv(a[0])] + [0] * (U - 1)
    for n in range(1, U):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) * r[0] % p
    return r
def uadd(a, b, c=1): return [(x + c * y) % p for x, y in zip(a, b)]

def rows_colour0(U, K, free, alpha, beta, gamma=Fr(7, 10)):
    A, _, _ = flat_recursion(build_P(), K, free)
    f0, f1, f2 = frob(U + 1)
    vp1 = [a + gamma * b for a, b in zip(f1, f2)]; vp0 = [a + alpha * v + beta * b for a, v, b in zip(f0, vp1, f2)]
    C0l = [fr(x) for x in vp0][:U + 1]; tc = umul([fr(x) for x in vp1], uinv(C0l, U + 1), U + 1)   # one extra order
    dtc = [(k + 1) * tc[k + 1] % p for k in range(U)]
    C1 = [(-dtc[k] + (dtc[k - 1] if k else 0)) % p for k in range(U)]          # D t_c = -(1-u) t_c'
    C0 = C0l[:U]
    Phi = [[fr(A[n].get(4 * n + 1 - 5 * k, Fr(0))) for k in range(U)] for n in range(K)]
    kh = [uinv(C1, U), umul(umul(C0, C0, U), umul(C1, C1, U), U), uinv(C1, U), uinv(C0, U)]
    e = [x - 1 for x in NU]
    Q = [[umul(Phi[n], uinv(C0, U), U) for n in range(K)]]
    for a in range(4):
        F = Q[-1]; G = [list(x) for x in F]
        for k in range(K - 1):
            c = (4 * k - e[a]) * inv(5) % p
            t = [(c * F[k][m] - m * F[k][m]) % p for m in range(U)]
            G[k + 1] = uadd(G[k + 1], uadd(t, [0] + t[:-1], -1))
        Q.append([umul(kh[a], g, U) for g in G])
    F = Q[4]; G = [list(x) for x in F]                                       # closing relation Qh_0 = (1/C0) A_0 Qh_4
    for k in range(K - 1):
        t = [(4 * k * inv(5) * F[k][m] - m * F[k][m]) % p for m in range(U)]
        G[k + 1] = uadd(G[k + 1], uadd(t, [0] + t[:-1], -1))
    closing = all(umul(uinv(C0, U), G[k], U) == Q[0][k] for k in range(K - 1))
    return Q, C0, C1, closing

def model_rows(K):                                      # T-system, colour 0: T_k = -1 + (t/5)(k - 4 theta_t)
    def T(k, f):
        r = [-x for x in f] + [Fr(0)]
        for n, x in enumerate(f):
            r[n + 1] += Fr(k - 4 * n, 5) * x
        return r
    def chain(f0):
        f1 = T(-1, f0); f2 = [-x for x in T(-2, f1)]; f3 = T(2, f2); f4 = [-x for x in T(1, f3)]
        return [f0, f1, f2, f3, f4], [-x for x in T(0, f4)]
    f0 = [Fr(1)] + [Fr(0)] * K
    for n in range(1, K):                               # [t^(n+1)] (closing - Fs_0) = 0 fixes [t^n] Fs_0 (pivot -4n)
        f0[n] = Fr(0); r0 = chain(f0)[1][n + 1] - f0[n + 1]
        f0[n] = Fr(1); r1 = chain(f0)[1][n + 1] - f0[n + 1]
        f0[n] = -r0 / (r1 - r0)
    f0 = f0[:K]
    rows, back = chain(f0)
    return [r[:K] for r in rows], all(back[n] == f0[n] for n in range(K))
