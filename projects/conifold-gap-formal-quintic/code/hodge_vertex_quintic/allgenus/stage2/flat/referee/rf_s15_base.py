# -*- coding: utf-8 -*-
"""rf_s15_base.py -- referee of stage2/flat (s.1.5 at genus 3).  OWN code, stdlib only, no project import.

Contents: arithmetic mod p (p = 2^89-1 or 2^61-1, both = 1 mod 5, so a 5th root of unity ZETA exists in F_p); truncated power series
in u (plain lists of length N); Witten-Kontsevich numbers by the DVV recursion (exact over Q); enumeration of stable leg-free graphs
of genus g with canonical forms and automorphism counts (also for graphs whose edges carry two colours, inner/outer); the FROZEN
model (LG model at v = 1) in its idempotent basis, derived only from the Frobenius algebra along the v-line:
    flat basis H^a (a in Z/5), H * H^b = beta_b H^(b+1), beta = (1, -1/v, 1, -1, -1), eta(H^a, H^b) = 5 [a + b = 3 mod 5];
    g^a := N_a v^(-c_a/5) H^a (N = (1,-1,-1,1,1), c = (0,-1,3,2,1)) satisfy g^a g^b = g^(a+b) (asserted below from beta);
    idempotents eps_j = (1/5) Sum_a zeta^(-ja) g^a,  H^a = N_a v^(c_a/5) Sum_j zeta^(ja) eps_j,  Delta_j = 1/eta(eps_j, eps_j) = zeta^(3j) v^(2/5),
    H * eps_j = U'_j eps_j,  U'_j = -zeta^j v^(-1/5).
R-matrix (idempotent frame): fundamental solution Psi R(z) e^(U/z) of z d/dv J = H * J, homogeneous (R_k ~ v^(-4k/5)):
    (A - s) R_k + R_k' = [U', R_(k+1)],  A = Psi^-1 Psi',  s = diag part of A (scalar, = -1/(5v));  solved at v = 1 (off-diagonal
    from the commutator, diagonal from the next order and homogeneity).  The sign of z is a parameter (zsign) fixed by calibration.
"""
import itertools, math
from fractions import Fraction as Fr
from functools import lru_cache

P61 = (1 << 61) - 1
P89 = (1 << 89) - 1
CEXP = (0, -1, 3, 2, 1)
NSG = (1, -1, -1, 1, 1)
BETA_C = (1, -1, 1, -1, -1)          # beta_b = BETA_C[b] v^(-[b == 1])


def inv(a, p): return pow(a % p, p - 2, p)
def fr(x, p): x = Fr(x); return x.numerator % p * inv(x.denominator, p) % p


def zeta5(p):
    assert (p - 1) % 5 == 0
    for g in range(2, 200):
        z = pow(g, (p - 1) // 5, p)
        if z != 1: return z
    raise ValueError


def ratrec(a, p):
    """rational reconstruction (display only)"""
    a %= p; r0, r1, s0, s1 = p, a, 0, 1; B = math.isqrt(p // 2)
    while r1 > B: q = r0 // r1; r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


# ---------------- truncated power series in u (lists of length N, mod p) ----------------
def smul(a, b, p, N):
    out = [0] * N
    for i in range(N):
        x = a[i]
        if x:
            for j in range(N - i):
                y = b[j]
                if y: out[i + j] += x * y
    return [x % p for x in out]


def sinv(a, p, N):
    assert a[0] % p
    i0 = inv(a[0], p); out = [i0] + [0] * (N - 1)
    for n in range(1, N): out[n] = -sum(a[k] * out[n - k] for k in range(1, n + 1)) * i0 % p
    return out


def spow_q(a, r, p, N):
    """a^r, a[0] = 1, r rational:  n b_n = Sum_{k=1..n} (r k - (n - k)) a_k b_{n-k}"""
    assert a[0] % p == 1
    rr = fr(r, p); b = [1] + [0] * (N - 1)
    for n in range(1, N): b[n] = sum((rr * k - (n - k)) * a[k] * b[n - k] for k in range(1, n + 1)) * inv(n, p) % p
    return b


def shift(a, k, N):
    """u^k a (k >= 0), truncated"""
    return ([0] * k + list(a))[:N]


# ---------------- Witten-Kontsevich numbers (DVV recursion, exact) ----------------
def dfact(n):           # (2m+1)!! style double factorial with (-1)!! = 1
    out = 1
    while n > 1: out *= n; n -= 2
    return out


@lru_cache(maxsize=None)
def wk(g, ds):
    """<tau_d1 .. tau_dn>_g for a sorted tuple ds (exact Fraction)"""
    ds = tuple(sorted(ds)); n = len(ds)
    if g < 0 or 2 * g - 2 + n <= 0 or any(d < 0 for d in ds): return Fr(0)
    if sum(ds) != 3 * g - 3 + n: return Fr(0)
    if g == 0 and ds == (0, 0, 0): return Fr(1)
    if g == 1 and ds == (1,): return Fr(1, 24)
    if ds[0] == 0:                                   # string equation
        rest = list(ds[1:]); tot = Fr(0)
        for j in range(len(rest)):
            if rest[j] > 0:
                r2 = rest[:]; r2[j] -= 1; tot += wk(g, tuple(sorted(r2)))
        return tot
    k = ds[-1] - 1; rest = list(ds[:-1]); tot = Fr(0)   # DVV on tau_(k+1)
    for j in range(len(rest)):
        dj = rest[j]; r2 = rest[:]; r2[j] = dj + k
        tot += Fr(dfact(2 * k + 2 * dj + 1), dfact(2 * dj - 1)) * wk(g, tuple(sorted(r2)))
    for r in range(0, k):
        s = k - 1 - r; c = Fr(dfact(2 * r + 1) * dfact(2 * s + 1), 2)
        tot += c * wk(g - 1, tuple(sorted(rest + [r, s])))
        m = len(rest)
        for mask in range(1 << m):
            I = [rest[i] for i in range(m) if mask >> i & 1]; J = [rest[i] for i in range(m) if not mask >> i & 1]
            for g1 in range(0, g + 1):
                a = wk(g1, tuple(sorted(I + [r])))
                if a: tot += c * a * wk(g - g1, tuple(sorted(J + [s])))
    return tot / dfact(2 * k + 3)


# ---------------- graphs ----------------
def canon2(h, ein, eout):
    """canonical form of a graph with genera h and two edge multisets (pairs a <= b)"""
    V = len(h); best = None
    for pi in itertools.permutations(range(V)):
        hh = tuple(h[pi.index(i)] for i in range(V))
        f = lambda E: tuple(sorted(tuple(sorted((pi[a], pi[b]))) for a, b in E))
        key = (hh, f(ein), f(eout))
        if best is None or key < best: best = key
    return best


def aut2(h, ein, eout):
    V = len(h); cnt = 0
    si = sorted(tuple(sorted(e)) for e in ein); so = sorted(tuple(sorted(e)) for e in eout)
    for pi in itertools.permutations(range(V)):
        if any(h[pi[i]] != h[i] for i in range(V)): continue
        if sorted(tuple(sorted((pi[a], pi[b]))) for a, b in si) != si: continue
        if sorted(tuple(sorted((pi[a], pi[b]))) for a, b in so) != so: continue
        cnt += 1
    loc = 1
    for E in (si, so):
        for e in set(E):
            m = E.count(e); loc *= math.factorial(m) * (2 ** m if e[0] == e[1] else 1)
    return cnt * loc


def connected(V, E):
    par = list(range(V))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for a, b in E: par[f(a)] = f(b)
    return len(set(f(x) for x in range(V))) == 1


def stable_graphs(g):
    """all stable connected leg-free graphs of genus g: list of (h, edges), edges sorted pairs (a <= b)"""
    out = {}
    for V in range(1, 2 * g - 1):
        pairs = [(a, b) for a in range(V) for b in range(a, V)]
        for E in range(V - 1, 3 * g - 2):
            sh = g - 1 - E + V
            if sh < 0: continue
            for ed in itertools.combinations_with_replacement(pairs, E):
                if not connected(V, ed): continue
                deg = [0] * V
                for a, b in ed: deg[a] += 1; deg[b] += 1
                for h in itertools.product(range(sh + 1), repeat=V):
                    if sum(h) != sh or any(2 * h[v] - 2 + deg[v] <= 0 for v in range(V)): continue
                    key = canon2(h, ed, ())
                    if key not in out: out[key] = (tuple(h), list(ed))
    return list(out.values())


def contract(h, ed, inner):
    """contract the inner edges (index set): returns (outer genera, outer edges as pairs of outer vertices)"""
    V = len(h); par = list(range(V))
    def f(x):
        while par[x] != x: x = par[x]
        return x
    for e in inner: a, b = ed[e]; par[f(a)] = f(b)
    roots = sorted(set(f(x) for x in range(V))); idx = {r: i for i, r in enumerate(roots)}
    hv = [0] * len(roots); nv = [0] * len(roots); ne = [0] * len(roots)
    for x in range(V): hv[idx[f(x)]] += h[x]; nv[idx[f(x)]] += 1
    for e in inner: ne[idx[f(ed[e][0])]] += 1
    ho = tuple(hv[i] + ne[i] - nv[i] + 1 for i in range(len(roots)))
    eo = [tuple(sorted((idx[f(ed[e][0])], idx[f(ed[e][1])]))) for e in range(len(ed)) if e not in inner]
    return ho, eo


# ---------------- the frozen model in its idempotent basis (v = 1) ----------------
class Frozen(object):
    def __init__(self, p, K, zsign=1):
        self.p = p; z = zeta5(p); self.z = z; Z = [pow(z, k, p) for k in range(5)]; self.Z = Z
        for a in range(5):            # g^1 g^a = g^(a+1) from beta: N_1 N_a beta_a v^((1-c_a)/5) = N_(a+1) v^(-c_(a+1)/5)
            ex = Fr(1 - CEXP[a], 5) - (1 if a == 1 else 0)
            assert NSG[1] * NSG[a] * BETA_C[a] == NSG[(a + 1) % 5] and ex == Fr(-CEXP[(a + 1) % 5], 5)
        for a in range(5):
            for b in range(5):
                if (a + b) % 5 == 3: assert NSG[a] * NSG[b] == 1 and CEXP[a] + CEXP[b] == 2
        self.Delta = [Z[(3 * j) % 5] for j in range(5)]
        self.Up = [(-Z[j]) % p for j in range(5)]
        i25 = inv(25, p)
        A = [[(-i25 * sum(CEXP[a] * Z[((l - j) * a) % 5] for a in range(5))) % p for j in range(5)] for l in range(5)]
        s = A[0][0]
        assert all(A[j][j] == s for j in range(5)) and s == fr(Fr(-1, 5), p)
        self.A = A
        R = [[[int(i == j) for j in range(5)] for i in range(5)]]
        for k in range(K):
            Rk = R[-1]; c = fr(Fr(4 * k, 5), p)
            Mk = [[(sum((A[l][m] if l != m else 0) * Rk[m][j] for m in range(5)) - c * Rk[l][j]) % p for j in range(5)] for l in range(5)]
            Rn = [[0] * 5 for _ in range(5)]
            for l in range(5):
                for j in range(5):
                    if l != j: Rn[l][j] = Mk[l][j] * inv(self.Up[l] - self.Up[j], p) % p
            cc = fr(Fr(5, 4 * (k + 1)), p)
            for j in range(5): Rn[j][j] = cc * sum(A[j][l] * Rn[l][j] for l in range(5) if l != j) % p
            R.append(Rn)
        if zsign == -1: R = [[[x * (-1) ** k % p for x in row] for row in Rk] for k, Rk in enumerate(R)]
        self.R = R; self.K = K
        self.Ri = self.zinv(R)
        self.Psi_inv = [[NSG[x] * Z[(l * x) % 5] % p for x in range(5)] for l in range(5)]   # [l][x]: eps_l-component of H^x

    def zinv(self, R):
        p = self.p; out = [[[int(i == j) for j in range(5)] for i in range(5)]]
        for n in range(1, self.K + 1):
            acc = [[0] * 5 for _ in range(5)]
            for k in range(1, n + 1):
                for i in range(5):
                    for j in range(5): acc[i][j] -= sum(R[k][i][m] * out[n - k][m][j] for m in range(5))
            out.append([[x % p for x in row] for row in acc])
        return out

    def symplectic_defect(self):
        """R(-z)^T eta_c R(z) - eta_c, eta_c = diag(1/Delta): list of (order, max abs residue != 0)"""
        p = self.p; eta = [inv(d, p) for d in self.Delta]; bad = []
        for n in range(1, self.K + 1):
            for i in range(5):
                for j in range(5):
                    s = sum((-1) ** a * self.R[a][m][i] * eta[m] * self.R[n - a][m][j] for a in range(n + 1) for m in range(5)) % p
                    if s: bad.append((n, i, j))
        return bad

    def unit_row(self, i, nmax):
        """eta(1, Psi R(z) eps_i)/eta(1, eps_i): coefficients of z^0..z^nmax (compare F_0 of UNIFORM_REDUCTION s.2 at colour 0)"""
        p = self.p; th = [inv(d, p) for d in self.Delta]
        return [sum(self.R[k][l][i] * th[l] for l in range(5)) * inv(th[i], p) % p for k in range(nmax + 1)]
