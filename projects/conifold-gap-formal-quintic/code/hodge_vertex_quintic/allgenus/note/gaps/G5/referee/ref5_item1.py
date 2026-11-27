# ref5_item1.py -- referee of gaps/G5, item 1 (thm:con-Pa: flat identity, comparison, polynomiality).  Evidence only.
# Exact over F_p, p = 2^61 - 1, at random points L0, with JETS L = L0 + eps (mod eps^K) to differentiate in L.
# The full-Stirling saddle series c^[g](z; L) = SPn_L[S; g A]/A(L,0), S'' = -5/(sigma(sigma^5-1)), A = half-density x Stirling,
# is computed by the MORSE-COORDINATE method (y with S(L+t) - S(L) = (S''(L)/2) y^2, Lagrange inversion for t(y)), not by Wick
# contraction of e^(S_3/z) as in g5_flat_hatted.py: an independent algorithm (equal by lem:fsp-cov).
# Checks: (a) c_n(L0) = L0^(4n) ph_n(L0^-5), ph_n the HPF polynomials (ref5_item4 recursion) with free coefficients [z^5m] Ginf;
#   (b) (1.1) (zD + L)(L c^[g]) = L c^[sigma g] for g = sigma^m (D L = (L^6 - L)/5); (c) PF5^L (L c) = 0.
# Controls (must fail): (a) with f_1 shifted, with L0^(4n-1); (b) for g = q(L) (an amplitude with a NON-constant coefficient, which
#   GAPS_G5 s.1.3 Step 2 admits but for which (1.1) is false); (c) with the half-density dropped.
# Usage: py -I -B ref5_item1.py [N] [npts]
import os, sys, random
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref5_item4 import build_P, flat_recursion

p = (1 << 61) - 1
def inv(a): return pow(a % p, p - 2, p)
def fr(x): return x.numerator % p * inv(x.denominator) % p

# ---- jets: lists of length K (Taylor coefficients in eps)
def jmul(a, b):
    K = min(len(a), len(b)); r = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i):
                r[i + j] = (r[i + j] + a[i] * b[j]) % p
    return r
def jadd(a, b, c=1):
    K = min(len(a), len(b)); return [(a[i] + c * b[i]) % p for i in range(K)]
def jsc(a, c): return [x * c % p for x in a]
def jinv(a):
    K = len(a); r = [0] * K; r[0] = inv(a[0])
    for n in range(1, K):
        r[n] = -sum(a[k] * r[n - k] for k in range(1, n + 1)) * r[0] % p
    return r
def jcon(c, K): return [c % p] + [0] * (K - 1)
def jder(a): return [(i + 1) * a[i + 1] % p for i in range(len(a) - 1)]      # d/dL of the jet (length drops by one)

# ---- series in t with jet coefficients, truncated at t^T
def smul(a, b, T):
    K = len(a[0]); r = [[0] * K for _ in range(T + 1)]
    for i in range(min(len(a), T + 1)):
        if any(a[i]):
            for j in range(min(len(b), T + 1 - i)):
                r[i + j] = jadd(r[i + j], jmul(a[i], b[j]))
    return r
def sinv(a, T):
    r = [jinv(a[0])]
    for n in range(1, T + 1):
        acc = [0] * len(a[0])
        for k in range(1, min(n, len(a) - 1) + 1):
            acc = jadd(acc, jmul(a[k], r[n - k]))
        r.append(jsc(jmul(acc, r[0]), -1))
    return r
def ssqrt1(a, T):                                       # sqrt of a series with constant term 1
    r = [jcon(1, len(a[0]))]; h = inv(2)
    for n in range(1, T + 1):
        acc = a[n] if n < len(a) else [0] * len(a[0])
        for i in range(1, n):
            acc = jadd(acc, jmul(r[i], r[n - i]), -1)
        r.append(jsc(acc, h))
    return r
def sexp0(a, T):                                        # exp of a series with zero constant term (coefficients jets)
    K = len(a[0]); r = [jcon(1, K)]
    for n in range(1, T + 1):
        acc = [0] * K
        for k in range(1, min(n, len(a) - 1) + 1):
            acc = jadd(acc, jsc(jmul(a[k], r[n - k]), k))
        r.append(jsc(acc, inv(n)))
    return r
def compose(F, ty, T):                                  # F(t(y)), t(y) = O(y)
    K = len(F[0]); R = [list(F[min(T, len(F) - 1)])] + [[0] * K for _ in range(T)]
    for k in range(min(T, len(F) - 1) - 1, -1, -1):
        R = smul(R, ty, T); R[0] = jadd(R[0], F[k])
    return R

def bernoulli(n):
    B = [Fr(1)] + [Fr(0)] * n
    for m in range(1, n + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B

def saddle(Lj, N, gpolys, variant=''):
    """c^[g]_n (n <= N) at L = Lj (a jet) for each g in gpolys (polynomials in sigma with jet coefficients, low degree first)."""
    K = len(Lj); T = 2 * N + 2; TT = T + N + 4
    sig = [Lj, jcon(1, K)] + [[0] * K for _ in range(TT - 1)]
    s5 = [jcon(1, K)] + [[0] * K for _ in range(TT)]
    for _ in range(5):
        s5 = smul(s5, sig, TT)
    s5m1 = [jadd(s5[0], jcon(1, K), -1)] + s5[1:]
    h = [jsc(x, -5) for x in sinv(smul(sig, s5m1, TT), TT)]             # S''(L + t)
    H = h[0]; C = jsc(jinv(H), -1)
    # w(t) = sqrt(1 + (2/H) sum_{k>=3} s_k t^(k-2)),  s_k = h[k-2]/(k(k-1))
    q = [jcon(1, K)] + [jmul(jsc(h[k - 2], 2 * inv(k * (k - 1))), jinv(H)) for k in range(3, T + 3)]
    w = ssqrt1(q, T + 1); phi = sinv(w, T + 1)
    pw = [jcon(1, K)] + [[0] * K for _ in range(T + 1)]; ty = [[0] * K for _ in range(T + 2)]
    for n in range(1, T + 2):                                           # Lagrange: [y^n] t = (1/n) [t^(n-1)] phi^n
        pw = smul(pw, phi, T + 1); ty[n] = jsc(pw[n - 1], inv(n))
    ty = ty[:T + 1]
    dty = [jsc(ty[i + 1], i + 1) for i in range(T)] + [[0] * K]
    # amplitude: half-density sqrt((sigma/L)(L^5-1)/(sigma^5-1)) and Stirling exponent
    rr = smul([jmul(x, jinv(Lj)) for x in sig], [jmul(x, s5m1[0]) for x in sinv(s5m1, TT)], TT)
    hd = ssqrt1(rr, TT) if variant != 'nohd' else [jcon(1, K)] + [[0] * K for _ in range(TT)]
    g5 = smul([jsc(x, 5) for x in smul(smul(sig, sig, TT), smul(sig, sig, TT), TT)], sinv(s5m1, TT), TT)
    inv5s = sinv([jsc(x, 5) for x in sig], TT)
    B = bernoulli(2 * N + 2); E = {}
    for k in range(1, N + 1):
        i = 2 * k - 1
        if i > N:
            break
        bk = fr(B[2 * k] / (2 * k * (2 * k - 1)))
        powi = [jcon(1, K)] + [[0] * K for _ in range(TT)]
        for _ in range(i):
            powi = smul(powi, inv5s, TT)
        sj = [jsc(g5[m + i - 1], (-1) ** (i - 1) * comb(m + i - 1, i - 1)) for m in range(TT - i)]   # sum_j (sigma-zeta^j)^-i
        E[i] = [jsc(jadd(powi[m], sj[m], -1), bk) for m in range(TT - i)]
    Ar = [[jcon(1, K)] + [[0] * K for _ in range(TT - N - 1)]]          # exp(sum_i z^i E_i): r A_r = sum_i i E_i A_(r-i)
    for r in range(1, N + 1):
        acc = [[0] * K for _ in range(TT - N)]
        for i, Ei in E.items():
            if i <= r:
                acc = [jadd(x, jsc(y, i)) for x, y in zip(acc, smul(Ei, Ar[r - i], TT - N - 1))]
        Ar.append([jsc(x, inv(r)) for x in acc])
    Ar = [smul(hd, a, TT - N - 1) for a in Ar]
    Ary = [smul(compose(a, ty, T), dty, T) for a in Ar]                  # A_r(t(y)) t'(y)
    out = []
    for g in gpolys:
        gs = [[0] * K for _ in range(T + 1)]
        for d, c in enumerate(g):                                       # g(L + t(y)) = sum_d c_d (L + t)^d
            pw = [jcon(1, K)] + [[0] * K for _ in range(T)]
            for _ in range(d):
                pw = smul(pw, [Lj] + ty[1:], T)
            gs = [jadd(x, jmul(c, y)) for x, y in zip(gs, pw)]
        cn = []
        for n in range(N + 1):
            tot = [0] * K
            for r in range(n + 1):
                m = n - r; integ = smul(gs, Ary[r], 2 * m)[2 * m]
                df = 1
                for o in range(1, 2 * m, 2):
                    df = df * o % p
                Cm = jcon(1, K)
                for _ in range(m):
                    Cm = jmul(Cm, C)
                tot = jadd(tot, jsc(jmul(integ, Cm), df))
            cn.append(tot)
        out.append(cn)
    return out


def ginf(nmax):                                         # Ginf(z) = exp(sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))
    B = bernoulli(5 * (nmax // 5) + 2)
    def bpoly(k, x): return sum(comb(k, j) * B[j] * x ** (k - j) for j in range(k + 1))
    lg = [Fr(0)] * (nmax + 1)
    for l in range(1, nmax // 5 + 1):
        lg[5 * l] = Fr((-1) ** (l + 1)) * bpoly(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1))
    g = [Fr(1)] + [Fr(0)] * nmax                        # exp by n g_n = sum k lg_k g_(n-k)
    for n in range(1, nmax + 1):
        g[n] = sum(k * lg[k] * g[n - k] for k in range(1, n + 1)) / n
    return g

def evalA(An, L0):                                      # A_n(L0) mod p
    return sum(fr(c) * pow(L0, m, p) if m >= 0 else fr(c) * pow(inv(L0), -m, p) for m, c in An.items()) % p

def opD(F, Lj):                                         # D on a jet: ((L^6 - L)/5) d/dL
    L6 = Lj
    for _ in range(5):
        L6 = jmul(L6, Lj)
    DL = jsc(jadd(L6, Lj, -1), inv(5))
    return jmul(DL[:len(F) - 1], jder(F))

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    npts = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    rng = random.Random(20261115)
    G = ginf(N); free = {5 * m: G[5 * m] for m in range(1, N // 5 + 1)}
    print('ref5_item1: N = %d, %d random points, F_p with p = 2^61 - 1; free coefficients [z^5m] Ginf =' % (N, npts),
          [str(free[k]) for k in sorted(free)])
    A, _, ok = flat_recursion(build_P(), N, free)
    Abad, _, _ = flat_recursion(build_P(), N, {k: v + Fr(1, 7) for k, v in free.items()})
    pts = []
    while len(pts) < npts:
        L0 = rng.randrange(2, p - 1)
        if L0 * (pow(L0, 5, p) - 1) % p:
            pts.append(L0)
    ra = rb = rc = True; c1 = c2 = c3 = False
    for L0 in pts:
        c = saddle([L0], N, [[[1]]])[0]
        ra &= all(c[n][0] == evalA(A[n], L0) * inv(L0) % p for n in range(N + 1))
        c1 |= any(c[n][0] != evalA(Abad[n], L0) * inv(L0) % p for n in range(N + 1))
        c2 |= any(c[n][0] != evalA(A[n], L0) * inv(L0) ** 2 % p for n in range(1, N + 1))
    print('(a) c_n(L0) = L0^(4n) ph_n(L0^-5) for n <= %d at all points: %s' % (N, ra))
    print('    controls: f_m shifted by 1/7 -> mismatch: %s; L0^(4n-1) -> mismatch: %s (both must be True)' % (c1, c2))
    q = lambda Lj: jsc(jadd(jcon(1, len(Lj)), jmul(jinv(Lj), jmul(jinv(Lj), jmul(jinv(Lj), jmul(jinv(Lj), jinv(Lj))))), -1), inv(3125))
    N2 = min(N, 8); okb = True; failq = True
    for L0 in pts[:6]:
        Lj = [L0, 1]
        gl = [[[0, 0]] * m + [[1, 0]] for m in range(6)] + [[q(Lj)], [[0, 0], q(Lj)]]
        cs = saddle(Lj, N2, gl)
        for m in range(5):                               # (zD + L)(L c^[sigma^m]) = L c^[sigma^(m+1)]
            for n in range(N2 + 1):
                lhs = jmul(jmul(Lj, Lj), cs[m][n])[0]
                if n >= 1:
                    lhs = (lhs + opD(jmul(Lj, cs[m][n - 1]), Lj)[0]) % p
                okb &= lhs == jmul(Lj, cs[m + 1][n])[0]
        bad = 0
        for n in range(1, N2 + 1):                      # control g = q(L): (zD + L)(L c^[q]) - L c^[sigma q]
            lhs = (jmul(jmul(Lj, Lj), cs[6][n])[0] + opD(jmul(Lj, cs[6][n - 1]), Lj)[0]) % p
            bad += lhs != jmul(Lj, cs[7][n])[0]
        failq &= bad > 0
    print('(b) (1.1) for g = sigma^m, m = 0..4, n <= %d, %d points: %s;  control g = q(L): (1.1) fails: %s' % (N2, len(pts[:6]), okb, failq))
    N3 = min(N, 7); okc = True; failc = True
    for L0 in pts[:3]:
        Lj = [L0, 1, 0, 0, 0, 0]
        for variant in ('', 'nohd'):
            c = saddle(Lj, N3, [[[1, 0, 0, 0, 0, 0]]], variant)[0]
            F = [jmul(Lj, x) for x in c]
            def step(F, a, b):                          # a (zD + L) F + b z F
                return [jadd(jsc(jmul(Lj[:len(F[n]) - 1], F[n][:len(F[n]) - 1]), a), (jadd(jsc(opD(F[n - 1], Lj), a), jsc(F[n - 1][:len(F[n]) - 1], b)) if n else [0] * (len(F[n]) - 1)))
                        for n in range(len(F))]
            X = F
            for _ in range(5):
                X = step(X, 1, 0)
            Y = F
            for k in range(1, 6):
                Y = step(Y, 5, k)
            qq = q(Lj)[:1]
            res = [(X[n][0] - F[n][0] - qq[0] * Y[n][0]) % p for n in range(N3 + 1)]   # q(L) includes 1/3125
            if variant == '':
                okc &= all(r == 0 for r in res)
            else:
                failc &= any(r != 0 for r in res)
    print('(c) PF5^L (L c) = 0 through z^%d, %d points: %s;  control without half-density: PF5 fails: %s' % (N3, len(pts[:3]), okc, failc))

if __name__ == '__main__':
    main()
