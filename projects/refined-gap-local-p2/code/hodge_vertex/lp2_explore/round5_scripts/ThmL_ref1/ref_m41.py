# Independent referee test of (4.1), written from the report's text only.  Exact mod P = 2^61-1.
# usage: py ref_m41.py tau_num tau_den G LMAX WMAX
import sys, itertools, functools, math
from fractions import Fraction as Fr
sys.setrecursionlimit(100000)
P = (1 << 61) - 1
def inv(a):
    return pow(a % P, P - 2, P)
tn, td, GG = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
tau = tn * inv(td) % P
sigma = tau * (1 + tau) % P
LMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 70
WMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 19
fact = [1] * 400
for i in range(1, 400):
    fact[i] = fact[i - 1] * i % P
def bern(n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1)
    for m in range(1, n + 1):
        B[m] = -sum(Fr(math.comb(m + 1, k)) * B[k] for k in range(m)) / (m + 1)
    return B
BB = bern(LMAX + 2)
def fm(fr):
    return fr.numerator % P * inv(fr.denominator) % P
# log A(z), A=(z/2)/sinh(z/2): log A = -Sum_{j>=1} B_{2j} z^{2j}/(2j (2j)!)
loga = {j: fm(-BB[2 * j] / (2 * j * math.factorial(2 * j))) for j in range(1, LMAX // 2 + 2)}
def smul(a, b, L):
    r = [0] * L
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j >= L:
                    break
                r[i + j] += x * y
    return [v % P for v in r]
def sexp(a, L):
    e = [0] * L
    e[0] = 1
    for n in range(1, L):
        s = 0
        for k in range(1, n + 1):
            if k < len(a):
                s += k * a[k] % P * e[n - k]
        e[n] = s % P * inv(n) % P
    return e
def partitions(n, m=None):
    m = n if m is None else m
    if n == 0:
        yield ()
        return
    for k in range(min(n, m), 0, -1):
        for r in partitions(n - k, k):
            yield (k,) + r
@functools.lru_cache(None)
def chi(nu, mu):
    if not mu:
        return 1 if not nu else 0
    k = mu[0]
    l = len(nu)
    beta = [nu[i] + l - 1 - i for i in range(l)]
    S = set(beta)
    tot = 0
    for b in beta:
        c = b - k
        if c >= 0 and c not in S:
            h = sum(1 for x in beta if c < x < b)
            nb = sorted((S - {b}) | {c}, reverse=True)
            nn = tuple(x for x in (nb[i] - (l - 1 - i) for i in range(l)) if x > 0)
            tot += (-1) ** h * chi(nn, mu[1:])
    return tot
def hooks(nu):
    conj = [sum(1 for x in nu if x > j) for j in range(nu[0])] if nu else []
    return [nu[i] - j + conj[j] - i - 1 for i in range(len(nu)) for j in range(nu[i])]
def Snu(nu, L):
    hs = hooks(nu)
    f2 = sum(j - i for i in range(len(nu)) for j in range(nu[i]))
    a = [0] * L
    if L > 1:
        a[1] = (tau + inv(2)) * f2 % P
    for j in range(1, L):
        if 2 * j < L:
            a[2 * j] = loga[j] * (sum(pow(h, 2 * j, P) for h in hs) % P) % P
    e = sexp(a, L)
    c = 1
    for h in hs:
        c = c * h % P
    c = inv(c)
    return [x * c % P for x in e]
SN = {}
L = LMAX
def getS(nu):
    if nu not in SN:
        SN[nu] = Snu(nu, L)
    return SN[nu]
def zmu(mu):
    z = 1
    for k in set(mu):
        m = mu.count(k)
        z *= k ** m * math.factorial(m)
    return z
Zc = {}
Fc = {}
def Z(mu):
    if mu in Zc:
        return Zc[mu]
    n = sum(mu)
    zi = inv(zmu(mu))
    r = [0] * L
    for nu in partitions(n):
        ch = chi(nu, mu) % P
        if ch:
            s = getS(nu)
            for i in range(L):
                r[i] += ch * s[i]
    r = [v % P * zi % P for v in r]
    Zc[mu] = r
    return r
def submultisets(mu):
    cnt = {}
    for k in mu:
        cnt[k] = cnt.get(k, 0) + 1
    keys = sorted(cnt, reverse=True)
    for ch in itertools.product(*[range(cnt[k] + 1) for k in keys]):
        yield tuple(sorted(sum(([k] * c for k, c in zip(keys, ch)), []), reverse=True))
def minus(mu, rho):
    m = list(mu)
    for k in rho:
        m.remove(k)
    return tuple(m)
def F(mu):
    # coefficient of p_mu in log Z, equal to x^{-|mu|} F(mu)(x)
    if mu in Fc:
        return Fc[mu]
    i = mu[0]
    rest = mu[1:]
    mi = mu.count(i)
    r = [mi * v % P for v in Z(mu)]
    for rho in submultisets(rest):
        if rho == rest:
            continue
        rp = tuple(sorted(rho + (i,), reverse=True))
        co = rp.count(i)
        o = minus(rest, rho)
        s = smul(Z(o), F(rp), L)
        for t in range(L):
            r[t] = (r[t] - co * s[t]) % P
    im = inv(mi)
    r = [v * im % P for v in r]
    Fc[mu] = r
    return r
def Pi(k):
    r = 1
    for a in range(1, k):
        r = r * (k * tau + a) % P
    return r * inv(fact[k - 1]) % P
def aut(mu):
    a = 1
    for k in set(mu):
        a *= math.factorial(mu.count(k))
    return a
def M(g, mu):
    l = len(mu)
    o = 2 * g - 2 + l + sum(mu)
    if o >= L:
        raise OverflowError
    c = F(tuple(sorted(mu, reverse=True)))[o]
    d = pow(sigma, g + l - 1, P)
    for k in mu:
        d = d * Pi(k) % P
    return c * aut(mu) % P * inv(d) % P
class TooBig(Exception):
    pass
@functools.lru_cache(None)
def X(g, U):
    r = len(U)
    J = sum(U)
    if g == 0:
        return fact[J - 1] * inv(fact[J - r + 2]) % P
    N = 3 * g - 3 + r
    if J + N > WMAX:
        raise TooBig
    Ls = N + 1
    un = [0] * Ls
    for i in range(1, Ls):
        un[i] = (-1) ** (i - 1) % P
    a = [(-v) % P for v in un]
    eu = sexp(a, Ls)
    y = smul(un, eu, Ls)
    Gs = [0] * Ls
    yp = [1] + [0] * (Ls - 1)
    for ai in range(0, N + 1):
        m = M(g, tuple(sorted(U + (1,) * ai, reverse=True))) * inv(fact[ai]) % P
        for i in range(Ls):
            Gs[i] = (Gs[i] + m * yp[i]) % P
        yp = smul(yp, y, Ls)
    eJ = sexp([J * v % P for v in a], Ls)
    Q = smul(Gs, eJ, Ls)
    e = 2 * g - 2 + r
    b = [1] + [0] * (Ls - 1)
    for i in range(1, Ls):
        b[i] = b[i - 1] * (-(e + i - 1)) % P * inv(i) % P
    Q = smul(Q, b, Ls)
    return sum(((-1) ** m) * q for m, q in enumerate(Q)) % P
def ff(d, k):
    r = 1
    for i in range(k):
        r *= (d - i)
    return r
def delta(g, U):
    return sum(U) - len(U) - 2 * g + 2
COPT = {'codd': -1}
def cc(b, k):
    if b % 2:
        return 0 if k == 0 else COPT['codd']
    return 2 if k == 0 else 1
def block(B, TB, n):
    U = tuple(sorted(TB + B, reverse=True))
    bb = len(B)
    out = {}
    w = 1
    for i in B:
        w = w * i % P * Pi(i) % P
    for g in range(0, n // 2 + 1):
        for k in range(0, n + 3):
            e = 2 * g - 2 + bb + k
            if e > n:
                break
            c = cc(bb, k)
            if c == 0:
                continue
            f = ff(delta(g, U), k)
            if f == 0:
                continue
            val = c * w % P * pow(sigma, bb + g - 1, P) % P * pow(tau, k, P) % P * inv(fact[k]) % P * (f % P) % P * X(g, U) % P
            out[e] = (out.get(e, 0) + val) % P
    return out
def setparts(lst):
    if not lst:
        yield []
        return
    h = lst[0]
    for p in setparts(lst[1:]):
        yield [[h]] + p
        for i in range(len(p)):
            yield p[:i] + [[h] + p[i]] + p[i + 1:]
def calpha(alpha, k0):
    c = Fr(1)
    for i in set(alpha):
        m = alpha.count(i)
        c *= Fr(-1, i) ** m / math.factorial(m)
    if alpha == (k0,):
        c += 1
    return fm(c)
def rel(k0, R, n):
    tot = 0
    for alpha in partitions(k0):
        ca = calpha(alpha, k0)
        for pi in setparts(list(range(len(alpha)))):
            Bs = [tuple(sorted((alpha[i] for i in b), reverse=True)) for b in pi]
            for assign in itertools.product(range(len(Bs)), repeat=len(R)):
                Ts = [tuple(R[j] for j in range(len(R)) if assign[j] == b) for b in range(len(Bs))]
                cur = {0: 1}
                for B, TB in zip(Bs, Ts):
                    d = block(B, TB, n)
                    new = {}
                    for e1, v1 in cur.items():
                        for e2, v2 in d.items():
                            if e1 + e2 <= n:
                                new[e1 + e2] = (new.get(e1 + e2, 0) + v1 * v2) % P
                    cur = new
                tot = (tot + ca * cur.get(n, 0)) % P
    return tot
def irr(maxw, mp=None):
    mp = maxw if mp is None else mp
    yield ()
    for k in range(min(maxw, mp), 1, -1):
        for r in irr(maxw - k, k):
            yield (k,) + r
if __name__ == "__main__":
    okc = fail = skip = 0
    for mu in [(1,), (2,), (3, 1), (2, 2, 1), (3, 2, 1, 1)]:
        ex = pow(sum(mu), len(mu) - 3, P) if len(mu) >= 3 else inv(pow(sum(mu), 3 - len(mu), P))
        assert M(0, mu) == ex, mu
    print("genus0 values ok", flush=True)
    res = []
    for G in range(1, GG + 1):
        for k0 in range(2, 6):
            for R in irr(4):
                for n in (2 * G, 2 * G + 1):
                    D = sum(R) - len(R) + k0 - 1
                    try:
                        r = rel(k0, R, n)
                    except (TooBig, OverflowError):
                        skip += 1
                        continue
                    if r == 0:
                        okc += 1
                    else:
                        fail += 1
                        print("FAIL", G, k0, R, n, flush=True)
                    if n == 2 * G and G >= 2:
                        T = tuple(sorted(R + (k0,), reverse=True))
                        try:
                            x = X(G, T)
                            res.append((G, T, D, x == 0, D <= 2 * G - 3))
                        except (TooBig, OverflowError):
                            pass
        print("G", G, "relations ok", okc, "fail", fail, "skipped", skip, flush=True)
    nr = sum(1 for t in res if t[4])
    bad = sum(1 for t in res if t[4] and not t[3])
    off = sum(1 for t in res if not t[4])
    offz = sum(1 for t in res if not t[4] and t[3])
    print("in-range points", nr, "nonzero in range", bad, "off-range", off, "off-range zero", offz)
    print("off-range zero list", [t[:3] for t in res if (not t[4]) and t[3]])
    COPT['codd'] = 1
    b = n_c = 0
    X.cache_clear
    for k0, R, G in [(4, (), 2), (5, (), 2), (3, (2,), 2), (2, (3,), 2), (2, (4,), 2), (2, (2, 2), 2)]:
        try:
            v = rel(k0, R, 2 * G)
            n_c += 1
            b += v != 0
        except (TooBig, OverflowError):
            pass
    print("control flipped sign: nonzero", b, "of", n_c)
