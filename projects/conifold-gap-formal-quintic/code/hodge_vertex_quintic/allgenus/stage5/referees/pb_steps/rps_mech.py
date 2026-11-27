"""rps_mech.py -- referee pb_steps: CFU map, reduction modulo exact forms, the residue side of P(b) step 6.

Own code (stdlib only, no project module).  Polynomials in T = X^5 (or Y^5) are lists of Fractions; u-series of them are
lists over the u-order.  CFU: Y = X y(T), T = X^5, with  V/y + T y^4/4 = Sum_m (u^(m-1)/m)(1 + T^m/(5m-1)),  solved order by
order: (T - 1) y_k + V_k = R_k(T), V_k = R_k(1).  Inverse X = Y x(Y^5).  beta(Y) dY = B(X(Y)) X'(Y) dY.  Reduction:
Y^n == (n-3) w Y^(n-4) + V Y^(n-5) (n >= 4), from  w d(Y^m e^(-f/w)) = (m w Y^(m-1) + V Y^(m-2) - Y^(m+3)) e^(-f/w) dY.
"""
from fractions import Fraction as Fr
from rps_core import add, fact

def padd(a, b, c=1):
    n = max(len(a), len(b))
    out = [Fr(0)] * n
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += c * x
    return out

def pmul(a, b):
    if not a or not b:
        return []
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return out

def smul(a, b, N):
    """u-series of T-polynomials, product mod u^(N+1)."""
    out = [[] for _ in range(N + 1)]
    for i in range(min(len(a), N + 1)):
        if a[i]:
            for j in range(0, min(len(b), N + 1 - i)):
                if b[j]:
                    out[i + j] = padd(out[i + j], pmul(a[i], b[j]))
    return out

def sinv(a, N):
    """1/a for a u-series of T-polynomials with a[0] = [1]."""
    assert a[0] == [Fr(1)] or (len(a[0]) >= 1 and a[0][0] == 1 and all(x == 0 for x in a[0][1:]))
    inv = [[Fr(1)]] + [[] for _ in range(N)]
    for n in range(1, N + 1):
        acc = []
        for i in range(1, n + 1):
            if i < len(a) and a[i]:
                acc = padd(acc, pmul(a[i], inv[n - i]), -1)
        inv[n] = acc
    return inv

def spow(a, e, N):
    out = [[Fr(1)]] + [[] for _ in range(N)]
    for _ in range(e):
        out = smul(out, a, N)
    return out

def cfu(N):
    """returns (y, V): y u-series of T-polys (y[0] = [1]), V list of Fractions (V[0] = 1)."""
    y = [[Fr(1)]] + [[] for _ in range(N)]
    V = [Fr(1)] + [Fr(0)] * N
    for k in range(1, N + 1):
        inv = sinv(y, k)
        y4 = spow(y, 4, k)
        known = []
        for l in range(0, k + 1):
            if V[l] and inv[k - l]:
                known = padd(known, [V[l] * x for x in inv[k - l]])
        known = padd(known, [Fr(0)] + [x / 4 for x in y4[k]])
        Wk = [Fr(1, k + 1)] + [Fr(0)] * k + [Fr(1, (k + 1) * (5 * k + 4))]
        R = padd(Wk, known, -1)
        Vk = sum(R, Fr(0))
        V[k] = Vk
        R[0] -= Vk
        q = [Fr(0)] * (len(R) - 1)          # (R(T) - R(1))/(T - 1), synthetic division
        acc = Fr(0)
        for i in range(len(R) - 1, 0, -1):
            acc += R[i]
            q[i - 1] = acc
        y[k] = q
    return y, V

def compose(p, t, N):
    """p(t) for a T-polynomial p and a u-series t of T-polynomials, mod u^(N+1)."""
    out = [[] for _ in range(N + 1)]
    pw = [[Fr(1)]] + [[] for _ in range(N)]
    for d, c in enumerate(p):
        if d > 0:
            pw = smul(pw, t, N)
        if c:
            for i in range(N + 1):
                if pw[i]:
                    out[i] = padd(out[i], [c * x for x in pw[i]])
    return out

def scompose(y, t, N):
    """y(t) = Sum_k u^k y_k(t) for a u-series y of T-polys and a u-series t of T-polys."""
    out = [[] for _ in range(N + 1)]
    for k in range(N + 1):
        if y[k]:
            c = compose(y[k], t, N - k)
            for i in range(N - k + 1):
                if c[i]:
                    out[i + k] = padd(out[i + k], c[i])
    return out

def cfu_inverse(y, N):
    """x(T) with X = Y x(Y^5): T_X = T y(T_X)^(-5) by u-adic iteration, then x = 1/y(T_X)."""
    TX = [[Fr(0), Fr(1)]] + [[] for _ in range(N)]
    for _ in range(N + 1):
        yt = scompose(y, TX, N)
        TX = smul([[Fr(0), Fr(1)]], sinv(spow(yt, 5, N), N), N)
    return sinv(scompose(y, TX, N), N)

def sderiv(a):
    return [[i * c for i, c in enumerate(p)][1:] for p in a]

def beta_Y(B, x, N, noJac=False):
    """beta(Y) = B(X(Y)) X'(Y) as dict {(i, n, j): coeff of u^i Y^n w^j}; X = Y x(Y^5), X' = x(T) + 5 T x'(T)."""
    dX = [padd(x[i], [Fr(0)] + [5 * c for c in sderiv(x)[i]]) for i in range(N + 1)]
    if noJac:
        dX = [[Fr(1)]] + [[] for _ in range(N)]
    pmax = max((p for n in range(N + 1) for (p, q) in B[n]), default=0)
    xp = [[[Fr(1)]] + [[] for _ in range(N)]]
    for _ in range(pmax):
        xp.append(smul(xp[-1], x, N))
    out = {}
    cache = {}
    for i in range(N + 1):
        for (p, q), c in B[i].items():
            if p not in cache:
                cache[p] = smul(xp[p], dX, N)
            s = cache[p]
            for i2 in range(0, N - i + 1):
                for d, e in enumerate(s[i2]):
                    if e:
                        add(out, (i + i2, p + 5 * d, q), c * e)
    return out

def reduce_forms(beta, V, N):
    """Y^n -> (n-3) w Y^(n-4) + V Y^(n-5) for n >= 4; returns pi[m] = {(i, j): coeff}, m = -1..3."""
    work = dict(beta)
    while True:
        top = max((n for (i, n, j) in work), default=-1)
        if top < 4:
            break
        for (i, n, j) in [k for k in work if k[1] == top]:
            c = work.pop((i, n, j))
            add(work, (i, n - 4, j + 1), (n - 3) * c)
            for l in range(0, N - i + 1):
                if V[l]:
                    add(work, (i + l, n - 5, j), c * V[l])
    pi = {m: {} for m in range(-1, 4)}
    for (i, n, j), c in work.items():
        pi[n][(i, j)] = c
    return pi

def vpowers(V, emax, N):
    """[V^e mod u^(N+1)] for e = 0..emax."""
    out = [[Fr(1)] + [Fr(0)] * N]
    for _ in range(emax):
        prev = out[-1]
        out.append([sum(prev[a] * V[n - a] for a in range(n + 1)) for n in range(N + 1)])
    return out

def resid_orig(pi, Vp, E, beta, imax):
    """[u^E z^beta] of Sum_m pi_m s^(-m) I^C_m(v(u), z) using pi only at u-orders <= imax (key lemma of PB s.2.7).
    term: p_(i,j) u^i w^j s^(-m) (-1/(4z))^i'/i'! (-v/z)^e/e!,  e = m + 4i' + 1,  [u^l] V^e,  l = E - i + (4 beta - 1)/5."""
    tot = Fr(0)
    if (4 * beta - 1) % 5:
        return tot
    sh = (4 * beta - 1) // 5
    for m, d in pi.items():
        for (i, j), c in d.items():
            if i > imax:
                continue
            l = E - i + sh
            r = j - m - 1 - beta
            if l < 0 or r < 0 or r % 5:
                continue
            ip = r // 5
            e = m + 4 * ip + 1
            if e < 0 or l >= len(Vp[0]):
                continue
            tot += c * Fr((-1) ** (ip + e), 4 ** ip * fact(ip) * fact(e)) * Vp[e][l]
    return tot

def resid_hat(pi, Vp, N, jmin, jmax):
    """hatted Sum_m pi_m Res_Y[Y^m e^(-f/w) dY] at u^N, w^j in [jmin, jmax];
    Res_Y[Y^m e^(-f/w) dY] = Sum_l (-V/w)^e/e! (-1/(4w))^l/l!, e = 4l + m + 1."""
    out = {}
    for m, d in pi.items():
        for (i, j1), c in d.items():
            if i > N:
                continue
            n2 = N - i
            l = 0
            while True:
                e = 4 * l + m + 1
                j = j1 - e - l
                if j < jmin:
                    break
                if e >= 0 and j <= jmax:
                    add(out, j, c * Fr((-1) ** (e + l), fact(e) * 4 ** l * fact(l)) * Vp[e][n2])
                l += 1
    return out
