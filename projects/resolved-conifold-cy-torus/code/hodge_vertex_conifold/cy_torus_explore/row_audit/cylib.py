# Generic-weight rows of the resolved conifold near Q = 1, t = log Q. Fresh code, evidence only.
# Arithmetic in F_p. Laurent series in t: (val, list) with fixed length L.
from fractions import Fraction as Fr
from math import factorial, comb
P = 2**127 - 1


def inv(a):
    return pow(a % P, P - 2, P)


def fp(x):
    x = Fr(x)
    return x.numerator % P * inv(x.denominator) % P


def bern(N):
    B = [Fr(0)] * (N + 1)
    B[0] = Fr(1)
    for m in range(1, N + 1):
        B[m] = -sum(comb(m + 1, k) * B[k] for k in range(m)) / (m + 1)
    return B


BERN = bern(120)


def bk(k):
    return fp(BERN[2 * k] / (2 * k * (2 * k - 1)))


class Ser:
    L = 16


def mul(a, b, L=None):
    L = L or Ser.L
    r = [0] * L
    for i, x in enumerate(a[:L]):
        if x:
            for j, y in enumerate(b[:L - i]):
                r[i + j] += x * y
    return [v % P for v in r]


def add(a, b):
    L = max(len(a), len(b))
    a = a + [0] * (L - len(a))
    b = b + [0] * (L - len(b))
    return [(x + y) % P for x, y in zip(a, b)]


def sub(a, b):
    L = max(len(a), len(b))
    a = a + [0] * (L - len(a))
    b = b + [0] * (L - len(b))
    return [(x - y) % P for x, y in zip(a, b)]


def scal(c, a):
    return [c * x % P for x in a]


def sinv(a, L=None):
    L = L or Ser.L
    i0 = inv(a[0])
    r = [0] * L
    r[0] = i0
    for n in range(1, L):
        s = 0
        for k in range(1, n + 1):
            if k < len(a):
                s += a[k] * r[n - k]
        r[n] = (-s * i0) % P
    return r


def const(c, L=None):
    L = L or Ser.L
    r = [0] * L
    r[0] = fp(c) if isinstance(c, Fr) else c % P
    return r


def deriv(a):
    return [(k + 1) * a[k + 1] % P for k in range(len(a) - 1)] + [0]


def spow(a, n, L=None):
    L = L or Ser.L
    if n < 0:
        return spow(sinv(a, L), -n, L)
    r = const(1, L)
    for _ in range(n):
        r = mul(r, a, L)
    return r


def ssqrt(a, root0, L=None):
    L = L or Ser.L
    r = [0] * L
    r[0] = root0
    i2 = inv(2 * root0)
    for n in range(1, L):
        s = a[n] if n < len(a) else 0
        for k in range(1, n):
            s -= r[k] * r[n - k]
        r[n] = s * i2 % P
    return r


def sexp(a, L=None):
    L = L or Ser.L
    assert a[0] == 0
    r = [0] * L
    r[0] = 1
    for n in range(1, L):
        s = 0
        for k in range(1, n + 1):
            if k < len(a):
                s += k * a[k] * r[n - k]
        r[n] = s * inv(n) % P
    return r


def sint(a, L=None):
    L = L or Ser.L
    r = [0] * L
    for k in range(L - 1):
        r[k + 1] = a[k] * inv(k + 1) % P
    return r


def tau_series(L=None):
    L = L or Ser.L
    r = [0] * L
    f = 1
    for k in range(1, L):
        f = f * k % P
        r[k] = (-(-1) ** k * inv(f)) % P
    return r


# Laurent: (val, ser)
def lmul(A, B, L=None):
    return (A[0] + B[0], mul(A[1], B[1], L))


def ladd(A, B, L=None):
    L = L or Ser.L
    v = min(A[0], B[0])
    a = [0] * (A[0] - v) + A[1]
    b = [0] * (B[0] - v) + B[1]
    a = (a + [0] * L)[:L]
    b = (b + [0] * L)[:L]
    return (v, [(x + y) % P for x, y in zip(a, b)])


def lscal(c, A):
    return (A[0], scal(c, A[1]))


def lneg(A):
    return lscal(P - 1, A)


def lcoef(A, k):
    i = k - A[0]
    if i < 0 or i >= len(A[1]):
        return 0
    return A[1][i]


def lder(A):
    v, s = A
    out = [(v + i) * s[i] % P for i in range(len(s))]
    return (v - 1, out)


class W:
    pass


def weights(a0, a1, a2, e1, e2):
    w = W()
    w.raw = tuple(Fr(x) for x in (a0, a1, a2, e1, e2))
    w.a0, w.a1, w.a2, w.e1, w.e2 = [fp(x) for x in (a0, a1, a2, e1, e2)]
    w.A = (w.a0 + w.a2) % P  # alpha0+alpha2
    w.B = (w.a0 + w.a1) % P  # alpha0+alpha1
    w.Ea = (w.a0 + w.a1 + w.a2) % P  # sum of the X weights
    w.Ee = (w.e1 + w.e2) % P
    w.cy = (w.Ea == w.Ee)
    w.ei = (w.e1 * w.e2 % P * w.a0 % P * w.a1 % P * w.a2 % P,
            w.e1 * w.e2 % P * w.a1 % P * w.a2 % P,
            (-w.e1 * w.e2 % P * w.Ea) % P)
    return w


def pol_table(K):
    pols = [{1: 1}]
    for n in range(1, K):
        d = {}
        for i, c in pols[-1].items():
            d[i + 1] = d.get(i + 1, 0) + c * i
            d[i] = d.get(i, 0) - c * i
        pols.append({i: c for i, c in d.items() if c})
    return pols


def star_series(w, tau, L=None):
    """s*(tau): (a0+a1) tau s^2 - (Ea + a0 tau) s + A = 0, s*(0)=A/Ea, as series in t"""
    L = L or Ser.L
    den = add(const(w.Ea, L), scal(w.a0, tau))
    idn = sinv(den, L)
    s = [0] * L
    s[0] = w.A * inv(w.Ea) % P
    for _ in range(L + 1):
        s2 = mul(s, s, L)
        s = mul(idn, add(const(w.A, L), scal(w.B, mul(tau, s2, L))), L)
    return s


def wick(wpow, wk, N, L):
    """wk: k -> series (k >= 3); wpow[n] = series (-1/w2)^n. Returns r[m], m = 0..N (C7)."""
    fact = [1] * (6 * N + 8)
    for i in range(1, len(fact)):
        fact[i] = fact[i - 1] * i % P
    ifact = [inv(f) for f in fact]
    c = {(0, 0): const(1, L)}
    D = 6 * N
    for d in range(3, D + 1):
        idinv = inv(d)
        jmax = d // 3
        jmin = max(1, (d + 1) // 2 - N)
        for j in range(jmin, jmax + 1):
            s = [0] * L
            for k in range(3, d + 1):
                if (d - k, j - 1) in c:
                    t = mul(wk[k], c[(d - k, j - 1)], L)
                    co = k * ifact[k] % P
                    for i in range(L):
                        s[i] += co * t[i]
            s = [x * idinv % P for x in s]
            if any(s):
                c[(d, j)] = s
    res = []
    dfs = [1]
    for n in range(1, 3 * N + 2):
        dfs.append(dfs[-1] * (2 * n - 1) % P)
    for m in range(N + 1):
        s = [0] * L
        for n in range(m, 3 * m + 1):
            j = n - m
            if (2 * n, j) in c:
                t = mul(wpow[n], c[(2 * n, j)], L)
                for i in range(L):
                    s[i] += dfs[n] * t[i]
        res.append([x % P for x in s])
    return res


def h_exponent(vs, N):
    """exponent of h_v(z): -Sum_k bk z^(2k-1) Sum_a v_a^(1-2k); index = z power"""
    ex = [0] * (N + 1)
    for k in range(1, N + 2):
        if 2 * k - 1 > N:
            break
        s = 0
        for v in vs:
            s += inv(pow(v, 2 * k - 1, P))
        ex[2 * k - 1] = (-bk(k) * s) % P
    return ex


def zexp(ex, N):
    r = [0] * (N + 1)
    r[0] = 1
    for n in range(1, N + 1):
        s = 0
        for k in range(1, n + 1):
            if ex[k]:
                s += k * ex[k] * r[n - k]
        r[n] = s * inv(n) % P
    return r


def h_series(vs, N):
    return zexp(h_exponent(vs, N), N)


class Rows:
    pass


def build_rows(w, N, L):
    """point-frame Wick rows (no constants) at both colours as t-series; L_i, Delta_i Laurent."""
    Ser.L = L
    R = Rows()
    R.w = w
    R.N = N
    R.L = L
    tau = tau_series(L)
    R.tau = tau
    itau_reg = sinv(tau[1:] + [0], L)  # t/tau
    R.itau_reg = itau_reg
    onemt = sub(const(1, L), tau)  # 1-tau = 1/Q
    s = star_series(w, tau, L)
    R.s = s
    U1 = sub(const(1, L), mul(tau, s, L))
    R.U1 = U1
    U2 = scal((-w.a2 * inv(w.B)) % P, mul(onemt, sinv(U1, L), L))
    R.U2 = U2
    KK = 2 * N + 3
    pols = pol_table(KK + 1)
    tX = sinv(s, L)
    tY = mul(onemt, sinv(sub(s, const(1, L)), L), L)
    powX = [const(1, L)]
    powY = [const(1, L)]
    for i in range(KK + 1):
        powX.append(mul(powX[-1], tX, L))
        powY.append(mul(powY[-1], tY, L))
    tpows = [const(1, L)]
    for i in range(KK + 1):
        tpows.append(mul(tpows[-1], tau, L))

    def w_rescaled(k):
        out = [0] * L
        for i, co in pols[k - 1].items():
            co = co % P
            t = mul(tpows[k - i], add(scal(w.a1, powY[i]), scal(w.A, powX[i])), L)
            for n in range(L):
                out[n] -= co * t[n]
        return [x % P for x in out]
    wk1 = {k: w_rescaled(k) for k in range(2, KK + 1)}
    R.wk1 = wk1
    i2 = sinv(scal(P - 1, wk1[2]), L)
    p = [const(1, L)]
    for n in range(1, 3 * N + 3):
        p.append(mul(p[-1], i2, L))
    R.r0pt1 = wick(p, wk1, N, L)
    X2 = sinv(sub(const(1, L), U2), L)
    Y2 = sinv(sub(const(1, L), mul(sinv(onemt, L), U2, L)), L)
    powX = [const(1, L)]
    powY = [const(1, L)]
    for i in range(KK + 1):
        powX.append(mul(powX[-1], X2, L))
        powY.append(mul(powY[-1], Y2, L))

    def w2(k):
        out = [0] * L
        for i, co in pols[k - 1].items():
            co = co % P
            t = add(scal(w.a1, powY[i]), scal(w.A, powX[i]))
            for n in range(L):
                out[n] -= co * t[n]
        return [x % P for x in out]
    wk2 = {k: w2(k) for k in range(2, KK + 1)}
    R.wk2 = wk2
    i2 = sinv(scal(P - 1, wk2[2]), L)
    p = [const(1, L)]
    for n in range(1, 3 * N + 3):
        p.append(mul(p[-1], i2, L))
    R.r0pt2 = wick(p, wk2, N, L)
    num1 = add(scal(w.B, U1), const(w.a2, L))
    R.L1 = (-1, mul(num1, itau_reg, L))
    num2 = add(scal(w.B, U2), const(w.a2, L))
    assert num2[0] == 0
    R.L2 = (0, mul(num2[1:] + [0], itau_reg, L))
    kap = w.a1 * w.A % P * w.e1 % P * w.e2 % P
    R.D1 = lmul((-2, spow(itau_reg, 2, L)), (0, scal(kap, wk1[2])), L)
    R.D2 = (0, scal(kap, wk2[2]))
    return R


def pair(w, a, b, L):
    """<a,b> for Laurent covectors, eta^-1 of Prop 13.2 at generic weights"""
    e00, e01, e11 = w.ei
    t = lscal(e00, lmul(a[0], b[0], L))
    t = ladd(t, lscal(e01, ladd(lmul(a[0], b[1], L), lmul(a[1], b[0], L), L)), L)
    t = ladd(t, lscal(e11, lmul(a[1], b[1], L)), L)
    return t


def framed_rows(R, cls):
    """rows v_i[m] = (r0, r1) in the frame with vertex class multiset cls:
       r0 = h_nw h_cls^-1 Wick. cls = None: h_nw only (point frame, class 1)."""
    w = R.w
    N = R.N
    L = R.L
    nw = [w.a1, w.A, (-w.e1) % P, (-w.e2) % P]
    ex = h_exponent(nw, N)
    if cls is not None:
        ex2 = h_exponent(cls, N)
        ex = [(x - y) % P for x, y in zip(ex, ex2)]
    fac = zexp(ex, N)
    out = {}
    for col, base, Lc, Dc in ((1, R.r0pt1, R.L1, R.D1), (2, R.r0pt2, R.L2, R.D2)):
        r0 = []
        for m in range(N + 1):
            s = [0] * L
            for i in range(m + 1):
                if fac[i]:
                    for n in range(L):
                        s[n] += fac[i] * base[m - i][n]
            r0.append((0, [x % P for x in s]))
        vD, sD = Dc
        dl = lmul(lder(Dc), (-vD, sinv(sD, L)), L)
        dl = lscal(inv(2), dl)
        r1 = []
        for m in range(N + 1):
            t = lmul(Lc, r0[m], L)
            if m >= 1:
                d = lder(r0[m - 1])
                d = ladd(d, lneg(lmul(dl, r0[m - 1], L)), L)
                t = ladd(t, d, L)
            r1.append(t)
        out[col] = (r0, r1)
    return out, fac


def edge_table(R, rows, i, j, M):
    """E_ij(m,n) for m+n <= M, Laurent (C11)"""
    w = R.w
    L = R.L
    v = rows[i]
    u = rows[j]
    Pq = {}
    for p in range(M + 2):
        for q in range(M + 2 - p):
            s = pair(w, (v[0][p], v[1][p]), (u[0][q], u[1][q]), L)
            Pq[(p, q)] = lscal((-(-1) ** (p + q)) % P, s)
    E = {}
    for m in range(M + 1):
        for n in range(M + 1 - m):
            t = (0, [0] * L)
            for k in range(n + 1):
                t = ladd(t, lscal((-1) ** k % P, Pq[(m + 1 + k, n - k)]), L)
            E[(m, n)] = t
    return E
