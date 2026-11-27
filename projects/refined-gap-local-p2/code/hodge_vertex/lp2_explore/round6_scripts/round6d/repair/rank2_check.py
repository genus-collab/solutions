# Evidence only, exact Fractions, truncated in Q and z. Rank-2 data of Proposition 13.2(a).
from fractions import Fraction as Fr
from math import factorial
from collections import Counter

N = 7
ZM = 3


def run(a0, a1, a2, e1, e2):
    a0, a1, a2, e1, e2 = map(Fr, (a0, a1, a2, e1, e2))

    def S(c=0):
        return [Fr(c)] + [Fr(0)] * N

    def add(x, y):
        return [p + q for p, q in zip(x, y)]

    def sub(x, y):
        return [p - q for p, q in zip(x, y)]

    def mul(x, y):
        r = [Fr(0)] * (N + 1)
        for i, p in enumerate(x):
            if p:
                for j in range(N + 1 - i):
                    r[i + j] += p * y[j]
        return r

    def scal(c, x):
        return [c * p for p in x]

    def inv(x):
        r = [Fr(0)] * (N + 1)
        r[0] = 1 / x[0]
        for n in range(1, N + 1):
            r[n] = -sum(x[k] * r[n - k] for k in range(1, n + 1)) / x[0]
        return r

    def th(x):
        return [k * p for k, p in enumerate(x)]

    def thinv(x, c0):
        assert x[0] == 0
        return [c0] + [x[k] / k for k in range(1, N + 1)]

    geo = [Fr(1)] * (N + 1)
    s1 = scal(Fr(-1), mul(add(S(a0), [Fr(0), a1 + a2] + [Fr(0)] * (N - 1)), geo))
    s2 = scal(-a1 * a2, mul([Fr(0), Fr(1)] + [Fr(0)] * (N - 1), geo))

    def root(l0):
        L = S(l0)
        for _ in range(N + 2):
            f = add(sub(mul(L, L), mul(s1, L)), s2)
            d = sub(scal(Fr(2), L), s1)
            L = sub(L, mul(f, inv(d)))
        return L

    Lp = root(0)
    Lm = root(-a0)
    assert add(Lp, Lm) == s1 and mul(Lp, Lm) == s2
    d = sub(Lp, Lm)
    assert d[0] == a0  # L_+ - L_- is a unit
    sig = a0 + a1 + a2
    ee = e1 * e2
    den = (a0 + a1) * (a0 + a2) * ee
    eta = [[sig / (a1 * a2 * den), 1 / den], [1 / den, -a0 / den]]
    Ep = a0 * a1 * a2 * ee
    Em = -a0 * (a0 + a1) * (a0 + a2) * ee
    assert eta[0][0] == 1 / Ep + 1 / Em and eta[0][1] == -a0 / Em and eta[1][1] == a0 * a0 / Em

    def Delta(Li, Lj):
        x = sub(S(eta[0][1]), scal(eta[0][0], Lj))  # eta(H-L_j,1)
        x = mul(x, inv(sub(Li, Lj)))  # eta(e_i,1)
        return inv(x)

    Dp = Delta(Lp, Lm)
    Dm = Delta(Lm, Lp)

    def form132b(L):
        return scal(ee, sub(add(S(a0 * a1 * a2), scal(2 * a1 * a2, L)), scal(sig, mul(L, L))))

    def form132a(Li, Lj):
        return mul(scal(a1 * a2 * (a0 + a1) * (a0 + a2) * ee, sub(Li, Lj)), inv(sub(S(a1 * a2), scal(sig, Lj))))

    def W2(L):
        A = scal(1 / (a0 + a2), sub(L, S(a2)))
        B = scal(-1 / a1, L)
        return sub(scal(-a1, mul(B, add(S(1), B))), scal(a0 + a2, mul(A, add(S(1), A))))

    res = {}
    for Li, Lj, D_ in ((Lp, Lm, Dp), (Lm, Lp, Dm)):
        assert D_ == form132b(Li) and D_ == form132a(Li, Lj)
        assert D_ == scal(a1 * (a0 + a2) * ee, W2(Li))
    assert Dp[0] == Ep and Dm[0] == Em
    ratio = Dp[0] / Dm[0]
    assert ratio == -a1 * a2 / ((a0 + a1) * (a0 + a2))
    res['(13.2) three forms, Delta_i(0)=e(wt_i), ratio'] = True

    def bern(m):
        B = [Fr(0)] * (m + 1)
        B[0] = Fr(1)
        for n in range(1, m + 1):
            B[n] = -sum(Fr(factorial(n + 1), factorial(k) * factorial(n + 1 - k)) * B[k] for k in range(n)) / (n + 1)
        return B

    B = bern(2 * ZM + 2)

    def h(v):
        ex = [Fr(0)] * (ZM + 1)
        for k in range(1, ZM + 2):
            if 2 * k - 1 <= ZM:
                ex[2 * k - 1] = -B[2 * k] / (2 * k * (2 * k - 1)) * sum(x ** (1 - 2 * k) for x in v)
        r = [Fr(0)] * (ZM + 1)
        r[0] = Fr(1)
        for n in range(1, ZM + 1):
            r[n] = sum(k * ex[k] * r[n - k] for k in range(1, n + 1)) / n
        return r

    wtp = (a0, a1, a2, -e1, -e2)
    wtm = (-a0, a0 + a1, a0 + a2, -e1, -e2)
    nw = (a1, a0 + a2, -e1, -e2)

    # rows from the QDE: theta c_n = -(D^2 c_{n-1})/(L_i-L_j), c_n(0) = [z^n] h_{wt_i}
    def rows(Li, Lj, D_, wt):
        g = scal(Fr(1, 2), mul(th(D_), inv(D_)))
        assert mul(g, sub(Li, Lj)) == th(Li)  # order z^0
        Dop = lambda x: sub(th(x), mul(g, x))
        hh = h(wt)
        c = [S(1)]
        iv = inv(sub(Li, Lj))
        for n in range(1, ZM + 1):
            rhs = scal(Fr(-1), mul(Dop(Dop(c[n - 1])), iv))
            c.append(thinv(rhs, hh[n]))
        r1 = []
        for n in range(ZM + 1):
            x = mul(Li, c[n])
            if n > 0:
                x = add(x, Dop(c[n - 1]))
            r1.append(x)
        return c, r1, g

    cp, r1p, gp = rows(Lp, Lm, Dp, wtp)
    cm, r1m, gm = rows(Lm, Lp, Dm, wtm)

    # Laplace rows: r_0 = h_nw * Lnorm[W], (C7), W = a1 log(1-QU) - a2 log U + (a0+a2) log(1-U)
    # W_k as polynomials in A = U/(1-U), B = QU/(1-QU); theta A = A(1+A), theta B = B(1+B)
    def dpoly(p):
        r = {}
        for (i, j), c_ in p.items():
            if i:
                r[(i, j)] = r.get((i, j), 0) + c_ * i
                r[(i + 1, j)] = r.get((i + 1, j), 0) + c_ * i
            if j:
                r[(i, j)] = r.get((i, j), 0) + c_ * j
                r[(i, j + 1)] = r.get((i, j + 1), 0) + c_ * j
        return r

    W1 = {(0, 0): -a2, (1, 0): -(a0 + a2), (0, 1): -a1}
    Wp = {1: W1}
    for k in range(2, 2 * ZM + 3):
        Wp[k] = dpoly(Wp[k - 1])

    def evalp(p, A, Bq):
        r = S(0)
        for (i, j), c_ in p.items():
            t = S(1)
            for _ in range(i):
                t = mul(t, A)
            for _ in range(j):
                t = mul(t, Bq)
            r = add(r, scal(c_, t))
        return r

    def partitions(n, minp=3):
        if n == 0:
            yield []
            return
        for p in range(minp, n + 1):
            for rest in partitions(n - p, p):
                yield [p] + rest

    def lnorm(L):
        A = scal(1 / (a0 + a2), sub(L, S(a2)))
        Bq = scal(-1 / a1, L)
        Wk = {k: evalp(Wp[k], A, Bq) for k in Wp}
        assert Wk[1] == S(0)
        iW2 = inv(Wk[2])
        out = [S(1)]
        for m in range(1, ZM + 1):
            tot = S(0)
            for n in range(m, 3 * m + 1):
                for p in partitions(2 * n):
                    if n - len(p) != m:
                        continue
                    c_ = Fr(factorial(2 * n) * (-1) ** n, factorial(n) * 2 ** n)
                    t = S(c_)
                    for _ in range(n):
                        t = mul(t, iW2)
                    for j, cc in Counter(p).items():
                        t = scal(Fr(1, factorial(j) ** cc * factorial(cc)), t)
                        for _ in range(cc):
                            t = mul(t, Wk[j])
                    tot = add(tot, t)
            out.append(tot)
        return out

    hn = h(nw)

    def laprow(L):
        ln = lnorm(L)
        out = []
        for m in range(ZM + 1):
            t = S(0)
            for k in range(m + 1):
                t = add(t, scal(hn[k], ln[m - k]))
            out.append(t)
        return out

    lp = laprow(Lp)
    lm = laprow(Lm)
    res['QDE row = Laplace row, colour +'] = all(lp[m] == cp[m] for m in range(ZM + 1))
    res['QDE row = Laplace row, colour -'] = all(lm[m] == cm[m] for m in range(ZM + 1))
    res['QDE row at -z = Laplace row, colour + (should be False)'] = all(
        lp[m] == [(-1) ** m * q for q in cp[m]] for m in range(ZM + 1))

    # R from rows: v_i(x) = Delta_i eta(R e_i, x), so R e_i = eta^-1 v_i / Delta_i
    detinv = 1 / (eta[0][0] * eta[1][1] - eta[0][1] ** 2)
    einv = [[eta[1][1] * detinv, -eta[0][1] * detinv], [-eta[0][1] * detinv, eta[0][0] * detinv]]

    def col(c, r1, D_):
        iD = inv(D_)
        out = []
        for n in range(ZM + 1):
            out.append([mul(iD, add(scal(einv[a][0], c[n]), scal(einv[a][1], r1[n]))) for a in (0, 1)])
        return out

    Cp = col(cp, r1p, Dp)
    Cm = col(cm, r1m, Dm)

    # H* in the basis {1,H}: H*1 = H, H*H = s1 H - s2 1
    def Hstar(x):
        return [scal(Fr(-1), mul(s2, x[1])), add(x[0], mul(s1, x[1]))]

    # flatness on columns: z(g C + theta C) = (H* - L_i) C
    ok = True
    for C, L, g in ((Cp, Lp, gp), (Cm, Lm, gm)):
        for n in range(ZM):
            lhs = [add(mul(g, C[n][a]), th(C[n][a])) for a in (0, 1)]
            H1 = Hstar(C[n + 1])
            rhs = [sub(H1[a], mul(L, C[n + 1][a])) for a in (0, 1)]
            ok = ok and lhs == rhs
        hs = Hstar(C[0])
        ok = ok and all(hs[a] == mul(L, C[0][a]) for a in (0, 1))
    res['flatness z theta(R Psi) = [H*,R] Psi from the rows'] = ok

    # R = Sum_i (R e_i) (x) e_i^*, e_i^*(1) = 1, e_i^*(H) = L_i
    R = []
    for n in range(ZM + 1):
        R.append([[add(mul(Cp[n][a], [S(1), Lp][b]), mul(Cm[n][a], [S(1), Lm][b])) for b in (0, 1)] for a in (0, 1)])
    res['R = 1 mod z'] = all(R[0][a][b] == S(1 if a == b else 0) for a in (0, 1) for b in (0, 1))

    # v_i o R(-z) = e_i^*, i.e. v_i = e_i^* o R^-1(-z)
    ok = True
    for c, r1, L in ((cp, r1p, Lp), (cm, r1m, Lm)):
        for n in range(ZM + 1):
            for b in (0, 1):
                t = S(0)
                for k in range(n + 1):
                    t = add(t, scal(Fr((-1) ** (n - k)), add(mul(c[k], R[n - k][0][b]), mul(r1[k], R[n - k][1][b]))))
                exp = ([S(1), L][b] if n == 0 else S(0))
                ok = ok and t == exp
    res['v_i = e_i^* o R^-1(-z)'] = ok

    def pair(u, v):
        return add(add(scal(einv[0][0], mul(u[0], v[0])), scal(einv[0][1], add(mul(u[0], v[1]), mul(u[1], v[0])))),
                   scal(einv[1][1], mul(u[1], v[1])))

    ok = True
    rowsall = ((cp, r1p, Dp), (cm, r1m, Dm))
    for i in (0, 1):
        for j in (0, 1):
            for n in range(ZM + 1):
                t = S(0)
                for k in range(n + 1):
                    t = add(t, scal(Fr((-1) ** (n - k)), pair((rowsall[i][0][k], rowsall[i][1][k]),
                                                              (rowsall[j][0][n - k], rowsall[j][1][n - k]))))
                exp = (rowsall[i][2] if (i == j and n == 0) else S(0))
                ok = ok and t == exp
    res['<v_i(z), v_j(-z)> = delta_ij Delta_i'] = ok
    return ratio, res


for w in ((2, 3, 5, 7, 11), (Fr(1, 3), Fr(-2, 5), 7, Fr(3, 2), -1), (5, 1, Fr(1, 2), 2, 3)):
    ratio, res = run(*w)
    print('weights', [str(x) for x in w], 'ratio at Q=0 =', ratio)
    for k, v in res.items():
        print('  ', k, v)
