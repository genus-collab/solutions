# Refutation test for A2.8: genus-1 divisor identity through R_2.
# theta Int_{Mbar_{1,1}} Omega(x) = Int_{Mbar_{1,2}} Omega(x,H),
# graph sum from PPZ Defs 2.2, 2.13, written here from scratch.
# Everything in the idempotent frame e_i, exact truncated series over Q.
# Finite order: can refute only.
from fractions import Fraction as Fr
from functools import reduce
from math import factorial as fa
N = 6


def sadd(a, b): return [x + y for x, y in zip(a, b)]
def sneg(a): return [-x for x in a]
def ssub(a, b): return [x - y for x, y in zip(a, b)]


def smul(a, b):
    r = [Fr(0)] * N
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j >= N:
                break
            r[i + j] += x * y
    return r


def sinv(a):
    r = [Fr(0)] * N
    r[0] = 1 / a[0]
    for n in range(1, N):
        s = sum((a[k] * r[n - k] for k in range(1, n + 1)), Fr(0))
        r[n] = -s / a[0]
    return r


def const(c):
    r = [Fr(0)] * N
    r[0] = Fr(c)
    return r


ZERO = const(0)
ONE = const(1)


def scal(c, a): return [Fr(c) * x for x in a]
def thy(a): return [Fr(k) * a[k] for k in range(N)]


def integ(a):
    assert a[0] == 0, a[0]
    r = [Fr(0)] * N
    for k in range(1, N):
        r[k] = a[k] / k
    return r


def mmul(A, B):
    n = len(A)
    m = len(B[0])
    l = len(B)
    R = [[ZERO for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = ZERO
            for k in range(l):
                s = sadd(s, smul(A[i][k], B[k][j]))
            R[i][j] = s
    return R


def madd(A, B): return [[sadd(a, b) for a, b in zip(r, s)] for r, s in zip(A, B)]
def msub(A, B): return [[ssub(a, b) for a, b in zip(r, s)] for r, s in zip(A, B)]
def mscal(c, A): return [[scal(c, a) for a in r] for r in A]
def mid(n): return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def minv(A):
    n = len(A)
    M = [list(A[i]) + [ONE if i == j else ZERO for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c][0] != 0)
        M[c], M[p] = M[p], M[c]
        iv = sinv(M[c][c])
        M[c] = [smul(iv, x) for x in M[c]]
        for i in range(n):
            if i != c:
                f = M[i][c]
                M[i] = [ssub(x, smul(f, y)) for x, y in zip(M[i], M[c])]
    return [r[n:] for r in M]


def transpose(A): return [list(r) for r in zip(*A)]


class Theory:
    pass


def p1theory(a0, a1):
    # (H+a0)(H+a1)=Q, H|P_i=-a_i, eta(1,1)=0, eta(1,H)=1
    T = Theory()
    T.r = 2
    a0, a1 = Fr(a0), Fr(a1)
    Q = [Fr(0)] * N
    Q[1] = Fr(1)

    def root(l0):
        l = const(l0)
        for _ in range(N + 1):
            lp = sadd(l, const(a0))
            lq = sadd(l, const(a1))
            f = ssub(smul(lp, lq), Q)
            df = sadd(lp, lq)
            l = ssub(l, smul(f, sinv(df)))
        return l
    lam = [root(-a0), root(-a1)]
    T.lam = lam
    T.Delta = [ssub(lam[0], lam[1]), ssub(lam[1], lam[0])]
    T.X = [[ONE, ONE], [lam[0], lam[1]]]
    T.theta = thy
    T.wts = [[a1 - a0], [a0 - a1]]
    return T


def lp2theory(al, e1, e2):
    T = Theory()
    T.r = 3
    al = [Fr(a) for a in al]
    e1, e2 = Fr(e1), Fr(e2)
    s1 = sum(al)
    s2 = al[0] * al[1] + al[0] * al[2] + al[1] * al[2]
    s3 = al[0] * al[1] * al[2]
    A = [Fr(fa(3 * d), fa(d) ** 3) * (-1) ** d for d in range(N)]
    Ainv = sinv(A)
    y = [Fr(0)] * N
    y[1] = Fr(1)

    def Lhat(Y):
        P = ONE
        for a in al:
            P = smul(P, sadd(Y, const(a)))
        return sadd(P, scal(27, smul(y, smul(Y, smul(Y, Y)))))

    def dLhat(Y):
        Y2 = smul(Y, Y)
        d = sadd(sadd(scal(3, Y2), scal(2 * s1, Y)), const(s2))
        return sadd(d, scal(81, smul(y, Y2)))
    Ys = []
    for a in al:
        Y = const(-a)
        for _ in range(N + 1):
            Y = ssub(Y, smul(Lhat(Y), sinv(dLhat(Y))))
        Ys.append(Y)

    def W(Y): return sadd(sadd(scal(s1, smul(Y, Y)), scal(2 * s2, Y)), const(3 * s3))
    T.Delta = [scal(3 * e1 * e2, W(Y)) for Y in Ys]
    lam = [smul(Y, Ainv) for Y in Ys]
    T.lam = lam
    W0 = []
    for i in range(3):
        j, k = [m for m in range(3) if m != i]
        W0.append(3 * al[i] * (al[j] - al[i]) * (al[k] - al[i]) * e1 * e2)
    eta = [[sum((-al[i]) ** (a + b) / W0[i] for i in range(3)) for b in range(3)] for a in range(3)]
    Dinv = [sinv(d) for d in T.Delta]

    def pw(s, k):
        r = ONE
        for _ in range(k):
            r = smul(r, s)
        return r
    S = [None] * 5
    for k in range(5):
        t = ZERO
        for i in range(3):
            t = sadd(t, smul(pw(lam[i], k), Dinv[i]))
        S[k] = t
    for k, (a, b) in zip(range(3), [(0, 0), (0, 1), (1, 1)]):
        assert S[k] == const(eta[a][b]), (k, S[k], eta[a][b])
    # S1=S2=0 identically; C2 = S3/m3, C1 = (S4 - m4 C2^2)/(2 S3)
    m3, m4 = eta[1][2], eta[2][2]
    assert S[1] == ZERO and S[2] == ZERO
    C2 = scal(1 / m3, S[3])
    C1 = smul(ssub(S[4], scal(m4, smul(C2, C2))), sinv(scal(2, S[3])))
    # check A^3 C2 = X = 1/(1+27y), that is A^2 D = X with C2 = D/A
    xi = [Fr(0)] * N
    xi[0] = Fr(1)
    xi[1] = Fr(27)
    assert smul(smul(A, smul(A, A)), C2) == sinv(xi), "A^2 D = X"
    C2inv = sinv(C2)
    T.X = [[ONE] * 3, [lam[0], lam[1], lam[2]],
           [smul(ssub(smul(l, l), smul(C1, l)), C2inv) for l in lam]]
    T.theta = lambda a: smul(Ainv, thy(a))
    T.ginv = A
    T.wts = []
    for i in range(3):
        j, k = [m for m in range(3) if m != i]
        T.wts.append([3 * al[i], al[j] - al[i], al[k] - al[i], -e1, -e2])
    return T


def rmatrix(T, nmax=2):
    """solve z(V rho + theta rho + rho G) = [Lambda,rho] in the e-frame, rho(Q=0)=1"""
    r = T.r
    th = T.theta
    X = T.X
    Xinv = minv(X)
    thX = [[th(a) for a in row] for row in X]
    V = transpose(mscal(-1, mmul(Xinv, thX)))
    G = [scal(Fr(1, 2), smul(th(T.Delta[i]), sinv(T.Delta[i]))) for i in range(r)]
    for i in range(r):
        assert sadd(V[i][i], G[i]) == ZERO, "V_ii+G_i"
    T.V = V
    T.G = G
    rho = [mid(r)]
    dinv = [[None if i == j else sinv(ssub(T.lam[i], T.lam[j])) for j in range(r)] for i in range(r)]
    ginv = getattr(T, 'ginv', ONE)  # theta_y = ginv * theta
    for n in range(nmax):
        cur = rho[n]
        Vr = mmul(V, cur)
        thr = [[th(a) for a in row] for row in cur]
        rG = [[smul(cur[i][j], G[j]) for j in range(r)] for i in range(r)]
        Lhs = madd(madd(Vr, thr), rG)
        nxt = [[ZERO] * r for _ in range(r)]
        for i in range(r):
            for j in range(r):
                if i != j:
                    nxt[i][j] = smul(Lhs[i][j], dinv[i][j])
        for i in range(r):
            t = ZERO
            for k in range(r):
                if k != i:
                    t = ssub(t, smul(V[i][k], nxt[k][i]))
            nxt[i][i] = integ(smul(ginv, t))
        rho.append(nxt)
    return rho


def hser(wts):
    # h_v(z)=exp(-Sum B_2k/(2k(2k-1)) p_{2k-1}(1/v) z^{2k-1}), through z^2
    a1 = -Fr(1, 6) / 2 * sum(1 / w for w in wts)
    return [Fr(1), a1, a1 * a1 / 2]


def apply_const(rho, T, mode='h', extra=None):
    # R = R~ diag(const series)
    r = T.r
    if mode == 'h':
        hs = [hser(w) for w in T.wts]
    elif mode == 'none':
        hs = [[Fr(1), Fr(0), Fr(0)] for w in T.wts]
    elif mode == 'arb':
        hs = [[Fr(1), Fr(i + 3, 7), Fr(i + 3, 7) ** 2 / 2] for i in range(r)]
    out = []
    for n in range(3):
        M = [[ZERO] * r for _ in range(r)]
        for i in range(r):
            for j in range(r):
                t = ZERO
                for k in range(n + 1):
                    t = sadd(t, scal(hs[j][n - k], rho[k][i][j]))
                M[i][j] = t
        out.append(M)
    if mode == 'Qdep':
        pass
    return out


def Qdep_const(R, T):
    # multiply on the right by diag exp(a_i(Q) z), a_i = Q (a Q-dependent 'constant'): must fail
    r = T.r
    a = [[Fr(0)] * N for _ in range(r)]
    for i in range(r):
        a[i][1] = Fr(i + 1)
    hs = [[ONE, a[i], scal(Fr(1, 2), smul(a[i], a[i]))] for i in range(r)]
    out = []
    for n in range(3):
        M = [[ZERO] * r for _ in range(r)]
        for i in range(r):
            for j in range(r):
                t = ZERO
                for k in range(n + 1):
                    t = sadd(t, smul(hs[j][n - k], R[k][i][j]))
                M[i][j] = t
        out.append(M)
    return out


def inv_series(rho):
    b1 = mscal(-1, rho[1])
    b2 = msub(mmul(rho[1], rho[1]), rho[2])
    return [rho[0], b1, b2]


def variant(rho, kind):
    if kind == 'R':
        return rho
    if kind == 'R(-z)':
        return [rho[0], mscal(-1, rho[1]), rho[2]]
    if kind == 'R*':   # R^*(z)=R^{-1}(-z)
        b = inv_series(rho)
        return [b[0], mscal(-1, b[1]), b[2]]
    if kind == 'R^-1(z)':
        return inv_series(rho)
    raise ValueError


def integrals(T, R):
    r = T.r
    D = T.Delta
    Di = [sinv(d) for d in D]
    b = inv_series(R)
    b1, b2 = b[1], b[2]
    I_ = [ONE] * r

    def ap(M, x):
        return [reduce(sadd, [smul(M[i][j], x[j]) for j in range(r)], ZERO) for i in range(r)]
    T2 = [sneg(v) for v in ap(b1, I_)]
    T3 = [sneg(v) for v in ap(b2, I_)]

    def om(g, *xs):
        t = ZERO
        for i in range(r):
            p = ONE
            for x in xs:
                p = smul(p, x[i])
            if g == 0:
                p = smul(p, Di[i])
            t = sadd(t, p)
        return t

    def e(s): return [ONE if i == s else ZERO for i in range(r)]
    E0 = [[sneg(smul(b1[s][t], D[t])) for t in range(r)] for s in range(r)]
    E1 = [[sneg(smul(b2[s][t], D[t])) for t in range(r)] for s in range(r)]
    sym = all(E0[s][t] == E0[t][s] for s in range(r) for t in range(r))

    def I1(x):
        t = scal(Fr(1, 24), sadd(om(1, ap(b1, x)), om(1, x, T2)))
        for s in range(r):
            for u in range(r):
                t = sadd(t, scal(Fr(1, 2), smul(om(0, x, e(s), e(u)), E0[s][u])))
        return t

    def J2(x, y):
        tot = ZERO
        bx, by = ap(b1, x), ap(b1, y)
        tot = sadd(tot, scal(Fr(1, 24), sadd(sadd(om(1, ap(b2, x), y), om(1, x, ap(b2, y))), om(1, bx, by))))
        tot = sadd(tot, scal(Fr(1, 12), sadd(om(1, bx, y, T2), om(1, x, by, T2))))
        tot = sadd(tot, scal(Fr(1, 24), om(1, x, y, T3)))
        tot = sadd(tot, scal(Fr(1, 12), om(1, x, y, T2, T2)))
        for s in range(r):
            for t_ in range(r):
                # (b) genus-1 vertex - (0,3){x,y}
                w = sadd(smul(om(1, e(s), T2), E0[s][t_]), smul(om(1, e(s)), E1[s][t_]))
                tot = sadd(tot, scal(Fr(1, 24), smul(om(0, x, y, e(t_)), w)))
                # (c) loop on (0,4)
                w = smul(E0[s][t_], sadd(sadd(om(0, bx, y, e(s), e(t_)), om(0, x, by, e(s), e(t_))), om(0, x, y, e(s), e(t_), T2)))
                w = sadd(w, smul(om(0, x, y, e(s), e(t_)), sadd(E1[s][t_], E1[t_][s])))
                tot = sadd(tot, scal(Fr(1, 2), w))
                # (e) (0,3){x,y} - (0,3) with loop
                for u in range(r):
                    for v in range(r):
                        tot = sadd(tot, scal(Fr(1, 2), smul(smul(om(0, x, y, e(s)), E0[s][t_]), smul(om(0, e(t_), e(u), e(v)), E0[u][v]))))
        # (d) banana
        for s1 in range(r):
            for s2 in range(r):
                ox = om(0, x, e(s1), e(s2))
                for t1 in range(r):
                    for t2 in range(r):
                        tot = sadd(tot, scal(Fr(1, 2), smul(smul(ox, om(0, y, e(t1), e(t2))), smul(E0[s1][t1], E0[s2][t2]))))
        return tot
    return I1, J2, sym


def first_nz(s, M):
    return next((k for k, v in enumerate(s[:M]) if v != 0), None)


def run(name, T, M):
    rho = rmatrix(T)
    print('==', name)
    for mode in ('h', 'none', 'arb'):
        R = apply_const(rho, T, mode)
        for kind in ('R', 'R(-z)', 'R*', 'R^-1(z)'):
            Rv = variant(R, kind)
            I1, J2, sym = integrals(T, Rv)
            xs = [T.X[a] for a in range(T.r)]
            H = T.X[1]
            one = T.X[0]
            res = [first_nz(ssub(T.theta(I1(x)), J2(x, H)), M) for x in xs]
            uni = [first_nz(J2(x, one), M) for x in xs]
            print('  const=%-5s %-8s R1 sym=%s; first order of failure theta I(x_a)-J(x_a,H): %s; unit J(x_a,1): %s'
                  % (mode, kind, sym, res, uni))
    R = apply_const(rho, T, 'h')
    Rq = Qdep_const(R, T)
    I1, J2, sym = integrals(T, Rq)
    xs = [T.X[a] for a in range(T.r)]
    print('  Q-dependent exp(a(Q) z):', [first_nz(ssub(T.theta(I1(x)), J2(x, T.X[1])), M) for x in xs],
          'unit', [first_nz(J2(x, T.X[0]), M) for x in xs])
    I1, J2, _ = integrals(T, R)
    print('  I(H)  =', [str(c) for c in I1(T.X[1])[:M]])
    print('  J(H,H)=', [str(c) for c in J2(T.X[1], T.X[1])[:M]])
    return rho


if __name__ == '__main__':
    M = N - 2
    run('equivariant P1 (2,-5)', p1theory(2, -5), M)
    run('equivariant P1 (1/3,7)', p1theory(Fr(1, 3), 7), M)
    run('local P2 x C2 alpha=(1,2,5) eps=(3,7)', lp2theory((1, 2, 5), 3, 7), M)
