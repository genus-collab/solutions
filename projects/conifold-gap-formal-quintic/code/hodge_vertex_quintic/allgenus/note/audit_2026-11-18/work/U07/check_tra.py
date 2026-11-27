# Independent sanity checks for Section 4 (transport) of the note.
# Exact arithmetic (fractions).  Nothing from the project tree is imported or run.
from fractions import Fraction as Fr
from math import factorial
import random

random.seed(12345)
NU = [0, -1, 3, 2, 1]
EPS = [1, -1, -1, 1, 1]
V = 7          # work modulo v^V


# ---------- truncated polynomials in v (lists of Fractions) ----------
def vz():
    return [Fr(0)] * V


def vmono(c, e):
    p = vz()
    if 0 <= e < V:
        p[e] = Fr(c)
    return p


def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vscal(c, a):
    return [c * x for x in a]


def vmul(a, b):
    r = vz()
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i + j < V:
                    r[i + j] += x * y
    return r


def viszero(a):
    return all(x == 0 for x in a)


# ---------- the vanishing solution ----------
def YC_coef(n, b):
    """[z^{-n}] Y^C_b as (coef, v-exponent) from the residue definition (def:gau-van)."""
    m = NU[b]
    # I^C_m = sum_i (-1/(4z))^i/i! * (-v/z)^{m+4i+1}/(m+4i+1)!,  z-power -(m+5i+1)
    tot = Fr(0)
    ve = None
    for i in range(0, n + 2):
        j = m + 4 * i + 1
        if j < 0:
            continue
        if m + 5 * i + 1 == n:
            tot += Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** j, factorial(j))
            ve = j
    return (EPS[b] * tot, ve)


def closed_form(n, b):
    r = n // 5
    if (n - NU[b] - 1) % 5 == 0:
        return (Fr(EPS[b] * (-1) ** n, 4 ** r * factorial(r) * factorial(n - r)), n - r)
    return (Fr(0), None)


for n in range(0, 40):
    for b in range(5):
        a = YC_coef(n, b)
        c = closed_form(n, b)
        assert a[0] == c[0] and (a[0] == 0 or a[1] == c[1]), (n, b, a, c)
print("closed form of Y^C (lem:gau-vanshape) agrees with the residue definition, n<40")


def Pi0(m, b):
    """Pi_0(z^m H^b) as a truncated polynomial in v:  (-1)^{m+1} [z^{-m}] Y^C_b."""
    c, e = YC_coef(m, b)
    if c == 0:
        return vz()
    return vmono((-1) ** (m + 1) * c, e)


# eq:tra-values
for m in range(0, 30):
    for b in range(5):
        r = m // 5
        if (b + m) % 5 == 1:
            exp = vmono(Fr(-EPS[b], 4 ** r * factorial(r) * factorial(m - r)), m - r)
        else:
            exp = vz()
        assert Pi0(m, b) == exp
print("eq:tra-values ok; Pi_0(H)=", Pi0(0, 1)[:2], " Pi_0(z1)=", Pi0(1, 0)[:3], " Pi_0(z^6 1) coefficient of v^5 =", Pi0(6, 0)[5])


def Pi(i, q):
    """q: dict (k,b) -> vpoly.  Pi_i(q) = sum_{k>=i} q_(k,b) Pi_0(z^{k-i}H^b)."""
    r = vz()
    for (k, b), c in q.items():
        if k >= i:
            r = vadd(r, vmul(c, Pi0(k - i, b)))
    return r


# ---------- frames via the recursion of lem:tra-frames ----------
def frame(k, a):
    """returns dict (j,b)->vpoly for frame(k,a), computed with exact polynomials (no truncation issue
    as long as exponents < V)."""
    f = {}
    for i in range(k - 1, -1, -1):
        val = vscal(Fr(-1), Pi0(k - i, a))
        r = 1
        while i + 5 * r < k:
            pir = Fr(1, 4 ** r * factorial(r) * factorial(4 * r))
            val = vadd(val, vscal(-pir, vmul(vmono(1, 4 * r), f[i + 5 * r])))
            r += 1
        f[i] = val
    q = {(k, a): vmono(1, 0)}
    for j, c in f.items():
        if not viszero(c):
            q[(j, 1)] = c
    return q


def show(q):
    out = []
    for (k, b), c in sorted(q.items(), reverse=True):
        terms = ["%s v^%d" % (x, e) for e, x in enumerate(c) if x]
        out.append("z^%d H^%d * (%s)" % (k, b, " + ".join(terms)))
    return " + ".join(out)


for (k, a) in [(1, 0), (2, 4), (3, 3), (4, 2), (6, 0)]:
    print("frame(%d,%d) =" % (k, a), show(frame(k, a)))

# frames lie in K_0 (as far as the truncation v^V allows: exponents (4(k-j)+nu_a+1)/5 < V)
triv = []
for a in [0, 2, 3, 4]:
    for k in range(0, 7):
        fr = frame(k, a)
        for i in range(0, k + 2):
            assert viszero(Pi(i, fr)), (k, a, i)
        if len(fr) == 1:
            triv.append((k, a))
        # exponent pattern
        for (j, b), c in fr.items():
            if b == 1:
                e = Fr(4 * (k - j) + NU[a] + 1, 5)
                assert e.denominator == 1 and e >= 1
                assert all(x == 0 for d, x in enumerate(c) if d != e)
print("frames in K_0 for k<=6; trivial frames:", sorted(triv))

# ---------- matrices / series in z over Q ----------
def mzero():
    return [[Fr(0)] * 5 for _ in range(5)]


def mid():
    m = mzero()
    for i in range(5):
        m[i][i] = Fr(1)
    return m


def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]


def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(5)] for i in range(5)]


def mscal(c, A):
    return [[c * A[i][j] for j in range(5)] for i in range(5)]


def mT(A):
    return [[A[j][i] for j in range(5)] for i in range(5)]


ETA = mzero()
for a in range(5):
    ETA[a][(3 - a) % 5] = Fr(5)
ETAI = mscal(Fr(1, 25), ETA)
assert mmul(ETA, ETAI) == mid()


def smul(A, B, K):
    """product of matrix series (lists indexed by z-power 0..K)."""
    C = [mzero() for _ in range(K + 1)]
    for i in range(K + 1):
        for j in range(K + 1 - i):
            C[i + j] = madd(C[i + j], mmul(A[i], B[j]))
    return C


def sexp(r, K):
    """exp of a matrix series with r[0]=0, mod z^{K+1}."""
    res = [mid()] + [mzero() for _ in range(K)]
    term = [mid()] + [mzero() for _ in range(K)]
    for n in range(1, K + 1):
        term = smul(term, r, K)
        term = [mscal(Fr(1, n), t) for t in term]
        res = [madd(x, y) for x, y in zip(res, term)]
    return res


def sinv(R, K):
    """inverse of a series with R[0] invertible diagonal or identity; here R[0]=1."""
    assert R[0] == mid()
    S = [mid()] + [mzero() for _ in range(K)]
    for n in range(1, K + 1):
        acc = mzero()
        for j in range(1, n + 1):
            acc = madd(acc, mmul(R[j], S[n - j]))
        S[n] = mscal(Fr(-1), acc)
    return S


def rnd():
    return Fr(random.randint(-5, 5), random.randint(1, 4))


def random_sympl_log(K, graded):
    """r(z)=sum_{k>=1} r_k z^k with r(-z)^T eta + eta r(z) = 0 (vector form)."""
    r = [mzero()]
    for k in range(1, K + 1):
        X = mzero()
        if graded:
            for b in range(5):
                X[(b - k) % 5][b] = rnd()   # vector form of a graded lowered matrix: shift -k
        else:
            X = [[rnd() for _ in range(5)] for _ in range(5)]
        # project: r_k := X - (-1)^k eta^{-1} X^T eta   solves (-1)^k r_k^T eta + eta r_k = 0
        Y = mmul(ETAI, mmul(mT(X), ETA))
        r.append(madd(X, mscal(Fr(-((-1) ** k)), Y)))
    return r


def check_sympl(R, K):
    # R(-z)^T eta R(z) = eta mod z^{K+1}
    Rm = [mscal(Fr((-1) ** k), mT(R[k])) for k in range(K + 1)]
    S = smul(Rm, [ETA] + [mzero()] * K, K)
    S = smul(S, R, K)
    assert S[0] == ETA and all(S[k] == mzero() for k in range(1, K + 1))


def edge_coeffs(S, kmax, lmax):
    """E_{k,l} = [psi^k psi'^l] E_R, with S=R^{-1}:  N_{a,b} = delta eta^{-1} - S_a eta^{-1} S_b^T,
    E_{k,l} = sum_{j=0}^{l} (-1)^j N_{k+1+j, l-j}."""
    def Ncoef(a, b):
        M = mscal(Fr(-1), mmul(S[a], mmul(ETAI, mT(S[b]))))
        if a == 0 and b == 0:
            M = madd(M, ETAI)
        return M
    E = {}
    for k in range(kmax + 1):
        for l in range(lmax + 1):
            acc = mzero()
            for j in range(l + 1):
                acc = madd(acc, mscal(Fr((-1) ** j), Ncoef(k + 1 + j, l - j)))
            E[(k, l)] = acc
    return E


# ---------- Test A: Lemma E over Q ----------
for graded in (False, True):
    M, Lmax = 5, 4
    K = M + Lmax + 2
    r = random_sympl_log(K, graded)
    R = sexp(r, K)
    check_sympl(R, K)
    S = sinv(R, K)
    E = edge_coeffs(S, M, Lmax)
    # divisibility: N(psi,-psi)=0, i.e. sum_{a+b=n} (-1)^b N_{a,b} = 0
    for n in range(0, K + 1):
        acc = mzero()
        for a in range(n + 1):
            b = n - a
            Mx = mscal(Fr(-1), mmul(S[a], mmul(ETAI, mT(S[b]))))
            if a == 0 and b == 0:
                Mx = madd(Mx, ETAI)
            acc = madd(acc, mscal(Fr((-1) ** b), Mx))
        assert acc == mzero()
    # symmetry E(psi,psi')^T = E(psi',psi)
    for k in range(min(M, Lmax) + 1):
        for l in range(min(M, Lmax) + 1):
            assert mT(E[(k, l)]) == E[(l, k)]
    # (psi+psi') E = N
    f = {m: [rnd() for _ in range(5)] for m in range(1, M + 1)}   # f = sum_m f_m z^{-m}
    # alpha_k = eta * [z^{-k-1}] f(-z) = eta * (-1)^{k+1} f_{k+1}
    C = []
    for l in range(Lmax + 1):
        vec = [Fr(0)] * 5
        for k in range(M):
            fk = [Fr((-1) ** (k + 1)) * x for x in f[k + 1]]
            al = [sum(ETA[a][c] * fk[c] for c in range(5)) for a in range(5)]
            for b in range(5):
                vec[b] += sum(al[a] * E[(k, l)][a][b] for a in range(5))
        C.append(vec)
    # RHS: -S(psi') (Rf)_{>=0}(psi')
    Rf = []
    for j in range(Lmax + 1):
        vec = [Fr(0)] * 5
        for m in range(1, M + 1):
            for a in range(5):
                vec[a] += sum(R[j + m][a][c] * f[m][c] for c in range(5))
        Rf.append(vec)
    RHS = []
    for l in range(Lmax + 1):
        vec = [Fr(0)] * 5
        for j in range(l + 1):
            for a in range(5):
                vec[a] -= sum(S[l - j][a][c] * Rf[j][c] for c in range(5))
        RHS.append(vec)
    assert C == RHS, (graded, C, RHS)
    print("Lemma E (eq:tra-lemE) verified, graded =", graded, "; nonzero sample:", C[0][:2])

# ---------- Test B: Lemma T (b), (d) and Cor E' for a graded dressing with (Sp), constant in u ----------
imax, Lmax = 3, 3
mmax = 7           # Pi_0(z^m H^b) = O(v^{m - floor(m/5)}), vanishes mod v^V for m >= mmax+1 ?
assert (mmax + 1) - (mmax + 1) // 5 >= V
K = imax + mmax + Lmax + 4
r = random_sympl_log(K, True)          # vector form, graded: X[(b+k)%5][b]
Rvec = sexp(r, K)
check_sympl(Rvec, K)
# P0vec diagonal with P0^T eta P0 = eta :  entries x_a x_{3-a} = 1
x0, x1 = Fr(3, 2), Fr(-2, 5)
P0vec = mzero()
for a, val in enumerate([x0, x1, 1 / x1, 1 / x0, Fr(-1)]):
    P0vec[a][a] = val
assert mmul(mT(P0vec), mmul(ETA, P0vec)) == ETA
Pvec = [mmul(P0vec, Rk) for Rk in Rvec]
# lowered form P = eta Pvec eta^{-1}
P = [mmul(ETA, mmul(Pk, ETAI)) for Pk in Pvec]
# graded in lowered form?  P_k = sum_b p_{k,b} E_{b+k,b}
for k in range(K + 1):
    for a in range(5):
        for b in range(5):
            if (a - b - k) % 5 != 0:
                assert P[k][a][b] == 0
d = P[0][3][3]
assert d == P0vec[0][0]
p00 = P[0][0][0]
assert d * p00 == 1


def Yprime_coef(j, a):
    """[z^j] (P Y^C)_a mod v^V as vpoly:  sum_n P_{j+n}[a][b] [z^{-n}] Y^C_b."""
    res = vz()
    for n in range(0, mmax + 1):
        if j + n < 0 or j + n > K:
            continue
        for b in range(5):
            c, e = YC_coef(n, b)
            if c != 0 and P[j + n][a][b] != 0:
                res = vadd(res, vmono(P[j + n][a][b] * c, e))
    return res


# q_T = (Pvec)^{-1}(-z 1) = -(z/d) Rvec^{-1}(z) 1
Svec = sinv(Rvec, K)
qT = {}
for k in range(1, K + 1):
    for b in range(5):
        c = -Svec[k - 1][b][0] / d
        if c != 0:
            qT[(k, b)] = vmono(c, 0)
TP = vscal(Fr(-1), Yprime_coef(-1, 0))
for i in range(0, imax + 1):
    lhs = Pi(i, qT)
    rhs = vscal(Fr((-1) ** (i + 1)), Yprime_coef(i - 1, 0))
    assert lhs == rhs, (i, lhs, rhs)
print("Lemma T(b): Pi_i(q_T) = (-1)^{i+1}[z^{i-1}]Y'_0 verified for i<=%d (mod v^%d); Pi_1(q_T)=0:" % (imax, V),
      viszero(Pi(1, qT)), "; Pi_2(q_T) nonzero without (A):", not viszero(Pi(2, qT)))
# (d): T_P - p00 v in v^2
diff = vadd(TP, vmono(-p00, 1))
assert diff[0] == 0 and diff[1] == 0
print("Lemma T(d): T_P = p_{0,0} v + O(v^2) ok; T_P =", TP[:4])

# Cor E': (Pi_i (x) 1) E_P  ==  (-1)^{i} * ( -Rvec^{-1}(psi') * [ z^{-(i+1)} P0vec^{-1} eta^{-1} (Y')_{>= i+1} ](psi') )
E = edge_coeffs(Svec, imax + mmax, Lmax)
P0vecI = mzero()
for a in range(5):
    P0vecI[a][a] = 1 / P0vec[a][a]
for i in range(0, imax + 1):
    k = i + 1
    # LHS
    LHS = []
    for l in range(Lmax + 1):
        vec = []
        for b in range(5):
            q = {}
            for kk in range(0, imax + mmax + 1):
                for a in range(5):
                    c = E[(kk, l)][a][b]
                    if c != 0:
                        q[(kk, a)] = vmono(c, 0)
            vec.append(Pi(i, q))
        LHS.append(vec)
    # g_j := [z^j] of z^{-k} P0vec^{-1} eta^{-1} (Y')_{>=k},  j>=0
    g = []
    for j in range(Lmax + 1):
        Yp = [Yprime_coef(j + k, a) for a in range(5)]
        w = []
        for a in range(5):
            acc = vz()
            for c in range(5):
                coef = sum(P0vecI[a][e] * ETAI[e][c] for e in range(5))
                if coef != 0:
                    acc = vadd(acc, vscal(coef, Yp[c]))
            w.append(acc)
        g.append(w)
    RHS = []
    sign = Fr((-1) ** (k + 1))          # Omega_s(z^{-k}J^C, .) = (-1)^{k+1} Pi_{k-1}
    for l in range(Lmax + 1):
        vec = []
        for a in range(5):
            acc = vz()
            for j in range(l + 1):
                for c in range(5):
                    coef = -Svec[l - j][a][c]
                    if coef != 0:
                        acc = vadd(acc, vscal(coef, g[j][c]))
            vec.append(vscal(sign, acc))      # divide by (-1)^{k+1} == multiply
        RHS.append(vec)
    assert LHS == RHS, (i,)
print("Cor E' formula: (Pi_i x 1)E_P = (-1)^i * ( -R^{-1} z^{-i-1} P0^{-1} eta^{-1} (Y')_{>=i+1} ) verified, i<=%d" % imax)

# Lemma theta: Theta = 1+N, inverse sum (-N)^r on finitely supported sequences (graded setting)
def Theta(chi):
    L = len(chi)
    return [Pi(i, {(j, 1): chi[j] for j in range(L)}) for i in range(L)]


def Nop(chi):
    L = len(chi)
    out = []
    for i in range(L):
        acc = vz()
        r = 1
        while i + 5 * r < L:
            pir = Fr(1, 4 ** r * factorial(r) * factorial(4 * r))
            acc = vadd(acc, vscal(pir, vmul(vmono(1, 4 * r), chi[i + 5 * r])))
            r += 1
        out.append(acc)
    return out


chi = [vmono(rnd(), random.randint(0, 2)) for _ in range(12)]
th = Theta(chi)
assert th == [vadd(a, b) for a, b in zip(chi, Nop(chi))]
inv = th
acc = th
term = th
for r_ in range(1, 4):
    term = [vscal(Fr(-1), t) for t in Nop(term)]
    acc = [vadd(a, b) for a, b in zip(acc, term)]
assert acc == chi
print("lem:tra-theta: Theta = 1+N and Theta^{-1} = sum (-N)^r verified on a random sequence (mod v^%d)" % V)
print("ALL CHECKS PASSED")
