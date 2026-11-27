# Referee U05_B: independent low-order test of Section 3 parts 5-6:
#   lem:gau-R1, lem:gau-diag, def:gau-van, lem:gau-vanmod, lem:gau-vanshape, lem:gau-vanq,
#   thm:gau-G (c), prop:gau-symp, on explicit cyclic systems with one simple pole.
# The graded gauge is built from the gauge equation z P' = A P - P At directly
# (component recursion p_{k,c}' = b_{c+k} p_{k+1,c} - bt_{c-1} p_{k+1,c-1}), not from the note's proof.
import sys
from fractions import Fraction as Fr
from math import factorial
from ser import S

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30   # u-adic working precision
K = int(sys.argv[2]) if len(sys.argv) > 2 else 11   # z-order of the gauge
MODE = sys.argv[3] if len(sys.argv) > 3 else "quintic"

u = S.u(N)
one = S.const(1, N)
nu = [0, -1, 3, 2, 1]
eps = [1, -1, -1, 1, 1]


def poly(cs):
    return S({i: Fr(c) for i, c in enumerate(cs)}, N)


# ---------------------------------------------------------------- the two systems
tau = S({n: Fr(4, 5 * n + 4) for n in range(N)}, N).inv()
if MODE == "quintic":
    # periods of PF4 = u D^4 - (1-u)(2D^3 + 7/5 D^2 + 2/5 D + 24/625), D = -(1-u) d/du
    def Dq(f):
        return -(one - u) * f.d()

    def PF4(f):
        d1 = Dq(f); d2 = Dq(d1); d3 = Dq(d2); d4 = Dq(d3)
        return u * d4 - (one - u) * (2 * d3 + Fr(7, 5) * d2 + Fr(2, 5) * d1 + Fr(24, 625) * f)

    def frob(k):
        # solution f = u^k + O(u^3) (k = 0,1,2), coefficient of u^m fixed by the u^{m-3} coefficient
        c = {k: Fr(1)}
        for m in range(3, N):
            f = S(c, N)
            r = PF4(f).c.get(m - 3, Fr(0))
            c[m] = -r / (m * (m - 1) ** 2 * (m - 2))
        f = S(c, N)
        res = PF4(f)
        assert all(e >= N - 5 for e in res.c), ("PF4 residual", res)
        return f

    f0, f1, f2 = frob(0), frob(1), frob(2)
    al, be = Fr(1, 3), Fr(-2, 7)          # a conifold frame (alpha, beta)
    w1 = f1 + Fr(7, 10) * f2
    w0 = f0 + al * w1 + be * f2
    tc = w1 / w0
    v = u * tau.powfr(Fr(-5, 4))
    b = [tc.d(), -(u * w0 * w0 * (one - u) * (one - u) * (one - u) * tc.d() * tc.d()).inv(), tc.d(),
         -w0 / (one - u), -w0 / (one - u)]
elif MODE == "selfdual":
    # generic self-dual pair with r != rt and b_a(0) != bt_a(0)
    v = poly([0, Fr(3, 2), Fr(1, 3), Fr(-2, 5), Fr(1, 7), 2, -1, Fr(1, 2)])
    f = poly([Fr(2), Fr(1, 2), -1, Fr(3, 4), Fr(1, 5), 1])
    g = poly([Fr(-3, 2), 1, Fr(2, 3), -1, Fr(1, 6)])
    Bp = -(v.d() * v.d() * v.d() * v.d() * v.d()) / v
    b = [f, Bp / (f * f * g * g), f, g, g]
else:
    # generic pair, not self-dual
    v = poly([0, Fr(1), Fr(1, 3), Fr(-2, 5), Fr(1, 7), 2, -1, Fr(1, 2)])
    bs = [poly([1, Fr(1, 2), -1, Fr(3, 4)]), None, poly([Fr(-2), 1, Fr(2, 3), -1]),
          poly([Fr(1, 3), 2, -1, Fr(1, 5)]), poly([Fr(5, 2), -1, 3, Fr(1, 2)])]
    Bp = -(v.d() * v.d() * v.d() * v.d() * v.d()) / v
    bs[1] = Bp / (bs[0] * bs[2] * bs[3] * bs[4])
    b = bs

vp = v.d()
bt = [vp, -vp / v, vp, -vp, -vp]
prodb = b[0] * b[1] * b[2] * b[3] * b[4]
prodbt = bt[0] * bt[1] * bt[2] * bt[3] * bt[4]
dif = prodb - prodbt
print("mode", MODE, " (EP) residual lowest exponent:", dif.val(), " precision", dif.prec)
r = (u * b[1]).coef(0)
rt = (u * bt[1]).coef(0)
print("residues r, rt:", r, rt, "  b_a(0), a=0,2,3,4:", [b[a].coef(0) for a in (0, 2, 3, 4)],
      " bt:", [bt[a].coef(0) for a in (0, 2, 3, 4)])
if MODE == "quintic":
    chk = prodb + (u * (one - u) * (one - u) * (one - u) * (one - u) * (one - u)).inv()
    print("prod b = -1/(u(1-u)^5):", chk.iszero(), " tc =", tc.trunc(4))


# ---------------------------------------------------------------- vanishing solutions
def vanishing(bb, nmax):
    ys = [[S({}, N) for _ in range(5)] for _ in range(nmax + 1)]
    ys[0][1] = S.const(-1, N)
    for n in range(nmax):
        for a in range(5):
            integrand = bb[a] * ys[n][(a + 1) % 5]
            assert integrand.val() >= 0, ("integrand has a pole", n, a)   # lem:gau-vanq (b)
            ys[n + 1][a] = integrand.integ().trunc(N)
    return ys


NMAX = (5 * N) // 4 + 2
yA = vanishing(b, NMAX)
yT = vanishing(bt, NMAX)

# lem:gau-vanq (a): support and exact order; leading coefficients
okv = True
for n in range(NMAX + 1):
    for a in range(5):
        e = Fr(4 * n + nu[a] + 1, 5)
        if n % 5 == (1 - a) % 5:
            if e < N:
                okv &= (yA[n][a].val() == e) and (yT[n][a].val() == e)
        else:
            okv &= yA[n][a].iszero() and yT[n][a].iszero()
print("lem:gau-vanq (a) support/orders:", okv)

# lem:gau-vanshape closed form at v = v(u), compared with Y^C_{At} (uniqueness, used in thm:gau-G)
vpow = [S.const(1, N)]
for i in range(1, NMAX + 2):
    vpow.append((vpow[-1] * v).trunc(N))
oks = True
for n in range(NMAX + 1):
    m = n // 5
    for a in range(5):
        if n % 5 == (nu[a] + 1) % 5:
            cf = vpow[n - m] * Fr(eps[a] * (-1) ** n, 4 ** m * factorial(m) * factorial(n - m))
        else:
            cf = S({}, N)
        dd = cf - yT[n][a]
        oks &= dd.iszero()
print("Y^C(v(u),z) (closed form of lem:gau-vanshape) == Y^C_{At} (integral recursion):", oks)

# direct residue definition of I^C_m against the closed form, and (i),(ii) of lem:gau-vanmod, in Q[1/z][[v]]
def IC(m, vmax=14):
    d = {}
    i = 0
    while m + 4 * i + 1 <= vmax:
        j = m + 4 * i + 1
        if j >= 0:
            # (-1/(4z))^i/i! * (-v/z)^j/j!  -> key (v-exp, zinv-exp)
            d[(j, i + j)] = Fr((-1) ** (i + j), 4 ** i * factorial(i) * factorial(j))
        i += 1
    return d


def dv(d):
    return {(j - 1, n): j * x for (j, n), x in d.items() if j > 0}


def sh(d, dvexp, dz):  # multiply by v^dvexp z^{dz}  (z^{dz}: zinv-exp decreases by dz)
    return {(j + dvexp, n - dz): x for (j, n), x in d.items()}


def comb(terms, vmax):
    out = {}
    for c, d in terms:
        for k, x in d.items():
            if k[0] <= vmax:
                out[k] = out.get(k, 0) + c * x
    return {k: x for k, x in out.items() if x != 0}


okm = True
for m in range(-6, 8):
    okm &= not comb([(1, sh(dv(IC(m)), 0, 1)), (1, IC(m - 1))], 11)                       # z d_v I_m = -I_{m-1}
    okm &= not comb([(m, sh(IC(m - 1), 0, 1)), (-1, IC(m + 3)), (1, sh(IC(m - 2), 1, 0))], 11)  # (ii)
print("lem:gau-vanmod (i),(ii) from the explicit series:", okm)

# ---------------------------------------------------------------- the graded gauge
def hom(k):
    # h_c with b_{c+k-1} h_c = bt_{c-1} h_{c-1}, h_0 = 1
    h = [one]
    for c in range(1, 5):
        h.append(h[c - 1] * bt[c - 1] / b[(c + k - 1) % 5])
    # closing relation (needs (EP))
    cl = b[(k - 1) % 5] * h[0] - bt[4] * h[4]
    assert cl.iszero(), "closing relation fails"
    return h


free = {0: Fr(1), 5: Fr(2, 3), 10: Fr(-5, 7), 15: Fr(1, 2)}   # free constants (arbitrary)
P = []      # P[k][c] = p_{k,c}
g = [S({}, N) for _ in range(5)]
hk = hom(0)
for k in range(K + 1):
    hk1 = hom(k + 1)
    # xi_c = phi + Xi_c
    sig = [g[c] / (b[(c + k - 1) % 5] * hk[c]) for c in range(5)]
    tot = sig[0] + sig[1] + sig[2] + sig[3] + sig[4]
    assert tot.iszero(), ("solvability fails at level", k, tot)
    Xi = [S({}, N)]
    for c in range(1, 5):
        Xi.append(Xi[-1] + sig[c])
    w = [(b[(c + k) % 5] * hk1[c]).inv() for c in range(5)]
    alpha = sum((w[c] * hk[c] for c in range(5)), S({}, N))
    beta = sum((w[c] * hk[c].d() for c in range(5)), S({}, N))
    gamma = -sum((w[c] * (hk[c] * Xi[c]).d() for c in range(5)), S({}, N))
    B = beta / alpha
    G = gamma / alpha
    assert B.val() >= -1, ("irregular singular point", k)
    rho = B.c.get(-1, Fr(0))
    n0 = G.val() + 1 if G.c else G.prec + 1
    nstar = None
    if rho.denominator == 1:
        nstar = -int(rho)
        n0 = min(n0, nstar)
    pr = min(B.prec + min(n0, 0), G.prec) + 1
    phi = {}
    for n in range(n0, pr):
        rhs = G.c.get(n - 1, Fr(0)) - sum(B.c.get(n - 1 - l, Fr(0)) * x for l, x in phi.items())
        piv = n + rho
        if piv == 0:
            assert rhs == 0, ("resonance obstruction at level", k, rhs)
            phi[n] = free.get(k, Fr(0))
        else:
            phi[n] = rhs / piv
    phi = S(phi, pr)
    pk = [hk[c] * (phi + Xi[c]) for c in range(5)]
    for c in range(5):
        assert pk[c].val() >= 0, ("gauge coefficient not regular", k, c, pk[c])
    P.append(pk)
    g = [pk[c].d() for c in range(5)]
    hk = hk1

minprec = min(P[k][c].prec for k in range(K + 1) for c in range(5))
print("graded gauge built to z^%d; regular; minimal u-precision %d" % (K, minprec))

# direct check of the gauge equation componentwise: p_{k,c}' = b_{c+k} p_{k+1,c} - bt_{c-1} p_{k+1,c-1}
okg = True
for k in range(-1, K):
    for c in range(5):
        lhs = P[k][c].d() if k >= 0 else S({}, N)
        rhs = b[(c + k) % 5] * P[k + 1][c] - bt[(c - 1) % 5] * P[k + 1][(c - 1) % 5]
        okg &= (lhs - rhs).iszero()
print("gauge equation z P' = A P - P At holds (to working precision):", okg)

# ---------------------------------------------------------------- lem:gau-R1
okr = True
c11 = {}
for k in range(K + 1):
    # (P_k)_{a,1} = p_{k,1} at a = 1+k ; (P_k)_{2,b} = p_{k,2-k} at b = 2-k
    col = P[k][1].coef(0)
    row = P[k][(2 - k) % 5].coef(0)
    if k % 5 == 0:
        c11[k] = col
        okr &= (r * row == rt * col)
    else:
        okr &= (col == 0 and row == 0)
print("lem:gau-R1 (column 1, row 2, r P22 = rt P11 at u = 0):", okr, "  P_11(0,z) coefficients:", c11)

# ---------------------------------------------------------------- lem:gau-diag
d = [one]
for a in range(4):
    d.append(d[a] * bt[a] / b[a])
prodp0 = P[0][0] * P[0][1] * P[0][2] * P[0][3] * P[0][4]
const = prodp0.coef(0)
okd = (prodp0 - const).iszero()
rho = [d[a].coef(0) / d[1].coef(0) for a in range(5)]
pred = P[0][1].coef(0) ** 5 * rho[0] * rho[1] * rho[2] * rho[3] * rho[4]
okd &= (const == pred)
okd &= all(P[0][a].coef(0) == P[0][1].coef(0) * rho[a] for a in range(5))
print("lem:gau-diag (prod p_{0,a} constant = p_{0,1}(0)^5 prod rho_a; p_{0,a}(0) = p_{0,1}(0) rho_a):", okd)
if MODE == "quintic":
    # normalise to p_{0,1}(0) = 1 and compare p_{0,0} with tau^{-1/4}/w0
    p00 = P[0][0] * (Fr(1) / P[0][1].coef(0))
    target = tau.powfr(Fr(-1, 4)) / w0
    print("  p_{0,a}(0):", [P[0][a].coef(0) for a in range(5)],
          "  p_{0,0} == tau^{-1/4}/varpi_0:", (p00 - target).iszero(), " precision", (p00 - target).prec)

# ---------------------------------------------------------------- thm:gau-G (c)
# coefficient of z^j of (P Yt)_a = sum_{k-n=j} p_{k,a-k} yt^{(n)}_{a-k};  of P11(0,z) Y_A: sum_{5m-n=j} c11[5m] y^{(n)}_a
okG = True
nchecked = 0
for j in range(-NMAX, K + 1):
    for a in range(5):
        lhs = S({}, N)
        for k in range(0, K + 1):
            n = k - j
            if 0 <= n <= NMAX:
                lhs = lhs + P[k][(a - k) % 5] * yT[n][(a - k) % 5]
        rhs = S({}, N)
        for k in range(0, K + 1, 5):
            n = k - j
            if 0 <= n <= NMAX:
                rhs = rhs + c11[k] * yA[n][a]
        # terms with k > K have n > K - j, hence u-order >= 4(K-j+1)/5
        bound = min((4 * (K - j + 1) + 4) // 5, lhs.prec, rhs.prec, minprec)
        dd = (lhs - rhs).trunc(bound)
        nchecked += max(bound, 0)
        if not dd.iszero():
            okG = False
            print("  (G) mismatch at z^%d, row %d:" % (j, a), dd)
print("thm:gau-G (c)  P Y^C_{At} = P_11(0,z) Y^C_A  (", nchecked, "coefficients compared):", okG)

# ---------------------------------------------------------------- prop:gau-symp
selfdual = (b[0] - b[2]).iszero() and (b[3] - b[4]).iszero()
# S_m = sum_{k+l=m} (-1)^k P_k^T eta P_l ; entry (c,d): sum_a (P_k)_{a,c} eta_{a,a'} (P_l)_{a',d},
# (P_k)_{a,c} = p_{k,c} delta_{a,c+k};  eta_{a,a'} = 5 delta_{a+a'=3}  =>  need c+k + d+l = 3 mod 5
cc = {}
for m in range(K + 1):
    cc[m] = sum(((-1) ** k) * c11.get(k, Fr(0)) * c11.get(m - k, Fr(0)) for k in range(m + 1))
oksym = True
for m in range(K + 1):
    for c in range(5):
        for dd_ in range(5):
            tot = S({}, N)
            for k in range(m + 1):
                l = m - k
                if (c + k + dd_ + l) % 5 == 3:
                    tot = tot + ((-1) ** k) * 5 * P[k][c] * P[l][dd_]
            target = (rt / r) * cc[m] * (5 if (c + dd_) % 5 == 3 else 0)
            if not (tot - target).iszero():
                oksym = False
print("self-dual:", selfdual, "  P(-z)^T eta P(z) = (rt/r) c(z)c(-z) eta through z^%d:" % K, oksym)
