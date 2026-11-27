
# Finite-order checks for RIGIDITY_PROOF.md (evidence only). Exact rational arithmetic, no sympy.
from fractions import Fraction as Fr
import random, math
random.seed(20261103)

def rnd():
    while True:
        v = Fr(random.randint(-30, 30), random.randint(1, 12))
        if v != 0:
            return v

# ---------- series helpers (lists of Fractions, index = power) ----------
def smul(a, b, N):
    c = [Fr(0)] * (N + 1)
    for i, x in enumerate(a[:N + 1]):
        if x == 0: continue
        for j, y in enumerate(b[:N + 1 - i]):
            c[i + j] += x * y
    return c

def sinv(a, N):
    # a[0] != 0
    b = [Fr(0)] * (N + 1)
    b[0] = 1 / a[0]
    for n in range(1, N + 1):
        s = sum(a[k] * b[n - k] for k in range(1, n + 1))
        b[n] = -s / a[0]
    return b

# ---------- A. beta_n, P_g ----------
N = 12
# S(y) = y/(2 sinh(y/2)) = 1 / sum_k y^{2k}/(4^k (2k+1)!)
den = [Fr(0)] * (N + 1)
for k in range(0, N // 2 + 1):
    den[2 * k] = Fr(1, 4 ** k * math.factorial(2 * k + 1))
S = sinv(den, N)
s = [S[2 * k] for k in range(N // 2 + 1)]
print("A. S(y) = y/(2 sinh(y/2)) even coefficients s_k:", [str(v) for v in s[:5]])
# P_g = [t^{2g}] S(e1 t) S(e2 t) = sum_{i+j=g} s_i s_j e1^{2i} e2^{2j}
def P_poly(g):
    return {(2 * i, 2 * (g - i)): s[i] * s[g - i] for i in range(g + 1)}
for g in (2, 3):
    P = P_poly(g)
    fac = math.factorial(2 * g - 3)
    print(f"   P_{g} =", {k: str(v) for k, v in P.items()})
    print(f"   beta_{2*g-2} = (2g-3)! P_g/(e1 e2) has e1 e2 * beta =", {k: str(fac * v) for k, v in P.items()})
print("   BS: beta_2 = (7 e1^4 + 10 e1^2 e2^2 + 7 e2^4)/(5760 e1 e2)  ->", [str(Fr(7, 5760)), str(Fr(10, 5760))])
print("   BS: beta_4 = -(31 e1^6 + 49 e1^4 e2^2 + 49 e1^2 e2^4 + 31 e2^6)/(161280 e1 e2) ->", [str(Fr(-31, 161280)), str(Fr(-49, 161280))])
# T_F specialisation (x,1): gap_2 = (7 + 10 x^2 + 7 x^4)/(5760 x)
x = rnd()
P2 = P_poly(2)
val = sum(v * x ** k[0] for k, v in P2.items()) / x
print("   gap_2 at (x,1) vs (7+10x^2+7x^4)/(5760x):", val == (7 + 10 * x ** 2 + 7 * x ** 4) / (5760 * x))

# ---------- B. Disc identities ----------
def Disc(Q, a0, a1, a2):
    return a0 ** 2 + (2 * a0 * (a1 + a2) + 4 * a1 * a2) * Q + (a1 - a2) ** 2 * Q ** 2
okB = True
for _ in range(50):
    a0, a1, a2, Q = rnd(), rnd(), rnd(), rnd()
    if Q == 1: continue
    E, A, B = a0 + a1 + a2, a0 + a2, a0 + a1
    tau = 1 - 1 / Q
    # (L1-L2)^2 tau^2 = (E - a0 tau)^2 - 4 a1 a2 tau
    lhs = Q ** 2 * ((E - a0 * tau) ** 2 - 4 * a1 * a2 * tau)
    okB &= (lhs == Disc(Q, a0, a1, a2))
    # discriminant of Disc in Q
    okB &= ((2 * a0 * (a1 + a2) + 4 * a1 * a2) ** 2 - 4 * (a1 - a2) ** 2 * a0 ** 2 == 16 * a1 * a2 * A * B)
    # Disc(1) = E^2
    okB &= (Disc(Fr(1), a0, a1, a2) == E ** 2)
    # (a1a2 - E L1)(a1a2 - E L2) = a1 a2 A B via Vieta: L1+L2 = E/tau - a0, L1 L2 = a1a2/tau
    s1, p = E / tau - a0, a1 * a2 / tau
    okB &= (a1 ** 2 * a2 ** 2 - E * a1 * a2 * s1 + E ** 2 * p == a1 * a2 * A * B)
print("B. Disc = Q^2 tau^2 (L1-L2)^2, disc(Disc) = 16 a1 a2 A B, Disc(1) = E^2, prod(a1a2 - E L_i) = a1 a2 A B:", okB)

# ---------- C. critical polynomial, L(U), Delta(U) mod P(U) ----------
def padd(a, b):
    n = max(len(a), len(b)); r = [Fr(0)] * n
    for i, v in enumerate(a): r[i] += v
    for i, v in enumerate(b): r[i] += v
    return r
def pmul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            r[i + j] += v * w
    return r
def pscale(a, c): return [c * v for v in a]
def pmod(a, m):
    a = a[:]
    while len(a) >= len(m):
        if a[-1] == 0:
            a.pop(); continue
        c = a[-1] / m[-1]; d = len(a) - len(m)
        for i, v in enumerate(m): a[d + i] -= c * v
        a.pop()
    return a
def iszero(a): return all(v == 0 for v in a)
okC = True
for _ in range(30):
    a0, a1, a2, Q, e1 = rnd(), rnd(), rnd(), rnd(), rnd()
    if Q == 1 or a1 == a2: continue
    E, A, B = a0 + a1 + a2, a0 + a2, a0 + a1
    e2 = E - e1
    tau = 1 - 1 / Q
    P = [-a2, -(a0 + Q * (a1 - a2)), B * Q]          # B Q U^2 - (a0 + Q(a1-a2)) U - a2
    # U W_U (1-U)(1-QU) with W = -a2 log U + A log(1-U) + a1 log(1-QU):
    one_U = [Fr(1), Fr(-1)]; one_QU = [Fr(1), -Q]
    expr = padd(padd(pscale(pmul(one_U, one_QU), -a2), pscale(pmul([Fr(0), Fr(1)], one_QU), -A)), pscale(pmul([Fr(0), Fr(1)], one_U), -a1 * Q))
    okC &= (padd(expr, pscale(P, -1)) == [Fr(0)] * 3 or iszero(padd(expr, pscale(P, -1))))
    # discriminant of P equals Disc
    okC &= ((a0 + Q * (a1 - a2)) ** 2 + 4 * B * Q * a2 == Disc(Q, a0, a1, a2))
    # L(U) = -a1 Q U/(1-QU) satisfies tau L^2 - (E - a0 tau) L + a1 a2 = 0 mod P: numerator
    L_num = [Fr(0), -a1 * Q]                         # L = L_num/(1-QU)
    num = padd(padd(pscale(pmul(L_num, L_num), tau), pscale(pmul(L_num, one_QU), -(E - a0 * tau))), pscale(pmul(one_QU, one_QU), a1 * a2))
    okC &= iszero(pmod(num, P))
    # Delta(U) = e1 e2 (2L - (E/tau - a0))(a1a2 - E L) == a1 A e1 e2 W_2(U) mod P, with
    # W_2 = -A U/(1-U)^2 - a1 Q U/(1-QU)^2 ; multiply everything by (1-U)^2 (1-QU)^2
    twoL = padd(pscale(L_num, 2), pscale(one_QU, -(E / tau - a0)))           # /(1-QU)
    a1a2EL = padd(pscale(one_QU, a1 * a2), pscale(L_num, -E))                 # /(1-QU)
    lhs = pscale(pmul(pmul(twoL, a1a2EL), pmul(one_U, one_U)), e1 * e2)       # * (1-U)^2 (1-QU)^2
    W2num = padd(pscale(pmul([Fr(0), Fr(1)], pmul(one_QU, one_QU)), -A), pscale(pmul([Fr(0), Fr(1)], pmul(one_U, one_U)), -a1 * Q))
    rhs = pscale(W2num, a1 * A * e1 * e2)
    okC &= iszero(pmod(padd(lhs, pscale(rhs, -1)), P))
    # eta^-1 consistency: Delta = e1 e2 (a0 a1 a2 + 2 a1 a2 L - E L^2) mod P
    alt = padd(padd(pscale(pmul(one_QU, one_QU), a0 * a1 * a2), pscale(pmul(L_num, one_QU), 2 * a1 * a2)), pscale(pmul(L_num, L_num), -E))
    alt = pscale(pmul(alt, pmul(one_U, one_U)), e1 * e2)
    okC &= iszero(pmod(padd(lhs, pscale(alt, -1)), P))
    # U = 1 and QU = 1 in P
    okC &= (sum(c * Fr(1) ** i for i, c in enumerate(P)) == A * (Q - 1))
    okC &= (sum(c * (1 / Q) ** i for i, c in enumerate(P)) == a1 * (1 - Q) / Q)
    # a0 = 0: P = a2 * (x Q U^2 + Q(1-x) U - 1), x = a1/a2
    P0 = [-a2, -(Q * (a1 - a2)), a1 * Q]; xx = a1 / a2
    okC &= (P0 == pscale([Fr(-1), Q * (1 - xx), xx * Q], a2))
print("C. P(U) from the phase, disc P = Disc, L(U) in the QDE quadratic, Delta = a1 A e1 e2 W_2 = e1e2(a0a1a2+2a1a2L-EL^2) mod P, P(1) = A(Q-1), P(1/Q) = a1(1-Q)/Q, a0=0 rescaling:", okC)

# ---------- C'. Q = oo limits (floats, evidence) ----------
a0, a1, a2, e1 = 1.3, 2.1, -0.7, 0.45
E = a0 + a1 + a2; A = a0 + a2; B = a0 + a1; e2 = E - e1
for Q in (1e6, 1e9):
    disc = (a0 + Q * (a1 - a2)) ** 2 + 4 * B * Q * a2
    U1 = ((a0 + Q * (a1 - a2)) + math.sqrt(disc)) / (2 * B * Q)
    U2 = ((a0 + Q * (a1 - a2)) - math.sqrt(disc)) / (2 * B * Q)
    def Lf(U): return -a1 * Q * U / (1 - Q * U)
    def W2(U): return -A * U / (1 - U) ** 2 - a1 * Q * U / (1 - Q * U) ** 2
    print(f"C'. Q={Q:.0e}: U1={U1:.6f} (-> (a1-a2)/B={(a1-a2)/B:.6f}), Q*U2={Q*U2:.6f} (-> -a2/(a1-a2)={-a2/(a1-a2):.6f}); L1={Lf(U1):.6f}, L2={Lf(U2):.6f} (-> a1,a2 = {a1},{a2})")
    print(f"    Delta1 = a1 A e1 e2 W2(U1) = {a1*A*e1*e2*W2(U1):.6f} vs -a1 B (a1-a2) e1 e2 = {-a1*B*(a1-a2)*e1*e2:.6f};  Delta2 = {a1*A*e1*e2*W2(U2):.6f} vs a2 A (a1-a2) e1 e2 = {a2*A*(a1-a2)*e1*e2:.6f}")

# ---------- D. Hensel root Q_1 of Disc in K[[a0]] ----------
# series in a0 with coefficients in Q (a1, a2 fixed random rationals)
Nn = 7
def ser_add(a, b): return [x + y for x, y in zip(a, b)]
def ser_scal(a, c): return [c * x for x in a]
okD = True
for _ in range(5):
    a1, a2 = rnd(), rnd()
    if a1 == a2: continue
    a0s = [Fr(0), Fr(1)] + [Fr(0)] * (Nn - 1)       # the series "a0"
    one = [Fr(1)] + [Fr(0)] * Nn
    def DiscS(Qs):
        c1 = ser_add(ser_scal(a0s, 2 * (a1 + a2)), ser_scal(one, 4 * a1 * a2))
        return ser_add(ser_add(smul(a0s, a0s, Nn), smul(c1, Qs, Nn)), ser_scal(smul(Qs, Qs, Nn), (a1 - a2) ** 2))
    def DiscpS(Qs):
        c1 = ser_add(ser_scal(a0s, 2 * (a1 + a2)), ser_scal(one, 4 * a1 * a2))
        return ser_add(c1, ser_scal(Qs, 2 * (a1 - a2) ** 2))
    Qs = [Fr(0)] * (Nn + 1)
    for it in range(8):
        Qs = ser_add(Qs, ser_scal(smul(DiscS(Qs), sinv(DiscpS(Qs), Nn), Nn), -1))
    okD &= iszero(DiscS(Qs))
    okD &= (Qs[0] == 0 and Qs[1] == 0 and Qs[2] == -1 / (4 * a1 * a2) and Qs[3] == (a1 + a2) / (8 * a1 ** 2 * a2 ** 2))
    c4 = Qs[4]
    pred4 = -(5 * a1 ** 2 + 6 * a1 * a2 + 5 * a2 ** 2) / (64 * a1 ** 3 * a2 ** 3)
    okD &= (c4 == pred4)
print("D. Hensel: Disc(Q_1) = 0 to O(a0^8), Q_1 = -a0^2/(4a1a2) + (a1+a2)a0^3/(8a1^2a2^2) - (5a1^2+6a1a2+5a2^2)a0^4/(64a1^3a2^3) + ... :", okD, " e.g. Q_1 coefficients:", [str(v) for v in Qs[:6]])

# ---------- E. Li_{-n} ----------
def Li_poly(n):
    # Li_{-n}(Q) = A_n(Q)/(1-Q)^{n+1}; A_0 = Q; A_{n+1} = Q[(1-Q)A_n' + (n+1)A_n]
    Apoly = [Fr(0), Fr(1)]
    for m in range(n):
        d = [i * Apoly[i] for i in range(1, len(Apoly))]
        t1 = pmul([Fr(1), Fr(-1)], d) if d else [Fr(0)]
        t2 = pscale(Apoly, m + 1)
        Apoly = pmul([Fr(0), Fr(1)], padd(t1, t2))
    return Apoly
def Li_rat(n, Q):
    Apoly = Li_poly(n)
    return sum(c * Q ** i for i, c in enumerate(Apoly)) / (1 - Q) ** (n + 1)
okE = True
for n in range(1, 8):
    Apoly = Li_poly(n)
    while Apoly and Apoly[-1] == 0: Apoly.pop()
    okE &= (Apoly[0] == 0 and len(Apoly) - 1 == n)
    for _ in range(10):
        Q = rnd()
        if Q == 1: continue
        okE &= (Li_rat(n, 1 / Q) == (-1) ** (n + 1) * Li_rat(n, Q))
    # series check at Q = 1/3
    Q = Fr(1, 3); ser = sum(Fr(d) ** n * Q ** d for d in range(1, 200))
    dif = ser - Li_rat(n, Q)
    okE &= (abs(dif) < Fr(1, 10**30))
    if not okE: print('   E fail at n', n, float(dif))
# n = 0 behaves differently at oo
okE &= (Li_rat(0, Fr(5)) == -1 - Li_rat(0, Fr(1, 5)))
print("E. Li_{-n} = A_n/(1-Q)^{n+1}, A_n(0)=0, deg A_n = n, Li_{-n}(1/Q) = (-1)^{n+1} Li_{-n}(Q) (n>=1), series agreement:", okE)
# Laurent expansion of Li_{-n}(e^t) at t = 0: polar part n! t^{-n-1} only.
M = 14
expt = [Fr(1, math.factorial(k)) for k in range(M + 1)]
g = [-Fr(1, math.factorial(k + 1)) for k in range(M + 1)]        # (1 - e^t)/t = -sum t^k/(k+1)!
f = smul(expt, sinv(g, M), M)                                      # e^t / ((1-e^t)/t) ; then Li_0(e^t) = f/t
# Laurent series as dict power -> coeff, starting at t^{-1}
lau = {k - 1: f[k] for k in range(M + 1)}
okE2 = True
for n in range(1, 7):
    cur = dict(lau)
    for _ in range(n):
        cur = {k - 1: k * v for k, v in cur.items() if k != 0}
    polar = {k: v for k, v in cur.items() if k < 0 and v != 0}
    okE2 &= (polar == {-(n + 1): Fr((-1) ** (n + 1) * math.factorial(n))})
    if n in (1, 3):
        print(f"   Li_{{-{n}}}(e^t) Laurent: polar part {dict((k,str(v)) for k,v in polar.items())}, t^0 coeff {cur[0]}, t^1 coeff {cur[1]}")
print("E'. polar part of Li_{-n}(e^t) is exactly (-1)^{n+1} n! t^{-n-1} (= +n! t^{-n-1} for odd n = 2g-3; no lower-order polar terms), n = 1..6:", okE2)

# ---------- F. h-constant cancellation: p_{2k-1}(1/v) sums ----------
okF = True
for _ in range(20):
    a0, a1, a2, e1 = rnd(), rnd(), rnd(), rnd()
    E, A, B = a0 + a1 + a2, a0 + a2, a0 + a1; e2 = E - e1
    if 0 in (E, A, B, e2): continue
    for k in (1, 2, 3):
        pw = 1 - 2 * k
        nw = [a1, A, -e1, -e2]; beta_c = [E, -a1, -A]; w = [E, -e1, -e2]
        okF &= (sum(v ** pw for v in nw) + sum(v ** pw for v in beta_c) == sum(v ** pw for v in w))
        # a0 = 0 three-torus: nw - w = (a1, a2, -E)
        nw0 = [a1, a2, -e1, -e2]; w0 = [a1 + a2, -e1, -e2]; three = [a1, a2, -(a1 + a2)]
        okF &= (sum(v ** pw for v in nw0) - sum(v ** pw for v in w0) == sum(v ** pw for v in three))
        # flop weights at Q = oo: h_{(a2-a1,-A,B)} h_nw = h_{(a2-a1, a1, B, -e1, -e2)} ; h_{(a2,-a1,a1-a2)} h_nw = h_{(a1-a2,a2,A,-e1,-e2)}
        if a1 != a2:
            okF &= (sum(v ** pw for v in [a2 - a1, -A, B]) + sum(v ** pw for v in nw) == sum(v ** pw for v in [a2 - a1, a1, B, -e1, -e2]))
            okF &= (sum(v ** pw for v in [a2, -a1, a1 - a2]) + sum(v ** pw for v in nw) == sum(v ** pw for v in [a1 - a2, a2, A, -e1, -e2]))
print("F. h_nw h_{(E,-a1,-A)} = h_{(E,-e1,-e2)} (power sums k=1..3), a0=0: h_nw h_w^-1 = h_{(a1,a2,-E)}, Q=oo flop constants:", okF)

# ---------- G. constant-map term F_g(0) alpha-dependence at g = 2 (evidence) ----------
# F_2(0) = sum_p Int_{M_2} prod_a Lambda^vee(w_a)/w_a, using Int lambda_1^3 = 1/2880... we avoid Hodge integrals here; skipped.
print("G. (F_g(0) alpha-dependence: not recomputed here; cited from cy_torus_explore critic flop_const.log as evidence.)")
