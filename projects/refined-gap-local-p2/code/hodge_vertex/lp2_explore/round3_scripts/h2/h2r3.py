# h2r3.py -- round 3, task H2: independent exact checks (Python 3, fractions only).
# Rows built literally from Definition 7 of ROUND2_PROOFS_rows_on_TF.md, starting from the closed form
#   P0^B = (v/K) S(v,z),  S = Sum_l r_l (vz)^l / prod_{j<=l+1} (1 + j z/E)    (Prop. 8; the constant K is dropped)
# with r_l computed here from BS's (s1,s2) Wick sum (own implementation).
# Variables: v = 1/q = (1+3t)/(E t); theta = y d/dy = theta(v) d/dv, theta(v) = 1/E - v/3;
# Delta_c = 1+27y = 1 - (1 - E v/3)^-3; Y = -E/3 + 1/v; c1 = -wc(Delta_c) (on shell) or -wc(1+Q Delta_c) (off shell).
# usage: py h2r3.py e1 e2 NV NZ [Q]
import sys
from fractions import Fraction as F

e1, e2 = F(sys.argv[1]), F(sys.argv[2]); E = e1 + e2
NV = int(sys.argv[3]); NZ = int(sys.argv[4]); Q = F(sys.argv[5]) if len(sys.argv) > 5 else F(0)
print(f"e1={e1} e2={e2} E={E} NV={NV} NZ={NZ} Q={Q}  ({'ON SHELL' if Q == 0 else 'OFF SHELL c1=-wc(1+Q Del)'})")

# ---------------- bivariate truncated Laurent series in x (v or t_c), z-polynomials to z^NZ -------------
class S2:
    def __init__(s, val, co, N):  # co[i] = list of NZ+1 coeffs of x^(val+i); known for exponents < N
        s.val, s.N = val, N
        s.co = [list(c) for c in co[:max(0, N - val)]]
    def get(s, e):
        i = e - s.val
        if e >= s.N: raise ValueError("precision")
        return s.co[i] if 0 <= i < len(s.co) else [F(0)] * (NZ + 1)
def zero(N=10**6): return S2(0, [], N)
def zpoly(c):  # z-free constant / z-polynomial
    return S2(0, [list(c) + [F(0)] * (NZ + 1 - len(c))], 10**6)
def const(a): return zpoly([F(a)])
def uni(val, cs, N):  # z-free series
    return S2(val, [[F(c)] + [F(0)] * NZ for c in cs], N)
def add(f, g):
    v = min(f.val, g.val); N = min(f.N, g.N)
    top = N if N < 10**5 else max(f.val + len(f.co), g.val + len(g.co))
    return S2(v, [[a + b for a, b in zip(f.get(e) if e < f.N else [F(0)]*(NZ+1), g.get(e) if e < g.N else [F(0)]*(NZ+1))] for e in range(v, top)], N)
def scal(a, f): return S2(f.val, [[a * x for x in c] for c in f.co], f.N)
def sub(f, g): return add(f, scal(F(-1), g))
def zmul(c1, c2):
    r = [F(0)] * (NZ + 1)
    for i, a in enumerate(c1):
        if a == 0: continue
        for j in range(NZ + 1 - i):
            b = c2[j]
            if b: r[i + j] += a * b
    return r
def mul(f, g):
    v = f.val + g.val; N = min(f.N + g.val, g.N + f.val, CAP + v if False else 10**6)
    N = min(N, NMAX)
    out = []
    for e in range(v, N):
        acc = [F(0)] * (NZ + 1)
        for i in range(len(f.co)):
            ef = f.val + i; eg = e - ef
            j = eg - g.val
            if j < 0 or j >= len(g.co): continue
            if not any(f.co[i]): continue
            p = zmul(f.co[i], g.co[j])
            acc = [a + b for a, b in zip(acc, p)]
        out.append(acc)
    return S2(v, out, N)
def strip(f):
    co = list(f.co); v = f.val
    while co and not any(co[0]): co.pop(0); v += 1
    return S2(v, co, f.N)
def inv_uni(f):  # inverse of z-free series
    g = strip(f); a0 = g.co[0][0]; R = g.N - g.val
    b = [1 / a0]
    for n in range(1, R):
        acc = F(0)
        for j in range(1, n + 1):
            if j < len(g.co): acc += g.co[j][0] * b[n - j]
        b.append(-acc / a0)
    return uni(-g.val, b, -g.val + R)
def der(f):  # d/dx
    return S2(f.val - 1, [[(f.val + i) * x for x in c] for i, c in enumerate(f.co)], f.N - 1 if f.N < 10**5 else 10**6)
def zmulz(f):  # multiply by z
    return S2(f.val, [[F(0)] + c[:NZ] for c in f.co], f.N)
def zflip(f):  # z -> -z
    return S2(f.val, [[x if n % 2 == 0 else -x for n, x in enumerate(c)] for c in f.co], f.N)
def compose(f, V, Vpows):  # f(x) with x = V(t), V of valuation 1; Vpows[k] = V^k for k in range
    acc = zero(NMAX)
    for i, c in enumerate(f.co):
        k = f.val + i
        if not any(c): continue
        acc = add(acc, mul(zpoly(c), Vpows[k]))
    return S2(acc.val, acc.co, min(acc.N, f.N))

NMAX = NV + 4
# ---------------- r_l from BS Wick sum ----------------
LMAX = NV + 4
D = (3 * LMAX) // 2 + 2
u = {(2, 2): F(1, 27), (3, 0): F(-4, 27), (1, 1): F(18, 27), (0, 3): F(-4, 27)}
G = {}; P = {(0, 0): F(1)}; ck = F(1)
for k in range(0, 2 * D + 2):
    for key, val in P.items(): G[key] = G.get(key, F(0)) + ck * val
    Qd = {}
    for (a, b), val in P.items():
        for (da, db), uv in u.items():
            if a + da > D or b + db > D: continue
            Qd[(a + da, b + db)] = Qd.get((a + da, b + db), F(0)) + val * uv
    P = Qd; ck = ck * (2 * k + 1) / (2 * k + 2)
fact = [F(1)]
for i in range(1, D + 1): fact.append(fact[-1] * i)
r = [F(0)] * (LMAX + 1)
for m in range(D + 1):
    for j in range(m // 3 + 1):
        l = m - j
        if l > LMAX or m - 3 * j > D: continue
        r[l] += G.get((m - 3 * j, m), F(0)) * fact[m] / fact[j] / F(3) ** m
print("r_l:", [str(x) for x in r[:8]])

# ---------------- closed form S(v,z), S1, S2 ----------------
def invprod(L):  # z-series of 1/prod_{j=1}^{L}(1+jz/E)
    c = [F(1)] + [F(0)] * NZ
    for j in range(1, L + 1):
        g = [(-F(j) / E) ** n for n in range(NZ + 1)]
        c = zmul(c, g)
    return c
Scl = S2(0, [[r[l] * x for x in ([F(0)] * l + invprod(l + 1))[:NZ + 1]] for l in range(NV)], NV)

# ---------------- z-free ingredients in v ----------------
thv = uni(0, [1 / E, F(-1, 3)], 10**6)
def theta(f): return mul(thv, der(f))
# Delta_c(v) = 1 - (1 - E v/3)^-3 = -Sum_{n>=1} binom(n+2,2) (E/3)^n v^n
Del = uni(0, [F(0)] + [-F((n + 1) * (n + 2), 2) * (E / 3) ** n for n in range(1, NV + 4)], NV + 4)
wcoef = [F(1)]
for n in range(1, NV + 6): wcoef.append(wcoef[-1] * (F(1, 3) + n - 1) * (F(2, 3) + n - 1) / (n * n))
def unicompose(coefs, X, N):  # Sum coefs[n] X^n, X z-free with val>=1
    acc = zero(N); Pw = const(1)
    for n, c in enumerate(coefs):
        if n > N: break
        acc = add(acc, scal(c, Pw)); Pw = mul(Pw, X)
    return S2(acc.val, acc.co, min(N, acc.N))
wc = unicompose(wcoef, Del, NV + 2)
# t_c(Delta) = Int_0^Delta wc(x)/(1-x) dx
gco = [sum(wcoef[:n + 1]) for n in range(NV + 5)]
tcco = [F(0)] + [gco[n - 1] / n for n in range(1, NV + 5)]
tc = unicompose(tcco, Del, NV + 2)
C1cf = scal(F(-1), wc)
c1 = mul(C1cf, add(const(1), scal(Q, Del)))
C1i = inv_uni(c1)
# checks
chk = strip(sub(theta(tc), C1cf)); print("check theta t_c = C1cf = -wc :", "OK" if not chk.co else "NO")
chk = strip(sub(theta(Del), sub(Del, const(1)))); print("check theta Delta_c = Delta_c - 1 :", "OK" if not chk.co else "NO")
Yv = uni(-1, [F(1), -E / 3], 10**6)
d2 = mul(scal(E / 2, sub(mul(C1i, C1i), const(1))), inv_uni(Del))
# ---------------- rows, Definition 7 ----------------
vser = uni(1, [F(1)], 10**6)
P0 = mul(vser, Scl)
def YzT(f): return add(mul(Yv, f), zmulz(theta(f)))
P1 = mul(C1i, YzT(P0))
P2 = mul(mul(Del, mul(c1, c1)), sub(YzT(P1), mul(d2, P1)))
iv = uni(-1, [F(1)], 10**6)
Ph = [mul(iv, P0), mul(iv, P1), mul(iv, P2)]   # Phat = Delta^(1/2) P (up to the dropped constant)
rho = mul(scal(-E, vser), inv_uni(tc))
sig_v = [mul(rho, x) for x in Ph]
# ---------------- pass to t_c ----------------
# invert t_c(v): v = V(t)
tcs = strip(tc); a = [tcs.get(e)[0] for e in range(1, NV + 1)]  # t_c = a1 v + a2 v^2 + ...
# reversion
Vc = [F(0), 1 / a[0]]
for n in range(2, NV + 1):
    # find Vc[n] s.t. t_c(V(t)) = t + O(t^(n+1))
    Vs = uni(1, Vc[1:] + [F(0)], n + 1)
    comp = zero(n + 1); Pw = const(1)
    for k in range(1, n + 1):
        Pw = mul(Pw, Vs); comp = add(comp, scal(a[k - 1], Pw))
    Vc.append(-comp.get(n)[0] / a[0])
V = uni(1, Vc[1:], NV + 1)
Vinv = inv_uni(V)
Vp = {0: const(1)}
for k in range(1, NV + 3): Vp[k] = mul(Vp[k - 1], V)
for k in range(1, 4): Vp[-k] = mul(Vp[-k + 1], Vinv)
def tov(f): return compose(f, V, Vp)
chk = strip(sub(tov(tc), uni(1, [F(1)], 10**6))); print("check t_c(V(t_c)) = t_c :", "OK" if not chk.co else "NO", "; v = ", [str(x) for x in Vc[1:5]])
sig = [tov(x) for x in sig_v]
NT = min(x.N for x in sig); print("sigma known in t_c to order <", NT)

# ---------------- leaf fits (psi = -z, l_j(psi) = 1/(E - j psi) = 1/(E + j z)) ----------------
def leaffit(c):
    for J in range(0, NZ - 3):
        Pz = [F(1)] + [F(0)] * NZ
        for j in range(1, J + 1): Pz = zmul(Pz, [E, F(j)] + [F(0)] * (NZ - 1))
        Nn = zmul(c, Pz)
        if any(Nn[J + 2:]): continue
        spare = NZ - (J + 2) + 1
        # polynomial part of N/Pz (degree <= 1): division
        num = Nn[:J + 2]
        # leading: N = (a1 z + a0) Pz + R, deg R < J
        a1 = num[J + 1] / Pz[J] if J + 1 <= J + 1 and Pz[J] != 0 else F(0)
        rem = [num[i] - a1 * (Pz[i - 1] if i >= 1 else F(0)) for i in range(J + 2)]
        a0 = rem[J] / Pz[J]
        rem = [rem[i] - a0 * Pz[i] for i in range(J + 2)]
        beta = []
        for j in range(1, J + 1):
            z0 = -E / j
            val = sum(rem[i] * z0 ** i for i in range(J + 2))
            den = F(1)
            for i in range(1, J + 1):
                if i != j: den *= (E + i * z0)
            beta.append(val / den)
        # psi-linear coefficient: a1 z = -a1 psi
        return dict(J=J, a_psi=-a1, a0=a0, beta=beta, spare=spare)
    return None
print("\n=== leaf form of sigma_a^(k) (t_c-expansion): sigma = a_psi psi + a0 + Sum_j beta_j l_j(psi) ===")
fits = {}
viol = []
for a_ in range(3):
    for k in range(-1, NT):
        c = sig[a_].get(k); f = leaffit(c); fits[(a_, k)] = f
        if f is None: print(f"  sigma_{a_}^({k}): NO FIT"); viol.append((a_, k)); continue
        top = max([j + 1 for j, b in enumerate(f['beta']) if b != 0], default=0)
        bnd = (k + 1) if a_ < 2 else k
        ok = (f['a_psi'] == 0) and (top <= max(bnd, 0)) and ((f['a0'] == 0) if k >= 0 else True)
        if not ok: viol.append((a_, k))
        bs = ", ".join(str(b) for b in f['beta'][:3]) + (" ..." if len(f['beta']) > 3 else "")
        print(f"  sigma_{a_}^({k}): a_psi={f['a_psi']}  a0={f['a0']}  top leaf={top} (bound {max(bnd,0)})  spare z-orders={f['spare']}  beta=[{bs}]")
print("H2-form violations (psi-linear, constants at k>=0, leaf bound sigma_0,1: k+1, sigma_2: k):", viol if viol else "none")
# Z-infinity: Sum_j beta_j/j of sigma_0^(k) = E delta_k0
zi = [sum(b / (j + 1) for j, b in enumerate(fits[(0, k)]['beta'])) for k in range(0, NT) if fits[(0, k)]]
print("Z-infinity: Sum_j beta_j/j [sigma_0^(k)], k=0..:", [str(x) for x in zi])

# ---------------- (I1), (I2) in t_c ----------------
tT = uni(1, [F(1)], 10**6); itT = uni(-1, [F(1)], 10**6)
L1 = tov(mul(C1i, Yv))
def Dc(f): return add(mul(add(L1, zmulz(itT)), f), zmulz(der(f)))
def iszero(f, upto):
    g = f; return all(not any(g.get(e)) for e in range(g.val, min(upto, g.N)))
I1 = sub(sig[1], Dc(sig[0]))
kap = tov(mul(mul(mul(c1, c1), c1), Del)); mu = tov(scal(E / 2, sub(mul(c1, c1), const(1))))
I2 = sub(sig[2], add(mul(kap, Dc(sig[1])), mul(mu, sig[1])))
print("\n(I1) sigma_1 = Dc sigma_0, Dc = (L_1 + z/t_c) + z d/dt_c, L_1 = C1i Y :", "exact to t_c^%d" % (I1.N - 1) if iszero(I1, I1.N) else "FAILS: valuation %s" % strip(I1).val)
print("(I2) sigma_2 = kappa Dc sigma_1 + mu sigma_1, kappa = c1^3 Del, mu = (E/2)(c1^2-1) :", "exact to t_c^%d" % (I2.N - 1) if iszero(I2, I2.N) else "FAILS: lowest t_c-order %s" % strip(I2).val)
# general (G) form valid for any section: sigma_1 = C1i [Y + z g (1/t_c + d/dt_c)] sigma_0, sigma_2 = Del c1^2 [Y + z g(1/t_c + d/dt_c) - d2] sigma_1, g = theta t_c
gT = tov(C1cf); YT = tov(Yv)
def Dg(f): return add(mul(YT, f), zmulz(mul(gT, add(mul(itT, f), der(f)))))
G1 = sub(sig[1], mul(tov(C1i), Dg(sig[0])))
G2 = sub(sig[2], mul(tov(mul(Del, mul(c1, c1))), sub(Dg(sig[1]), mul(tov(d2), sig[1]))))
print("(G1),(G2) general-section forms:", "OK" if iszero(G1, G1.N) and iszero(G2, G2.N) else "FAIL", " (to t_c^%d)" % (min(G1.N, G2.N) - 1))
# L_1 alternatives
L1b = tov(mul(scal(E, theta(vser)), mul(iv, C1i)))   # E C1i theta log v
s13 = tov(sub(const(1), inv_uni(mul(scal(E/3, vser), uni(0,[F(1)],10**6)))))  # placeholder unused
print("L_1 = E C1i theta(log v) = E d log v/dt_c :", "OK" if iszero(sub(L1, L1b), L1.N) else "NO", "; L_1 coefficients t_c^-1..:", [str(L1.get(e)[0]) for e in range(-1, 3)])
print("kappa:", [str(kap.get(e)[0]) for e in range(0, 4)], " mu:", [str(mu.get(e)[0]) for e in range(0, 4)])

# ---------------- symplecticity in sigma form ----------------
et00 = e1 * e2 * E ** 3 / 9; et11 = -3 * E * e1 * e2; et12 = -3 * e1 * e2
def pair(A, B):
    return add(add(scal(et00, mul(A[0], B[0])), scal(et11, mul(A[1], B[1]))), scal(et12, add(mul(A[1], B[2]), mul(A[2], B[1]))))
sigm = [zflip(x) for x in sig]
Kp = pair(sig, sigm)
target = uni(-2, [3 * e1 * e2 * E ** 3], 10**6)
print("\nsymplecticity <sigma(psi), sigma(-psi)> = 3 e1 e2 E^3 / t_c^2 :", "exact to t_c^%d, all z^0..z^%d" % (Kp.N - 1, NZ) if iszero(sub(Kp, target), Kp.N) else "FAILS")

# ---------------- residue lemma for E_00 ----------------
# Y = t_c sigma; Nb(psi,psi') = 3e1e2E^3 - <Y(psi), Y(psi')>; Ebar = Nb/(psi+psi') ; t_c^2 rho^2 Ehat_00 = Ebar
# sigma(psi) coefficients: [psi^m] = (-1)^m [z^m]
def psi_coeffs(c): return [x if m % 2 == 0 else -x for m, x in enumerate(c)]
Ysig = [[psi_coeffs(sig[a_].get(k)) for k in range(-1, NT)] for a_ in range(3)]  # Y^(k+1) = sigma^(k)
NE = NZ - 2
etab = {(0, 0): et00, (1, 1): et11, (1, 2): et12, (2, 1): et12}
def Nb_coef(n):  # [t_c^n] Nb as matrix [m][m'] (psi^m psi'^m')
    M = [[F(0)] * (NZ + 1) for _ in range(NZ + 1)]
    if n == 0: M[0][0] += 3 * e1 * e2 * E ** 3
    for aa in range(0, n + 1):
        bb = n - aa
        if aa > NT or bb > NT: continue
        for (x, y), et in etab.items():
            A = Ysig[x][aa]; B = Ysig[y][bb]
            for m in range(NZ + 1):
                if A[m] == 0: continue
                for mp in range(NZ + 1 - m):
                    M[m][mp] -= et * A[m] * B[mp]
    return M
def Ebar_direct(n):
    M = Nb_coef(n)
    Eb = [[F(0)] * (NE + 1) for _ in range(NE + 1)]
    for m in range(NE + 1):
        for k2 in range(NE + 1 - m):
            Eb[m][k2] = sum(((-1) ** k) * M[m + 1 + k][k2 - k] for k in range(k2 + 1))
    # check divisibility: Nb(-psi',psi') = 0 coefficientwise
    okdiv = all(sum(((-1) ** m) * M[m][d - m] for m in range(d + 1)) == 0 for d in range(NZ + 1))
    return Eb, okdiv
def beta_of(a_, kk):  # leaf coefficients of sigma_a^(kk) in l_j basis, kk>=0
    f = fits[(a_, kk)]; return f['beta']
def Ebar_pred(n):  # Sum_{i,l>=1} il/(E(i+l)) C_il l_i l_l',  C_il = -Sum_{k1+k2=n-2} [l_i l_l'] <sigma^(k1), sigma^(k2)>
    C = {}
    for k1 in range(0, n - 1):
        k2 = n - 2 - k1
        if k1 >= NT or k2 >= NT: return None, None
        for (x, y), et in etab.items():
            bx = beta_of(x, k1); by = beta_of(y, k2)
            for i, bi in enumerate(bx):
                for l, bl in enumerate(by):
                    C[(i + 1, l + 1)] = C.get((i + 1, l + 1), F(0)) - et * bi * bl
    Eb = [[F(0)] * (NE + 1) for _ in range(NE + 1)]
    for (i, l), cil in C.items():
        if cil == 0: continue
        coef = F(i * l) / (E * (i + l)) * cil
        for m in range(NE + 1):
            for mp in range(NE + 1 - m):
                Eb[m][mp] += coef * F(i) ** m * F(l) ** mp / E ** (m + mp + 2)
    return Eb, C
print("\n=== residue lemma: Ebar = t_c^2 rho^2 Ehat_00 ===")
for n in range(0, NT + 1):
    Ed, okdiv = Ebar_direct(n)
    Ep, C = Ebar_pred(n)
    if Ep is None: break
    eq = all(Ed[m][mp] == Ep[m][mp] for m in range(NE + 1) for mp in range(NE + 1 - m))
    supp = all((i - 1) + (l - 1) <= n - 2 for (i, l), c in C.items() if c != 0)
    zeroE = all(Ed[m][mp] == 0 for m in range(NE + 1) for mp in range(NE + 1 - m))
    msg = f"  t_c^{n}: Nb(-psi',psi')=0: {okdiv}; direct == residue formula ({(NE+1)*(NE+2)//2} entries): {eq}; support (i-1)+(l-1)<={n-2}: {supp}"
    if n <= 1: msg += f"; Ebar coefficient identically 0: {zeroE}"
    if n == 2: msg += f"; C_11/(2E) = {C.get((1,1),0)/(2*E)} vs -e1 e2 E^4/27 = {-e1*e2*E**4/27}; nonzero C: {[k for k,c in C.items() if c!=0]}"
    print(msg)
# also: Ehat_00 in the code's s-variable at s^0, s^1 (for comparison with edges_sing_out): s = E v/(E v - 3)
print("done")
# off-shell prediction of the psi-linear coefficient of sigma_2: a_psi = E Del wc^2 theta(f) / t_c, f = 1 + Q Del
fpred = tov(mul(mul(scal(E, Del), mul(wc, wc)), mul(theta(add(const(1), scal(Q, Del))), inv_uni(tc))))
cmp = [(str(fits[(2,k)]['a_psi']) if fits[(2,k)] else 'nofit', str(fpred.get(k)[0])) for k in range(-1, NT)]
print("sigma_2 psi-linear coefficient: fitted vs predicted E Del wc^2 theta(f)/t_c:", cmp, "ALL EQUAL" if all(a == b for a, b in cmp) else "DIFFER")
# ---- comparison with the code's edges_sing a_il (Ehat_00 = Ebar/(t_c^2 rho^2), expanded in s = 1+3t, s = E v/(E v - 3)) ----
KM = NT - 1
Acoef = {}
for k in range(0, KM):
    _, C = Ebar_pred(k + 2)
    Acoef[k] = {il: F(il[0] * il[1]) / (E * sum(il)) * c for il, c in C.items() if c != 0}
rho_t = tov(rho); w = inv_uni(mul(rho_t, rho_t))   # rho^-2 in t_c
B = {}
for k in range(0, KM):
    acc = {}
    for j in range(0, k + 1):
        wj = w.get(j)[0]
        for il, c in Acoef[k - j].items(): acc[il] = acc.get(il, F(0)) + wj * c
    B[k] = acc            # [t_c^k] Ehat_00 in l_i l_l'
# t_c as a series in s: v(s) = -(3/E) s/(1-s)
vs = uni(1, [-F(3) / E] * (NV + 3), NV + 4)
Vsp = {0: const(1)}
for k in range(1, NV + 3): Vsp[k] = mul(Vsp[k - 1], vs)
Ts = compose(tc, vs, Vsp)
Tsp = {0: const(1)}
for k in range(1, KM + 1): Tsp[k] = mul(Tsp[k - 1], Ts)
for n in range(0, min(KM, 5)):
    acc = {}
    for k in range(0, n + 1):
        tk = Tsp[k].get(n)[0]
        for il, c in B[k].items(): acc[il] = acc.get(il, F(0)) + tk * c
    print(f"[s^{n}] Ehat_00 a_il:", {il: str(c) for il, c in sorted(acc.items()) if c != 0})
# ---- closed forms of Phat1, Phat2 (Proof B) ----
def Sk(shift):  # Sum_l r_l (vz)^l / prod_{j<=l+1-shift}(1+jz/E)
    return S2(0, [[r[l] * x for x in ([F(0)] * l + invprod(max(l + 1 - shift, 0)))[:NZ + 1]] for l in range(NV)], NV)
S1c, S2c = Sk(1), Sk(2)
Ph1c = mul(scal(E, mul(C1i, mul(thv, iv))), S1c)
lg = mul(theta(mul(C1i, thv)), inv_uni(mul(C1i, thv)))   # theta log(C1i theta v)
br = add(mul(scal(E, mul(thv, mul(iv, iv))), S2c), mul(sub(zmulz(lg), d2), mul(iv, S1c)))
Ph2c = mul(scal(E, mul(mul(Del, c1), thv)), br)
d1 = strip(sub(Ph[1], Ph1c)); d2c = strip(sub(Ph[2], Ph2c))
print("closed form Phat1 = E C1i (theta v / v) S^(1):", "OK" if not any(any(c) for c in d1.co) else "NO", "(to v^%d)" % (d1.N - 1))
print("closed form Phat2 = (Del/C1i) E theta v [E theta v S^(2)/v^2 + (z theta log(C1i theta v) - d2) S^(1)/v]:", "OK" if not any(any(c) for c in d2c.co) else "NO", "(to v^%d)" % (d2c.N - 1))
