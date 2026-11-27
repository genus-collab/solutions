# Independent sanity checks for unit U14 (lem:con-cfu, lem:con-red, lem:con-S2, lem:con-id(b), prop:con-res first half)
# Exact arithmetic with fractions. Sparse polynomials: dict {exponent tuple: Fraction}.
from fractions import Fraction as Fr
from math import comb, factorial
import sys

NU = int(sys.argv[1]) if len(sys.argv) > 1 else 5   # truncation order in u (keep u^0..u^NU)

def add(a, b, c=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + c * v
        if r[k] == 0:
            del r[k]
    return r

def scal(a, c):
    return {k: v * c for k, v in a.items() if v * c != 0}

def mul(a, b, nvar, trunc_idx=0, trunc=None):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(ka[i] + kb[i] for i in range(nvar))
            if trunc is not None and k[trunc_idx] > trunc:
                continue
            r[k] = r.get(k, 0) + va * vb
    return {k: v for k, v in r.items() if v != 0}

def one(nvar):
    return {(0,) * nvar: Fr(1)}

def exp_series(a, nvar, trunc_idx, trunc, maxterms):
    # a has no constant term, positive order in variable trunc_idx (or nilpotent under truncation)
    r = one(nvar)
    term = one(nvar)
    for n in range(1, maxterms + 1):
        term = scal(mul(term, a, nvar, trunc_idx, trunc), Fr(1, n))
        if not term:
            break
        r = add(r, term)
    return r

# ---------------- univariate u-series helpers (lists) ----------------
def useries_mul(a, b):
    r = [Fr(0)] * (NU + 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if i + j > NU:
                break
            r[i + j] += x * y
    return r

def useries_pow_binom(a, alpha):
    # a = 1 + w, w = O(u); (1+w)^alpha
    w = [Fr(0)] + a[1:]
    r = [Fr(0)] * (NU + 1); r[0] = Fr(1)
    term = [Fr(0)] * (NU + 1); term[0] = Fr(1)
    for n in range(1, NU + 1):
        term = useries_mul(term, w)
        coef = Fr(1)
        for i in range(n):
            coef *= (alpha - i)
        coef /= factorial(n)
        r = [r[i] + coef * term[i] for i in range(NU + 1)]
    return r

def useries_inv(a):
    return useries_pow_binom(a, -1)

tau_inv = [Fr(4, 5 * n + 4) for n in range(NU + 1)]
tau = useries_inv(tau_inv)
V_expected = useries_pow_binom(tau_inv, Fr(5, 4))      # tau^{-5/4}

# ---------------- lem:con-cfu: solve for y(T) and V ----------------
# polys in (u, T): keys (i, d)
NV2 = 2
def RHS_cfu():
    r = {}
    for m in range(1, NU + 2):
        if m - 1 > NU:
            break
        r = add(r, {(m - 1, 0): Fr(1, m)})
        r = add(r, {(m - 1, m): Fr(1, m * (5 * m - 1))})
    return r

def lhs_cfu(y, V):
    # V/y + T y^4/4, truncated at u^NU; y, V in (u,T) polys
    # 1/y: y = 1 + w
    w = add(y, one(2), -1)
    inv = one(2); term = one(2)
    for n in range(1, NU + 1):
        term = scal(mul(term, w, 2, 0, NU), -1)
        inv = add(inv, term)
    y2 = mul(y, y, 2, 0, NU); y4 = mul(y2, y2, 2, 0, NU)
    return add(mul(V, inv, 2, 0, NU), mul({(0, 1): Fr(1, 4)}, y4, 2, 0, NU))

y = one(2)
V = one(2)
Yk = {}
Vk = {0: Fr(1)}
rhs = RHS_cfu()
for k in range(1, NU + 1):
    l = lhs_cfu(y, V)
    Rk = {}
    for (i, d), c in add(rhs, l, -1).items():
        if i == k:
            Rk[d] = Rk.get(d, 0) + c
        elif i < k and c != 0:
            raise Exception("lower order mismatch at u^%d T^%d: %s" % (i, d, c))
    R1 = sum(Rk.values())
    Vk[k] = R1
    # Y_k = (R_k(T) - R_k(1))/(T-1): synthetic division
    deg = max(Rk.keys()) if Rk else 0
    coeffs = [Rk.get(d, Fr(0)) for d in range(deg + 1)]
    coeffs[0] -= R1
    # divide by (T - 1)
    q = [Fr(0)] * deg
    rem = Fr(0)
    acc = Fr(0)
    for d in range(deg, 0, -1):
        acc = coeffs[d] + acc
        q[d - 1] = acc
    # check remainder
    if sum(coeffs) != 0:
        raise Exception("division remainder nonzero")
    Yk[k] = q
    for d, c in enumerate(q):
        if c != 0:
            y = add(y, {(k, d): c})
    V = add(V, {(k, 0): R1})

print("lem:con-cfu checks (NU=%d)" % NU)
print("  Y_1 =", Yk[1], " (expected [1/18, 1/18])")
print("  deg Y_k <= k:", all(len(Yk[k]) - 1 <= k for k in Yk), [len(Yk[k]) - 1 for k in Yk])
print("  V == tau^{-5/4}:", [Vk[k] for k in range(NU + 1)] == V_expected)
# Y(1)^5 == V: Y(1) = y(T=1)
y1 = [Fr(0)] * (NU + 1)
for (i, d), c in y.items():
    y1[i] += c
y1_5 = useries_mul(useries_mul(useries_mul(useries_mul(y1, y1), y1), y1), y1)
print("  Y(1)^5 == V:", y1_5 == V_expected)
print("  Y(1) == tau^{-1/4}:", y1 == useries_pow_binom(tau_inv, Fr(1, 4)))

# ---------------- inverse map X(Y) = Y g(Y^5) ----------------
def compose_y(gpoly):
    # y(T * g^5) where y in (u,T), g in (u,T); returns (u,T) poly
    g2 = mul(gpoly, gpoly, 2, 0, NU); g4 = mul(g2, g2, 2, 0, NU); g5 = mul(g4, gpoly, 2, 0, NU)
    Tg5 = mul({(0, 1): Fr(1)}, g5, 2, 0, NU)
    # powers of Tg5
    maxd = max(d for (i, d) in y.keys())
    pw = [one(2)]
    for d in range(1, maxd + 1):
        pw.append(mul(pw[-1], Tg5, 2, 0, NU))
    r = {}
    for (i, d), c in y.items():
        r = add(r, mul({(i, 0): c}, pw[d], 2, 0, NU))
    return r

def inv_uT(a):
    w = add(a, one(2), -1)
    inv = one(2); term = one(2)
    for n in range(1, NU + 1):
        term = scal(mul(term, w, 2, 0, NU), -1)
        inv = add(inv, term)
    return inv

g = one(2)
for it in range(NU + 1):
    g = inv_uT(compose_y(g))
# check: g(T) y(T g^5) == 1
chk = mul(g, compose_y(g), 2, 0, NU)
print("  inverse map consistent:", chk == one(2))
print("  deg X_k <= k:", all(d <= i for (i, d) in g.keys()))

# ---------------- hatted amplitude ahat(X; zh) in (u, X, zh) ----------------
B2k = {1: Fr(1, 6), 2: Fr(-1, 30), 3: Fr(1, 42), 4: Fr(-1, 30), 5: Fr(5, 66), 6: Fr(-691, 2730), 7: Fr(7, 6)}
def bk(k):
    return B2k[k] / (2 * k * (2 * k - 1))

loga = {}
for n in range(1, NU + 1):
    loga = add(loga, {(n, 5 * n, 0): Fr(1, 2 * n)})
for k in range(1, NU + 1):
    if 2 * k - 1 > NU:
        break
    base = (2 * k - 1, 2 * k - 1, 2 * k - 1)
    loga = add(loga, {base: bk(k) * Fr(5) ** (1 - 2 * k)})
    for P in range(0, NU + 1):
        if 2 * k - 1 + P > NU:
            break
        loga = add(loga, {(2 * k - 1 + P, 2 * k - 1 + 5 * P, 2 * k - 1): -5 * bk(k) * comb(2 * k - 2 + 5 * P, 5 * P)})
ahat = exp_series(loga, 3, 0, NU, NU + 1)

# ---------------- fa_Y(Y; zh) = ahat(X(Y)) X'(Y), variables (u, Y, zh) ----------------
def uT_to_uY(p):
    return {(i, 5 * d, 0): c for (i, d), c in p.items()}
gY = uT_to_uY(g)
XofY = mul({(0, 1, 0): Fr(1)}, gY, 3, 0, NU)
# X'(Y) = d/dY (Y g(Y^5))
dXdY = {}
for (i, a, j), c in XofY.items():
    if a != 0:
        dXdY = add(dXdY, {(i, a - 1, j): c * a})
maxa = max(a for (i, a, j) in ahat.keys())
Xpow = [one(3)]
for a in range(1, maxa + 1):
    Xpow.append(mul(Xpow[-1], XofY, 3, 0, NU))
faY = {}
for (i, a, j), c in ahat.items():
    faY = add(faY, mul({(i, 0, j): c}, Xpow[a], 3, 0, NU))
faY = mul(faY, dXdY, 3, 0, NU)
print("  fa_Y has only Y^a zh^b with a = b mod 5:", all((a - j) % 5 == 0 for (i, a, j) in faY.keys()))

# ---------------- reduction modulo exact forms ----------------
Vpoly = {(k, 0, 0): Vk[k] for k in range(NU + 1) if Vk[k] != 0}
red = dict(faY)
while True:
    big = [k for k in red if k[1] >= 4]
    if not big:
        break
    nmax = max(k[1] for k in big)
    for key in [k for k in big if k[1] == nmax]:
        c = red.pop(key)
        i, n, j = key
        red = add(red, {(i, n - 4, j + 1): c * (n - 3)})
        red = add(red, mul({(i, n - 5, j): c}, Vpoly, 3, 0, NU))
pi = {m: {} for m in range(-1, 4)}
for (i, n, j), c in red.items():
    assert -1 <= n <= 3
    pi[n][(i, j)] = pi[n].get((i, j), 0) + c
print("  pi_m grading k = m mod 5:", all((j - m) % 5 == 0 for m in pi for (i, j) in pi[m]))
print("  pi_0 at u^0:", {k: v for k, v in pi[0].items() if k[0] == 0})

# ---------------- T-system: Fs_b(t) up to t^K ----------------
KT = 4
def T_op(k, F):
    # T_k F = -F + (t/5)(k F - 4 t F'), F list of coefficients length KT+2
    n = len(F)
    r = [-F[i] for i in range(n)]
    for i in range(n - 1):
        r[i + 1] += Fr(k - 4 * i, 5) * F[i]
    return r
def closing(F0):
    F1 = T_op(-1, F0)
    F2 = [-x for x in T_op(-2, F1)]
    F3 = T_op(2, F2)
    F4 = [-x for x in T_op(1, F3)]
    F0b = [-x for x in T_op(0, F4)]
    return [F0, F1, F2, F3, F4], F0b
L = KT + 2
F0 = [Fr(0)] * L; F0[0] = Fr(1)
for n in range(1, KT + 1):
    # solve for F0[n] from t^{n+1} coefficient of closing relation (linear in F0[n])
    F0[n] = Fr(0)
    _, r0 = closing(F0)
    res0 = r0[n + 1] - F0[n + 1]
    F0[n] = Fr(1)
    _, r1 = closing(F0)
    res1 = r1[n + 1] - F0[n + 1]
    F0[n] = -res0 / (res1 - res0)
Fs, F0b = closing(F0)
print("T-system: Fs_0 =", F0[:4], " closing ok up to t^%d:" % KT, all(F0b[i] == F0[i] for i in range(KT + 1)))
eps = [1, -1, -1, 1, 1]
nu = [0, -1, 3, 2, 1]
print("  Fs_b(0) =", [Fs[b][0] for b in range(5)])

# ---------------- generic 1-variable Wick SPn with u-series coefficients ----------------
# phase F(c+t) = F(c) + H t^2/2 + F3(t); C = -1/H.  Data: polys in (u, t) [keys (i, d)], C as u-series list
def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r
def spn(F3, amp_by_z, Cser, K, nu_trunc):
    # F3: dict (i,d) -> c, d>=3 ; amp_by_z: list over zh-power m of dict (i,d)->c ; returns list over k of u-series
    # moments: <t^{2j}> = (2j-1)!! C^j
    Cpow = [[Fr(1)] + [Fr(0)] * nu_trunc]
    def umul(a, b):
        r = [Fr(0)] * (nu_trunc + 1)
        for i, x in enumerate(a):
            if x == 0: continue
            for j2, y2 in enumerate(b):
                if i + j2 > nu_trunc: break
                r[i + j2] += x * y2
        return r
    for j in range(1, 3 * K + 2):
        Cpow.append(umul(Cpow[-1], Cser))
    F3pow = [one(2)]
    for r_ in range(1, 2 * K + 1):
        F3pow.append(mul(F3pow[-1], F3, 2, 0, nu_trunc))
    out = []
    for k in range(K + 1):
        acc = [Fr(0)] * (nu_trunc + 1)
        for m in range(0, k + 1):
            if m >= len(amp_by_z): break
            for r_ in range(0, 2 * (k - m) + 1):
                prod = mul(F3pow[r_], amp_by_z[m], 2, 0, nu_trunc)
                deg = 2 * (k - m + r_)
                for (i, d), c in prod.items():
                    if d == deg:
                        j = d // 2
                        mom = Fr(dfact(2 * j - 1))
                        for ii in range(nu_trunc + 1 - i):
                            acc[i + ii] += c * mom * Cpow[j][ii] / factorial(r_)
        out.append(acc)
    return out

# lem:con-S2 check: F = -W_1(1+t), amplitude (1+t)^{nu_b}; no u (use nu_trunc = 0)
TD = 3 * KT + 4
W1 = {}
for d in range(0, TD + 1):
    c = Fr(comb(4, d), 4) if d <= 4 else Fr(0)
    c += Fr((-1) ** d)
    if d >= 3 and c != 0:
        W1[(0, d)] = -c
H1 = -(Fr(comb(4, 2), 4) + 1)  # -W1''(1)/1 ... F''/2 coefficient: F = -W1 => t^2 coeff = -(6/4+1) = -5/2, H = -5
C1 = [Fr(1, 5)]
def binom_gen(alpha, d):
    r = Fr(1)
    for i in range(d):
        r *= (alpha - i)
    return r / factorial(d)
ok_S2 = True
for b in range(5):
    amp = [{(0, d): binom_gen(nu[b], d) for d in range(0, TD + 1) if binom_gen(nu[b], d) != 0}]
    res = spn(W1, amp, C1, KT, 0)
    lhs = [res[k][0] for k in range(KT + 1)]
    rhs = [eps[b] * Fs[b][k] for k in range(KT + 1)]
    if lhs != rhs:
        ok_S2 = False
        print("  b=%d mismatch:" % b, lhs, rhs)
print("lem:con-S2 (SPn_1[-W_1; p^nu_b] = eps_b Fs_b) up to t^%d:" % KT, ok_S2)

# ---------------- PhiSP by Wick at X = 1 ----------------
KP = 3
TDP = 6 * KP + 2
# phase F = -hatW(1+t): w_k = sum_m u^{m-1} [(-1)^k + binom(5m-1,k)/(5m-1)]/m
F3h = {}
for d in range(3, TDP + 1):
    for m in range(1, NU + 2):
        c = (Fr((-1) ** d) + Fr(comb(5 * m - 1, d), 5 * m - 1)) / m
        if c != 0:
            F3h = add(F3h, {(m - 1, d): -c})
Ch = [Fr(1, 5), Fr(-1, 5)] + [Fr(0)] * (NU - 1)   # C = (1-u)/5
# amplitude ahat(1+t; zh)/ahat(1;0), split by zh power
inv_a10 = useries_pow_binom([Fr(1), Fr(-1)] + [Fr(0)] * (NU - 1), Fr(1, 2))  # (1-u)^{1/2}
amp_by_z = []
for jz in range(KP + 1):
    p = {}
    for (i, a, j), c in ahat.items():
        if j != jz: continue
        # X^a = (1+t)^a
        for d in range(0, min(a, TDP) + 1):
            p = add(p, {(i, d): c * comb(a, d)})
    p = mul(p, {(i, 0): inv_a10[i] for i in range(NU + 1) if inv_a10[i] != 0}, 2, 0, NU)
    amp_by_z.append(p)
PhiSP = spn(F3h, amp_by_z, Ch, KP, NU)
print("PhiSP coefficients ph_n(u) (u^0..u^%d):" % NU)
for n in range(KP + 1):
    print("  ph_%d =" % n, PhiSP[n])

# RHS of lem:con-id(b): sum_b eps_b pi_{nu_b} tau^{(1-nu_b)/4} Fs_b(tau zh)
def useries_from(pdict_u):
    r = [Fr(0)] * (NU + 1)
    for i, c in pdict_u.items():
        if i <= NU: r[i] += c
    return r
rhs = [[Fr(0)] * (NU + 1) for _ in range(KP + 1)]
taupow = {}
for b in range(5):
    tp = useries_pow_binom(tau_inv, Fr(-(1 - nu[b]), 4))  # tau^{(1-nu)/4} = (tau^{-1})^{-(1-nu)/4}
    pim = pi[nu[b]]
    for n in range(KP + 1):
        # coefficient of zh^n: sum_{j + k = n} [zh^j]pi * Fs_b[k] tau^k
        acc = [Fr(0)] * (NU + 1)
        for k in range(0, n + 1):
            j = n - k
            pj = useries_from({i: c for (i, jj), c in pim.items() if jj == j})
            tk = [Fr(1)] + [Fr(0)] * NU
            for _ in range(k): tk = useries_mul(tk, tau)
            term = useries_mul(useries_mul(pj, tk), tp)
            acc = [acc[i] + eps[b] * Fs[b][k] * term[i] for i in range(NU + 1)]
        rhs[n] = [rhs[n][i] + acc[i] for i in range(NU + 1)]
print("lem:con-id(b) check up to zh^%d, u^%d:" % (KP, NU), all(rhs[n] == PhiSP[n] for n in range(KP + 1)))
for n in range(KP + 1):
    if rhs[n] != PhiSP[n]:
        print("  mismatch n=%d:" % n, rhs[n], PhiSP[n])
# polynomiality / HPF facts: ph_1 = -3/20
print("  ph_1 == -3/20 const:", PhiSP[1] == [Fr(-3, 20)] + [Fr(0)] * NU)

# ---------------- residue: cV (def:con-PhiSP) vs Sigma (lem:con-key) ----------------
EM = int(sys.argv[2]) if len(sys.argv) > 2 else 3
# flat amplitude and exp(-G/z) fa in variables (x, z) [z may be negative]; need x-degree <= EM-1
XD = EM
logfa = {}
for n in range(1, XD // 5 + 2):
    if 5 * n <= XD: logfa = add(logfa, {(5 * n, 0): Fr(1, 2 * n)})
for k in range(1, XD + 1):
    if 2 * k - 1 > XD: break
    logfa = add(logfa, {(2 * k - 1, 2 * k - 1): bk(k) * Fr(5) ** (1 - 2 * k)})
    for P in range(0, XD + 1):
        if 2 * k - 1 + 5 * P > XD: break
        logfa = add(logfa, {(2 * k - 1 + 5 * P, 2 * k - 1): -5 * bk(k) * comb(2 * k - 2 + 5 * P, 5 * P)})
for m in range(1, XD + 1):
    if 5 * m - 1 > XD: break
    logfa = add(logfa, {(5 * m - 1, -1): Fr(-1, m * (5 * m - 1))})
E1 = exp_series(logfa, 2, 0, XD, XD + 1)
Xi = [Fr(0)] + [Fr(1, n) for n in range(1, EM + 1)]
cV = {}
Xik = [Fr(1)] + [Fr(0)] * EM
for k in range(1, EM + 1):
    Xik = [sum(Xik[i] * Xi[j - i] for i in range(j + 1)) for j in range(EM + 1)]
    coef = Fr((-1) ** k, factorial(k))
    for (a, zz), c in E1.items():
        if a != k - 1: continue
        for E in range(EM + 1):
            if Xik[E] != 0:
                key = (E, zz - k)
                cV[key] = cV.get(key, 0) + coef * Xik[E] * c
cV = {k: v for k, v in cV.items() if v != 0}

# Sigma
Sig = {}
for b in range(5):
    m = nu[b]
    for (i, j), c in pi[m].items():
        ip = 0
        while True:
            e = m + 4 * ip + 1
            beta = j - 5 * ip - m - 1
            # V^e expansion
            Ve = useries_pow_binom(Vk_list := [Vk[t] for t in range(NU + 1)], e) if e > 0 else [Fr(1)] + [Fr(0)] * NU
            base = c * Fr((-1) ** ip, 4 ** ip * factorial(ip)) * Fr((-1) ** e, factorial(e))
            num = 4 * beta - 1
            assert num % 5 == 0
            for l in range(NU + 1):
                if Ve[l] == 0: continue
                E = i + l - num // 5
                if E > EM: continue
                Sig[(E, beta)] = Sig.get((E, beta), 0) + base * Ve[l]
            ip += 1
            # stop when E lower bound exceeds EM: E >= i - (4*beta-1)/5 increases by 4 per ip
            if i - (4 * (j - 5 * ip - m - 1) - 1) // 5 > EM:
                break
# completeness: coefficient (E,beta) complete iff E + (4beta-1)/5 <= NU
ok = True; ncheck = 0
keys = set(cV) | set(Sig)
for (E, beta) in sorted(keys):
    if E < 0 or E > EM: 
        if Sig.get((E, beta), 0) != 0 and E + (4 * beta - 1) // 5 <= NU:
            print("  nonzero Sigma at E<0:", E, beta, Sig[(E, beta)]); ok = False
        continue
    if E + (4 * beta - 1) // 5 > NU: continue
    ncheck += 1
    if cV.get((E, beta), 0) != Sig.get((E, beta), 0):
        ok = False
        print("  mismatch at u^%d z^%d: cV=%s Sigma=%s" % (E, beta, cV.get((E, beta), 0), Sig.get((E, beta), 0)))
print("prop:con-res first half (cV = Sigma), %d complete coefficients with E<=%d checked:" % (ncheck, EM), ok)
print("  cV:", sorted(cV.items())[:8])
checked = []
for E in range(0, EM + 1):
    for beta in range(-5 * EM - 6, 30):
        if (4 * beta - 1) % 5 != 0: continue
        if E + (4 * beta - 1) // 5 > NU: continue
        if (E, beta) in cV or (E, beta) in Sig:
            checked.append(((E, beta), cV.get((E, beta), 0), Sig.get((E, beta), 0)))
print("  nonzero complete coefficients:", checked)
