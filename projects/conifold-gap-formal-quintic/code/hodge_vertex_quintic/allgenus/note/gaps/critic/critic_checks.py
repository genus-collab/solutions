"""critic_checks.py -- gap-phase critic (key gaps/critic): four cheap exact checks, each with failing controls.

Evidence only (VERIFIED, exact over Q or Q[c]/(c^5 - 1)); no proof rests on it. Routes differ from the provers' and referees':
  A  Givental 1996: Theorem 11.6 rewritten for s_i gives the factor f exp{[g_1 Sum lambda_alpha - Sum lambda'_a g_(l_a)]/(hbar f)}
     (Calabi-Yau case); controls: non-CY degrees, and Cor. 11.7 read literally (Sum_alpha l_alpha = n + 1).
  B  trivial frames of lem:tra-frames: by the support rule of Pi_0 alone (combinatorics), ten frames; control: the rule k <= nu_a + 1.
  C  colour invariance of PF5: [z^j] PF5^(u_i) f = c^-j [z^j] PF5^(u_0) f (c = zeta^i formal, c^5 = 1), z^0 part 0, z^1 part on L g
     equal to 5 c^4 D g (the colour-i transport operator); control: PF5 + q z^2.
  D  prop:con-deg (L0)-(L4) by applying P to monomials x^m (not as normal-ordered symbols); controls: P + z^2 x^3, P + z x^6.
Run: py -I -B critic_checks.py   (prints to stdout; no paths are printed)
"""
import random
from fractions import Fraction as Fr

FAILS = []
def report(tag, ok, expect=True):
    good = (ok == expect)
    if not good:
        FAILS.append(tag)
    print(("PASS  " if good else "FAIL  ") + tag + ("" if expect else "   [control: must be False; is %s]" % ok))

# ---------------- A: Givental Thm 11.6 -> corrected Cor. 11.7 -----------------------------------------------------------
def givental(n, ls, rng, literal=False):
    lam = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(n + 1)]
    lamp = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in ls]
    f, hb = Fr(rng.randint(1, 9), rng.randint(1, 5)), Fr(rng.randint(1, 9), rng.randint(1, 5))
    g = {l: Fr(rng.randint(-9, 9), rng.randint(1, 7)) for l in set(ls) | {1}}
    ok = True
    for i in range(n + 1):
        Xi = (sum((l * lam[i] - lp) * g[l] for l, lp in zip(ls, lamp)) - sum(lam[i] - la for la in lam) * g[1]) / (f * hb)
        TmT = sum(l * (g[l] - g[1]) for l in ls) / f                       # T - t (operation (1))
        Yi = Xi - lam[i] * TmT / hb                                         # exponent in s*_i = f e^(Yi) s_i(T)
        s = sum(ls) if literal else sum(lam)                                # literal reading: Sum_alpha l_alpha = Sum_a l_a
        Zi = (g[1] * s - sum(lp * g[l] for l, lp in zip(ls, lamp))) / (hb * f)
        ok = ok and (Yi == Zi)
    return ok

rng = random.Random(20261115)
report("A  CY case n=4, l=(5) (formal quintic): exponent = [g_1 Sum lambda - lambda' g_5]/(hbar f), 20 random instances",
       all(givental(4, [5], rng) for _ in range(20)))
report("A  CY case n=4, l=(2,3): same identity, 20 random instances", all(givental(4, [2, 3], rng) for _ in range(20)))
report("A  control: non-CY n=4, l=(4)", all(givental(4, [4], rng) for _ in range(20)), expect=False)
report("A  control: Cor. 11.7 read literally (Sum_alpha l_alpha)", all(givental(4, [5], rng, literal=True) for _ in range(20)),
       expect=False)

# ---------------- B: trivial frames by the support rule --------------------------------------------------------------------
nu = {0: 0, 1: -1, 2: 3, 3: 2, 4: 1}
def pi0_nonzero(m, b):          # lem:gau-vanshape / FL Step 0: Pi_0(z^m H^b) != 0 iff b + m = 1 mod 5
    return (b + m) % 5 == 1
def trivial(k, a):              # frame(k,a) = z^k H^a iff Pi_i(z^k H^a) = Pi_0(z^(k-i) H^a) = 0 for all 0 <= i <= k
    return all(not pi0_nonzero(k - i, a) for i in range(k + 1))
triv = sorted((k, a) for k in range(60) for a in (0, 2, 3, 4) if trivial(k, a))
print("   trivial frames:", triv)
report("B  trivial frames = {(k,a): k < nu_a + 1}, ten in all", triv == sorted((k, a) for a in (0, 2, 3, 4) for k in range(nu[a] + 1))
       and len(triv) == 10)
report("B  (2,2), (3,2) are trivial (missing from the sources' list of eight)", (2, 2) in triv and (3, 2) in triv)
report("B  control: rule k <= nu_a + 1 (14 frames)", triv == sorted((k, a) for a in (0, 2, 3, 4) for k in range(nu[a] + 2)),
       expect=False)

# ---------------- C: colour invariance of PF5 and the transport operator ----------------------------------------------------
# elements: dict {(zpow, Lpow, cpow): Fr}, cpow mod 5; D L = (L^6 - L)/5, q = (1 - L^-5)/3125, D acts on L-powers only.
def add(*xs):
    out = {}
    for x in xs:
        for k, v in x.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
def scal(x, s):
    return {k: v * s for k, v in x.items() if v * s != 0}
def mulmono(x, zp, Lp, cp, coef=Fr(1)):
    return {(k[0] + zp, k[1] + Lp, (k[2] + cp) % 5): v * coef for k, v in x.items()}
def Dop(x):                     # D on coefficients: D(L^m) = m (L^(m+5) - L^m)/5
    out = {}
    for (zp, Lp, cp), v in x.items():
        if Lp:
            for dl, s in ((5, Fr(Lp, 5)), (0, -Fr(Lp, 5))):
                key = (zp, Lp + dl, cp)
                out[key] = out.get(key, 0) + v * s
    return {k: v for k, v in out.items() if v != 0}
def theta_plus(x, cpow):        # (zD + c^cpow L) x
    return add(mulmono(Dop(x), 1, 0, 0), mulmono(x, 0, 1, cpow))
def qmul(x):                    # q = (1 - L^-5)/3125
    return add(scal(x, Fr(1, 3125)), mulmono(x, 0, -5, 0, Fr(-1, 3125)))
def pf5(f, cpow, extra=False):
    y = f
    for _ in range(5):
        y = theta_plus(y, cpow)
    w = f
    for k in range(1, 6):        # prod_k (5 (zD + cL) + k z), factors commute
        w = add(scal(theta_plus(w, cpow), 5), mulmono(w, 1, 0, 0, Fr(k)))
    out = add(y, scal(f, -1), scal(qmul(w), -1))
    if extra:
        out = add(out, mulmono(qmul(f), 2, 0, 0))
    return out
def zpart(x, j):
    return {k: v for k, v in x.items() if k[0] == j}
def cshift(x, j):               # multiply by c^-j
    return {(k[0], k[1], (k[2] - j) % 5): v for k, v in x.items()}
okC, okT, okctl = True, True, True
for m in range(-3, 6):
    f = {(0, m, 0): Fr(1)}
    base, base_x = pf5(f, 0), pf5(f, 0, extra=True)
    for i in range(1, 5):
        col, col_x = pf5(f, i), pf5(f, i, extra=True)
        for j in range(6):       # c = zeta^i, so c^-j = zeta^(-ij)
            okC = okC and zpart(col, j) == cshift(zpart(base, j), i * j)
        okctl = okctl and all(zpart(col_x, j) == cshift(zpart(base_x, j), i * j) for j in range(6))
    for i in range(5):           # transport operator: f = L g, g = L^(m-1)
        col = pf5(f, i)
        g = {(0, m - 1, 0): Fr(1)}
        okT = okT and not zpart(col, 0) and zpart(col, 1) == mulmono(scal(Dop(g), 5), 1, 0, 4 * i)
report("C  [z^j] PF5^(u_i) L^m = zeta^(-ij) [z^j] PF5^(u_0) L^m, i = 1..4, j = 0..5, m = -3..5", okC)
report("C  [z^0] PF5^(u_i)(L g) = 0 and [z^1] PF5^(u_i)(L g) = 5 zeta^(4i) D g, g = L^(m-1)", okT)
report("C  control: PF5 + q z^2 is colour covariant", okctl, expect=False)

# ---------------- D: prop:con-deg (L0)-(L4) on monomials ----------------------------------------------------------------------
# elements: dict {(zpow, xpow): Fr}; M = x - (z/5)(1 - x^5) theta, theta = x d/dx; P = M^5 - 1 + (x^-5 - 1) Pi,
# Pi = 5^-5 prod_k (5M + kz).
def Mop(y):
    out = {}
    for (zp, xp), v in y.items():
        for key, s in (((zp, xp + 1), v), ((zp + 1, xp), -v * Fr(xp, 5)), ((zp + 1, xp + 5), v * Fr(xp, 5))):
            if s:
                out[key] = out.get(key, 0) + s
    return {k: v for k, v in out.items() if v != 0}
def Pop(m, ctl=None):
    f = {(0, m): Fr(1)}
    y = f
    for _ in range(5):
        y = Mop(y)
    w = f
    for k in range(1, 6):
        w5 = {kk: 5 * v for kk, v in Mop(w).items()}
        for (zp, xp), v in w.items():
            w5[(zp + 1, xp)] = w5.get((zp + 1, xp), 0) + k * v
        w = {kk: v for kk, v in w5.items() if v != 0}
    w = {kk: v / 3125 for kk, v in w.items()}
    out = dict(y)
    out[(0, m)] = out.get((0, m), 0) - 1
    for (zp, xp), v in w.items():
        out[(zp, xp - 5)] = out.get((zp, xp - 5), 0) + v
        out[(zp, xp)] = out.get((zp, xp), 0) - v
    if ctl == "z2x3":
        out[(2, m + 3)] = out.get((2, m + 3), 0) + 1
    if ctl == "zx6":
        out[(1, m + 6)] = out.get((1, m + 6), 0) + 1
    return {k: v for k, v in out.items() if v != 0}
def divisible_1mx5(coeffs):     # Laurent polynomial {xp: v} divisible by 1 - x^5 iff each residue class mod 5 sums to 0
    s = {}
    for xp, v in coeffs.items():
        s[xp % 5] = s.get(xp % 5, 0) + v
    return all(v == 0 for v in s.values())
def lcheck(ctl=None, mrange=range(-6, 26)):
    L0 = L1 = L2 = L3 = L4 = True
    for m in mrange:
        P = Pop(m, ctl)
        parts = {}
        for (zp, xp), v in P.items():
            parts.setdefault(zp, {})[xp] = v
        L0 = L0 and not parts.get(0)
        for j, c in parts.items():
            if j >= 1:
                L1 = L1 and divisible_1mx5(c)
                L4 = L4 and max(c) <= m + 4 * j
                if m >= 1:
                    L2 = L2 and min(c) >= 1
        want = {} if m == 1 else {m - 1: Fr(1 - m), m + 4: Fr(m - 1)}
        L3 = L3 and parts.get(1, {}) == want
    return L0, L1, L2, L3, L4
L0, L1, L2, L3, L4 = lcheck()
report("D  (L0) P_0 x^m = 0, m = -6..25", L0)
report("D  (L1) 1 - x^5 divides P_j x^m, j >= 1", L1)
report("D  (L2) P_j(x^m) in x Q[x] for m >= 1", L2)
report("D  (L3) P_1 x^m = (1 - m)(x^(m-1) - x^(m+4))", L3)
report("D  (L4) deg P_j(x^m) <= m + 4j (so deg RHS_n <= 4n + 5)", L4)
report("D  control P + z^2 x^3: (L1)", lcheck("z2x3")[1], expect=False)
report("D  control P + z x^6: (L4)", lcheck("zx6")[4], expect=False)

print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "FAILED: %s" % FAILS)
