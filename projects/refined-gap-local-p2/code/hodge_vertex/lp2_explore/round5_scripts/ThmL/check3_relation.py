# Check 3: the dilaton-point relation M[k0,n,R] of the note, in the uniform form (parts >= 1, no empty-profile table),
# with values X_g(U) from the MV side, exact at rational tau.
#   Sum_alpha c_alpha Sum_{pi} Sum_{assign R} Sum_{(g_B,k_B)} prod_B phi(B,T_B,g_B,k_B) = 0,  budget Sum e_B = n, e_B = 2g-2+|B|+k
#   phi = c(|B|,k) prod_{i in B}(i Pi_i) sigma^{|B|+g-1} tau^k/k! (delta_g(U))_k X_g(U),  U = T_B + B
#   c: |B| odd: -1 (k>=1); |B| even: 2 (k=0), 1 (k>=1)
# Controls: perturb target; beta -> -beta (tau_{-1}); rising instead of falling; wrong c_(k0); wrong sign of even blocks at k=0
import sys, itertools
from fractions import Fraction as Fr
from math import factorial
from tl_lib import *
from mvx import *

tau = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(2, 7)
gmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
maxdeg = int(sys.argv[3]) if len(sys.argv) > 3 else 16
kmaxmem = int(sys.argv[4]) if len(sys.argv) > 4 else 5
mv = MV(tau, L=2 * maxdeg + 2 * gmax + 4)
sigma = mv.sigma

class TooBig(Exception): pass

def Xval(g, U):
    U = tuple(sorted(U, reverse=True))
    if g == 0:
        return x0(U)
    # strip units by string where possible, to keep |mu| small
    core = tuple(u for u in U if u >= 2)
    v = len(U) - len(core)
    if not core:
        if g == 1:
            return factorial(v - 1) * (-(1 + 1 / sigma) / 24)
        return Fr(0)
    N = 3 * g - 3 + len(core)
    if sum(core) + N > maxdeg: raise TooBig((g, core))
    val = mv.X(g, core)
    J = sum(core)
    for i in range(v): val *= (J + i)
    return val

def relation(k0, n, R, mode="plain", target=None):
    P = member_Pk(k0)
    if mode == "wrong_ck":
        P = dict(P); P[(k0,)] = P[(k0,)] + 1
    total = Fr(0); nterms = 0
    Rl = list(R); r = len(Rl)
    for alpha, ca in P.items():
        al = list(alpha)
        for pi in set_partitions(list(range(len(al)))):
            m = len(pi)
            # assignments of labelled R to blocks
            for assign in itertools.product(range(m), repeat=r):
                blocks = []
                for bi, b in enumerate(pi):
                    B = [al[i] for i in b]
                    T = [Rl[j] for j in range(r) if assign[j] == bi]
                    blocks.append((B, T))
                # DP over budget
                dp = {0: Fr(1)}
                for (B, T) in blocks:
                    U = tuple(sorted(T + B, reverse=True))
                    sz = len(B)
                    pref = Fr(1)
                    for i in B: pref *= i * Pi(i, tau)
                    opts = {}
                    for g in range(0, n // 2 + 2):
                        for k in range(0, n + 3):
                            e = 2 * g - 2 + sz + k
                            if e < 0 or e > n: continue
                            if sz % 2 == 1:
                                if k == 0: continue
                                c = Fr(-1)
                                if mode == "beta_flip": c = Fr(-1) * (-1)**k * (-1) if False else Fr(-1)
                            else:
                                c = Fr(2) if k == 0 else Fr(1)
                            if mode == "beta_flip":
                                # tau_1 -> tau_{-1}: beta -> -beta
                                c = c * (-1)**k if k >= 1 else c
                            if mode == "even0" and sz % 2 == 0 and k == 0:
                                c = Fr(1)
                            d = delta(g, U)
                            if mode == "rising":
                                ff = 1
                                for t in range(k): ff *= (d + t)
                            else:
                                ff = falling(d, k)
                            if ff == 0: continue
                            xv = Xval(g, U)
                            if target is not None and g == target[0] and U == target[1]:
                                xv = xv + target[2]
                            if xv == 0: continue
                            val = c * pref * sigma**(sz + g - 1) * tau**k / factorial(k) * ff * xv
                            opts[e] = opts.get(e, 0) + val
                    nd = {}
                    for e0, v0 in dp.items():
                        for e1, v1 in opts.items():
                            if e0 + e1 <= n:
                                nd[e0 + e1] = nd.get(e0 + e1, 0) + v0 * v1
                    dp = nd
                total += ca * dp.get(n, 0)
                nterms += 1
    return total, nterms

def multisets_upto(Jmax, rmax):
    res = [()]
    for r in range(1, rmax + 1):
        for lam in partitions_with(Jmax, r):
            res.append(lam)
    return res
def partitions_with(Jmax, r):
    for J in range(2 * r, Jmax + 1):
        for lam in partitions(J):
            if len(lam) == r and lam[-1] >= 2:
                yield lam

stats = {}
fails = 0
ctrl_ok = {"perturb": [0, 0], "beta_flip": [0, 0], "rising": [0, 0], "wrong_ck": [0, 0], "even0": [0, 0]}
ctrl_tc = [0, 0]
Rs = multisets_upto(int(sys.argv[5]) if len(sys.argv) > 5 else 8, 3)
for k0 in range(2, kmaxmem + 1):
    for G in range(0, gmax + 1):
        for eps in (0, 1):
            n = 2 * G + eps
            for R in Rs:
                try:
                    val, nt = relation(k0, n, R)
                except TooBig:
                    stats[(k0, "skip")] = stats.get((k0, "skip"), 0) + 1
                    continue
                key = (k0, G, eps)
                s = stats.setdefault(key, [0, 0])
                s[0] += 1
                if val != 0:
                    fails += 1
                    print("FAIL", k0, G, eps, R, val)
                # target test & controls only at even n, target coefficient nonzero
                if eps == 0 and G >= 1 and len(R) <= 2:
                    U = tuple(sorted(R + (k0,), reverse=True))
                    d = delta(G, U)
                    if d != 0:
                        try:
                            v2, _ = relation(k0, n, R, target=(G, U, Fr(1)))
                            expected = -(k0 - 1) * Pi(k0, tau) * sigma**G * tau * d
                            ctrl_ok["perturb"][0] += 1
                            # perturbation of X_G(U) changes the relation by the target coefficient plus products with lower blocks
                            ctrl_ok["perturb"][1] += (v2 != 0)
                            if k0 not in R:
                                ctrl_tc[0] += 1
                                ctrl_tc[1] += (v2 - val == expected)
                        except TooBig:
                            pass
                    for mode in ("beta_flip", "rising", "wrong_ck", "even0"):
                        try:
                            v3, _ = relation(k0, n, R, mode=mode)
                            ctrl_ok[mode][0] += 1
                            ctrl_ok[mode][1] += (v3 != 0)
                        except TooBig:
                            pass
print("tau", tau, "gmax", gmax, "maxdeg", maxdeg)
tot = 0
for key in sorted(k for k in stats if k[1] != "skip"):
    print("  k0=%d G=%d eps=%d: relations %d" % (key[0], key[1], key[2], stats[key][0]))
    tot += stats[key][0]
print("total relations", tot, "failures", fails, "skips", {k: v for k, v in stats.items() if k[1] == "skip"})
print("controls (tested, nonzero):", ctrl_ok)
print("target coefficient -(k0-1) Pi_k0 sigma^G tau delta = response of the relation to X_G(R+k0), k0 not in R (tested, equal):", ctrl_tc)
