# -*- coding: utf-8 -*-
"""rv_compose.py -- independent check of the composition law used by fq4_split.Split at the genus-4 parameters (KZ = 9, K = 6).

From the hatted rows Q (rebuilt with the same project calls as Split.__init__) and an INDEPENDENT construction of
M^-1(zhat; u) = M_0^-1(tau(u) zhat)  (tau from the hypergeometric series by the Pochhammer ratio, M_0 = colour matrices of the symbol),
with own series arithmetic and own (recursive) division by (psi + psi'):
  (1) M_0 M_0^-1 = 1 and M^-1(z) M^-1(-z)^T = 1 (symplectic) mod zhat^(KZ+1), all u-orders < K;
  (2) E_R(a, b) recomputed from Q equals Split.ER;  E_M recomputed from M^-1 equals Split.EM;
  (3) R' := R M^-1 (rows Q' = Q M^-1);  E_R' recomputed from Q';  M^-1 E_R' M^-T (transport) equals Split.EO = E_R - E_M;
  (4) leaves: T_R = T_M + M^-1 T_R' with T_M = -M^-1_(d-1) 1,  1_j = Qhat_(0,0)^(j);  T_R equals Split.TR;
  (5) [u^0] E_R' = 0, [u^0] T_R' = 0 (R' = Psi(u)(1 + O(u)));  val_u of every E_out entry >= 1;
  (6) natural-form support: G = Psi(u)^-1 R' has [u^k zhat^n] G = 0 for n > floor((5k+2)/4), k = 1..K-1, n <= KZ (NATURAL_FORM.md s.1.1);
      the same test for the mutated tau of the mutants m04 / m05 (must FAIL: the support rule sees tau_4, tau_5).
Usage: py -I rv_compose.py [89|61] [alpha,beta]     -> out/compose_<prime>[_ab].log
"""
import os, sys, time
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq_fp import FpRing, P1, P2

prime = sys.argv[1] if len(sys.argv) > 1 else "89"; abs_ = sys.argv[2] if len(sys.argv) > 2 else "0,0"
al, be = [Fr(x) for x in abs_.split(",")]
g = 4; KZ = 3 * g - 3; K = 2 * g - 2
p = P2 if prime == "89" else P1; R = FpRing(p)
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "compose_%s%s.log" % (prime, "" if abs_ == "0,0" else "_ab")), "w", encoding="utf-8")
FAILED = []
def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

# ---------------------------------------------------------------- own series arithmetic: lists of K residues
Z = [0] * K; ONE = [1] + [0] * (K - 1)
def smul(a, b):
    out = [0] * K
    for i in range(K):
        if a[i]:
            for j in range(K - i):
                if b[j]: out[i + j] += a[i] * b[j]
    return [c % p for c in out]
def sadd(a, b): return [(x + y) % p for x, y in zip(a, b)]
def ssub(a, b): return [(x - y) % p for x, y in zip(a, b)]
def sneg(a): return [(-x) % p for x in a]
def norm(a): return Z if a is None else list(a)

t0 = time.time()
from fq_rows import Rows
from fq_conifold import Conifold
from fq_sser import SRows
from fq_symbol import Symbol
from fq3_msum import HatRowsU
rows = Rows(R, Nz=KZ); C = Conifold(R, Nu=K + 4, alpha=al, beta=be); S = SRows(rows, C, prec_abs=5 * K + 4)
sym = Symbol(S); H = HatRowsU(S, K)
Q = [[[[H.Q[a][i][n].coef(k) % p for k in range(K)] for n in range(KZ + 1)] for i in range(5)] for a in range(5)]
eta = dict(S.eta_inv)
log("rv_compose: prime 2^%s-1, frame (%s, %s), KZ = %d, K = %d; rows built in %.0f s" % (prime, al, be, KZ, K, time.time() - t0))

# ---------------------------------------------------------------- tau, independently: 2F1(1, 4/5; 9/5; u) by the Pochhammer ratio
def tau_series(shift=None):
    f = [Fr(1)]
    for k in range(1, K): f.append(f[-1] * Fr(k) * (Fr(4, 5) + k - 1) / ((Fr(9, 5) + k - 1) * k))     # (1)_k (4/5)_k / ((9/5)_k k!)
    if shift is not None: f[shift] = Fr(4, 5 * shift + 5)                                             # the mutants m04 / m05
    t = [Fr(1)]
    for k in range(1, K): t.append(-sum(f[j] * t[k - j] for j in range(1, k + 1)))
    return [R.of(x) for x in t], t

def build_Minv(tau):
    Cm = sym.colour_matrix(); tp = [ONE]
    for k in range(KZ): tp.append(smul(tp[-1], tau))
    Minv = [[[[x * (((-1) ** k) * Cm[k][j][i]) % p for x in tp[k]] for j in range(5)] for i in range(5)] for k in range(KZ + 1)]
    Mfwd = [[[[x * Cm[k][i][j] % p for x in tp[k]] for j in range(5)] for i in range(5)] for k in range(KZ + 1)]
    return Minv, Mfwd

def bil(X, Y, sign):
    """B[(M, N)][i][j] = delta - sign^(M+N) Sum_ab eta^ab X[a][i][M] Y[b][j][N]   (X, Y rows indexed [a][i][n])"""
    B = {}
    for M_ in range(KZ + 1):
        for N_ in range(KZ + 1 - M_):
            tab = [[None] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    acc = [0] * K
                    for (a, b), e in eta.items(): acc = sadd(acc, [x * e % p for x in smul(X[a][i][M_], Y[b][j][N_])])
                    if (M_ + N_) % 2 and sign < 0: acc = sneg(acc)
                    acc = sneg(acc)
                    if i == j and M_ == 0 and N_ == 0: acc = sadd(acc, ONE)
                    tab[i][j] = acc
            B[(M_, N_)] = tab
    return B

def divide(B):
    """E(m, n) with B(M, N) = E(M-1, N) + E(M, N-1): E(m, n) = B(m+1, n) - E(m+1, n-1); returns (E, remainder test B(0, N) = E(0, N-1))"""
    E = {}
    for tot in range(KZ):
        for n in range(tot + 1):
            m = tot - n
            E[(m, n)] = B[(m + 1, n)] if n == 0 else [[ssub(B[(m + 1, n)][i][j], E[(m + 1, n - 1)][i][j]) for j in range(5)] for i in range(5)]
    ok = all(B[(0, N_)][i][j] == (E[(0, N_ - 1)][i][j] if N_ else Z) for N_ in range(KZ + 1) for i in range(5) for j in range(5))
    return E, ok

def bilM(Minv):
    B = {}
    for M_ in range(KZ + 1):
        for N_ in range(KZ + 1 - M_):
            tab = [[None] * 5 for _ in range(5)]
            for i in range(5):
                for j in range(5):
                    acc = [0] * K
                    for k in range(5): acc = sadd(acc, smul(Minv[M_][i][k], Minv[N_][j][k]))
                    acc = sneg(acc)
                    if i == j and M_ == 0 and N_ == 0: acc = sadd(acc, ONE)
                    tab[i][j] = acc
            B[(M_, N_)] = tab
    return B

def primed_rows(Minv):
    """Q'[a][i][n] = Sum_{k <= n} Sum_j Q[a][j][n-k] Minv_k[j][i]"""
    Qp = [[[None] * (KZ + 1) for i in range(5)] for a in range(5)]
    for a in range(5):
        for i in range(5):
            for n in range(KZ + 1):
                acc = [0] * K
                for k in range(n + 1):
                    for j in range(5): acc = sadd(acc, smul(Q[a][j][n - k], Minv[k][j][i]))
                Qp[a][i][n] = acc
    return Qp

def support_violations(Qp):
    """G[n][j][i] = Sum_ab Psi^-1 ... = Sum_ab Q'[b][j][0] eta^ba Q'[a][i][n] (colour-colour form); list of (k, n) with [u^k] G_n != 0, n > floor((5k+2)/4)"""
    bad = set(); nz = {}
    for n in range(KZ + 1):
        for j in range(5):
            for i in range(5):
                acc = [0] * K
                for (a, b), e in eta.items(): acc = sadd(acc, [x * e % p for x in smul(Qp[b][j][0], Qp[a][i][n])])
                for k in range(K):
                    if acc[k]:
                        nz[k] = max(nz.get(k, 0), n)
                        if n > (5 * k + 2) // 4 or (k == 0 and n > 0): bad.add((k, n))
    return sorted(bad), nz

tau, tauQ = tau_series()
log("   tau = %s" % ", ".join(str(x) for x in tauQ))
from fq4_split import Split, tau_pred
check("tau: own Pochhammer series = fq4_split.tau_pred", tauQ == tau_pred(K))
Minv, Mfwd = build_Minv(tau)
# (1) inverse and symplecticity
ok1 = ok2 = True
for k in range(KZ + 1):
    for i in range(5):
        for j in range(5):
            a1 = [0] * K; a2 = [0] * K
            for m in range(k + 1):
                for l in range(5):
                    a1 = sadd(a1, smul(Mfwd[m][i][l], Minv[k - m][l][j]))
                    t = smul(Minv[m][i][l], Minv[k - m][j][l]); a2 = sadd(a2, t if (k - m) % 2 == 0 else sneg(t))
            tgt = ONE if (k == 0 and i == j) else Z
            ok1 = ok1 and a1 == tgt; ok2 = ok2 and a2 == tgt
check("(1) M(z) M^-1(z) = 1 mod zhat^%d, all u-orders < %d" % (KZ + 1, K), ok1)
check("(1) M^-1(z) M^-1(-z)^T = 1 (M symplectic) mod zhat^%d" % (KZ + 1), ok2)

sp = Split(p, g, (al, be), K, 'nat')
keys = sorted(sp.ER)
# (2) E_R and E_M
ER, okdiv = divide(bil(Q, Q, -1))
check("(2) symplectic remainder of R: B(0, N) = E_R(0, N-1)", okdiv)
check("(2) E_R recomputed from the rows = Split.ER (%d keys x 25 colours x %d u-orders)" % (len(keys), K),
      all(ER[key][i][j] == norm(sp.ER[key][i][j]) for key in keys for i in range(5) for j in range(5)))
EM, okdivM = divide(bilM(Minv))
check("(2) symplectic remainder of M", okdivM)
check("(2) E_M recomputed from M^-1 = M_0^-1(tau zhat) = Split.EM", all(EM[key][i][j] == norm(sp.EM[key][i][j]) for key in keys for i in range(5) for j in range(5)))
# (3) R' = R M^-1 and the transport
Qp = primed_rows(Minv); ERp, okdivp = divide(bil(Qp, Qp, -1))
check("(3) symplectic remainder of R' = R M^-1", okdivp)
bad = []
for (a, b) in keys:
    for i in range(5):
        for j in range(5):
            acc = [0] * K
            for k1 in range(a + 1):
                for k2 in range(b + 1):
                    Ep = ERp[(a - k1, b - k2)]
                    for i2 in range(5):
                        for j2 in range(5):
                            acc = sadd(acc, smul(smul(Minv[k1][i][i2], Ep[i2][j2]), Minv[k2][j][j2]))
            if acc != norm(sp.EO[(a, b)][i][j]): bad.append((a, b, i, j))
check("(3) M^-1 E_R' M^-T = Split.EO = E_R - E_M (all keys, colours, u-orders)", not bad, bad[:5])
check("(5) [u^0] E_R' = 0 for all keys", all(ERp[key][i][j][0] == 0 for key in keys for i in range(5) for j in range(5)))
check("(5) every entry of Split.EO has u-valuation >= 1 (needed for the restriction popcount < K)", min(sp.valO.values()) >= 1, sorted(set(sp.valO.values())))
check("(5) Split.frozen_ok", sp.frozen_ok is True)
# (4) leaves
unit = [Q[0][j][0] for j in range(5)]; okT = okTR = okT0 = True
TRp = {d: [[x * ((-1) ** d) % p for x in Qp[0][j][d - 1]] for j in range(5)] for d in range(1, KZ + 2)}
check("(4) T_R'(1) = 0 (R' has the same unit)", all(sadd(TRp[1][j], unit[j]) == Z for j in range(5)))
for d in range(2, KZ + 2):
    for i in range(5):
        TR = [x * ((-1) ** d) % p for x in Q[0][i][d - 1]]
        okTR = okTR and TR == norm(sp.TR[d][i])
        acc = [0] * K
        for j in range(5): acc = ssub(acc, smul(Minv[d - 1][i][j], unit[j]))                        # T_M(d)
        for dd in range(2, d + 1):
            for j in range(5): acc = sadd(acc, smul(Minv[d - dd][i][j], TRp[dd][j]))                # M^-1 T_R'
        okT = okT and acc == TR; okT0 = okT0 and TRp[d][i][0] == 0
check("(4) T_R recomputed = Split.TR", okTR)
check("(4) T_R = T_M + M^-1 T_R' (d = 2..%d)" % (KZ + 1), okT)
check("(5) [u^0] T_R' = 0", okT0)
# (6) natural-form support rule, and its sensitivity to tau_4, tau_5
badS, nz = support_violations(Qp)
check("(6) support of G = Psi^-1 R': [u^k zhat^n] G = 0 for n > floor((5k+2)/4), k <= %d, n <= %d" % (K - 1, KZ), not badS, badS)
log("      highest non-zero zhat-order of [u^k] G: %s" % {k: nz[k] for k in sorted(nz)})
for sh in (4, 5):
    tm, _ = tau_series(sh); Mi, _ = build_Minv(tm); badm, nzm = support_violations(primed_rows(Mi))
    check("(6) control: the support rule FAILS for the mutated tau (1/tau coefficient of u^%d -> 4/%d)" % (sh, 5 * sh + 5), bool(badm), badm[:8])
log("RV_COMPOSE (2^%s-1, frame %s) %s  (%.0f s)" % (prime, abs_, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
