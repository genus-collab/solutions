# -*- coding: utf-8 -*-
"""rf_res.py -- REFEREE: (a) the formal residue at x = 0 of the Mellin-Barnes/Stirling integrand versus C0 (Y^C_Q)_0 (POLY s.5,
Theorem P steps 6-7), in BOTH conifold frames; (b) the CFU map of POLY s.5, exactly over Q.

(a) Original variables x = 1/sigma: integrand exp(-(Lam/x + G(x))/z) at(x,z) dx, Lam = -log(1-u).  Here G is NOT taken from the
closed form: G(x) := Sum_j (1/x - zeta^j) log(1 - zeta^j x) is expanded in F_p with zeta_5 in F_p and compared with
Sum_l x^(5l-1)/(l(5l-1)); log at = -(1/2) log(1-x^5) + Sum_k b_k (xz)^(2k-1) [5^(1-2k) - Sum_j (1-zeta^j x)^(1-2k)], the zeta-sum
again in F_p (rf_sp.stirling_sigma).  V = Sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)] (e^(-G/z) at).  Claim: V = C0 (Y^C_Q)_0 (kappa = C0,
s = 1).  The referee's bookkeeping of Theorem P steps 4-6 gives kappa(u) = C0 exactly (REFEREE.md s.3), in every frame.
Controls: (i) Stirling terms removed (b_k = 0), (ii) gamma = 1/2 instead of 7/10 in the vanishing period.
(b) CFU: Y = X(1 + Sum_k u^k Y_k(T)), T = X^5, V(u) with W^(X;u) = V/Y + Y^4/4; order u^k: (T-1) Y_k + V_k = R_k(T).
Checks: Y_k in Q[T] (exact division), V = tau^(-5/4), Y_1 = (1+T)/18, and Y(X=1)^5 = V (critical point to critical point).
Run: py -I -B rf_res.py NU KCFU
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import P61, P89, md, Ser, bernoulli, zeta5, frame_data
from rf_vanish import quintic_vanishing
from rf_sp import stirling_sigma


def ladd(A, B, c=1, p=None):
    out = dict(A)
    for k, v in B.items(): out[k] = (out.get(k, 0) + c * v) % p
    return {k: v for k, v in out.items() if v}


def lmul(A, B, p):
    out = {}
    for i, x in A.items():
        for j, y in B.items(): out[i + j] = (out.get(i + j, 0) + x * y) % p
    return {k: v for k, v in out.items() if v}


def residue_V(p, NU, stirling=True):
    zt = zeta5(p); D = NU + 1; Bn = bernoulli(2 * NU + 4)
    Gx = [0] * (D + 1)                                     # Sum_j (1/x - zeta^j) log(1 - zeta^j x), log(1-a x) = -Sum a^n x^n/n
    for j in range(5):
        a = pow(zt, j, p)
        for n in range(1, D + 2):
            c = -pow(a, n, p) * pow(n, -1, p) % p
            if n - 1 <= D: Gx[n - 1] = (Gx[n - 1] + c) % p
            if n <= D: Gx[n] = (Gx[n] - a * c) % p
    closed = [md(Fr(1, ((m + 1) // 5) * m), p) if (m + 1) % 5 == 0 and m > 0 else 0 for m in range(D + 1)]
    assert Gx == closed, "cyclotomic identity for G fails"
    Lg = [dict() for _ in range(D)]                        # log(e^(-G/z) at) as series in x, Laurent in z
    for m in range(1, D):
        if Gx[m]: Lg[m] = ladd(Lg[m], {-1: (-Gx[m]) % p}, 1, p)
        if m % 5 == 0: Lg[m] = ladd(Lg[m], {0: md(Fr(1, 2 * (m // 5)), p)}, 1, p)
    if stirling:
        sig = stirling_sigma(p, NU, NU // 5 + 1)
        for k in range(1, NU + 1):
            e = 2 * k - 1; bk = md(Bn[2 * k] / (2 * k * e), p)
            for P in range(0, NU // 5 + 2):
                if e + 5 * P >= D: break
                c = ((md(Fr(1, 5 ** e), p) if P == 0 else 0) - sig[k][P]) * bk % p
                Lg[e + 5 * P] = ladd(Lg[e + 5 * P], {e: c}, 1, p)
    E = [dict() for _ in range(D)]; E[0] = {0: 1}
    for n in range(1, D):
        acc = {}
        for k in range(1, n + 1):
            if Lg[k] and E[n - k]: acc = ladd(acc, lmul(Lg[k], E[n - k], p), k, p)
        E[n] = {a: b * pow(n, -1, p) % p for a, b in acc.items()}
    S = Ser(p, D + 1); Lam = [0] + [pow(m, -1, p) for m in range(1, D + 1)]
    V = [dict() for _ in range(D + 1)]; Lk = S.one(); fact = 1
    for k in range(1, D + 1):
        Lk = S.mul(Lk, Lam); fact = fact * k % p
        cf = (-1) ** k * pow(fact, -1, p) % p
        for i in range(D + 1):
            if Lk[i] and k - 1 < D:
                for zp, c in E[k - 1].items():
                    V[i][zp - k] = (V[i].get(zp - k, 0) + cf * Lk[i] * c) % p
    return [{a: b for a, b in d.items() if b} for d in V]


def vanishing_row0(p, NU, al, be, gamma=Fr(7, 10)):
    S = Ser(p, NU + 1); fd = frame_data(S, al, be, gamma); ys = quintic_vanishing(S, fd['b'], NU)
    out = [dict() for _ in range(NU + 1)]
    for n, y in enumerate(ys):
        t = S.mul(fd['C0'], y[0])
        for i in range(NU + 1):
            if t[i]: out[i][-n] = t[i]
    return out


def cfu(K):
    """exact CFU data: Y[k] (polys in T), V (list)"""
    def pm(a, b):
        out = [Fr(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): out[i + j] += x * y
        return out
    def bm(A, B):                                          # bivariate: list over u of T-polys, truncated at u^K
        out = [[Fr(0)] for _ in range(K)]
        for i in range(K):
            for j in range(K - i):
                pr = pm(A[i], B[j]); q = out[i + j]
                out[i + j] = [(q[t] if t < len(q) else 0) + (pr[t] if t < len(pr) else 0) for t in range(max(len(q), len(pr)))]
        return out
    XW = [[Fr(1, m + 1)] + [Fr(0)] * m + [Fr(1, (m + 1) * (5 * m + 4))] for m in range(K)]
    Y = [[Fr(0)] for _ in range(K)]; V = [Fr(1)] + [Fr(0)] * (K - 1)
    for k in range(1, K):
        S1 = [[Fr(1)]] + [Y[j] for j in range(1, K)]       # 1 + S with Y_k = 0, V_k = 0 at this stage
        inv = [[Fr(1)]] + [[Fr(0)] for _ in range(K - 1)]
        for n in range(1, K):
            acc = [Fr(0)]
            for j in range(1, n + 1):
                pr = pm(S1[j], inv[n - j]); acc = [(acc[t] if t < len(acc) else 0) - (pr[t] if t < len(pr) else 0) for t in range(max(len(acc), len(pr)))]
            inv[n] = acc
        s2 = bm(S1, S1); s4 = bm(s2, s2); Vs = [[V[j]] for j in range(K)]
        a = bm(Vs, inv)[k]; b4 = [Fr(0)] + [x / 4 for x in s4[k]]
        n = max(len(XW[k]), len(a), len(b4))
        Rk = [(XW[k][t] if t < len(XW[k]) else 0) - (a[t] if t < len(a) else 0) - (b4[t] if t < len(b4) else 0) for t in range(n)]
        Vk = sum(Rk); num = list(Rk); num[0] -= Vk
        q = [Fr(0)] * (len(num) - 1); r = list(num)
        for t in range(len(num) - 1, 0, -1): q[t - 1] = r[t]; r[t - 1] += r[t]; r[t] = Fr(0)
        assert r[0] == 0
        Y[k] = q; V[k] = Vk
    return Y, V


if __name__ == "__main__":
    NU = int(sys.argv[1]) if len(sys.argv) > 1 else 14; KC = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    t0 = time.time(); ok = True
    print("rf_res: residue V versus C0 (Y^C_Q)_0, u^1..u^%d, both frames, primes 2^61-1, 2^89-1; CFU exact, k < %d" % (NU, KC))
    for p in (P61, P89):
        V = residue_V(p, NU)
        pos = [i for i in range(1, NU + 1) if any(a > 0 for a in V[i])]
        for al, be in ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7))):
            W = vanishing_row0(p, NU, al, be)
            same = [i for i in range(1, NU + 1) if V[i] == W[i]]
            print("p=%d frame (%s,%s): V == C0 (Y^C_Q)_0 at u^1..u^%d: %d/%d orders; positive z-powers in V: %s"
                  % (p, al, be, NU, len(same), NU, pos if pos else "none"))
            ok = ok and len(same) == NU and not pos
        Vn = residue_V(p, NU, stirling=False); Wg = vanishing_row0(p, NU, Fr(0), Fr(0), Fr(1, 2))
        W0 = vanishing_row0(p, NU, Fr(0), Fr(0))
        c1 = [i for i in range(1, NU + 1) if Vn[i] != W0[i]]; c2 = [i for i in range(1, NU + 1) if V[i] != Wg[i]]
        print("   CONTROL no Stirling terms: first mismatch u^%s; CONTROL gamma = 1/2: first mismatch u^%s -> %s"
              % (c1[0] if c1 else None, c2[0] if c2 else None, "fail as they should" if (c1 and c2) else "UNEXPECTED"))
        ok = ok and bool(c1) and bool(c2)
    Y, Vc = cfu(KC); St = Ser(P89, KC)
    tau = St.inv([md(Fr(4, 5 * n + 4), P89) for n in range(KC)]); Vt = St.powr(tau, Fr(-5, 4))
    okV = [md(x, P89) for x in Vc] == Vt
    Y1 = Y[1] == [Fr(1, 18), Fr(1, 18)]
    y1 = [sum(Y[k]) for k in range(KC)]; y1[0] = Fr(1)                    # Y(X=1) = 1 + Sum u^k Y_k(1)
    y1m = [md(x, P89) for x in y1]; crit = St.mul(St.mul(St.mul(y1m, y1m), St.mul(y1m, y1m)), y1m) == [md(x, P89) for x in Vc]
    print("CFU (exact): Y_k in Q[T] for k < %d; V = tau^(-5/4) (mod 2^89-1 after exact computation): %s; Y_1 = (1+T)/18: %s; "
          "Y(1)^5 = V: %s; V = %s" % (KC, okV, Y1, crit, [str(x) for x in Vc[:4]]))
    ok = ok and okV and Y1 and crit
    print("RESULT rf_res: %s  [%.1f s]" % ("PASS" if ok else "FAIL", time.time() - t0))
