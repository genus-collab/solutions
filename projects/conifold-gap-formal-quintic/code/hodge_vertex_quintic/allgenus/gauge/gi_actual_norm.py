# -*- coding: utf-8 -*-
"""gi_actual_norm.py -- is the gauge of the ACTUAL R-matrix of the formal quintic the CANONICAL gauge (P_11(0, z) = 1)?
(GAUGE_IDENTITY.md s.6.)  Exact mod p.

Actual side (no gauge input): the referee's from-scratch hatted rows natural_form/referee/ref_hat.HatRows (QRR normalisation at q = 0,
zhat <= NZ, u < K), the symbol M_0 = R~(zhat; 0) from the same rows, tau = 1/2F1(1, 4/5; 9/5; u), and
        Ghat^act(u, zhat) := Qhat(u, zhat) M_0(tau zhat)^-1 Psi0^-1        (flat x flat; Psi0[b][i] = eps_b zeta^(i(b+1)))
which is the hatted form of the gauge of the actual R-matrix if (H) holds (UNIFORM_REDUCTION s.4: Qhat = Ghat R^(tau zhat)).
Canonical side (no R-matrix input): ur_gauge.gauge over Q from (varpi_0, t_c, tau), normalised by P_11(0, z) = 1, hatted:
        [zhat^k] Ghat_{b+k, b} = tau^((1 - c_b)/4) u^mu(b,k) p_{k,b}(u).
Checks: (0) grading and the valuations u^mu(b,k) of Ghat^act (these are (H) + Theorem U, not assumed);
        (1) Ghat^act = Ghat^can entry by entry, all k <= NZ, all u^N, N < K;   (C) control: the stored normalisation differs;
        (2) the normalisation slots: p^act_{5q,1}(0) = [u^4q zhat^5q] Ghat^act_11 for q = 1..4, i.e. c^act(z) := P^act_11(0, z) through z^20.
Since both G^act and the canonical gauge are symplectic, c^act(z) c^act(-z) = 1, i.e. c^act = exp(a_1 z^5 + a_3 z^15 + ...): the
z^10 and z^20 coefficients are FORCED by those at z^5, z^15, and only z^5 and z^15 are independent tests (printed as such).
Run: py -I -B gi_actual_norm.py [alpha beta] [NZ K]   (defaults 0 0 24 18; both primes)   log: gi_actual_norm_<frame>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import PROJ, Q, Ser, gauge_frame, p11_at_zero, inv_poly, normalise
sys.path.insert(0, os.path.join(PROJ, "natural_form", "referee"))
from fq_fp import FpRing, P1, P2
from ref_hat import HatRows, M0_of, tau_pred
CEXP = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)


def mu(b, k):
    v = 4 * k + CEXP[b] - CEXP[(b + k) % 5]; assert v % 5 == 0
    return v // 5


def smul(a, b, p, K):
    out = [0] * K
    for i, x in enumerate(a[:K]):
        if x:
            for j in range(min(len(b), K - i)): out[i + j] = (out[i + j] + x * b[j]) % p
    return out


def hatted_can(pg, tau, NZ, K):
    """{(k, b): list of K Fractions} = [u^0..u^(K-1)] of [zhat^k] Ghat_{b+k, b} for the gauge pg (over Q)"""
    lt = tau.log(); out = {}
    for b in range(5):
        f = lt.scale(Fr(1 - CEXP[b], 4)).exp().truncate(K)
        for k in range(min(NZ, len(pg) - 1) + 1):
            m = mu(b, k)
            if m >= K: out[(k, b)] = [Fr(0)] * K; continue
            s = f.mul(pg[k][b].truncate(K - m)).shift(m)
            assert s.prec >= K or not s.c, "canonical gauge precision too low at (k, b) = (%d, %d): O(u^%d)" % (k, b, s.prec)
            out[(k, b)] = [s.coef(N) if N >= m else Fr(0) for N in range(K)]
    return out


def actual(p, NZ, K, al, be, log):
    R = FpRing(p); z = R.zeta; i5 = R.inv(5)
    H = HatRows(p, NZ, K, al, be, zeta=z)
    log("   referee rows: own checks %s" % H.checks)
    M0, M0inv = M0_of(H.Rt, NZ, p)
    tq = [R.of(x) for x in tau_pred(K)]
    tp = [[1] + [0] * (K - 1)]
    for k in range(1, NZ + 1): tp.append(smul(tp[-1], tq, p, K))
    psi0_ok = all(H.Q[b][i][0][0] % p == EPS[b] * pow(z, (i * (b + 1)) % 5, p) % p for b in range(5) for i in range(5))
    pinv = [[EPS[b] * pow(z, (-i * (b + 1)) % 5, p) * i5 % p for b in range(5)] for i in range(5)]
    G = {}                                                    # G[(n, a, b)] = list of K residues: [zhat^n] Ghat^act_ab
    for a in range(5):
        for n in range(NZ + 1):
            Qp = []                                           # Q'[a][i][n], i = 0..4
            for i in range(5):
                acc = [0] * K
                for k in range(n + 1):
                    for j in range(5):
                        c = M0inv[k][j][i]
                        if c:
                            t = smul(H.Q[a][j][n - k], tp[k], p, K)
                            acc = [(x + c * y) % p for x, y in zip(acc, t)]
                Qp.append(acc)
            for b in range(5):
                G[(n, a, b)] = [sum(Qp[i][N] * pinv[i][b] for i in range(5)) % p for N in range(K)]
    return R, G, psi0_ok, H.checks


def main():
    al = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(0); be = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(0)
    NZ = int(sys.argv[3]) if len(sys.argv) > 3 else 24; K = int(sys.argv[4]) if len(sys.argv) > 4 else 18
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "gi_actual_norm_%s.log" % tag), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(name)

    t0 = time.time()
    log("gi_actual_norm: frame (%s, %s), referee rows zhat <= %d, u < %d; canonical gauge from the periods over Q" % (al, be, NZ, K))
    res = gauge_frame(NZ, 72, al, be); ps = res['p']; tau = res['tau']
    c11 = p11_at_zero(ps); pc = normalise(ps, inv_poly(c11, NZ))
    can = hatted_can(pc, tau, NZ, K); sto = hatted_can(ps, tau, NZ, K)
    log("   canonical gauge (k <= %d) in %.0f s; stored normalisation c(z) = P^s_11(0, z) = %s" % (NZ, time.time() - t0, " + ".join("(%s) z^%d" % (c11[k], k) for k in sorted(c11) if k <= 10) + " + ..."))
    for p, name in ((P2, "2^89-1"), (P1, "2^61-1")):
        t1 = time.time()
        R, G, psi0_ok, chk = actual(p, NZ, K, al, be, log)
        log("-- p = %s: actual Ghat from the referee rows in %.0f s" % (name, time.time() - t1))
        check("[%s] referee rows pass their own checks; Psi0 = Qhat(u = 0, zhat = 0) = eps_b zeta^(i(b+1))" % name, all(chk.values()) and psi0_ok)
        grad = all(not any(G[(n, a, b)]) for (n, a, b) in G if (a - b - n) % 5)
        vals = all(not any(G[(n, (b + n) % 5, b)][:min(mu(b, n), K)]) for n in range(NZ + 1) for b in range(5))
        check("[%s] (0) Ghat^act is Z/5-graded (a = b + k) and [zhat^k] Ghat^act_{b+k,b} = O(u^mu(b,k)) for all k <= %d" % (name, NZ), grad and vals)
        tested = 0; nonzero = 0; bad = []; diff_st = []
        for n in range(NZ + 1):
            for b in range(5):
                a = (b + n) % 5; g = G[(n, a, b)]
                for N in range(K):
                    if N < mu(b, n): continue
                    tested += 1; c = R.of(can[(n, b)][N]); s = R.of(sto[(n, b)][N])
                    if g[N]: nonzero += 1
                    if g[N] != c: bad.append((n, a, b, N))
                    if g[N] != s: diff_st.append((n, a, b, N))
        check("[%s] (1) Ghat^act = Ghat^can: all %d slots [zhat^k u^N] (k <= %d, N < %d, N >= mu(b,k); %d of them non-zero) agree"
              % (name, tested, NZ, K, nonzero), not bad, bad[:5])
        check("[%s] (C) control: the stored normalisation P^s differs from the actual gauge (%d slots), first at (k, a, b, N) = %s"
              % (name, len(diff_st), min(diff_st) if diff_st else None), bool(diff_st) and min(n for (n, a, b, N) in diff_st) == 5)
        rows = []
        for q in range(1, 5):
            k = 5 * q
            if 4 * q >= K or k > NZ: continue
            act11 = G[(k, 1, 1)][4 * q]; act22 = G[(k, 2, 2)][4 * q]
            rows.append((k, act11, act22, R.of(c11[k])))
        log("   [%s] normalisation slots [u^4q zhat^5q] of Ghat_11, Ghat_22 (= p_{5q,1}(0), p_{5q,2}(0)); stored p^s_{5q,1}(0) for contrast:" % name)
        for (k, x1, x2, s) in rows:
            log("      z^%-2d: actual (%s, %s)   stored %s   %s" % (k, x1, x2, s, "INDEPENDENT test (odd multiple of z^5)" if (k // 5) % 2 else "forced by symplecticity from the lower odd orders"))
        check("[%s] (2) c^act(z) = P^act_11(0, z) = P^act_22(0, z) = 1 through z^%d (z^5, z^15 independent; z^10, z^20 forced)" % (name, max(k for (k, *_ ) in rows)),
              all(x1 == 0 and x2 == 0 for (k, x1, x2, s) in rows) and G[(0, 1, 1)][0] == 1 and G[(0, 2, 2)][0] == 1)
    log("GI_ACTUAL_NORM %s: %s   (%.0f s)" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
