# -*- coding: utf-8 -*-
"""v2_run.py -- driver for v2_gauge.py (adversarial check of UNIFORM_REDUCTION.md s.3-5).

(a) own gauge (proof route) for k <= K: holomorphy asserted at every step; quoted numbers of the note re-derived:
    p_{1,b}(0), valuations, pi_3 = P_30 at z^3, (R1)-(R3), p_{5,b}(0);
(b) agreement with natural_form/theory/ur_gauge.gauge (a copy in this folder) coefficient by coefficient;
(c) OWN comparison with the project's R-matrix: Ghat := Qhat(u, zhat) . R^(tau zhat)^(-1), the inverse taken as a genuine matrix-series
    inverse (Neumann series; no use of the symplectic formula), zhat^0..6, u^0..5; frames (0, 0) AND (1/3, -2/7), both primes;
    the constant c_5 is fitted on the entry (0, 0) at u^4 zhat^5 and compared with -p_{5,1}(0) in EACH frame.
Run: py -I -B v2_run.py
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from v2_gauge import gauge, S
CEXP = (0, -1, 3, 2, 1)
OUT = open(os.path.join(HERE, "v2_run.log"), "w", encoding="utf-8")
BAD = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: BAD.append(name)


def mu(b, k):
    w = 4 * k + CEXP[b] - CEXP[(b + k) % 5]; assert w % 5 == 0
    return w // 5


# ---------------------------------------------------------------- bivariate truncated series mod p: dict (n, j) -> value, n <= KZ, j < KU
def bmul(A, B, p, KZ, KU):
    out = {}
    for (n1, j1), x in A.items():
        for (n2, j2), y in B.items():
            n = n1 + n2; j = j1 + j2
            if n > KZ or j >= KU: continue
            out[(n, j)] = (out.get((n, j), 0) + x * y) % p
    return {k: v for k, v in out.items() if v}


def mmul(A, B, p, KZ, KU):
    out = [[{} for _ in range(5)] for _ in range(5)]
    for a in range(5):
        for b in range(5):
            acc = {}
            for c in range(5):
                for key, val in bmul(A[a][c], B[c][b], p, KZ, KU).items(): acc[key] = (acc.get(key, 0) + val) % p
            out[a][b] = {k: v for k, v in acc.items() if v}
    return out


def inv5(M, p):
    n = 5; A = [list(r) + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] % p)
        A[c], A[piv] = A[piv], A[c]
        iv = pow(A[c][c], p - 2, p); A[c] = [x * iv % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] % p:
                f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]


def project_ghat(st, tau, KZ, KU):
    """Ghat[a][b] with Qhat[a][i] = Sum_b Ghat[a][b] R^[b][i](tau zhat)"""
    R = st['R']; p = R.p; H = st['H']; sym = st['sym']
    tp = [S(0, [Fr(1)] + [Fr(0)] * (KU - 1))]
    for n in range(KZ): tp.append(tp[-1] * tau)
    Q = [[{(n, j): H.Q[a][i][n].coef(j) % p for n in range(KZ + 1) for j in range(KU)} for i in range(5)] for a in range(5)]
    Rt = [[{(n, j): sym.Rh[b][i][n] * R.of(tp[n].coef(j)) % p for n in range(KZ + 1) for j in range(KU)} for i in range(5)] for b in range(5)]
    R0 = [[sym.Rh[b][i][0] % p for i in range(5)] for b in range(5)]
    R0i = inv5(R0, p)
    R0iM = [[{(0, 0): R0i[i][b]} for b in range(5)] for i in range(5)]
    Rhigh = [[{k: v for k, v in Rt[b][i].items() if k != (0, 0)} for i in range(5)] for b in range(5)]
    assert all(Rt[b][i].get((0, j), 0) == 0 for b in range(5) for i in range(5) for j in range(1, KU))
    Nm = mmul(R0iM, Rhigh, p, KZ, KU)                      # R0^-1 (Rt - R0), colour x colour, zhat-order >= 1
    Nm = [[{k: (-v) % p for k, v in Nm[i][j].items()} for j in range(5)] for i in range(5)]
    acc = [[({(0, 0): 1} if i == j else {}) for j in range(5)] for i in range(5)]; term = acc
    for _ in range(KZ):
        term = mmul(term, Nm, p, KZ, KU)
        acc = [[{k: (acc[i][j].get(k, 0) + term[i][j].get(k, 0)) % p for k in set(acc[i][j]) | set(term[i][j])} for j in range(5)] for i in range(5)]
    Rinv = mmul(acc, R0iM, p, KZ, KU)                      # (1 + R0^-1 Rhigh)^-1 R0^-1
    one = mmul(Rt, Rinv, p, KZ, KU)
    assert all(one[a][b] == ({(0, 0): 1} if a == b else {}) for a in range(5) for b in range(5)), "matrix inverse"
    return mmul(Q, Rinv, p, KZ, KU)


def hat_pred(res, KZ, KU):
    tau = res['tau']; out = {}
    for k in range(KZ + 1):
        for b in range(5):
            f = tau.powq(Fr(1 - CEXP[b], 4)) * res['p'][k][b]
            out[(k, b)] = [f.coef(n - mu(b, k)) if n >= mu(b, k) else Fr(0) for n in range(KU)]
    return out


def main():
    K = 8; N = 22; KZ = 6; KU = 6
    t0 = time.time(); res = gauge(K, N, 0, 0, log); p = res['p']; log("own gauge, frame (0,0): %.0f s" % (time.time() - t0))
    tc = res['tc']; log("t_c = " + ", ".join(str(tc.coef(n)) for n in range(1, 6)) + "   [literature: 1, 7/10, 41/75, 1133/2500, 6089/15625]")
    log("tau = " + ", ".join(str(res['tau'].coef(n)) for n in range(5)) + ";  v = " + ", ".join(str(res['v'].coef(n)) for n in range(1, 5)))
    c0 = lambda k, b, n=0: p[k][b].coef(n)
    check("p_{1,b}(0) = (8/45, 0, 8/45, -1/15, -1/15)", [c0(1, b) for b in range(5)] == [Fr(8, 45), 0, Fr(8, 45), Fr(-1, 15), Fr(-1, 15)])
    check("(R2) p_{1,1}, p_{6,1} = O(u^2), with [u^2] p_{1,1} = -127/3150", p[1][1].v == 2 and p[6][1].v == 2 and c0(1, 1, 2) == Fr(-127, 3150))
    check("(R1) column 1 / row 2 of P(0,z) vanish off the diagonal, k <= %d; P_11 = P_22 at z^5" % K,
          all(c0(k, 1) == 0 and c0(k, (2 - k) % 5) == 0 for k in range(1, K + 1) if k % 5) and c0(5, 1) == c0(5, 2))
    check("(R3) pi_3 = p_{3,0}(0) = 127/3150 = [u^1] p_{2,0} = -[u^1] p_{2,1} = -[u^2] p_{1,1}",
          c0(3, 0) == Fr(127, 3150) == c0(2, 0, 1) == -c0(2, 1, 1) == -c0(1, 1, 2))
    check("(R3) pi_8 = p_{8,0}(0) = -77141089/362679778125000 = [u^1] p_{7,0} = -[u^1] p_{7,1} = -[u^2] p_{6,1}",
          c0(8, 0) == Fr(-77141089, 362679778125000) == c0(7, 0, 1) == -c0(7, 1, 1) == -c0(6, 1, 2))
    check("p_{5,b}(0) = (0, -43397/807975000, -43397/807975000, 18629/403987500, -94033/807975000)",
          [c0(5, b) for b in range(5)] == [0, Fr(-43397, 807975000), Fr(-43397, 807975000), Fr(18629, 403987500), Fr(-94033, 807975000)])
    fo = {k: min(mu(b, k) + p[k][b].v for b in range(5)) for k in range(K + 1)}
    check("first u-order of zhat^k, k = 0..8, is 0,1,2,2,3,4,5,6,6 (degrees 0,1,3,4,5,6 for n = 0..5; zhat^8 at n = 6)", [fo[k] for k in range(K + 1)] == [0, 1, 2, 2, 3, 4, 5, 6, 6], fo)
    import ur_gauge
    theirs = ur_gauge.gauge(K, N + 6)['p']
    bad = [(k, b, n) for k in range(K + 1) for b in range(5) for n in range(min(p[k][b].prec, theirs[k][b].prec, 10)) if p[k][b].coef(n) != theirs[k][b].coef(n)]
    check("own gauge (proof route) = ur_gauge.gauge, k <= %d, all b, u^0..u^9" % K, not bad, bad[:6])
    from fq_fp import P1, P2
    from fq3_par import setup
    for (al, be) in (("0", "0"), ("1/3", "-2/7")):
        rs = res if al == "0" else gauge(6, N, Fr(al), Fr(be))
        pred = hat_pred(rs, KZ, KU); p51 = rs['p'][5][1].coef(0)
        log("-" * 100); log("frame (alpha, beta) = (%s, %s):  p_{1,b}(0) = %s;  -p_{5,1}(0) = %s" % (al, be, [str(rs['p'][1][b].coef(0)) for b in range(5)], -p51))
        for prime, tag in ((P2, "2^89-1"), (P1, "2^61-1")):
            st = setup((prime, al, be, 34, KU)); R = st['R']; pp = R.p
            G = project_ghat(st, rs['tau'], KZ, KU)
            off = [(a, b, key) for a in range(5) for b in range(5) for key in G[a][b] if (a - b - key[0]) % 5]
            check("[%s] Z/5 grading of Qhat R^(tau zhat)^-1" % tag, not off, off[:4])
            val = lambda k, b, n: G[(b + k) % 5][b].get((k, n), 0)
            bad = [(k, b, n) for k in range(5) for b in range(5) for n in range(KU) if val(k, b, n) != R.of(pred[(k, b)][n])]
            check("[%s] Qhat R^(tau zhat)^-1 = own gauge, zhat^0..4, all entries, u^0..u^5" % tag, not bad, bad[:6])
            c5 = (val(5, 0, 4) - R.of(pred[(5, 0)][4])) % pp
            check("[%s] fitted c_5 = -p_{5,1}(0) = -p_{5,2}(0) (P_11(0,z) = P_22(0,z) = 1 mod z^10) in THIS frame" % tag, (c5 + R.of(p51)) % pp == 0 and rs['p'][5][2].coef(0) == p51)
            bad = []
            for k in (5, 6):
                for b in range(5):
                    for n in range(KU):
                        corr = R.of(pred[(k - 5, b)][n - 4]) * c5 % pp if n >= 4 else 0
                        if (R.of(pred[(k, b)][n]) + corr) % pp != val(k, b, n): bad.append((k, b, n))
            check("[%s] zhat^5, zhat^6, all entries, u^0..u^5, with the one constant c_5" % tag, not bad, bad[:6])
    log("V2_RUN " + ("ALL PASS" if not BAD else "FAILURES: %s" % BAD)); OUT.close()


if __name__ == "__main__":
    main()
