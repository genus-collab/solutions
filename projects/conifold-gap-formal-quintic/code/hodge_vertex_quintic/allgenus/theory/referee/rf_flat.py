# -*- coding: utf-8 -*-
"""rf_flat.py -- referee check, in the FLAT frame and exactly over Q, of Corollary E' and Lemma T (REDUCTION_DRAFT_2.md s.3.1) for the
canonical Theorem-U gauge P of the formal quintic, with the stored normalisation of ur_gauge as control.  Helpers: rf_flat_core.py.

Conventions (REDUCTION_DRAFT.md s.1): P acts on LOWERED components; vector form Pv = eta^-1 P eta; V^(k)(q) = (-1)^k [z^(k-1)] <Y^C(-z), q(z)>.
Checks:  (M) Y^C solves the model system (monomials);  (Sp) P(-z)^T eta P(z) = eta;  (G) (P Y^C)_+ = 0, (P Y^C)_0 = -e_1;
(T) T_P := -[z^-1](P Y^C)_0 = t_c;  (LT) V^(1)(q_T) = T_P and V^(k)(q_T) = 0 (k = 2..KV) for q_T = Pv^-1(-z phi_0), by DIRECT pairing
(the symplectic identity of Lemma T is not used);  (D) d = (P_0)_33 = 1/(P_0)_00 = tau^(1/4) varpi_0, s = d T_P - v = O(u^2);
(E') the flat edge E_P (divisibility by psi + psi1 checked) contracted with V^(k) in slot 1 vanishes, k = 1..KV, on every complete
coefficient; control: the same contraction without the sign (-1)^n.  Only u-orders that are complete in the truncation are compared.
Run: py -I -B rf_flat.py [K NU a b]   Log: rf_flat.log (appended).
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat_core as C
from rf_flat_core import (s0, sc, sadd, sscale, smul, sinv, spow, sval, nz, mzero, madd, mT, mmul, mvec, mscale, zinv, eta_conj,
                          tau_v, yc_series, check_model_system, load_gauge, canonical)
ETA = 5


def sympl_orders(P, K):
    """orders n <= K where (P(-z)^T eta P(z))_n differs from eta delta_{n0} (lowered form)"""
    N = C.NU[0]; bad = []
    for n in range(K + 1):
        S = mzero()
        for k in range(n + 1):
            A = mscale(mT(P[k]), Fr((-1) ** k))
            # eta = 5 x (permutation a <-> 3 - a):  A eta B = 5 A[:, 3-.] B
            Bp = [[sscale(P[n - k][(3 - i) % 5][j], Fr(ETA)) for j in range(5)] for i in range(5)]
            S = madd(S, mmul(A, Bp))
        for i in range(5):
            for j in range(5):
                target = sc(ETA) if (n == 0 and (i + j) % 5 == 3) else s0()
                if S[i][j] != target: bad.append((n, i, j)); break
    return bad


def run(K, a, b, log):
    t0 = time.time(); N = C.NU[0]
    tau, v, itau = tau_v()
    nmax = (5 * N) // 4 + 2
    log("  (M) model system on monomials, z^-n, n <= %d: mismatches %d" % (nmax, check_model_system(nmax)))
    Ps, tc, var0, res = load_gauge(K, a, b)
    log("  v(u) own = ur_gauge v: %s;  t_c = %s ...; varpi_0 = %s ..." % (v == C.from_ser(res['v']), tc[:4], var0[:4]))
    Pc, cser = canonical(Ps, K)
    log("  stored P_11(0, z) = %s" % ", ".join("%s z^%d" % (cser[k], k) for k in sorted(cser)))
    Y = yc_series(v, nmax)
    out = {}
    for name, P in (("canonical", Pc), ("stored", Ps)):
        log("  --- %s gauge (K = %d, u < %d)" % (name, K, N))
        sb = sympl_orders(P, K); log("  (Sp) symplectic defects (first 3): %s" % sb[:3])
        # Y' = P Y^C : coefficient of z^j
        Yp = {}
        for k in range(K + 1):
            for n in range(nmax + 1):
                j = k - n; vec = mvec(P[k], Y[n])
                Yp[j] = [sadd(x, y) for x, y in zip(Yp.get(j, [s0() for _ in range(5)]), vec)]
        def complete(j): return (4 * (K - j + 1) - 1) // 5          # u-orders <= this are complete at z^j
        pos_bad = []; pos_n = 0
        for j in range(1, K + 1):
            cu = min(complete(j), N - 1)
            for bb in range(5):
                for Nn in range(cu + 1):
                    pos_n += 1
                    if Yp[j][bb][Nn] != 0: pos_bad.append((j, bb, Nn))
        z0 = Yp[0]; cu0 = min(complete(0), N - 1)
        z0ok = all(z0[bb][Nn] == (-1 if (bb == 1 and Nn == 0) else 0) for bb in range(5) for Nn in range(cu0 + 1))
        TP = sscale(Yp[-1][0], Fr(-1)); cum1 = min(complete(-1), N - 1)
        log("  (G) positive powers: %d complete coefficients, %d non-zero (first %s); z^0 part = -e_1: %s (u < %d)"
            % (pos_n, len(pos_bad), pos_bad[:3], z0ok, cu0 + 1))
        log("  (T) T_P = t_c through u^%d: %s;  T_P = %s ..." % (cum1, TP[:cum1 + 1] == tc[:cum1 + 1], TP[:4]))
        # leaf point q_T = Pv^-1(-z e_0), vector form; V^(k)(q_T) by direct pairing
        Pi = zinv(P, K); Pvi = {k: eta_conj(Pi[k]) for k in Pi}
        q = {k + 1: [sscale(Pvi[k][i][0], Fr(-1)) for i in range(5)] for k in Pvi}
        KV = 8; lt = []
        for kk in range(1, KV + 1):
            acc = s0(); cu = min((4 * (K + 3 - kk) - 1) // 5, N - 1)
            for n in range(0, nmax + 1):
                m = kk - 1 + n
                if m not in q: continue
                for aa in range(5):
                    if nz(Y[n][aa]) and nz(q[m][aa]): acc = sadd(acc, sscale(smul(Y[n][aa], q[m][aa]), Fr((-1) ** n)))
            Vk = sscale(acc, Fr((-1) ** kk))
            target = TP if kk == 1 else s0()
            lt.append((kk, cu, Vk[:cu + 1] == target[:cu + 1], sval(Vk)))
        log("  (LT) [k, complete u-orders <= , V^(k)(q_T) = (T_P if k = 1 else 0), u-valuation]: %s" % lt)
        # d and s
        d = P[0][3][3]; d0 = P[0][0][0]
        dd = smul(spow(itau, Fr(-1, 4)), var0)
        s = sadd(smul(d, TP), v, -1)
        log("  (D) d (P_0)_00 = 1: %s;  d = tau^(1/4) varpi_0: %s;  s = d T_P - v: valuation %s, s = %s ..."
            % (smul(d, d0) == sc(1), d == dd, sval(s[:cum1 + 1]), s[:4]))
        # flat edge in vector form and its contraction with V^(k)
        EI = mzero()
        for i in range(5): EI[i][(3 - i) % 5] = sc(Fr(1, 5))
        Num = {}
        for m in range(K + 1):
            for n in range(K + 1 - m):
                Num[(m, n)] = mscale(mmul(mmul(Pvi[m], EI), mT(Pvi[n])), Fr(-1))
        Num[(0, 0)] = madd(Num[(0, 0)], EI)
        E = {}
        for m in range(K):
            for n in range(K - m):
                acc = mzero()
                for t in range(n + 1): acc = madd(acc, Num[(m + 1 + t, n - t)], (-1) ** t)
                E[(m, n)] = acc
        divbad = [n for n in range(1, K) if madd(Num[(0, n)], E[(0, n - 1)], -1) != mzero()]
        if Num[(0, 0)] != mzero(): divbad = [0] + divbad
        ec = dict(n=0, bad=[], nz=0, ctrl=0)
        for kk in range(1, KV + 1):
            for n in range(0, K):
                cu = min((4 * (K - n - kk + 1) - 1) // 5, N - 1)
                if cu < 0: continue
                vec = [s0() for _ in range(5)]; ctl = [s0() for _ in range(5)]; terms = 0
                for m in range(kk - 1, K - n):
                    nn = m + 1 - kk
                    if nn > nmax: break
                    for aa in range(5):
                        if not nz(Y[nn][aa]): continue
                        for bb in range(5):
                            if nz(E[(m, n)][aa][bb]):
                                t = smul(Y[nn][aa], E[(m, n)][aa][bb]); terms += 1
                                vec[bb] = sadd(vec[bb], sscale(t, Fr((-1) ** (kk + nn))))
                                ctl[bb] = sadd(ctl[bb], sscale(t, Fr((-1) ** kk)))
                for bb in range(5):
                    ec['n'] += cu + 1
                    if any(vec[bb][:cu + 1]): ec['bad'].append((kk, n, bb, sval(vec[bb])))
                    if any(ctl[bb][:cu + 1]): ec['ctrl'] += 1
        log("  (E') edge divisibility defects %s; V^(k)-contractions: %d complete coefficients, %d non-zero (first %s); "
            "control without (-1)^n: %d of the (k, n, b) slots non-zero" % (divbad[:3], ec['n'], len(ec['bad']), ec['bad'][:4], ec['ctrl']))
        out[name] = dict(sb=sb, pos_bad=pos_bad, z0ok=z0ok, lt=lt, ebad=ec['bad'], divbad=divbad, ctrl=ec['ctrl'],
                         tp=(TP[:cum1 + 1] == tc[:cum1 + 1]), d_ok=(d == dd), s_val=sval(s[:cum1 + 1]))
    log("  time %.0f s" % (time.time() - t0))
    return out


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    C.NU[0] = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    a = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); b = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    LOG = open(os.path.join(HERE, "rf_flat.log"), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 110); log("rf_flat: K = %d, NU = %d, frame (%s, %s) [%s]" % (K, C.NU[0], a, b, time.strftime("%Y-%m-%dT%H:%M:%S")))
    o = run(K, a, b, log)
    c = o['canonical']; s = o['stored']
    okc = (not c['sb'] and not c['pos_bad'] and c['z0ok'] and c['tp'] and all(x[2] for x in c['lt']) and not c['ebad']
           and not c['divbad'] and c['ctrl'] > 0 and c['d_ok'] and (c['s_val'] is None or c['s_val'] >= 2))
    okctl = bool(s['pos_bad']) and bool(s['ebad']) and not all(x[2] for x in s['lt'])
    log("RF_FLAT: canonical %s; stored-normalisation control fails as expected: %s" % ("ALL PASS" if okc else "FAIL", okctl))
    LOG.close()


if __name__ == "__main__":
    main()
