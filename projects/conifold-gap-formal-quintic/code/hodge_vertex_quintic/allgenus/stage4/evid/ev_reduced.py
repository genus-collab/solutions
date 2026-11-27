# -*- coding: utf-8 -*-
"""ev_reduced.py -- Stage 4 evidence, part (2): the REDUCED criterion of the connection problem, colour 0, exact over Q or mod p
(EVID.md s.2).  Uses only ev_core (periods, QRR constants, model rows F_b, the colour-0 WKB polynomials ph_n, hatted rows).

Ghat_ac[zhat^N] := Sum_{n<=N} Qh_{a,n}(u) (-tau)^(N-n) F_{3-c,N-n}  (Qh = colour-0 hatted rows, F = model rows) is the actual hatted
gauge Qhat R^(tau zhat)^-1 with the five colours summed by covariance (EVID s.2.2; validated slot by slot in ev_validate.py).
For the polynomial-type WKB solution with ANY free coefficients f_5q = [u^4q] ph_5q it is the hatted P_can c(z), c in Q[[z^5]], so
   REG_q:  [u^j] Ghat_11[zhat^5q] = 0 for j < 4q        (automatic, Theorem U)
   K_q  :  [u^4q] Ghat_11[zhat^5q] = [z^5q] c(z)       (also = [u^4q] Ghat_22[zhat^5q], by (R1))
and K_q = f_5q + (ph_n, n < 5q; periods; model).  D = 1 iff (H) ph_n(1) = h_n (5 not | n) and (N) K_q = 0 (all q) in the QRR
normalisation f_5q = h_5q - (ph_5q(1) - f_5q).  Mode 'can' (no QRR input on the conifold side): f_5q is fixed by K_q = 0 and the
POLYNOMIAL CONNECTION IDENTITIES ph^can_n(1) = h_n are tested for EVERY n <= 5 QMAX + 4.
Run (cwd outside the project): py -I -B ev_reduced.py FIELD MODE QMAX [alpha beta QREG CTRL]
   FIELD Q | 89 | 61 | 127;  MODE qrr | can;  QREG: REG_q and Ghat_22 checks for q <= QREG (default QMAX);  CTRL 1: failing controls
   log: logs/ev_reduced_<FIELD>_<MODE>_<frame>_<5 QMAX>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ev_core as E
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1, "127": (1 << 127) - 1}
M5 = Fr(-43397, 807975000)          # P_stored = P_can (1 + m5 z^5 + ..)  (UNIFORM_REDUCTION s.4, gauge s.1)


def short(F, x):
    if not isinstance(F, E.QF): return str(x)
    s = str(x)
    return s if len(s) <= 48 else "%s..(%d/%d digits)" % (s[:14], len(str(x.numerator)), len(str(x.denominator)))


def kcoef(F, Qa, Fb, tp, N, J):
    """[u^J] Sum_{n<=N} Qa[n] (-tau)^(N-n) F_b[N-n]"""
    s = 0
    for n in range(N + 1):
        m = N - n; c = Fb[m] if m % 2 == 0 else -Fb[m]
        if not c: continue
        a = Qa[n]; t = tp[m]; acc = 0
        for j in range(min(J, len(a) - 1) + 1): acc += a[j] * t[J - j]
        s += c * acc
    return F.red(s)


def coef(c, k): return c[k] if k < len(c) else 0


def run(F, mode, Q, al, be, QREG, log, h=None, perturb=None, tau_one=False, quiet=False):
    t0 = time.time(); KU = 4 * Q + 1; NZ = 5 * Q + 4; QREG = min(QREG, Q)
    h = E.qrr_h(F, NZ) if h is None else h
    C0, C1, tc, g = E.periods(F, KU, al, be)
    tau = ([F.of(1)] + [0] * (KU - 1)) if tau_one else E.tau_series(F, KU)
    Fm = E.model_F(F, 5 * Q)
    tp = [[F.of(1)] + [0] * (KU - 1)]
    for m in range(1, 5 * Q + 1): tp.append(E.smul(F, tp[-1], tau, KU))
    KR = 4 * QREG + 1
    R1 = E.Rows(F, C0, C1, KU, amax=1); R2 = E.Rows(F, C0[:KR], C1[:KR], KR, amax=2)
    W = E.Wkb(F); R1.set(0, W.ph[0]); R2.set(0, W.ph[0])
    slope = F.red(Fm[2][0] * R1.kh[0][0] * R1.C0i[0])           # dK_q/df_5q = F_2(0)/(C1(0) C0(0)) = 1
    rule = E.qrr_free(F, h); res = dict(H=[], q=[], t_setup=time.time() - t0)
    if not quiet: log("   setup (periods, tau^m, model rows to t^%d, u^%d) in %.1f s; dK_q/df_5q = %s" % (5 * Q, KU - 1, res['t_setup'], slope))
    for n in range(1, NZ + 1):
        q, r = divmod(n, 5)
        W.step(rule if (r == 0 and mode == 'qrr') else None)
        if r == 0 and perturb and n in perturb: W.set_free(n, F.red(W.free[n] + F.of(perturb[n])))
        if n <= 5 * Q: R1.set(n, W.ph[n])
        if n <= 5 * QREG: R2.set(n, W.ph[n])
        if r == 0 and q <= Q:
            kf = 4 * q; Kq = kcoef(F, R1.Q[1], Fm[2], tp, n, kf); row = dict(q=q, K=Kq)
            if mode == 'can':
                W.set_free(n, F.red(-Kq * F.inv(slope))); R1.set(n, W.ph[n])
                if n <= 5 * QREG: R2.set(n, W.ph[n])
                row['K_after'] = kcoef(F, R1.Q[1], Fm[2], tp, n, kf)
            row['f'] = W.free.get(n, coef(W.ph[n], kf))
            row['f_qrr'] = F.red(h[n] - sum(W.ph[n]) + coef(W.ph[n], kf))
            if q <= QREG:
                G11 = E.pair(F, R1.Q[1], Fm[2], tp, n, kf + 1); G22 = E.pair(F, R2.Q[2], Fm[1], tp, n, kf + 1)
                row['reg'] = (all(F.red(x) == 0 for x in G11[:kf]), all(F.red(x) == 0 for x in G22[:kf]))
                row['k11_full'] = F.red(G11[kf] - (row.get('K_after', Kq) if mode == 'can' else Kq)) == 0
                row['k22'] = G22[kf]
            row['t'] = time.time() - t0; res['q'].append(row)
        res['H'].append((n, F.red(sum(W.ph[n]) - h[n]) == 0))
    res['W'] = W; res['t'] = time.time() - t0
    return res


def main():
    a = sys.argv[1:]
    fld = a[0] if len(a) > 0 else "Q"; mode = a[1] if len(a) > 1 else "qrr"; Q = int(a[2]) if len(a) > 2 else 4
    al = Fr(a[3]) if len(a) > 3 else Fr(0); be = Fr(a[4]) if len(a) > 4 else Fr(0)
    QREG = int(a[5]) if len(a) > 5 else Q; ctrl = len(a) > 6 and a[6] == "1"
    F = E.QF() if fld == "Q" else E.PF(PR[fld])
    fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    name = "ev_reduced_%s_%s_%s_%d%s.log" % (fld, mode, fr, 5 * Q, "_ctrl" if ctrl else "")
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)

    log("ev_reduced: field %s, mode %s, frame (%s, %s), q <= %d (zhat^%d, u^%d), (H) for n <= %d, REG/Ghat_22 for q <= %d"
        % (F.name, mode, al, be, Q, 5 * Q, 4 * Q, 5 * Q + 4, min(QREG, Q)))
    log("   log file logs/%s" % name)
    res = run(F, mode, Q, al, be, QREG, log)
    for row in res['q']:
        q = row['q']; n = 5 * q
        Hblk = [ok for (m, ok) in res['H'] if n - 4 <= m <= n - 1]
        s = "   z^%-3d: " % n
        if mode == 'qrr': s += "f^QRR = %s; K_q = %s" % (short(F, row['f']), short(F, row['K']))
        else: s += "f^can = %s; ph^can(1) - h = %s" % (short(F, row['f']), short(F, F.red(row['f'] - row['f_qrr'])))
        if 'reg' in row: s += "; REG %s/%s; Ghat_22 slot %s" % ("ok" if row['reg'][0] else "BAD", "ok" if row['reg'][1] else "BAD",
                                                               short(F, row['k22']))
        s += "; (H) n = %d..%d %s; %.0f s" % (n - 4, n - 1, "ok" if all(Hblk) else "FAIL", row['t'])
        log(s)
    Hbad = [m for (m, ok) in res['H'] if not ok and m % 5]
    check("(H): ph_n(1) = h_n for all %d n <= %d with 5 not | n (no free parameter: the WKB recursion alone)"
          % (sum(1 for (m, ok) in res['H'] if m % 5), 5 * Q + 4), not Hbad, Hbad[:5])
    rs = res['W'].resonance
    check("no resonance obstruction at the %d multiples of 5 <= %d (%d of them because deg rho_n < 4n/5): polynomial-type solution exists"
          % (len(rs), 5 * Q + 4, sum(1 for x in rs if x[2])), len(rs) == Q and all(x[1] for x in rs))
    if mode == 'qrr':
        bad = [r['q'] for r in res['q'] if r['K'] != 0]
        check("(N): K_q = [z^5q] c^act = 0 for q = 1..%d, i.e. a_1 = a_3 = .. = 0 through z^%d" % (Q, 5 * Q), not bad, bad[:5])
    else:
        bad = [r['q'] for r in res['q'] if F.red(r['f'] - r['f_qrr']) != 0 or r['K_after'] != 0]
        check("(N), connection form: ph^can_5q(1) = h_5q for q = 1..%d (f^can fixed at the conifold by K_q = 0, no QRR input)" % Q,
              not bad, bad[:5])
    regs = [r for r in res['q'] if 'reg' in r]
    check("REG and (R1): [u^j] Ghat_11 = [u^j] Ghat_22 = 0 (j < 4q) and Ghat_22 slot = Ghat_11 slot, q <= %d" % len(regs),
          all(r['reg'][0] and r['reg'][1] and r['k11_full'] and F.red(r['k22'] - (r['K_after'] if mode == 'can' else r['K'])) == 0
              for r in regs))
    if ctrl:
        Qc = min(Q, 4); log("-- failing controls (q <= %d); each must FAIL the test it targets" % Qc)
        cs = run(F, 'qrr', Qc, al, be, Qc, log, perturb={5: M5}, quiet=True)
        k1 = cs['q'][0]['K']
        check("control S (stored normalisation at z^5, f_5 -> f_5 + m5): K_1 = m5 = %s exactly, REG still passes (blind to it)" % M5,
              F.red(k1 - F.of(M5)) == 0 and all(r['reg'][0] and r['reg'][1] for r in cs['q']), short(F, k1))
        ct = run(F, 'qrr', min(Qc, 2), al, be, min(Qc, 2), log, tau_one=True, quiet=True)
        nb = sum(1 for r in ct['q'] if not (r['reg'][0] and r['reg'][1]))
        check("control T (tau -> 1, i.e. v = u, the wrong natural variable): REG fails for %d of %d orders" % (nb, len(ct['q'])), nb > 0)
        cw = run(F, 'qrr', Qc, al, be, Qc, log, h=E.qrr_h(F, 5 * Qc + 4, twist=False), quiet=True)
        hb = [m for (m, ok) in cw['H'] if not ok and m % 5]
        check("control W (QRR weights without -5 zeta^i): (H) fails at n = %s; K_1 = %s" % (hb[:4], short(F, cw['q'][0]['K'])), bool(hb))
    log("EV_REDUCED %s %s %s q<=%d: %s   (%.0f s)" % (F.name, mode, fr, Q, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED,
                                                   res['t']))
    OUT.close()


if __name__ == "__main__":
    main()
