# -*- coding: utf-8 -*-
"""rf_kq.py -- REFEREE: independent re-derivation of the reduced criterion of EVID.md s.2 (colour 0) and of its controls.
Uses only rf_core (no project or worker module).  For the WKB function Phi (rf_core.Wkb, built directly from the PF operator) with
free coefficients chosen by MODE (qrr: Phi(zhat, u = 1) = h at 5 | n;  can: the free coefficient fixed by K_q = 0):
  K11_q = [u^4q] Sum_(n<=5q) Qh_(1,n)(u) (-tau)^(5q-n) F_(2,5q-n),   K22_q = [u^4q] Sum Qh_(2,n) (-tau)^(5q-n) F_(1,5q-n)
(the colour-reduced (1,1), (2,2) slots of Qhat R^(tau zhat)^-1; the reduction itself is checked against a 5-colour matrix inverse in
rf_direct.py).  Checks: (H) Phi(1) = h at every 5 not | n; K11_q = 0 (qrr) or Phi^can(1) = h at 5 | n (can); K22_q = K11_q; REG:
lower u-orders of the (1,1) slot vanish; ph_n(0) = F_(0,n) (model); deg ph_n.  Controls (CTRL = 1): S f_5 -> f_5 + m5 (m5 from
rf_gauge, own stored normalisation), T tau = 1, W untwisted weights, B1 B_96 + 1, B2 B_8 + 1.
Run (cwd outside the project): py -I -B rf_kq.py FIELD MODE QMAX [alpha beta QREG Q22 CTRL]     FIELD Q | 89 | 61
log: logs/rf_kq_<FIELD>_<MODE>_<frame>_<5 QMAX>[_ctrl].log (relative paths only)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
PR = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def pair(F, Qa, Fb, tp, N, J, full=False):
    """[u^j] Sum_(n<=N) Qa[n] (-tau)^(N-n) F_b[N-n] for j = J (or all j <= J if full)"""
    out = [0] * (J + 1)
    for n in range(N + 1):
        m = N - n; c = Fb[m] if m % 2 == 0 else -Fb[m]
        if not c: continue
        a = Qa[n]; t = tp[m]
        for j in (range(J + 1) if full else (J,)):
            out[j] += c * sum(a[i] * t[j - i] for i in range(j + 1) if a[i])
    return [F.red(x) for x in out] if full else F.red(out[J])


def run(F, mode, Q, al=0, be=0, QREG=0, Q22=0, h=None, tau_one=False, shift=None, quiet=True, log=print):
    t0 = time.time(); KU = 4 * Q + 1; NZ = 5 * Q + 4
    h = C.qrr(F, NZ) if h is None else h
    C0, C1, tc, gam = C.periods(F, KU, al, be)
    tau = ([F.of(1)] + [0] * (KU - 1)) if tau_one else C.inv(F, C.itau(F, KU), KU)
    Fm = C.model(F, 5 * Q); iC0 = C.inv(F, C0, KU); iC1 = C.inv(F, C1, KU); k2 = C.mul(F, C.mul(F, C0, C0, KU), C.mul(F, C1, C1, KU), KU)
    tp = [[F.of(1)] + [0] * (KU - 1)]
    for m in range(1, 5 * Q + 1): tp.append(C.mul(F, tp[-1], tau, KU))
    if not quiet: log("   setup %.1f s; gamma (u^2 of the vanishing period) = %s, a_0 = %s" % (time.time() - t0, gam[1], gam[0]))
    W = C.Wkb(F, F.of(1)); rule = C.qrr_rule(F, h)
    Q0, Q1, Q2 = [], [], []

    def rows(n):
        del Q0[n:], Q1[n:], Q2[n:]
        Q0.append(C.mul(F, W.ph(n), iC0, KU))
        X = C.A_at(F, -1, 1, Q0[n], Q0[n - 1] if n else None, n, KU)
        Q1.append(C.mul(F, X, iC1, KU))
        if n <= 5 * Q22:
            Y = C.A_at(F, -2, 1, Q1[n], Q1[n - 1] if n else None, n, KU)
            Q2.append(C.mul(F, Y, k2, KU))
    rows(0); res = dict(q=[], H=[], deg=[], zero=[], t=0)
    for n in range(1, NZ + 1):
        q, r = divmod(n, 5)
        W.step(rule if (r == 0 and mode == "qrr") else None)
        if r == 0 and shift and n in shift: W.psi[n][-1] = F.red(W.psi[n].get(-1, 0) + F.of(shift[n]))
        if n <= 5 * Q: rows(n)
        if r == 0 and q <= Q:
            K = pair(F, Q1, Fm[2], tp, n, 4 * q); row = dict(q=q, K=K)
            if mode == "can":
                W.psi[n][-1] = F.red(W.psi[n].get(-1, 0) - K); rows(n); row['K_after'] = pair(F, Q1, Fm[2], tp, n, 4 * q)
                row['conn'] = F.red(sum(W.psi[n].values()) - h[n])
            row['f'] = W.psi[n].get(-1, 0)
            if q <= QREG: row['reg'] = not any(pair(F, Q1, Fm[2], tp, n, 4 * q - 1, full=True))
            if q <= Q22: row['K22'] = pair(F, Q2, Fm[1], tp, n, 4 * q)
            row['t'] = time.time() - t0; res['q'].append(row)
            if not quiet and (q <= 8 or q % 10 == 0): log("   q = %3d done, %.0f s" % (q, row['t']))
        ph = W.ph(n); res['H'].append((n, F.red(sum(ph) - h[n]) == 0)); res['deg'].append((n, len(ph) - 1))
        res['zero'].append(F.red(ph[0] - (Fm[0][n] if n <= 5 * Q else ph[0])) == 0)
    res['W'] = W; res['t'] = time.time() - t0
    return res


def short(x):
    s = str(x)
    return s if len(s) <= 46 else "%s..(%d/%d digits)" % (s[:14], len(str(Fr(x).numerator)), len(str(Fr(x).denominator)))


def main():
    a = sys.argv[1:]
    fld = a[0]; mode = a[1]; Q = int(a[2]); al = Fr(a[3]) if len(a) > 3 else Fr(0); be = Fr(a[4]) if len(a) > 4 else Fr(0)
    QREG = int(a[5]) if len(a) > 5 else min(Q, 12); Q22 = int(a[6]) if len(a) > 6 else min(Q, 12); ctrl = len(a) > 7 and a[7] == "1"
    F = C.QF() if fld == "Q" else C.PF(PR[fld])
    fr = "%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    name = "rf_kq_%s_%s_%s_%d%s.log" % (fld, mode, fr, 5 * Q, "_ctrl" if ctrl else "")
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)
    log("rf_kq (referee, independent code): field %s, mode %s, frame (%s, %s), q <= %d (zhat^%d), (H) n <= %d, REG q <= %d, K22 q <= %d"
        % (F.name, mode, al, be, Q, 5 * Q, 5 * Q + 4, QREG, Q22))
    log("   log file logs/%s" % name)
    res = run(F, mode, Q, al, be, QREG, Q22, quiet=False, log=log)
    for row in res['q']:
        if row['q'] <= 8 or row['q'] % 5 == 0:
            s = "   z^%-3d f = %s; K11 = %s" % (5 * row['q'], short(row['f']), short(row['K']))
            if mode == "can": s += "; Phi^can(1) - h = %s" % short(row['conn'])
            if 'K22' in row: s += "; K22 = %s" % short(row['K22'])
            if 'reg' in row: s += "; REG %s" % ("ok" if row['reg'] else "BAD")
            log(s)
    Hbad = [n for (n, ok) in res['H'] if n % 5 and not ok]
    check("(H) Phi(zhat, 1) = h at all %d orders n <= %d with 5 not | n (WKB from the PF operator, no free parameter)"
          % (sum(1 for (n, _) in res['H'] if n % 5), 5 * Q + 4), not Hbad, Hbad[:5])
    if mode == "qrr": bad = [r['q'] for r in res['q'] if r['K'] != 0]; nm = "(N) K11_q = 0 for q = 1..%d" % Q
    else: bad = [r['q'] for r in res['q'] if r['conn'] != 0 or r['K_after'] != 0]; nm = "connection form Phi^can(1) = h at n = 5q, q <= %d" % Q
    check(nm, not bad, bad[:5])
    check("K22_q = K11_q (the (2,2) slot, rows to a = 2) for q <= %d" % Q22,
          all(r['K22'] == r.get('K_after', r['K']) for r in res['q'] if 'K22' in r))
    check("REG: lower u-orders of the (1,1) slot vanish, q <= %d" % QREG, all(r['reg'] for r in res['q'] if 'reg' in r))
    check("ph_n(0) = F_(0,n) of the model for all n <= %d" % (5 * Q), all(res['zero'][:5 * Q]))
    degbad = [(n, d) for (n, d) in res['deg'] if n >= 3 and d != 4 * n // 5] + [(n, d) for (n, d) in res['deg'] if n < 3 and d != 0]
    check("deg ph_n = floor(4n/5) for 3 <= n <= %d, ph_1, ph_2 constant (%s, %s)" % (5 * Q + 4, W1(res, 1), W1(res, 2)), not degbad, degbad[:5])
    nres = len(res['W'].res)
    check("resonances unobstructed (%d of the %d multiples of 5 had a non-trivial rho at u^(4n/5))" % (nres, (5 * Q + 4) // 5),
          all(ok for (_, ok) in res['W'].res))
    if ctrl:
        sys.path.insert(0, HERE); import rf_gauge as G
        Qc = min(Q, 4); m5 = G.stored_m5(); log("-- controls (q <= %d); m5 from rf_gauge (own Theorem-U solver, exact): %s" % (Qc, m5))
        cs = run(F, "qrr", Qc, al, be, Qc, Qc, shift={5: m5})
        check("control S: f_5 -> f_5 + m5 gives K11_1 = K22_1 = m5 exactly, REG still passes",
              cs['q'][0]['K'] == F.of(m5) and cs['q'][0]['K22'] == F.of(m5) and all(r['reg'] for r in cs['q']), short(cs['q'][0]['K']))
        ct = run(F, "qrr", 2, al, be, 2, 0, tau_one=True)
        check("control T: tau = 1 breaks REG (%d of 2 orders)" % sum(1 for r in ct['q'] if not r['reg']), not all(r['reg'] for r in ct['q']))
        cw = run(F, "qrr", Qc, al, be, 0, 0, h=C.qrr(F, 5 * Qc + 4, twist=False))
        hb = [n for (n, ok) in cw['H'] if n % 5 and not ok]; check("control W: untwisted weights break (H), first at n = %s" % hb[:3], hb[:1] == [1])
        for idx, label, Qb in ((96, "B1", 20), (8, "B2", 2)):
            B = C.bernoulli(F, 5 * Qb + 5); B[idx] = F.red(B[idx] + 1); hb_ = C.qrr(F, 5 * Qb + 4, B=B)
            cb = run(F, "qrr", Qb, al, be, 0, 0, h=hb_)
            Hb = [n for (n, ok) in cb['H'] if n % 5 and not ok]; Kb = [r['q'] for r in cb['q'] if r['K'] != 0]
            if label == "B1":
                pred = F.of(Fr(1, 96 * 95) * Fr(1, 5 ** 95))
                check("control B1 (B_96 + 1): (H) intact, first K != 0 at q = %s, K_19 = 5^-95/9120 exactly" % Kb[:1],
                      not Hb and Kb[:1] == [19] and cb['q'][18]['K'] == pred, short(cb['q'][18]['K']))
            else: check("control B2 (B_8 + 1): (H) fails first at n = %s" % Hb[:1], Hb[:1] == [7])
    log("RF_KQ %s %s %s q<=%d: %s   (%.0f s)" % (F.name, mode, fr, Q, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, res['t']))
    OUT.close()


def W1(res, n): return res['W'].ph(n)


if __name__ == "__main__":
    main()
