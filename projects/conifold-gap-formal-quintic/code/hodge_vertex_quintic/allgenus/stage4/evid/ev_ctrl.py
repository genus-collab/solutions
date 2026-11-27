# -*- coding: utf-8 -*-
"""ev_ctrl.py -- Stage 4 evidence: detection power of the reduced criterion at HIGH order (EVID.md s.5).  Uses ev_core, ev_reduced.
Perturb one Bernoulli number inside the QRR constants h_{w_0} and rerun the reduced test (frame (0,0), q <= 20):
  B1: B_96 -> B_96 + 1.  This multiplies h by exp(dl z^95), dl = -S_95/(96*95), S_m = Sum_{a in w_0} a^-m, i.e. it is a change of
      normalisation by a series in z^5.  Prediction: (H) unchanged at every n; K_q = 0 for q <= 18; K_19 = dl EXACTLY (slope 1).
      So the z^95 test (a_19, genus 33) sees B_96 and nothing else does.
  B2: B_8 -> B_8 + 1 (the z^7 term).  Prediction: (H) fails first at n = 7 (not a multiple of 5).
Run (cwd outside the project): py -I -B ev_ctrl.py [FIELD]   FIELD Q | 89 | 61 (default 89);  log logs/ev_ctrl_<FIELD>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ev_core as E
import ev_reduced as RD


def S(F, m):
    p = [F.of(4), F.of(2), F.of(0), F.of(-1), F.of(Fr(-4, 5))]
    for k in range(5, m + 1): p.append(F.red(2 * p[k - 1] - 2 * p[k - 2] + p[k - 3] - F.of(Fr(1, 5)) * p[k - 4]))
    return F.red(p[m] + F.of(Fr(-1, 5) ** m))


def main():
    fld = sys.argv[1] if len(sys.argv) > 1 else "89"
    F = E.QF() if fld == "Q" else E.PF(RD.PR[fld]); Q = 20
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "ev_ctrl_%s.log" % fld
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)

    t0 = time.time(); log("ev_ctrl: field %s, frame (0,0), q <= %d; log file logs/%s" % (F.name, Q, name))
    E.bern(5 * Q + 6)
    for (idx, label) in ((96, "B1"), (8, "B2")):
        old = E._B[idx]; E._B[idx] = old + 1
        try:
            h = E.qrr_h(F, 5 * Q + 4)
        finally:
            E._B[idx] = old
        res = RD.run(F, 'qrr', Q, 0, 0, 0, log, h=h, quiet=True)
        Hbad = [m for (m, ok) in res['H'] if not ok and m % 5]; Kbad = [r['q'] for r in res['q'] if r['K'] != 0]
        if label == "B1":
            dl = F.red(-S(F, 95) * F.of(Fr(1, 96 * 95)))
            k19 = res['q'][18]['K']
            log("   B1: first K_q != 0 at q = %s; K_19 = %s; predicted dl = -S_95/(96*95) = %s" % (Kbad[:1], RD.short(F, k19), RD.short(F, dl)))
            check("B1 (B_96 + 1): (H) unchanged at all n <= %d, K_q = 0 for q <= 18, K_19 = dl exactly (the a_19 test sees B_96)" % (5 * Q + 4),
                  not Hbad and Kbad[:1] == [19] and F.red(k19 - dl) == 0)
        else:
            check("B2 (B_8 + 1): (H) fails first at n = 7 (found %s)" % Hbad[:3], Hbad[:1] == [7])
    log("EV_CTRL %s: %s   (%.0f s)" % (F.name, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
