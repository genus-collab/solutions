# -*- coding: utf-8 -*-
"""rf_pfree.py -- REFEREE observation: the conifold PERIODS DROP OUT of the normalisation functional K_q of EVID.md s.2.3.
Claim (REFEREE.md s.3): for ANY polynomial-type colour-0 WKB function Phi (any free coefficients),
      K_q = -[u^4q] S_q,   S_q(u) := Sum_(n<=5q) (A_-1 Phi)_n(u) (-tau)^(5q-n) F_(2,5q-n),   and S_q = O(u^4q),
where (A_-1 Phi)_n = ph_n + (1-u)((4n-3)/5 - theta_u) ph_(n-1).  Only Phi (the PF operator + the free coefficients), tau and the model row
F_2 enter; C0, C1 (hence the conifold frame and the vanishing period) do not.  Proof: Qh_1 = (1/(C0 C1)) [A_-1 Phi + zhat (1-u)
(theta C0/C0) Phi]; the second part pairs to (theta C0/C0) C0 Ghat_01[zhat^(5q-1)] = O(u) O(u^4q) because p_(5q-1,1)(0) = [z^(5q-1)]
P_01(0, z) = 0 by (R1) (Lemma N: Ghat = hatted P_can c); and C0(0) C1(0) = -1 in every conifold frame.
Test: perturb EVERY free coefficient (f_5q -> f_5q + 1/(q+2)) so that K_q != 0, and compare -[u^4q] S_q with K_q of rf_kq (rows with
the periods) in two frames; also the unperturbed QRR values (both 0) and REG' (S_q = O(u^4q)).
Run (cwd outside the project): py -I -B rf_pfree.py FIELD QMAX      log logs/rf_pfree_<FIELD>_<5 QMAX>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_core as C
import rf_kq as RK


def pfree(F, Q, shift=None):
    KU = 4 * Q + 1; NZ = 5 * Q; h = C.qrr(F, NZ); Fm = C.model(F, NZ); tau = C.inv(F, C.itau(F, KU), KU)
    tp = [[F.of(1)] + [0] * (KU - 1)]
    for m in range(1, NZ + 1): tp.append(C.mul(F, tp[-1], tau, KU))
    W = C.Wkb(F, F.of(1)); rule = C.qrr_rule(F, h); A = [[F.of(1)] + [0] * (KU - 1)]
    for n in range(1, NZ + 1):
        W.step(rule if n % 5 == 0 else None)
        if shift and n in shift: W.psi[n][-1] = F.red(W.psi[n].get(-1, 0) + F.of(shift[n]))
        P = (W.ph(n) + [0] * KU)[:KU]; Pm = (W.ph(n - 1) + [0] * KU)[:KU]; A.append(C.A_at(F, -1, 1, P, Pm, n, KU))
    out = []
    for q in range(1, Q + 1):
        S = RK.pair(F, A, Fm[2], tp, 5 * q, 4 * q, full=True); out.append((not any(S[:4 * q]), F.red(-S[4 * q])))
    return out


def main():
    fld = sys.argv[1]; Q = int(sys.argv[2]); F = C.QF() if fld == "Q" else C.PF(RK.PR[fld])
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True); name = "rf_pfree_%s_%d.log" % (fld, 5 * Q)
    OUT = open(os.path.join(HERE, "logs", name), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(nm, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(nm)
    t0 = time.time(); log("rf_pfree (referee): field %s, q <= %d; log logs/%s" % (F.name, Q, name))
    sh = {5 * q: Fr(1, q + 2) for q in range(1, Q + 1)}
    for shift, label in ((None, "QRR"), (sh, "perturbed f_5q + 1/(q+2)")):
        pf = pfree(F, Q, shift)
        Ks = [[r['K'] for r in RK.run(F, "qrr", Q, al, be, 0, 0, shift=shift)['q']] for (al, be) in ((0, 0), (Fr(1, 3), Fr(-2, 7)))]
        log("   %s: -[u^4q] S_q for q = 1..3: %s" % (label, [RK.short(x) for (_, x) in pf[:3]]))
        check("%s: REG' S_q = O(u^4q) for q <= %d (no periods used)" % (label, Q), all(r for (r, _) in pf))
        check("%s: -[u^4q] S_q = K_q (rf_kq, with C0, C1) in frames (0,0) and (1/3,-2/7), q <= %d" % (label, Q),
              [x for (_, x) in pf] == Ks[0] == Ks[1])
        check("%s: values %s" % (label, "all zero (D = 1)" if shift is None else "non-zero (the comparison is not vacuous)"),
              all(x == 0 for (_, x) in pf) if shift is None else all(x != 0 for (_, x) in pf))
    log("RF_PFREE %s q<=%d: %s   (%.0f s)" % (F.name, Q, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
