# -*- coding: utf-8 -*-
"""rf_deep.py -- two targeted checks (referee; own M_0 from the LG recursion in MVertex6, third prime p3; no worker module).

(1) The only singular list within MVertex6's deficit range (<= 6) at genus <= 1 whose frame depends on the i >= 1 terms of the vanishing
    series Y^C (rf_mla_scope.log):  <frame(6,0) H^5>_1,  frame(6,0) = z^6 phi_0 + v z^5 H - (v^5/120) H  (e = -1).
    ML_Lambda(1) predicts 0, i.e.  <tau_6(1) H^5>_1 + <tau_5(H) H^5>_1 = (1/120) <H^6>_1 = (1/120) F_1^(6)(1) = 1/12  (F_1 = -(1/12) log v).
    The frame built from the WRONG series (control C4: sign of the (-1/(4z))^i factor flipped) has +(v^5/120) H and predicts 1/6 != 0.
    <H^6>_1 (deficit 6, six legs) is not computed with MVertex6 here; it is taken from F_1 = -(1/12) log v (LG_MODEL.md s.5.4).
(2) One data point of the model worker's refutation of a genus-independent singular structure: Hess_h = <tau_1(H) tau_2(H^3)>_h / F_h'(1)
    (both insertions have PiC = 0), claimed -23/12 (h = 2) and -13/48 (h = 3) in ../../model/MODEL_SINGULAR_body.md s.4.
Run: py -I -B rf_deep.py   Log: rf_deep.log (appended).
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_mla_core import Corr, frame, frames_upto, w
from rf_mla import p3
from rf_flat_core import yc_monomials


def main():
    LOG = open(os.path.join(HERE, "rf_deep.log"), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    t0 = time.time(); C = Corr(p3()); R = C.R; p = C.p
    log("=" * 100); log("rf_deep [%s], p3 = %d" % (time.strftime("%Y-%m-%dT%H:%M:%S"), p))
    f60 = frame(6, 0)
    Ybad = [(b, j, n, cf * (-1) ** (n - j)) for (b, j, n, cf) in yc_monomials(60)]
    f60bad = frames_upto(w(6, 0), Y=Ybad)[(6, 0)]
    log("  frame(6,0) = %s;  C4 frame = %s" % (f60, f60bad))
    a = C.flat(1, [(6, 0)] + [(0, 1)] * 5); b = C.flat(1, [(5, 1)] + [(0, 1)] * 5)
    log("  <tau_6(1) H^5>_1 = %s,  <tau_5(H) H^5>_1 = %s,  sum = %s  (ML_Lambda(1) predicts 1/12)" % (C.fmt(a), C.fmt(b), C.fmt(a + b)))
    H6 = Fr(10)                                   # F_1^(6)(1) for F_1 = -(1/12) log v
    for lab, fr in (("ML_Lambda list <frame(6,0) H^5>_1", f60), ("control C4 (wrong series)", f60bad)):
        tot = 0
        for (k, c), (x, e) in fr.items():
            if (k, c) == (0, 1): val = R.of(H6)
            else: val = C.flat(1, [(k, c)] + [(0, 1)] * 5)
            tot = (tot + R.of(x) * val) % p
        log("  %-36s = %s" % (lab, C.fmt(tot)))
    for h, claim in ((2, Fr(-23, 12)), (3, Fr(-13, 48))):
        val = C.flat(h, [(1, 1), (2, 3)])
        Fp = Fr(2 - 2 * h) * {2: Fr(1, 1200), 3: Fr(1, 25200)}[h]
        hess = R.ratrec(val * R.inv(R.of(Fp)) % p)
        log("  genus %d: <tau_1(H) tau_2(H^3)>_h = %s, F_h'(1) = %s, Hess_h = %s (model worker: %s) %s"
            % (h, C.fmt(val), Fp, hess, claim, "agree" if hess == claim else "DIFFER"))
    log("  (%.0f s)" % (time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
