# -*- coding: utf-8 -*-
"""rf_run_res.py -- REFEREE driver for Theorem A(i) and the vanishing-period closed form (own code: rf_ring, rf_res).
Criterion (no use of the project's Y^C_Q): psi := Res_{x=0} Omega; rows Y_a from psi with the frame periods; if Y solves the cyclic
system (closing relation), has no positive power of z and Y(0) = -e_1, then Y = Y^C_Q by uniqueness of regular solutions
(GAUGE_IDENTITY Lemma 1, a 3-line induction re-checked in REFEREE.md s.3), i.e. Res Omega = C0 Y^C_Q,0 with K = 1.
Run: py -I -B rf_run_res.py NU   (working directory outside the project)   log: logs/rf_run_res_<NU>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_ring import mul, add, scal, log1, exp0, stirling_b
from rf_res import periods, psi_res, rows, frame, F

FRAMES = ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5)))


def zero_mod(f, n): return not any(f[:n])


def analyse(psi, C0, C1, NU):
    """returns (closing ok, no positive z-powers, Y(0) == -e_1, details)"""
    Y, close = rows(psi, C0, C1, NU)
    prec = [NU, NU - 1, NU - 1, NU - 2, NU - 3]
    ok_close = all(zero_mod(f, NU - 4) for f in close.values())
    pos = [(a, e) for a in range(5) for e, f in Y[a].items() if e > 0 and not zero_mod(f, prec[a])]
    y0 = {(a, e): f[0] for a in range(5) for e, f in Y[a].items() if f and f[0] != 0}
    return ok_close, not pos, y0 == {(1, 0): Fr(-1)}, dict(pos=pos[:4], y0=y0, tc=[-x for x in Y[0].get(-1, [])[:5]])


def main():
    NU = int(sys.argv[1]) if len(sys.argv) > 1 else 26
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "rf_run_res_%d.log" % NU), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("rf_run_res: exact over Q, mod u^%d (closing relation mod u^%d), all powers of z at these orders" % (NU, NU - 4))
    f0, f1, f2, g, gamma = periods(NU + 2)
    check("(P0) own Frobenius for L4: gamma FORCED by the log partner = 7/10; g = u + 7/10 u^2 + 41/75 u^3 + ...",
          gamma == Fr(7, 10) and g[:4] == [0, 1, Fr(7, 10), Fr(41, 75)], [str(x) for x in g[:6]])
    lu = log1(F, [Fr(1), Fr(-1)] + [Fr(0)] * NU, NU)
    top = [Fr(0)] * NU
    for k in range(1, NU // 2 + 1):
        if 2 * k - 1 < NU: top[2 * k - 1] = stirling_b(k) * (Fr(1, 5 ** (2 * k - 1)) - 5)
    alpha = exp0(F, top, NU)                                              # A_top(y) = Stirling factor of Gamma(5w)/Gamma(w)^5
    form = [Fr(0)] * NU; lk = [Fr(1)] + [Fr(0)] * (NU - 1); fact = 1
    for k in range(1, NU):
        lk = mul(F, lk, lu, NU); fact *= k; form = add(F, form, scal(F, lk, -alpha[k - 1] / fact))
    check("(P1) vanishing period closed form g = -Sum_k alpha_(k-1) (log(1-u))^k/k!, alpha = [y^m] exp(Sum b_k (5^(1-2k) - 5) y^(2k-1)), mod u^%d" % NU,
          form == g[:NU], [str(x) for x in alpha[:4]])
    psi = psi_res(NU)
    check("(R0) psi = Res Omega: only z^-n with n = 1 mod 5 (Z/5 grading of row 0), no z^0 or positive power",
          all(e < 0 and (-e) % 5 == 1 for e in psi), sorted(psi)[-4:])
    check("(R1) [z^-1] psi = -g (= -C0 t_c in every frame) mod u^%d" % NU, psi[-1] == [-x for x in g[:NU]])
    for (al, be) in FRAMES:
        C0, C1 = frame(f0[:NU], f2[:NU], g[:NU], al, be, NU)
        assert C0[0] == 1 and C1[0] == -1
        okc, okp, oky, det = analyse(psi, C0, C1, NU)
        check("(R2) frame (%s, %s): rows of psi/C0 satisfy the cyclic system incl. the closing relation 4 -> 0" % (al, be), okc)
        check("(R2) frame (%s, %s): no positive power of z in any row (within precision)" % (al, be), okp, det['pos'])
        check("(R2) frame (%s, %s): Y(0) = -e_1 (row 1 = -1 at u = 0: K = 1), so Res Omega = C0 Y^C_Q,0 (Lemma 1)" % (al, be), oky, det['y0'])
    log("     t_c in the last frame (7/11, -3/5) = %s" % [str(x) for x in det['tc']])
    NG = 8; f0g, f1g, f2g, gg, _ = periods(NG + 2); psig = psi_res(NG); bad = []           # all frames: Y(0) has degree <= 2
    for al in (Fr(0), Fr(1, 3), Fr(-2)):                                                  # in (alpha, beta) (REFEREE.md s.3.6)
        for be in (Fr(0), Fr(-2, 7), Fr(5)):
            C0, C1 = frame(f0g[:NG], f2g[:NG], gg[:NG], al, be, NG)
            if not analyse(psig, C0, C1, NG)[2]: bad.append((al, be))
    check("(R3) Y(0) = -e_1 on the 3 x 3 grid alpha in {0, 1/3, -2}, beta in {0, -2/7, 5}: hence for EVERY conifold frame", not bad, bad)
    C0, C1 = frame(f0[:NU], f2[:NU], g[:NU], 0, 0, NU)
    psi0 = psi_res(NU, stirling=False); okc, okp, oky, det = analyse(psi0, C0, C1, NU)
    check("(C1) control: without the Stirling corrections the criterion FAILS", not (okc and okp and oky),
          "closing %s, no-positive %s, Y(0) %s" % (okc, okp, oky))
    C0b, C1b = frame(f0[:NU], f2[:NU], add(F, g[:NU], scal(F, f2[:NU], Fr(1, 7))), 0, 0, NU)   # wrong gamma: 7/10 + 1/7
    okc, okp, oky, det = analyse(psi, C0b, C1b, NU)
    check("(C2) control: periods with a wrong vanishing period (gamma + 1/7) -> criterion FAILS", not (okc and okp and oky),
          "closing %s, no-positive %s, Y(0) %s" % (okc, okp, oky))
    okc, okp, oky, det = analyse({e: scal(F, f, 2) for e, f in psi.items()}, C0, C1, NU)
    check("(C3) control: residue normalisation K = 2 -> Y(0) = -2 e_1 (criterion FAILS)", not oky, det['y0'])
    log("RF_RUN_RES %d: %s   (%.0f s)" % (NU, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
