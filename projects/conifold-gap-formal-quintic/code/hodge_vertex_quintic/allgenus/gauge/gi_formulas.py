# -*- coding: utf-8 -*-
"""gi_formulas.py -- exact checks over Q of the closed formulas used in the proofs of GAUGE_IDENTITY.md (s.2, s.5, s.7).

 (V)  valuations: ord_u [z^-n] Y^C_{Q,a}(u) = (4n + c_a + 1)/5 EXACTLY for every non-zero slot, and the slot is zero iff n != 1 - a mod 5
      (quintic, iterated integrals with zero constants), and the same for the model Y^C(v(u)); leading coefficients agree.
 (D)  the diagonal of the canonical gauge: p_{0,a}(u) = p_{0,0}(u) prod_{j<a} bt_j/b_j and p_{0,0}(u) = tau^(-1/4)/varpi_0
      (so Ghat(u, 0)_00 = 1/xhat and the hatted unit row is Q'^(i)_{0,0} = zeta^i/xhat).
 (F)  the vanishing covector in closed form:  (-1)^m Sum_b (Psi0^-1)_{ib} [t^-m] y_b(t) = (1/5) zeta^(i(m-2)) / (4^q q! (m-q)!),  q = floor(m/5),
      checked as: exactly one b contributes (b = 1 - m mod 5), Psi0^-1[i][b] = (1/5) eps_b zeta^(-i(b+1)) = (1/5) eps_b zeta^(i(m-2)), and
      (-1)^m eps_b [t^-m] y_b = 1/(4^q q! (m-q)!)  for m <= 80.
Run: py -I -B gi_formulas.py      log: gi_formulas.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import Q, Ser, gauge_frame, model_vanishing, quintic_vanishing, fact
from gi_hat_core import y_coeffs, CEXP, EPS
OUT = open(os.path.join(HERE, "gi_formulas.log"), "w", encoding="ascii"); FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
    if not ok: FAILED.append(name)


def main():
    t0 = time.time()
    for al, be in ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7))):
        NU = 40; NMAX = 45
        res = gauge_frame(0, NU + 6, al, be); bq = res['bq']; bt = res['bt']; v = res['v']; tau = res['tau']; C = res['C']
        YQ = quintic_vanishing(bq, NU, NMAX); YC = model_vanishing(v, NU, NMAX)
        bad = []; tested = 0; lead = []
        for n in range(NMAX + 1):
            for a in range(5):
                q_ = YQ.get(-n); m_ = YC.get(-n)
                sq = q_[a] if q_ is not None else Ser.zero(Q, NU); sm = m_[a] if m_ is not None else Ser.zero(Q, NU)
                w5 = 4 * n + CEXP[a] + 1
                graded = (n - (1 - a)) % 5 == 0
                if not graded:
                    if sq.c or sm.c: bad.append(("non-graded slot non-zero", n, a))
                    continue
                w = w5 // 5; assert w5 % 5 == 0
                if w >= NU - 1: continue
                tested += 1
                if not sq.c or sq.v != w: bad.append(("quintic valuation", n, a, sq.v if sq.c else None, w))
                if not sm.c or sm.v != w: bad.append(("model valuation", n, a, sm.v if sm.c else None, w))
                if sq.c and sm.c and sq.coef(w) != sm.coef(w): bad.append(("leading coefficient", n, a))
                if n <= 6: lead.append((n, a, w, str(sq.coef(w))))
        log("frame (%s, %s): leading terms [z^-n] Y^C_Q,a = coefficient u^w, (n, a, w, coefficient) for n <= 6: %s" % (al, be, lead))
        check("(V) frame (%s, %s): ord_u [z^-n] Y^C_{Q,a} = (4n + c_a + 1)/5 exactly, zero off the grading n = 1 - a mod 5, same valuation and leading "
              "coefficient for the model Y^C(v(u)) (%d graded slots, n <= %d)" % (al, be, tested, NMAX), not bad, bad[:4])
        p0 = res['p'][0]; P = NU
        okr = all(p0[a].truncate(P - 2).sub(p0[a - 1].mul(bt[a - 1]).mul(bq[a - 1].inv()).truncate(P - 2)).truncate(P - 2).iszero() for a in range(1, 5))
        lt = tau.log(); target = lt.scale(Fr(-1, 4)).exp().mul(C.varpi0.truncate(P).inv()).truncate(P - 2)
        check("(D) frame (%s, %s): p_{0,a} = p_{0,0} prod_{j<a} bt_j/b_j and p_{0,0} = tau^(-1/4)/varpi_0 through u^%d (P_11(0, 0) = %s)" % (al, be, P - 3, p0[1].coef(0)),
              okr and p0[0].truncate(P - 2).sub(target).truncate(P - 2).iszero() and p0[1].coef(0) == 1)
    y = y_coeffs(80); bad = []
    for m in range(81):
        bs = [b for b in range(5) if y[b].get(m, 0) != 0]
        q = m // 5
        if bs != [(1 - m) % 5]: bad.append(("support", m, bs)); continue
        b = bs[0]
        if ((b + 1) - (2 - m)) % 5: bad.append(("phase", m))
        val = (1 if m % 2 == 0 else -1) * EPS[b] * y[b][m]
        if val != Fr(1, 4 ** q * fact(q) * fact(m - q)): bad.append(("value", m, str(val)))
    check("(F) phi_{i,m} = (u tau^(-1/2)/5) tau^(-m) zeta^(i(m-2)) / (4^q q! (m-q)!), q = floor(m/5): support, phase and value for m <= 80", not bad, bad[:4])
    log("GI_FORMULAS: %s   (%.0f s)" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
