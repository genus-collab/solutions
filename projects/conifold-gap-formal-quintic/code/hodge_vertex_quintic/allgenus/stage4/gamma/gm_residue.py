# -*- coding: utf-8 -*-
"""gm_residue.py -- the residue at the pole of the full-Stirling 1-form is the quintic's vanishing solution (normalisation K = 1).

One-dimensional form of the full-Stirling Mellin-Barnes integrand (x = 1/sigma, GAMMA.md s.3):
   Omega = exp(W(x; u)/z) A(x; z) dx,   W = log(1-u)/x + G(x),  G(x) = Int_0^x log(1 - t^5) t^-2 dt = -Sum_l x^(5l+4)/((l+1)(5l+4)),
   A = (1 - x^5)^(-1/2) exp( Sum_k b_k z^(2k-1) x^(2k-1) [5^(1-2k) + Sum_j (zeta^j x - 1)^(1-2k)] ),   b_k = B_2k/(2k(2k-1)).
psi^MB := Res_{x=0} Omega = Sum_{k>=1} (log(1-u))^k/k! z^-k [x^(k-1)] (e^(G/z) A)   (exact over Q, order by order in u).
Claims tested (frame (0,0)): (R1) psi^MB has only negative powers of z; (R2) psi^MB = C_0 Y^C_Q,0 (identity row of the quintic's vanishing
solution, built by iterated integrals from the periods) for every z^-n, n <= NZ, mod u^NU; (R3) row 1: zD(psi^MB/C_0)/C_1 = Y^C_Q,1
(value -1 at u = 0: K = 1); (R4) [z^-1] psi^MB = -C_0 t_c = Sum_k alpha_(k-1) (log(1-u))^k/k!, alpha_m = [y^m] A_top,
A_top(y) = exp(Sum_k b_k (5^(1-2k) - 5) y^(2k-1)) (the Borel transform of the conifold Stirling factor, integrated);
control: dropping the Stirling corrections (A -> (1-x^5)^(-1/2)) fails (R2).
Read-only imports: natural_form/theory/ur_gauge.coefficients (b_a, varpi_0, t_c) and allgenus/gauge/gi_core.quintic_vanishing.
Run: py -I -B gm_residue.py [NU NZ [alpha beta]]   (default 14 10, frame (0,0))   log: logs/gm_residue_<NU>_<NZ>[_ab_<frame>].log
A frame (alpha, beta) != (0,0) is selected by rebinding ur_gauge.Conifold IN MEMORY, as allgenus/gauge/gi_core.gauge_frame does.
"""
import os, sys, time
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gm_series import Rq, bernoulli_list, smul, sexp0, spow
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
if not os.path.isdir(os.path.join(PROJ, "natural_form")):                 # scratch copy: fall back to the project folder
    PROJ = os.environ.get("FQ_PROJ", PROJ)
sys.path.insert(0, os.path.join(PROJ, "allgenus", "gauge"))
R = Rq()


def A_bivariate(NX, stirling=True):
    """A(x; z) as {power of z: x-series of length NX}"""
    B = bernoulli_list(2 * NX + 2)
    expo = {}                                                          # z^j -> x-series
    if stirling:
        for k in range(1, NX + 1):
            j = 2 * k - 1
            if j >= NX: break
            bk = B[2 * k] / (2 * k * (2 * k - 1)); ser = [Fr(0)] * NX
            for r in range(0, NX - j):
                c = (Fr(5) ** (1 - 2 * k) if r == 0 else Fr(0)) - (5 * comb(2 * k - 2 + r, r) if r % 5 == 0 else 0)
                ser[j + r] = bk * c                                    # x^(2k-1) x^r
            expo[j] = ser
    base = spow(R, [Fr(1)] + [Fr(0)] * 4 + [Fr(-1)] + [Fr(0)] * (NX - 6), Fr(-1, 2))[:NX]
    out = {0: base}                                                    # exp of Sum_j z^j expo[j] (a series in z)
    E = {0: [Fr(1)] + [Fr(0)] * (NX - 1)}
    for n in range(1, NX):
        acc = [Fr(0)] * NX
        for j, s in expo.items():
            if j <= n:
                prod = smul(R, s, E[n - j], NX); acc = [a + j * b for a, b in zip(acc, prod)]
        E[n] = [a / n for a in acc]
    return {n: smul(R, base, E[n], NX) for n in E if any(E[n])}


def psi_MB(NU, stirling=True):
    """{-n: list of NU Fractions (u-series)}"""
    NX = NU + 1
    Gs = [Fr(0)] * NX
    for l in range(0, NX):
        if 5 * l + 4 < NX: Gs[5 * l + 4] = Fr(-1, (l + 1) * (5 * l + 4))
    Gp = {0: [Fr(1)] + [Fr(0)] * (NX - 1)}                             # e^(G/z) = Sum_j G^j/j! z^-j
    for j in range(1, NX // 4 + 1): Gp[j] = [x / j for x in smul(R, Gp[j - 1], Gs, NX)]
    A = A_bivariate(NX, stirling)
    ell = [Fr(0)] + [Fr(-1, n) for n in range(1, NU)]                  # log(1-u)
    ellk = [Fr(1)] + [Fr(0)] * (NU - 1); out = {}
    for k in range(1, NU):
        ellk = [x / k for x in smul(R, ellk, ell, NU)]                 # ell^k / k!
        for j, g in Gp.items():
            for i, a in A.items():
                c = sum(g[r] * a[k - 1 - r] for r in range(k) if g[r] and a[k - 1 - r])
                if not c: continue
                zp = -k - j + i
                acc = out.setdefault(zp, [Fr(0)] * NU)
                for N in range(NU): acc[N] += c * ellk[N]
    return {m: s for m, s in out.items() if any(s)}


def main():
    NU = int(sys.argv[1]) if len(sys.argv) > 1 else 14; NZ = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    al = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0); be = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    tag = "" if (al, be) == (0, 0) else "_ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    OUT = open(os.path.join(HERE, "logs", "gm_residue_%d_%d%s.log" % (NU, NZ, tag)), "w", encoding="ascii"); FAIL = []
    def log(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(nm, ok, det=""):
        log(("PASS " if ok else "FAIL ") + nm + (": " + str(det) if det != "" else ""))
        if not ok: FAIL.append(nm)
    t0 = time.time(); log("gm_residue: mod u^%d, z^-n for n <= %d, exact over Q, frame (%s, %s)" % (NU, NZ, al, be))
    import gi_core, ur_gauge, fq_conifold
    orig = ur_gauge.Conifold
    ur_gauge.Conifold = lambda Rr, Nu=40: fq_conifold.Conifold(Rr, Nu=Nu, alpha=al, beta=be)
    try:
        bq, bt, tau, v, C = ur_gauge.coefficients(NU + 4)
    finally:
        ur_gauge.Conifold = orig
    assert C.alpha == al and C.beta == be
    YQ = gi_core.quintic_vanishing(bq, NU + 2, NZ + 2)
    C0 = [C.varpi0.coef(N) for N in range(NU)]; tc = [C.tc.coef(N) for N in range(NU + 2)]
    psi = psi_MB(NU)
    check("(R1) psi^MB has only negative powers of z", all(m < 0 for m in psi), sorted(psi)[-3:])
    bad = []
    for n in range(1, NZ + 1):
        y0 = [YQ[-n][0].coef(N) if -n in YQ else Fr(0) for N in range(NU)] if -n in YQ else [Fr(0)] * NU
        rhs = smul(R, C0, y0, NU); lhs = psi.get(-n, [Fr(0)] * NU)
        if lhs != rhs: bad.append((n, [N for N in range(NU) if lhs[N] != rhs[N]][:3]))
    check("(R2) Res_{x=0} Omega = C_0 Y^C_Q,0 (C_0 = varpi_0 of the frame) for z^-n, n = 1..%d, mod u^%d" % (NZ, NU), not bad, bad[:3])
    # row 1: Y_1 = zD(psi/C_0)/C_1, D = -(1-u) d/du, C_1 = -(1-u) t_c'
    C0i = [Fr(0)] * NU; C0i[0] = 1 / C0[0]
    for k in range(1, NU): C0i[k] = -sum(C0[j] * C0i[k - j] for j in range(1, k + 1)) / C0[0]
    tcp = [(N + 1) * tc[N + 1] for N in range(NU)]; C1 = [-(tcp[N] - (tcp[N - 1] if N else 0)) for N in range(NU)]
    C1i = [Fr(0)] * NU; C1i[0] = 1 / C1[0]
    for k in range(1, NU): C1i[k] = -sum(C1[j] * C1i[k - j] for j in range(1, k + 1)) / C1[0]
    bad1 = []
    for n in range(0, NZ):
        f = smul(R, psi.get(-(n + 1), [Fr(0)] * NU), C0i, NU)           # z D f z^-(n+1) -> z^-n
        Df = [-((N + 1) * f[N + 1] - N * f[N]) if N + 1 < NU else None for N in range(NU)]
        Df = [x for x in Df if x is not None]
        y1 = smul(R, Df, C1i[:len(Df)], len(Df))
        ref = [YQ[-n][1].coef(N) if -n in YQ else Fr(0) for N in range(len(y1))]
        if y1 != ref: bad1.append(n)
    check("(R3) row 1: zD(psi^MB/C_0)/C_1 = Y^C_Q,1 for z^-n, n = 0..%d (value -1 at u = 0, so K = 1)" % (NZ - 1), not bad1, bad1[:3])
    By = bernoulli_list(2 * NU + 2); ex = [Fr(0)] * NU
    for k in range(1, NU):
        if 2 * k - 1 < NU: ex[2 * k - 1] = By[2 * k] / (2 * k * (2 * k - 1)) * (Fr(5) ** (1 - 2 * k) - 5)
    At = sexp0(R, ex); ell = [Fr(0)] + [Fr(-1, n) for n in range(1, NU)]
    ellk = [Fr(1)] + [Fr(0)] * (NU - 1); integ = [Fr(0)] * NU
    for k in range(1, NU):
        ellk = [x / k for x in smul(R, ellk, ell, NU)]; integ = [a + At[k - 1] * b for a, b in zip(integ, ellk)]
    g = smul(R, C0, [YQ[-1][0].coef(N) for N in range(NU)], NU)
    check("(R4) vanishing period: -C_0 t_c = [z^-1] psi^MB = Sum_k [y^(k-1)]A_top (log(1-u))^k/k!, mod u^%d" % NU, integ == g and psi[-1] == g,
          [str(x) for x in g[:4]])
    ctl = psi_MB(NU, stirling=False)
    y0 = [YQ[-1][0].coef(N) for N in range(NU)]
    check("(C) control: without the Stirling corrections the residue is NOT the vanishing solution", ctl.get(-1) != smul(R, C0, y0, NU))
    log("GM_RESIDUE %d %d: %s   (%.0f s)" % (NU, NZ, "ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL, time.time() - t0))
    OUT.close()


if __name__ == "__main__":
    main()
