# -*- coding: utf-8 -*-
"""gi_check.py -- exact check over Q of the gauge identity (G) of GAUGE_IDENTITY.md s.2:

        P(u, z) Y^C(v(u), z) = P_11(0, z) Y^C_Q(u, z)      in C((z))[[u]]^5,

P the regular gauge of ur_gauge.gauge (z^0..z^K), Y^C the residue (vanishing-cycle) solution of the model, Y^C_Q the solution of
the quintic system with value -e_1 at u = 0 (iterated integrals, zero constants; a series in 1/z).
 (a) after P -> P c(z), c = 1/P_11(0, z): every coefficient of a POSITIVE power of z of P Y^C vanishes;
 (b) the non-positive part equals Y^C_Q, which satisfies the quintic system and has -t_c at z^-1 in component 0.
A coefficient [z^m u^n] is tested only if it is complete: m + floor(5n/4) <= K (gi_core.complete) and n < NUT.
Component a of P Y^C contains only powers z^m with m = a - 1 mod 5 (Z/5 grading); other slots are 0 = 0 and are NOT counted.
Run:  py -I -B gi_check.py [K] [NU] [alpha] [beta]        (defaults 16, 50, 0, 0)   log: gi_check_<frame>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import (Q, Ser, gauge_frame, model_vanishing, quintic_vanishing, apply_gauge, mul_scalar_z, p11_at_zero, inv_poly,
                     complete, coef_or_none, normalise)


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    NU = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    al = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0)
    be = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "gi_check_%s.log" % tag), "w", encoding="ascii"); FAILED = []

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
        if not ok: FAILED.append(name)

    t0 = time.time()
    log("gi_check: frame (alpha, beta) = (%s, %s), gauge to z^%d, working precision O(u^%d); exact rationals" % (al, be, K, NU))
    res = gauge_frame(K, NU, al, be); p = res['p']; bq = res['bq']; v = res['v']; C = res['C']
    NUT = min(min(x.prec for x in pk) for pk in p) - 1        # safe common u-precision of the gauge
    NUT = min(NUT, v.prec - 2)
    log("   gauge computed in %.0f s; common precision of the p_{k,b}: O(u^%d)" % (time.time() - t0, NUT))
    log("   varpi_0 = %s" % C.varpi0.fmt(3)); log("   t_c     = %s" % C.tc.fmt(4))
    hol = all((not x.c) or x.v >= 0 for pk in p for x in pk)
    check("the gauge is regular: every p_{k,b} is a power series in u (k <= %d)" % K, hol)
    c11 = p11_at_zero(p)
    log("   P_11(0, z) = " + " + ".join("(%s) z^%d" % (c11[k], k) for k in sorted(c11)))
    check("(R1) column 1 of P(0, z) is P_11(0, z) e_1 (p_{k,1}(0) = 0 unless 5 | k)", all(p[k][1].coef(0) == 0 for k in range(K + 1) if k % 5))
    stored_canonical = all(c11[k] == (1 if k == 0 else 0) for k in c11)
    log("   the STORED gauge %s in the canonical normalisation P_11(0, z) = 1" % ("IS" if stored_canonical else "is NOT"))
    cinv = inv_poly(c11, K)
    log("   c(z) = 1/P_11(0, z) = " + " + ".join("(%s) z^%d" % (cinv[k], k) for k in sorted(cinv)))
    # ---- Y^C(v(u)), P Y^C
    nmaxY = (5 * NUT) // 4 + 2
    YC = model_vanishing(v, NUT, nmaxY)
    mlo = -nmaxY
    raw = apply_gauge(p, YC, mlo, K, NUT)                       # P Y^C, stored normalisation
    W = mul_scalar_z(cinv, raw, mlo, K, NUT)                    # (P c) Y^C, canonical normalisation
    # ---- (a) positive powers
    for name, vec in (("STORED normalisation (control: must FAIL unless the stored gauge is canonical)", raw), ("CANONICAL normalisation P_11(0, z) = 1", W)):
        tested = 0; nonzero = []; trivial = 0
        for m in range(1, K + 1):
            for a in range(5):
                if (m - (a - 1)) % 5: continue
                for n in range(NUT):
                    if not complete(m, n, K): continue
                    tested += 1
                    val = coef_or_none(vec[m][a], n) if m in vec else Fr(0)
                    assert val is not None
                    if val != 0: nonzero.append((m, a, n, val))
                    # how many terms p_{k,b} * Y^C_b[z^(m-k)] can contribute (a = b + k): a sum with >= 2 terms is a cancellation
        if name.startswith("STORED"):
            log("   %s: %d graded complete coefficients [z^m u^n], m >= 1; non-zero: %d; first: %s" % (name, tested, len(nonzero), [(m, a, n, str(x)) for (m, a, n, x) in nonzero[:3]]))
            check("control: in the stored normalisation the positive powers vanish iff the stored gauge is canonical", (len(nonzero) == 0) == stored_canonical)
        else:
            check("(a) %s: ALL %d complete graded coefficients [z^m u^n] (m = 1..%d) of the positive part of P Y^C vanish" % (name, tested, K), not nonzero,
                  [(m, a, n, str(x)) for (m, a, n, x) in nonzero[:4]])
            log("       tested slots by z-power m: %s" % {m: sum(1 for a in range(5) if (m - a + 1) % 5 == 0 for n in range(NUT) if complete(m, n, K)) for m in range(1, K + 1)})
    # how many of the tested slots are genuine cancellations in the CANONICAL gauge pc = p c(z): number of non-zero products
    # pc_{k,b}[u^n1] * Y^C_b[z^(m-k)][u^n2] entering each slot (0: empty slot, 1 would be a contradiction, >= 2: a cancellation)
    pc = normalise(p, cinv)
    Wd = apply_gauge(pc, YC, 1, K, NUT)
    check("(a') the same from the normalised coefficients pc_{k,b} = Sum_j c_j p_{k-j,b} directly", all(
        (m not in Wd) or all(coef_or_none(Wd[m][a], n) == 0 for a in range(5) for n in range(NUT) if complete(m, n, K)) for m in range(1, K + 1)))
    genuine = 0; single = 0; empty = 0; maxterms = 0
    for m in range(1, K + 1):
        for a in range(5):
            if (m - (a - 1)) % 5: continue
            for n in range(NUT):
                if not complete(m, n, K): continue
                cnt = 0
                for k in range(K + 1):
                    src = YC.get(m - k)
                    if src is None: continue
                    b = (a - k) % 5
                    if not src[b].c or not pc[k][b].c: continue
                    for n1 in range(pc[k][b].v, n + 1):
                        if n - n1 < src[b].v: break
                        if pc[k][b].coef(n1) != 0 and src[b].coef(n - n1) != 0: cnt += 1
                maxterms = max(maxterms, cnt)
                if cnt >= 2: genuine += 1
                elif cnt == 1: single += 1
                else: empty += 1
    log("       canonical-gauge products entering the %d tested slots: %d slots are sums of >= 2 non-zero products (genuine cancellations, up to %d terms), %d have exactly one, %d are empty (0 = 0)" % (genuine + single + empty, genuine, maxterms, single, empty))
    check("no tested slot consists of a single non-zero product", single == 0)
    # ---- (b) the non-positive part = Y^C_Q
    YQ = quintic_vanishing(bq, NUT, nmaxY)
    tc = C.tc
    y1 = YQ[-1][0]
    check("(b0) Y^C_Q = -e_1 - t_c e_0/z + ...: the z^-1 coefficient is -t_c e_0 (iterated integral vs fq_conifold t_c, O(u^%d))" % min(NUT, tc.prec),
          y1.add(tc).truncate(min(NUT, tc.prec) - 1).iszero() and all(not YQ[-1][a].c for a in range(1, 5)))
    tested = 0; bad = []
    for m in range(mlo, 1):
        for a in range(5):
            if (m - (a - 1)) % 5: continue
            for n in range(NUT):
                if not complete(m, n, K): continue
                lhs = coef_or_none(W[m][a], n) if m in W else Fr(0)
                rhs = coef_or_none(YQ[m][a], n) if m in YQ else Fr(0)
                if lhs is None or rhs is None: continue
                tested += 1
                if lhs != rhs: bad.append((m, a, n, str(lhs), str(rhs)))
    check("(b) canonical P Y^C = Y^C_Q on all %d complete graded coefficients [z^m u^n], %d <= m <= 0, n < %d" % (tested, mlo, NUT), not bad, bad[:4])
    cnt_nz = sum(1 for m in YQ for a in range(5) if YQ[m][a].c)
    log("       Y^C_Q: %d non-zero (component, z-power) series; valuations %s" % (cnt_nz, {(m, a): YQ[m][a].v for m in sorted(YQ, reverse=True)[:8] for a in range(5) if YQ[m][a].c}))
    # direct check of the quintic system on W (the product, not the iterated integrals): d/du W_a[z^(m-1)] = b_a W_{a+1}[z^m]
    tested = 0; bad = []
    for m in range(mlo + 1, 1):
        for a in range(5):
            lhs = W[m - 1][a].deriv() if (m - 1) in W else Ser.zero(Q, NUT)
            rhs = bq[a].mul(W[m][(a + 1) % 5]) if m in W else Ser.zero(Q, NUT)
            for n in range(0, NUT - 2):
                if not (complete(m - 1, n + 1, K) and complete(m, n + 1, K)): continue
                l = coef_or_none(lhs, n); r = coef_or_none(rhs, n)
                if l is None or r is None: continue
                if l == 0 and r == 0: continue
                tested += 1
                if l != r: bad.append((m, a, n, str(l), str(r)))
    check("(b') the canonical product P Y^C satisfies the quintic system th Y_a = b_a Y_{a+1} (%d non-trivial coefficient identities)" % tested, not bad, bad[:4])
    # the order-1/z statement
    w1 = W[-1][0]
    nn = min(NUT, tc.prec) - 1
    okn = [n for n in range(nn) if complete(-1, n, K)]
    check("(b'') [z^-1] of (P Y^C)_0 = -t_c through u^%d" % (max(okn)), all(w1.coef(n) == -tc.coef(n) for n in okn))
    log("GI_CHECK %s: %s   (%.0f s)" % (tag, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
