# -*- coding: utf-8 -*-
"""tr_run.py -- driver of the allgenus/tr task (TR_MODEL.md).  py -I -B tr_run.py <block> [args]   (writes tr_run_<block>.log)
Blocks:
  airy    Witten-Kontsevich regression of the TR engine (Q(zeta), 2^89-1, 2^61-1)
  rmat    exact (Q(zeta)) check R_TR(z) = M_0(z/kappa) and the unit/dilaton-leaf relation (2-32) through z^8, kappa = 1, 4, -1;
          M_0 from lg_symbol_rows_Q.json against the project's dumped symbol (both primes)
  exact   exact Q(zeta) TR on x = W through (3,1): symmetry, remainders, dilaton, string, flat tables, F_2, F_3
  frozen  TR mod p (both primes, zeta of fq_fp) against the project's MVertex6(M_0) tensors, colour basis, all d with D <= 6
  fg      F_g for g = 2..G (arg G, default 5) mod both primes on x = W, against B_2g/(2g(2g-2)) and c_g = (-5)^(1-g) F_g
  v0      u -> 0 at fixed p: expansion of omega_{g,n} at p = infinity (both primes; exact Q(zeta) for (g,n) <= (2,1))
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from tr_field import QzField, make_fp, P89, P61
from tr_core import ModelCurve
from tr_rec import TR
from tr_xi import XiBasis, laplace_Rinv, laplace_unit, pairing
from tr_compare import correlators, frozen_compare, string_failures, dilaton_failures, bernoulli, Fg_predicted

LIST_EXACT = ((0, 3), (1, 1), (0, 4), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1))


class Log(object):
    def __init__(self, name):
        self.f = open(os.path.join(HERE, "tr_run_%s.log" % name), "w"); self.t0 = time.time(); self.fails = 0
    def __call__(self, s=""):
        line = "%s" % s; print(line, flush=True); self.f.write(line + "\n"); self.f.flush()
    def check(self, ok, what):
        self("  [%s] %s" % ("PASS" if ok else "FAIL", what))
        if not ok: self.fails += 1
    def close(self):
        self("== %s (%.1f s, %d failures)" % ("ALL CHECKS PASSED" if not self.fails else "SOME CHECKS FAILED", time.time() - self.t0, self.fails))
        self.f.close()


def block_airy(log):
    from tr_test_airy import run
    for F in (QzField(), make_fp(P89), make_fp(P61)):
        log.check(run(F, log), "Airy curve over %s" % F.name)


def block_rmat(log):
    from tr_frozen import M0_coeffs, dump_check
    for tag in ("89", "61"):
        ok, cnt = dump_check(tag); log.check(ok, "M_0 from the model rows == project symbol dump lg_symbol_%s.json (%d entries, zhat^0..6)" % (tag, cnt))
    F = QzField(); NM = 8; C = M0_coeffs(F, NM)
    for kap in (1, 4, -1):
        cur = ModelCurve(F, 2 * NM + 8, kap); X = laplace_Rinv(cur, NM); U = laplace_unit(cur, NM)
        ci = [F.one / cur.c[k] for k in range(5)]
        A = []; Bt = []; Uok = []; Uf = []
        for n in range(NM + 1):
            f = Fr(-1, kap) ** n                                         # M_0(z/kappa)^-1 = M_0(-z/kappa)^T
            A.append(all((X[n][j][i] - C[n][i][j] * f).iszero() for i in range(5) for j in range(5)))
            Bt.append(all((X[n][j][i] - C[n][j][i] * f).iszero() for i in range(5) for j in range(5)))
            Uok.append(all((U[n][j] - sum((X[n][j][k] * ci[k] for k in range(5)), F.zero)).iszero() for j in range(5)))
            Uf.append(all((U[n][j] - sum((C[n][k][j] * f * ci[k] for k in range(5)), F.zero)).iszero() for j in range(5)))
        log("kappa = %s (x = %s W, Lambda = %s):" % (kap, kap, -5 * kap))
        log.check(all(A), "(2-31) Laplace transform of B = M_0(z/kappa)^-1 (standard reading, row = branch point), z^0..z^%d, 225 entries" % NM)
        log("  transposed reading X = (M_0^-1)^T holds at z^n for n = %s only (expected: n <= 1)" % [n for n in range(NM + 1) if Bt[n]])
        log.check(all(Uok), "(2-32) Laplace transform of dy (omega_{0,1}) = R^-1(z) 1^ with R^-1 from (2-31), z^0..z^%d" % NM)
        log.check(all(Uf), "(2-32) = M_0(z/kappa)^-1 1^ of the frozen theory (unit 1^ = Sum_k et_k/c_k), z^0..z^%d" % NM)
        if kap == 1:
            log("  dilaton leaf T(z) = z (1^ - M_0^-1(z) 1^) = -Sum_{n>=1} z^(n+1) (Laplace dy)_n, colour 0, z^2..z^5: %s" % [str(-U[n][0]) for n in range(1, 5)])


def block_exact(log):
    F = QzField(); cur = ModelCurve(F, 38, 1); tr = TR(cur); XB = XiBasis(cur, 7); flat = {}
    for g, n in LIST_EXACT:
        t = time.time(); W = tr.omega(g, n)
        log("omega_{%d,%d}: %d coefficients (%.1f s)" % (g, n, len(W), time.time() - t))
        log.check(tr.symmetric_failures(g, n) == 0, "omega_{%d,%d} symmetric in its arguments" % (g, n))
        cc, nb1 = correlators(tr, XB, g, n, 'c'); ff, nb2 = correlators(tr, XB, g, n, 'f'); flat[(g, n)] = ff
        log.check(nb1 == 0 and nb2 == 0, "omega_{%d,%d} lies in the span of the d xi (no remainder, colour and flat)" % (g, n))
        log.check(all(v.is_rational() for v in ff.values()), "flat correlators of (%d,%d) are rational (%d non-zero)" % (g, n, len(ff)))
        if n >= 2 and 2 * g - 3 + n > 0:
            D = tr.dilaton(W); low = tr.omega(g, n - 1)
            ok = set(D) == set(low) and all((D[k] - low[k] * (3 - 2 * g - n)).iszero() for k in low)
            log.check(ok, "TR dilaton: Sum Res Phi omega_{%d,%d}(., J) = (2-2g-|J|) omega_{%d,%d}(J)" % (g, n, g, n - 1))
        if (g, n - 1) in flat:
            c1, b1 = string_failures(ff, g, n, flat[(g, n - 1)]); c2, b2 = dilaton_failures(ff, g, n, flat[(g, n - 1)])
            log.check(b1 == 0 and b2 == 0, "flat string and dilaton equations (%d, %d tuples)" % (c1, c2))
    okp = True
    for aa in range(5):
        for dd in range(4):
            V = XB.vec[('f', aa, dd)]; vy = pairing(cur, V, 'y'); vp = pairing(cur, V, 'Phi')
            okp = okp and (vy - (1 if (aa, dd) == (0, 0) else 0)).iszero() and (vp - (-1 if (aa, dd) == (0, 1) else 0)).iszero()
    log.check(okp, "pairings: Sum Res y dxif_{a,d} = delta_(a,d),(0,0); Sum Res Phi dxif_{a,d} = -delta_(a,d),(0,1) (a <= 4, d <= 3)")
    log("\nFlat correlators of the frozen theory from TR (x = W): <prod tau_d(H^a)>_g, H^a = N_a Sum zeta^(i(a+1)) et_i")
    for (g, n), ff in flat.items():
        items = sorted((k, v.rational()) for k, v in ff.items() if tuple(sorted(k)) == k)
        log("  (g,n) = (%d,%d), %d sorted entries: " % (g, n, len(items)) + "; ".join("<%s> = %s" % (" ".join("t%d(H%d)" % (d, a) for a, d in k), v) for k, v in items[:40]))
    prim = {(1, 1): Fr(1, 12), (1, 2): Fr(1, 12), (1, 3): Fr(1, 6), (2, 1): Fr(1, 600), (2, 2): Fr(1, 200)}
    for (g, n), val in prim.items():
        key = tuple((1, 0) for _ in range(n)); mine = flat[(g, n)].get(key)
        mine = mine.rational() if mine is not None else Fr(0)
        log.check(mine * (-1) ** n == val, "primary sector <phi^%d>_%d = %s (phi = -H = Sum zeta^(2i) et_i, the convention of OU_EXTENSION_genus3 s.6d)" % (n, g, val))
    for g in (2, 3):
        Fg = tr.Fg(g).rational()
        log.check(Fg == Fg_predicted(g), "F_%d(x = W) = %s = B_%d/(%d.%d); frozen F_%d = (-5)^(1-g) F_g = %s" % (g, Fg, 2 * g, 2 * g, 2 * g - 2, g, Fg * Fr(-5) ** (1 - g)))


def block_frozen(log):
    from tr_frozen import Frozen
    for p in (P89, P61):
        Fz = Frozen(p); F = Fz.F; cur = ModelCurve(F, 38, 1); tr = TR(cur); XB = XiBasis(cur, 7); tot = 0; totnz = 0
        for g, n in LIST_EXACT + ((0, 5),):
            cc, nb = correlators(tr, XB, g, n, 'c'); ncmp, nnz, bad = frozen_compare(cc, Fz, g, n, F); tot += ncmp; totnz += nnz
            log.check(not bad and nb == 0, "p = 2^%d-1: (g,n) = (%d,%d): TR x (-5)^(1-g) == MVertex6(M_0) on %d colour-basis correlators (D <= 6), %d of them non-zero" % (p.bit_length(), g, n, ncmp, nnz))
            if bad: log("     first mismatches: %s" % bad[:3])
        for g in (2, 3):
            fz = Fz.mv.tensor(g, ())[0]; mine = (tr.Fg(g) * F.of(Fr(-5) ** (1 - g))).v
            log.check(fz == mine, "p = 2^%d-1: F_%d frozen (MVertex6) = %s == (-5)^(1-g) F_g(TR)" % (p.bit_length(), g, F.fmt(F.cls(fz))))
        log("  p = 2^%d-1: %d correlators compared, %d non-zero" % (p.bit_length(), tot, totnz))


def block_fg(log, G):
    for p in (P89, P61):
        F = make_fp(p); kmx = 6 * G - 3; cur = ModelCurve(F, 2 * kmx + 8, 1); tr = TR(cur)
        for g in range(2, G + 1):
            t = time.time(); Fg = tr.Fg(g); pred = Fg_predicted(g)
            log.check((Fg - F.of(pred)).iszero(), "p = 2^%d-1: F_%d(x = W) = %s (B_%d/(2g(2g-2)) = %s); c_%d = (-5)^(1-g) F_g = %s (%.1f s)" % (
                p.bit_length(), g, F.fmt(Fg), 2 * g, pred, g, pred * Fr(-5) ** (1 - g), time.time() - t))


def block_v0(log):
    from tr_v0 import v0_profile
    cases = ((0, 3), (0, 4), (0, 5), (1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (4, 1))
    for F in (make_fp(P89), make_fp(P61)):
        cur = ModelCurve(F, 50, 1); tr = TR(cur)
        for g, n in cases:
            byK, Kmin, bad = v0_profile(F, cur.c, tr.omega(g, n), g, n); reg = 2 * g - 2 + n
            lead = byK[Kmin]
            txt = ", ".join("%s: %s" % (m, F.fmt(x)) for m, x in sorted(lead.items())) if len(lead) <= 6 else "%d monomials" % len(lead)
            log.check(Kmin == reg and bad == 0, "%s: omega_{%d,%d}: c_m = 0 for |m| < %d, so omega^(u) = u^0 Omega_%d(p) + O(u) (regular at u = 0); limit: %s" % (
                F.name[:14], g, n, 5 * reg, reg, txt))


if __name__ == "__main__":
    blk = sys.argv[1]; log = Log(blk if blk != 'fg' else 'fg')
    log("tr_run.py %s  (allgenus/tr, %s)" % (" ".join(sys.argv[1:]), time.strftime("%Y-%m-%d %H:%M")))
    if blk == 'airy': block_airy(log)
    elif blk == 'rmat': block_rmat(log)
    elif blk == 'exact': block_exact(log)
    elif blk == 'frozen': block_frozen(log)
    elif blk == 'fg': block_fg(log, int(sys.argv[2]) if len(sys.argv) > 2 else 5)
    elif blk == 'v0': block_v0(log)
    log.close()
