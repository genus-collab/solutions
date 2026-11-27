# -*- coding: utf-8 -*-
"""ref_run4.py -- REFEREE item 4 (dictionary Fourier component <-> H^a (x) H^b) and two side checks.

 (A) Dictionary.  e~_j = Sum_a v_j^a Hh^a with v_j^a = Sum_b eta^{ab} Qh_{b,0}^{(j)} = (1/5) Qh_{3-a,0}^{(j)} (Hh^a the hatted flat classes).
     A colour matrix X is the bivector Sum_ij X[i][j] e~_i (x) e~_j = Sum_ab X^{ab} Hh^a (x) Hh^b,  X^{ab} = Sum_ij X[i][j] v_i^a v_j^b.
     Claim to check: X^{ab} = psi_{3-a} psi_{3-b} Xhat(-a-1, -b-1) (psi_b := Qh_{b,0}^{(0)}(u = 0), all nonzero), i.e. the Fourier
     component (m, m') is the component H^(-m-1) (x) H^(-m'-1), and Sum_j zeta^(m j) e~_j is proportional to H^(m-1).
     Then: the flat components of G_1 = [u^1 zhat^1] G computed DIRECTLY (no Fourier transform).
 (B) The LG hand check: formal stationary-phase expansion of Int exp((x^4 + u/x)/z) dx at x_c (x_c^5 = u/4), eps = z/x_c^4, against the
     identity row f(zhat) of the symbol with eps = -4 zhat, to order zhat^NLG (exact rationals; f reconstructed at 2^521-1).
 (C) Cross-check of the referee rows against the project row machinery (fq_rows.Rows -> fq_conifold -> fq_sser.SRows -> fq3_msum.HatRowsU),
     Nz = 6, prec_abs = 34, K = 6, at the prime given (only with --project; imports project modules read-only).
Run: py -I ref_run4.py --prime 89 [--project]
"""
import os, sys, time, argparse, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, inv, fmt, ratrec, mzero
from ref_hat import HatRows, tau_pred, M0_of, G_of, fourier

PR = {'89': P89, '61': P61, '521': P521}


def flat_components(X, V, p):
    return {(a, b): sum(X[i][j] * V[i][a] * V[j][b] for i in range(5) for j in range(5)) % p for a in range(5) for b in range(5)}


def stationary_phase(N):
    """coefficients a_m (m = 0..N) of the normalised expansion in delta^2 = -eps/20 of Int exp(phi(1+t)/eps) dt, phi(y) = y^4 + 4/y"""
    kmax = 2 * N + 2
    A = {}
    for k in range(4, kmax + 1):
        phik = Fr(5) if k == 4 else Fr(4 * (-1) ** k)
        A[(k - 2, k)] = -phik / 20                                   # delta^(k-2) xi^k
    def mul(X, Y):
        out = {}
        for (d1, x1), c1 in X.items():
            for (d2, x2), c2 in Y.items():
                d = d1 + d2
                if d > 2 * N: continue
                out[(d, x1 + x2)] = out.get((d, x1 + x2), 0) + c1 * c2
        return out
    tot = {(0, 0): Fr(1)}; power = {(0, 0): Fr(1)}; f = 1
    for j in range(1, N + 1):
        power = mul(power, A); f *= j
        for key, c in power.items(): tot[key] = tot.get(key, 0) + c / f
    def dfact(n):
        r = 1
        while n > 1: r *= n; n -= 2
        return r
    a = [Fr(0)] * (N + 1)
    for (d, x), c in tot.items():
        if x % 2 == 0 and d % 2 == 0: a[d // 2] += c * dfact(x - 1)
    return a


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--project", action="store_true")
    ap.add_argument("--NLG", type=int, default=10); a = ap.parse_args(); p = PR[a.prime]; NZ = 8; K = 6
    H = HatRows(p, NZ, K); z = H.zeta; F = lambda x: fmt(x, p, z)
    print("=" * 100); print("ref_run4: p = 2^%s-1; checks %s" % (a.prime, H.checks))
    # ---------------------------------------------------------------- (A)
    psi = [H.Q[b][0][0][0] for b in range(5)]
    print("(A) psi_b = Qh_{b,0}^(0)(u=0), b = 0..4: %s   (all nonzero: %s)" % ([F(x) for x in psi], all(psi)))
    i5 = inv(5, p); V = [[H.Q[(3 - aa) % 5][j][0][0] * i5 % p for aa in range(5)] for j in range(5)]      # V[j][a] = v_j^a
    ok_vec = True
    for m in range(5):
        comp = [sum(pow(z, (m * j) % 5, p) * V[j][aa] for j in range(5)) % p for aa in range(5)]
        ok_vec &= all((comp[aa] != 0) == (aa == (m - 1) % 5) for aa in range(5))
    print("    Sum_j zeta^(m j) e~_j has the single flat component H^(m-1) for every m: %s" % ok_vec)
    rng = random.Random(7); ok = True
    for _ in range(5):
        X = [[rng.randrange(p) for _ in range(5)] for _ in range(5)]; fc = flat_components(X, V, p); fo = fourier(X, z, p)
        for aa in range(5):
            for bb in range(5):
                pred = psi[(3 - aa) % 5] * psi[(3 - bb) % 5] * fo.get(((-aa - 1) % 5, (-bb - 1) % 5), 0) % p
                ok &= (pred == fc[(aa, bb)])
    print("    X^{ab} = psi_{3-a} psi_{3-b} Xhat(-a-1, -b-1) for 5 random colour matrices (125 components): %s" % ok)
    M0, M0i = M0_of(H.Rt, NZ, p); tau = [of(x, p) for x in tau_pred(K)]; G = G_of(H.Rt, M0i, tau, p, NZ, K)
    for n, k in ((1, 1), (2, 1), (3, 1)):
        fc = flat_components(G[(k, n)], V, p)
        print("    flat components X^{ab} of [u^%d zhat^%d] G (direct): %s" % (n, k, {key: F(v) for key, v in sorted(fc.items()) if v}))
        print("       H (x) H component X^{11}: %s" % F(fc[(1, 1)]))
    # ---------------------------------------------------------------- (B)
    N = a.NLG; aa_ = stationary_phase(N); pred = [aa_[m] / Fr(5) ** m for m in range(N + 1)]
    Hb = HatRows(P521, N, 1); f = [ratrec(Hb.Q[0][0][n][0], P521) for n in range(N + 1)]
    print("(B) stationary phase, in eps: %s" % ", ".join(str(aa_[m] * Fr(-1, 20) ** m) for m in range(4)))
    print("    predicted f (eps = -4 zhat): %s" % ", ".join(str(x) for x in pred[:6]))
    print("    identity row f of the symbol: %s" % ", ".join(str(x) for x in f[:6]))
    print("    agreement through zhat^%d: %s" % (N, [pred[m] == f[m] for m in range(N + 1)]))
    # ---------------------------------------------------------------- (C)
    if a.project:
        sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
        from fq_fp import FpRing
        from fq_rows import Rows
        from fq_conifold import Conifold
        from fq_sser import SRows
        from fq3_msum import HatRowsU
        t0 = time.time(); R = FpRing(p); assert R.zeta == z
        rows = Rows(R, Nz=6); C = Conifold(R, Nu=10, alpha=Fr(0), beta=Fr(0)); S = SRows(rows, C, prec_abs=34); HP = HatRowsU(S, 6)
        print("(C) project rows built in %.0f s (Nz = 6, prec_abs = 34, K = 6)" % (time.time() - t0))
        Hr = HatRows(p, 6, 6); ncmp = 0; bad = 0; minprec = 99
        for aa in range(5):
            for i in range(5):
                for n in range(7):
                    q = HP.Q[aa][i][n]; minprec = min(minprec, q.prec)
                    for k in range(min(q.prec, 6)):
                        ncmp += 1
                        if q.coef(k) != Hr.Q[aa][i][n][k]: bad += 1
        print("    coefficients compared: %d (minimal u-precision of the project rows: %d); mismatches: %d" % (ncmp, minprec, bad))
        print("    project conifold data: gamma = %s (referee %s); f0[3] = %s (referee %s)" % (C.gamma, F(Hr.g[2]), C.f[0].coef(3), F(Hr.f[0][3])))


if __name__ == "__main__":
    main()
