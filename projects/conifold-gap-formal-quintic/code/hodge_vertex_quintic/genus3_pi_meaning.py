# -*- coding: utf-8 -*-
"""genus3_pi_meaning.py -- what the plane Pi of genus3_pi.py is: constant zhat-linear dressings on the LEFT of M_0 (genus 3, O(u^3)).

Claim.  Pi = { Delta = M_0^{-1} (zhat X) M_0 :  X a constant symmetric Z/5-covariant matrix whose (H, H) component vanishes }.
In the normalised canonical basis (eta = 1) put f_m = (zeta^{m j})_j; f_m is proportional to the flat class H^{m-1} (Delta_j^{-1/2} ~ zeta^j,
lambda_j = zeta^j; f_m . f_m' = 5 delta_{m+m' = 0 mod 5} reproduces eta(H^a, H^b) = 5 delta_{a+b = 3 mod 5}).  A symmetric order-1 circulant X has
Fourier components Xhat(m, m') = f_m^T X f_m' only for m + m' = 1 mod 5: (0, 1), (2, 4), (3, 3); as a bivector X = (1/25) Sum Xhat(m, m')
f_{-m} (x) f_{-m'}:
    Xhat(2, 4) <-> 1 (x) H^2 + H^2 (x) 1        (the unit on one side: string equation)
    Xhat(0, 1) <-> H^3 (x) H^4 + H^4 (x) H^3
    Xhat(3, 3) <-> H (x) H                       (excluded: this is the direction that the conditions see).
With M = M_0 exp(u Delta + ..) = exp(u zhat X) M_0 (..) the primed rows are R' = Psi (1 - u zhat X) + O(u^2): the edge factor at order u is the
CONSTANT bivector X (no psi classes).
Checks (per prime): (1) for the stored A, B of Pi: M_0 Delta M_0^{-1} has no component of zhat-orders 2..5 and its order-1 part has
Xhat(3, 3) = 0;  (2) conversely Delta_X = M_0^{-1} (zhat X) M_0 for the two generators and random combinations is circulant and satisfies EVERY
N_1-condition (4 linear rows, 5 quadrics, 9 cubics), while the generator Xhat(3, 3) violates them;  (3) X^2 = 0 for the string generator (so
exp(zhat X) = 1 + zhat X and the edge bivector of exp(-zhat X) is exactly X); X^3 = 0, X^2 != 0 for the other generator.
Run: py -I genus3_pi_meaning.py --prime 61|89
"""
import os, sys, json, argparse, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_mfit import fp_solve
from fq3_steps import State, frame_of
from fq3_rec import Rec
from fq3_m import sp_mats, sp_params, in_sp, NSP, KMAX
from fq_symbol import mat_series_mul
from fq3_par import setup
from fq3_circsys import CircSys


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61")
    a = ap.parse_args(); S = CircSys(a.prime); p = S.p; rec = Rec(S.R); z = S.R.zeta; cs = State(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_cubic_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("-" * 110); log("genus3_pi_meaning: Pi = Ad(M_0)^-1 (zhat X), p = 2^%s-1" % a.prime)
    st = setup(frame_of(p, (0, 0), 34, 4)); M0 = st['M0']; inv25 = pow(25, p - 2, p)
    Z = lambda: [[0] * 5 for _ in range(5)]
    def series(v15): return sp_mats([sum(S.cols[j][k] * v15[j] for j in range(15)) % p for k in range(NSP)], p)
    def fourier(X): return {(m, mp): sum(pow(z, m * i, p) * X[i][j] * pow(z, mp * j, p) for i in range(5) for j in range(5)) % p for m in range(5) for mp in range(5)}
    def from_fourier(Xh): return [[inv25 * sum(pow(z, (-m * i) % 5, p) * v * pow(z, (-mp * j) % 5, p) for (m, mp), v in Xh.items()) % p for j in range(5)] for i in range(5)]
    def circ_coords(N):
        if not in_sp(N, p): return None
        return fp_solve([[S.cols[j][v] for j in range(15)] for v in range(NSP)], sp_params(N, p), p)[0]
    # (1) the stored plane
    for nm in ("Pi_A", "Pi_B"):
        N = series(cs.d[nm]); D = mat_series_mul(mat_series_mul(M0.M, N, p, KMAX), M0.Minv, p, KMAX); Xh = fourier(D[1])
        sup = sorted(k for k, v in Xh.items() if v)
        check("(1) %s: M_0 Delta M_0^-1 vanishes in the zhat-orders 2, 3, 4, 5" % nm, all(not any(x % p for r in D[k] for x in r) for k in (2, 3, 4, 5)))
        check("(1) %s: the order-1 part X has Fourier support %s, i.e. Xhat(3, 3) = 0 (no H (x) H component)" % (nm, sup), set(sup) <= {(0, 1), (1, 0), (2, 4), (4, 2)})
        log("      Xhat(2, 4)/25 = %s  [1 (x) H^2 + H^2 (x) 1],  Xhat(0, 1)/25 = %s  [H^3 (x) H^4 + H^4 (x) H^3]" % (rec.fmt(Xh[(2, 4)] * inv25 % p), rec.fmt(Xh[(0, 1)] * inv25 % p)))
    # (2) the converse
    def delta_of(x24, x01, x33):
        Xh = {}
        if x24: Xh[(2, 4)] = x24; Xh[(4, 2)] = x24
        if x01: Xh[(0, 1)] = x01; Xh[(1, 0)] = x01
        if x33: Xh[(3, 3)] = x33
        X = from_fourier(Xh); zX = [Z() for _ in range(KMAX + 1)]; zX[1] = X
        return X, mat_series_mul(mat_series_mul(M0.Minv, zX, p, KMAX), M0.M, p, KMAX)
    rng = random.Random(8)
    for nm, co in (("string generator Xhat(2,4)", (1, 0, 0)), ("generator Xhat(0,1)", (0, 1, 0)), ("random combination", (rng.randrange(p), rng.randrange(p), 0)), ("random combination", (rng.randrange(p), rng.randrange(p), 0))):
        X, N = delta_of(*co); c = circ_coords(N)
        sym = all(X[i][j] == X[j][i] for i in range(5) for j in range(5))
        ok = c is not None and not any(S.residuals(c).values())
        check("(2) %s: X symmetric, Delta_X = M_0^-1 (zhat X) M_0 is in sp and circulant, and ALL N_1-conditions vanish (4 linear rows, 5 quadrics, 9 cubics)" % nm, sym and ok)
        if co[2] == 0 and (co[0] == 0 or co[1] == 0):
            X2 = [[sum(X[i][k] * X[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
            X3 = [[sum(X2[i][k] * X[k][j] for k in range(5)) % p for j in range(5)] for i in range(5)]
            if co[0]: check("(3) %s: X^2 = 0" % nm, not any(x for r in X2 for x in r))
            else: check("(3) %s: X^2 != 0, X^3 = 0" % nm, any(x for r in X2 for x in r) and not any(x for r in X3 for x in r))
            log("      circulant coordinates of Delta_X (orders 1, 2): %s" % ["%s: %s" % (S.names[j], rec.fmt(c[j])) for j in range(5)])
    X, N = delta_of(0, 0, 1); c = circ_coords(N); res = S.residuals(c) if c is not None else None
    bad = [k for k, v in res.items() if v] if res else None
    check("(2) the excluded generator Xhat(3, 3) (H (x) H): Delta_X is circulant but violates the conditions", c is not None and bool(bad), bad)
    # plane spanned by the two generators = stored Pi (compare the low parts)
    from fq3_quad import rank_fp
    g1 = circ_coords(delta_of(1, 0, 0)[1]); g2 = circ_coords(delta_of(0, 1, 0)[1])
    check("the low parts (orders 1..3) and the order-4 parts of the two generators span the same plane as the stored A, B", rank_fp([g1[:10], g2[:10], cs.d['Pi_A'][:10], cs.d['Pi_B'][:10]], p) == 2)
    cs.put('Pi_gen_string', g1); cs.put('Pi_gen_H3H4', g2); cs.put('pi_meaning_failed', FAILED)
    log("GENUS3_PI_MEANING ALL CHECKS PASSED" if not FAILED else "GENUS3_PI_MEANING FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
