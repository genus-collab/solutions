# run_tests.py -- referee tests (1)-(12) of LC_CY_PROOF.md (Theorem LC_CY, Q_CY, CF_CY, Cor. 12.1, Remark 0.1).
# usage: py -I run_tests.py N Z [pointset]
import sys, time, math
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else '.')
from reflib import *

N = int(sys.argv[1]) if len(sys.argv) > 1 else 11
Z = int(sys.argv[2]) if len(sys.argv) > 2 else 11
which = sys.argv[3] if len(sys.argv) > 3 else "all"
KMAX1 = int(sys.argv[5]) if len(sys.argv) > 5 else 8
KMAX23 = int(sys.argv[6]) if len(sys.argv) > 6 else 5

RESULTS = []
def rec(point, name, ok, detail=""):
    RESULTS.append((point.tag, name, ok, detail))
    print("  [%s] %-48s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)

def firstdiff(a, b):
    k, prec = s_cmp(a, b)
    return k, prec

def zdiff(ra, rb, zlist):
    """compare two z-dicts on the given z-powers; return (first (zpow, tpow, value) mismatch or None, min t-precision)"""
    minprec = None
    for m in zlist:
        pa, pb = zprec(ra), zprec(rb)
        x = ra.get(m, S(pa, []))
        y = rb.get(m, S(pb, []))
        k, prec = s_cmp(x, y)
        minprec = prec if minprec is None else min(minprec, prec)
        if k is not None:
            return (m, k, show((x - y).coef(k))), minprec
    return None, minprec

def run_point(pt, controls=False):
    print("\n=== %s  weights a0=%s a1=%s a2=%s e1=%s e2=%s  (E=%s, A=%s, B=%s, u=%s, CY=%s)  N=%d Z=%d  frame=%s" % (
        pt.tag, pt.fr['a0'], pt.fr['a1'], pt.fr['a2'], pt.fr['e1'], pt.fr['e2'], pt.fr['E'], pt.fr['A'], pt.fr['B'], pt.fr['u'], pt.isCY, pt.N, pt.Z, pt.frame), flush=True)
    E, a1, a2, A, B, a0, e1, e2, u = pt.E, pt.a1, pt.a2, pt.A, pt.B, pt.a0, pt.e1, pt.e2, pt.u
    prec = pt.rowprec
    # sanity: critical point, W_1 = 0; eigenvalue identities
    rec(pt, "sanity: W_1(U_i) = 0 (critical points)", all(pt.W[i][1].iszero() or pt.W[i][1].v >= pt.W[i][1].prec - 1 for i in (1, 2)),
        "prec %d" % pt.W[1][1].prec)
    s12 = pt.L[1] + pt.L[2] - (const(E, N) / pt.tau - a0)
    p12 = pt.L[1] * pt.L[2] - const(a1 * a2 % P, N) / pt.tau
    rec(pt, "sanity: L1+L2 = E/tau - a0, L1 L2 = a1 a2/tau", s12.iszero() and p12.iszero(), "prec %d,%d" % (s12.prec, p12.prec))
    rec(pt, "sanity: r^pt_{0,i}[0] = 1", all(s_cmp(pt.rpt[i][0], const(1, N))[0] is None for i in (1, 2)))
    # V1: Delta_i = a1 A e1 e2 W_2(U_i) and e1 e2 (L_i - L_j)(a1 a2 - E L_i)
    for i, j in ((1, 2), (2, 1)):
        d1 = pt.Delta[i] - pt.W[i][2].scale(a1 * A % P * e1 % P * e2 % P)
        d2 = pt.Delta[i] - ((pt.L[i] - pt.L[j]) * (a1 * a2 - E * pt.L[i])).scale(e1 * e2 % P)
        rec(pt, "sanity (V1): Delta_%d = a1 A e1 e2 W_2 = e1e2 (Li-Lj)(a1a2 - E Li)" % i, d1.iszero() and d2.iszero(), "prec %d" % min(d1.prec, d2.prec))

    # ---------------- (1) KEY FORMULA
    rhs = pt.key_rhs(Z)
    lhs = zscale(pt.r0[1], pt.rho)
    neg = [m for m in rhs if m < 0 and not (rhs[m].iszero() or rhs[m].v >= prec)]
    rec(pt, "(1) KEY RHS: negative z-powers cancel (mod t^%d)" % prec, not neg, "non-cancelling z-powers: %s (lowest z-power computed %d)" % (sorted(neg), min(rhs)))
    bad, mp = zdiff(lhs, rhs, range(0, Z + 1))
    rec(pt, "(1) KEY FORMULA rho r_{0,1} = (E/(E+z))(tau/t) F Ephi", bad is None, "z^0..z^%d, t^0..t^%d%s" % (Z, mp - 1, "" if bad is None else "; first diff (z^%d, t^%d) = %s" % bad))
    # t = 0 constant: r^pt_{0,1}(z;0) = h_{(E,-a1,-A)} E/(E+z)  (Prop 3.15_CY, Stirling)
    hc = h_series([E, (-a1) % P, (-A) % P], Z, N)
    rhs0 = zmul(hc, geom(E, 1, Z, N), Z)
    rhs0 = zscale(rhs0, const(E, N))
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(pt.rpt[1], m, prec).coef(0)
        y = zcoef(rhs0, m, N).coef(0)
        if x != y:
            ok = False; det = "first diff z^%d: %s vs %s" % (m, show(x), show(y)); break
    rec(pt, "(1) r^pt_{0,1}(z;0) = h_{(E,-a1,-A)} E/(E+z)", ok, det or "z^0..z^%d" % Z)

    # ---------------- (2) F1
    unit = pt.Delta[1] * pt.t * pt.t * inv((-e1 * e2 % P * pow(E, 3, P)) % P)
    ok = unit.v == 0 and unit.coef(0) == 1
    exp1 = (-u * inv(E * E % P)) % P
    rec(pt, "(2) F1: Delta_1 t^2/(-e1e2E^3) unit, = 1 - u t/E^2 + ...", ok and unit.coef(1) == exp1,
        "val %d, c0 = %s, c1 = %s (expected %s), prec %d" % (unit.v, show(unit.coef(0)), show(unit.coef(1)), show(exp1), unit.prec))
    rho1 = pt.rho.coef(1)
    rec(pt, "(2) rho_1 = u/(2E^2)", rho1 == u * inv(2 * E * E % P) % P, "rho_1 = %s" % show(rho1))
    sst = (1 - pt.U[1]) / pt.tau
    delta = const(E, N) + a0 * pt.tau - (2 * B) * pt.tau * sst
    alg2 = (delta * pt.U[1]).scale(a1 * A % P * pow(inv(E), 3, P) % P) / (sst * (1 - sst))
    alg = s_sqrt1(alg2)
    chk = pt.rho * pt.t * alg - pt.tau
    rec(pt, "(2) F1: rho = tau/(t alg), alg^2 = a1 A delta U_1/(E^3 s*(1-s*)), delta^2 = E^2-(E^2+u)tau+a0^2 tau^2",
        chk.iszero() and s_cmp(delta * delta, const(E * E % P, N) - ((E * E + u) % P) * pt.tau + (a0 * a0 % P) * pt.tau * pt.tau)[0] is None, "prec %d" % chk.prec)
    D2 = pt.Delta[2]
    d20 = (e1 * e2 % P * a1 % P * a2 % P * A % P * B % P * inv(E)) % P
    rec(pt, "(2) F1: Delta_2 regular, Delta_2(0) = e1e2 a1a2 A B/E, Delta_2'(0) = 0", D2.v >= 0 and D2.coef(0) == d20 and D2.coef(1) == 0,
        "val %d, Delta_2(0) = %s (exp %s), Delta_2'(0) = %s, prec %d" % (D2.v, show(D2.coef(0)), show(d20), show(D2.coef(1)), D2.prec))

    # ---------------- (3) F2
    g1 = geom(E, 1, Z, N)
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(pt.r0[1], m, prec).coef(0); y = zcoef(g1, m, N).coef(0) * E % P
        if x != y:
            ok = False; det = "first diff z^%d: %s vs %s" % (m, show(x), show(y)); break
    rec(pt, "(3) F2: r_{0,1}(z;0) = E/(E+z)", ok, det or "z^0..z^%d" % Z)
    rr1 = zscale(pt.r1[1], pt.rho)
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(rr1, m, prec)
        # expected: E/t (m = 0) + a1 a2 [z^m] 1/(E+z) + O(t)
        exp_m1 = E if m == 0 else 0
        exp_0 = a1 * a2 % P * zcoef(g1, m, N).coef(0) % P
        if x.coef(-1) != exp_m1 or x.coef(0) != exp_0 or (x.v < -1):
            ok = False; det = "z^%d: t^-1 %s (exp %s), t^0 %s (exp %s), val %d" % (m, show(x.coef(-1)), show(exp_m1), show(x.coef(0)), show(exp_0), x.v); break
    rec(pt, "(3) F2: rho r_{1,1} = E/t + a1a2/(E+z) + O(t)", ok, det or "z^0..z^%d, prec %d" % (Z, zprec(rr1)))

    # ---------------- (4) L1
    # rho q_t(psi) = -psi (rho r_{0,1})(-psi): [psi^a] = (-1)^a ... : q coefficient psi^(m+1) = (-1)^(m+1) r[m]
    rq = {m + 1: (lhs[m] if (m + 1) % 2 == 0 else -lhs[m]) for m in range(Z + 1)}  # psi^(m+1) coefficients, m+1 = 1..Z+1
    rq[0] = const(0, prec)
    rqprec = zprec(rq)
    kmaxL1 = min(KMAX1, rqprec - 1)
    dphi = pt.dphi()
    allok = True; details = []
    for k in range(kmaxL1 + 1):
        coeffs = [rq[a].coef(k) if a in rq else 0 for a in range(Z + 2)]
        # polynomial test: prod_{j=1}^{k+1} (E - j psi) * [t^k](rho q_t) has degree <= k, exactly k
        poly = [1]
        for j in range(1, k + 2):
            poly = polymul(poly, [E, (-j) % P])
        prod = polymul(poly, coeffs)[:Z + 2]  # exact up to psi^(Z+1)
        c, bad, M = leaf_fit(coeffs, k, E, 1)
        lam0 = c[0]
        if k == 0:
            # q^(0) = E - E^2 lam_1: constant E plus lam_1; (E - psi) q^(0) = -E psi
            hi_zero = all(x == 0 for x in prod[2:Z + 2]) and prod[0] == 0
            deg_k = prod[1] == (-E) % P
        else:
            hi_zero = all(x == 0 for x in prod[k + 1:Z + 2])
            Pk = 1
            for i in range(1, k + 1):
                Pk = Pk * (((k + 1) * a2 - i * E) % P) % P * (((k + 1) * A - i * E) % P) % P
            # sharpness (Cor. 11.3): top leaf lam_{k+1} present iff P_k != 0; a_{k+1}^(k) = (-1)^(k+1) C_{k+1}, a_k^(k) = (-1)^(k+1) k C_k u/(2E^2)
            Ck1 = Pk * inv(math.factorial(k + 1) * math.factorial(k) % P) % P * pow(inv(E), 2 * k - 2, P) % P
            top_exp = Ck1 if (k + 1) % 2 == 0 else (-Ck1) % P
            Pk_1 = 1
            for i in range(1, k):
                Pk_1 = Pk_1 * ((k * a2 - i * E) % P) % P * ((k * A - i * E) % P) % P
            Ck = Pk_1 * inv(math.factorial(k) * math.factorial(k - 1) % P) % P * (pow(inv(E), 2 * k - 4, P) if 2 * k >= 4 else pow(E, 4 - 2 * k, P)) % P
            sub_exp = k * Ck % P * u % P * inv(2 * E * E % P) % P
            if (k + 1) % 2: sub_exp = (-sub_exp) % P
            deg_k = ((c[k + 1] != 0) == (Pk != 0)) and c[k + 1] == top_exp and c[k] == sub_exp
            if Pk == 0:
                details.append("note k=%d: P_k = 0 at this point (E = (k+1)a2/i or (k+1)A/i), top leaf coefficient = %s (predicted 0)" % (k, show(c[k + 1])))
            if prod[k] == 0:
                details.append("note k=%d: polynomial degree < k (Sum_j a_j/j = 0)" % k)
        lam0_ok = (k == 0 and lam0 == E * E % P) or (k >= 1 and lam0 == 0)  # c_0 lam_0 = c_0/E ; constant E means c_0 = E^2
        ok = hi_zero and deg_k and bad is None and lam0_ok
        allok &= ok
        if k == 0:
            ok0 = c[1] == (-E * E) % P
            allok &= ok0
            details.append("k=0: const %s (exp E), lam_1 coeff %s (exp -E^2)" % (show(lam0 * inv(E) % P), show(c[1])))
        if not ok:
            details.append("k=%d: hi_zero=%s deg_k=%s fit_bad=%s lam0=%s" % (k, hi_zero, deg_k, bad, show(lam0)))
        # closed form A_j
        cf_ok = True
        for j in range(1, k + 2):
            Aj = pt.A_j_closed(j, dphi)
            if Aj.coef(k) != c[j]:
                cf_ok = False
                details.append("k=%d j=%d: a_j^(k) = %s, closed form %s" % (k, j, show(c[j]), show(Aj.coef(k))))
        allok &= cf_ok
    rec(pt, "(4) L1: [t^k](rho q_t) in Lam_{k+1}, no lam_0, top leaf iff P_k!=0, Cor11.3 (k<=%d, psi^<=%d)" % (kmaxL1, Z + 1), allok, "; ".join(details[:4]))
    # A_j also against the full series rho q_t = E + Sum_j A_j lam_j (t^0..t^(rqprec-1)) for psi^a, a <= Z+1
    ok = True; det = ""
    for a in range(Z + 2):
        tot = const(E if a == 0 else 0, rqprec)
        for j in range(1, rqprec + 1):
            Aj = pt.A_j_closed(j, dphi)
            tot = tot + Aj.scale(pow(j, a, P) * pow(inv(E), a + 1, P) % P)
        k, pr = s_cmp(tot, rq[a])
        if k is not None:
            ok = False; det = "first diff psi^%d t^%d: %s" % (a, k, show((tot - rq[a]).coef(k))); break
    rec(pt, "(4) CF: rho q_t = E + Sum_j A_j lam_j (all t^<%d, psi^<=%d)" % (rqprec, Z + 1), ok, det)
    # control: z -> -z in the leaves (lam_j = 1/(E + j psi))
    ctrl_fail = False
    for k in range(1, min(3, kmaxL1) + 1):
        coeffs = [rq[a].coef(k) if a in rq else 0 for a in range(Z + 2)]
        c, bad, M = leaf_fit(coeffs, k, E, -1)
        if bad is not None:
            ctrl_fail = True
    rec(pt, "(12d) control: leaves with z -> -z must FAIL L1", ctrl_fail, "fails as required" if ctrl_fail else "UNEXPECTEDLY PASSES")

    # ---------------- (7) F3 regularity of Dil, E_ij
    dmax = min(6, Z)
    ok = True; det = []
    for i in (1, 2):
        for d in range(1, dmax + 1):
            x = pt.dil(i, d)
            if x.v < 0:
                ok = False; det.append("Dil_%d(%d) val %d" % (i, d, x.v))
    # Dil_1(d;0) = -E^(1-d), Dil_1(1) = -1 at every t
    x = pt.dil(1, 1)
    if not (x.v == 0 and x.coef(0) == P - 1 and all(x.coef(k) == 0 for k in range(1, x.prec))):
        ok = False; det.append("Dil_1(1) != -1")
    for d in range(1, dmax + 1):
        if pt.dil(1, d).coef(0) != (-pow(inv(E), d - 1, P)) % P:
            ok = False; det.append("Dil_1(%d;0) wrong" % d)
    rec(pt, "(7) F3: Dil_i(d) regular, d<=%d; Dil_1(d;0) = -E^(1-d)" % dmax, ok, "; ".join(det[:3]))
    EE = {}
    ok = True; det = []
    mmax = min(6, Z - 1)
    for i in (1, 2):
        for j in (1, 2):
            for m in range(0, Z):
                for n in range(0, Z - m):
                    x = pt.edge(i, j, m, n)
                    EE[(i, j, m, n)] = x
                    if m <= 6 and n <= 6 and x.v < 0:
                        ok = False; det.append("E_%d%d(%d,%d) val %d" % (i, j, m, n, x.v))
    eprec = min(x.prec for x in EE.values())
    rec(pt, "(7) F3: E_ij(m,n) regular at t=0, m,n<=%d (m+n<=%d)" % (mmax, Z - 1), ok, "; ".join(det[:3]) or "t-prec %d" % eprec)
    # symmetry and transposition
    ok = all(s_cmp(EE[(1, 1, m, n)], EE[(1, 1, n, m)])[0] is None for m in range(Z) for n in range(Z - m))
    ok2 = all(s_cmp(EE[(2, 2, m, n)], EE[(2, 2, n, m)])[0] is None for m in range(Z) for n in range(Z - m))
    ok3 = all(s_cmp(EE[(1, 2, m, n)], EE[(2, 1, n, m)])[0] is None for m in range(Z) for n in range(Z - m))
    rec(pt, "(5) E_11, E_22 symmetric; E_21(b,a) = E_12(a,b)", ok and ok2 and ok3, "E11 %s E22 %s E12/E21 %s" % (ok, ok2, ok3))
    # t = 0 values: E_11 = -(kappa2/2) lam_1 lam_1 ; E_22 and E_12 (uses hreg)
    kap2 = e1 * e2 % P * a1 % P * a2 % P * A % P * B % P
    ok = True; det = ""
    for m in range(Z):
        for n in range(Z - m):
            exp = (-kap2 * inv(2) % P) * pow(inv(E), m + n + 2, P) % P
            if EE[(1, 1, m, n)].coef(0) != exp:
                ok = False; det = "first diff (%d,%d): %s vs %s" % (m, n, show(EE[(1, 1, m, n)].coef(0)), show(exp)); break
        if not ok: break
    rec(pt, "(7) E_11(z,w)|_{t=0} = -(kappa2/2) lam_1(z) lam_1(w)", ok, det)
    hreg = h_series([a1, a2, A, B, (-E) % P, (-E) % P], Z, N)
    hreg_inv = zinv_scalar(hreg, Z, N)
    # E_12(psi,w)|_{t=0} = lam_1(psi) chi12(w), chi12 = -kappa2 hreg(w)^-1/(E+w)
    chi = zscale(zmul(hreg_inv, geom(E, 1, Z, N), Z), const((-kap2) % P, N))
    ok = True; det = ""
    for m in range(Z):
        for n in range(Z - m):
            exp = pow(inv(E), m + 1, P) * zcoef(chi, n, N).coef(0) % P
            if EE[(1, 2, m, n)].coef(0) != exp:
                ok = False; det = "first diff (%d,%d): %s vs %s" % (m, n, show(EE[(1, 2, m, n)].coef(0)), show(exp)); break
        if not ok: break
    rec(pt, "(9) E_12|_{t=0} = lam_1(psi) (-kappa2 hreg(w)^-1/(E+w))", ok, det)
    # E_22|_{t=0} = (kappa2/E)[1 - hreg(z)^-1 hreg(w)^-1]/(z+w): coefficient [z^m w^n] via EdN with rows hreg(1, a1a2/E) (pairing)
    # direct: [1 - hi(z) hi(w)]/(z+w) where hi = hreg^-1: use the Lemma 1.1 extraction with P(z,w) = -hi(-z) hi(-w)*(-Delta2(0))?  simpler: expand
    hi = [zcoef(hreg_inv, k, N).coef(0) for k in range(Z + 1)]
    ok = True; det = ""
    for m in range(Z):
        for n in range(Z - m):
            # [z^m w^n] (1 - hi(z) hi(w))/(z + w): write Num = 1 - sum hi_p hi_q z^p w^q; division by (z+w): coefficient c_{m,n} satisfies
            # Num_{m+1, n} ... use recursion: Num(z,w) = (z + w) C(z,w)  ->  Num_{p,q} = C_{p-1,q} + C_{p,q-1}; solve C_{m,n} = Num_{m+1,n} - C_{m+1,n-1}
            def Numc(p, q):
                v = (1 if (p == 0 and q == 0) else 0) - hi[p] * hi[q]
                return v % P
            val = 0
            sgn = 1
            for k in range(n + 1):
                val = (val + sgn * Numc(m + 1 + k, n - k)) % P
                sgn = -sgn
            exp = kap2 * inv(E) % P * val % P
            if EE[(2, 2, m, n)].coef(0) != exp:
                ok = False; det = "first diff (%d,%d): %s vs %s" % (m, n, show(EE[(2, 2, m, n)].coef(0)), show(exp)); break
        if not ok: break
    rec(pt, "(9) E_22|_{t=0} = (kappa2/E)[1 - hreg(z)^-1 hreg(w)^-1]/(z+w)", ok, det)

    # ---------------- (5) L2: legs f_{j,b}(psi) = Sum_a E_1j(a,b) psi^a, k <= 5, b <= 3
    rhoinv = pt.rho.inverse()
    allok = True; det = []
    for j in (1, 2):
        for b in range(0, 4):
            amax = Z - 1 - b
            legs = [EE[(1, j, a, b)] for a in range(amax + 1)]
            rlegs = [x * pt.rho for x in legs]
            lp = min(x.prec for x in rlegs)
            for k in range(0, min(KMAX23, lp - 1) + 1):
                for nm, L in (("rho f", rlegs), ("f", legs)):
                    coeffs = [x.coef(k) for x in L]
                    if len(coeffs) < k + 3:
                        det.append("j=%d b=%d k=%d: too few psi-orders (%d) to test" % (j, b, k, len(coeffs)))
                        continue
                    c, bad, M = leaf_fit(coeffs, k, E, 1)
                    if bad is not None or c[0] != 0:
                        allok = False
                        det.append("FAIL j=%d b=%d k=%d %s: fit mismatch at psi^%s, lam0 = %s" % (j, b, k, nm, bad, show(c[0])))
            # t = 0 value of rho f_{1,b}: -(kappa2/2) E^(-b-1) lam_1
            if j == 1:
                for a in range(amax + 1):
                    exp = (-kap2 * inv(2) % P) * pow(inv(E), b + 1, P) % P * pow(inv(E), a + 1, P) % P
                    if rlegs[a].coef(0) != exp:
                        allok = False; det.append("rho f_{1,%d} t^0 psi^%d wrong" % (b, a)); break
    rec(pt, "(5) L2: [t^k](rho f_{j,b}), [t^k] f_{j,b} in Lam_{k+1}, j=1,2, b<=3, k<=%d" % KMAX23, allok, "; ".join(det[:4]) or "psi^<=%d (b=0) .. %d (b=3)" % (Z - 1, Z - 4))

    # ---------------- (6) L3 and the flow
    rho2 = pt.rho * pt.rho
    C = {}
    M2 = Z - 1
    allok = True; det = []
    for k in range(0, min(KMAX23, eprec - 1) + 1):
        # use the square block a, b <= M2 with a + b <= Z-1: fill a full square of side Mside where all entries exist: a, b <= (Z-1)//2
        Ms = (Z - 1) // 2
        Cm = [[(EE[(1, 1, a, b)] * rho2).coef(k) for b in range(Ms + 1)] for a in range(Ms + 1)]
        if Ms < k + 2:
            det.append("k=%d: square block too small (Ms=%d)" % (k, Ms)); continue
        c, bad = leaf_fit2(Cm, k, E, 1)
        # also verify predictions on the off-square entries a + b <= Z-1 with a or b > Ms
        iE = inv(E)
        def lamc(j, a):
            if j == 0: return iE if a == 0 else 0
            return pow(j, a, P) * pow(iE, a + 1, P) % P
        bad2 = None
        for a in range(Z):
            for b in range(Z - a):
                if a <= Ms and b <= Ms: continue
                pred = sum(c[i][l] * lamc(i, a) % P * lamc(l, b) for i in range(k + 2) for l in range(k + 2)) % P
                if pred != (EE[(1, 1, a, b)] * rho2).coef(k):
                    bad2 = (a, b); break
            if bad2: break
        supp_ok = all(c[i][l] == 0 for i in range(k + 2) for l in range(k + 2) if i == 0 or l == 0 or (i - 1) + (l - 1) > k)
        sym_ok = all(c[i][l] == c[l][i] for i in range(k + 2) for l in range(k + 2))
        # anti-diagonal (7.6): e^(k)_{k1+1,k2+1} = (-1)^(k+1) e1e2 a1a2 A B (k1+1)(k2+1) P_k1 P_k2 /((k+2) (k1!)^2 (k2!)^2 E^(2k))
        def Pk(kk):
            v = 1
            for i in range(1, kk + 1):
                v = v * (((kk + 1) * a2 - i * E) % P) % P * (((kk + 1) * A - i * E) % P) % P
            return v
        anti_ok = True
        for k1 in range(k + 1):
            k2 = k - k1
            exp = kap2 * (k1 + 1) * (k2 + 1) % P * Pk(k1) % P * Pk(k2) % P * inv((k + 2) * math.factorial(k1) ** 2 * math.factorial(k2) ** 2 % P) % P * pow(iE, 2 * k, P) % P
            if (k + 1) % 2: exp = (-exp) % P
            if c[k1 + 1][k2 + 1] != exp:
                anti_ok = False; det.append("k=%d anti-diag (%d,%d): %s vs (7.6) %s" % (k, k1 + 1, k2 + 1, show(c[k1 + 1][k2 + 1]), show(exp)))
            if exp == 0:
                det.append("note k=%d: P_%d P_%d = 0 at this point" % (k, k1, k2))
        ok = bad is None and bad2 is None and supp_ok and sym_ok and anti_ok
        allok &= ok
        if not ok:
            det.append("k=%d: fit=%s offsq=%s supp=%s sym=%s anti=%s" % (k, bad, bad2, supp_ok, sym_ok, anti_ok))
    rec(pt, "(6) L3: [t^k](rho^2 E_11) support (i-1)+(l-1)<=k, symmetric, anti-diagonal = (7.6), k<=%d" % min(KMAX23, eprec - 1), allok, "; ".join(det[:4]))
    # flow: d/dt [t^2 rho^2 E_11(psi,psi')] = Theta Dg(psi) Dg(psi'), Theta = -e1e2 B tau s* L_2, g(psi;t) = (rho r_{0,1})(-psi)
    sstar = (1 - pt.U[1]) / pt.tau
    Theta = (pt.tau * sstar * pt.L[2]).scale((-e1 * e2 % P * B) % P)
    Theta_b = pt.Delta[2] / (pt.L[2] - pt.L[1])
    Theta_c = (a1 * a2 - E * pt.L[2]).scale(e1 * e2 % P)
    rec(pt, "(6) Theta: three forms agree, Theta_0 = 0, Theta_1 = -kappa2/E^2", s_cmp(Theta, Theta_b)[0] is None and s_cmp(Theta, Theta_c)[0] is None and Theta.v >= 1 and Theta.coef(1) == (-kap2 * inv(E * E % P)) % P,
        "Theta_1 = %s, prec %d" % (show(Theta.coef(1)), Theta.prec))
    Dg = {}
    for m in range(Z + 1):
        gm = lhs[m] if m % 2 == 0 else -lhs[m]   # [psi^m] g = (-1)^m [z^m] G
        Dg[m] = (pt.t * gm).der()
    ok = True; det = ""; fprec = None
    for m in range(Z):
        for n in range(Z - m):
            lhsf = (pt.t * pt.t * rho2 * EE[(1, 1, m, n)]).der()
            rhsf = Theta * Dg[m] * Dg[n]
            k, pr = s_cmp(lhsf, rhsf)
            fprec = pr if fprec is None else min(fprec, pr)
            if k is not None:
                ok = False; det = "first diff (psi^%d psi'^%d, t^%d): %s" % (m, n, k, show((lhsf - rhsf).coef(k))); break
        if not ok: break
    rec(pt, "(6) flow d/dt[t^2 rho^2 E_11] = Theta Dg Dg (m+n<=%d)" % (Z - 1), ok, det or "t-prec %d" % fprec)

    # ---------------- (8) Theorem Q_CY: second QDE row at both colours; symplecticity
    for i in (1, 2):
        Li = pt.L[i]
        dd = pt.Delta[i].der() / pt.Delta[i] * inv(2)
        ok = True; det = ""; qp = None
        for m in range(Z + 1):
            r1m = zcoef(pt.r1[i], m, prec)
            x = Li * r1m
            if m >= 1:
                f = zcoef(pt.r1[i], m - 1, prec)
                x = x + f.der() - dd * f
            res = x + (const(a1 * a2 % P, N) / pt.tau) * zcoef(pt.r0[i], m, prec) - (const(E, N) / pt.tau - a0) * r1m
            qp = res.prec if qp is None else min(qp, res.prec)
            if not res.iszero():
                ok = False; det = "z^%d: residual val %d, coef %s" % (m, res.v, show(res.c[0])); break
        rec(pt, "(8) QDE row 2 at colour %d, z^0..z^%d" % (i, Z), ok, det or "t-prec %d" % qp)
    ok = True; det = ""
    for i in (1, 2):
        for j in (1, 2):
            for m in range(Z + 1):
                x = pt.symp(i, j, m)
                exp = pt.Delta[i] if (i == j and m == 0) else const(0, x.prec)
                k, pr = s_cmp(x, exp)
                if k is not None:
                    ok = False; det = "(%d,%d) w^%d t^%d: %s" % (i, j, m, k, show((x - exp).coef(k))); break
            if not ok: break
        if not ok: break
    rec(pt, "(8) <v_i(w), v_j(-w)> = delta_ij Delta_i, w^0..w^%d" % Z, ok, det)

    # ---------------- (9) regular colour at t = 0
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(pt.r1[2], m, prec).coef(0); y = a1 * a2 % P * inv(E) % P * zcoef(pt.r0[2], m, prec).coef(0) % P
        if x != y:
            ok = False; det = "z^%d: %s vs %s" % (m, show(x), show(y)); break
    rec(pt, "(9) v_2(z;0) proportional to (1, a1a2/E)", ok, det)
    ha = h_series([a2, (-E) % P, B], Z, N)
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(pt.rpt[2], m, prec).coef(0); y = zcoef(ha, m, N).coef(0)
        if x != y:
            ok = False; det = "z^%d: %s vs %s" % (m, show(x), show(y)); break
    rec(pt, "(9) r^pt_{0,2}(z;0) = h_{(a2,-E,B)}", ok, det)
    ok = True; det = ""
    for m in range(Z + 1):
        x = zcoef(pt.r0[2], m, prec).coef(0); y = zcoef(hreg, m, N).coef(0)
        if x != y:
            ok = False; det = "z^%d: %s vs %s" % (m, show(x), show(y)); break
    rec(pt, "(9) r_{0,2}(z;0) = hreg = h_{(a1,a2,A,B,-E,-E)}", ok, det)
    rec(pt, "(9) [t^1] r_{0,2} = 0", all(zcoef(pt.r0[2], m, prec).coef(1) == 0 for m in range(Z + 1)))

    # ---------------- (12c) control: A -> a2 in the KEY FORMULA (only meaningful when A != a2)
    if pt.fr['A'] != pt.fr['a2']:
        rhsc = pt.key_rhs(Z, Apar=a2)
        bad, mp = zdiff(lhs, rhsc, range(0, Z + 1))
        rec(pt, "(12c) control: KEY FORMULA with A -> a2 must FAIL", bad is not None, "first diff (z^%d, t^%d) = %s" % bad if bad else "UNEXPECTEDLY PASSES")
    return EE

def control_frames(a0, a1, a2, e1, e2, N, Z):
    """(12a) wrong class (E+1, -e1, -e2); (12b) omit h_nw.  Report KEY FORMULA, F2, L1 at the first few orders."""
    base = dict(a0=a0, a1=a1, a2=a2, e1=e1, e2=e2)
    E = Fr(a0) + Fr(a1) + Fr(a2)
    for tag, kw in (("ctrl-wrongclass", dict(frame=[E + 1, -Fr(e1), -Fr(e2)])), ("ctrl-no-hnw", dict(use_hnw=False))):
        pt = Point(N=N, Z=Z, tag=tag, **base, **kw)
        print("\n=== %s  (a0,a1,a2,e1,e2) = %s" % (tag, (a0, a1, a2, e1, e2)), flush=True)
        prec = pt.rowprec
        rhs = pt.key_rhs(Z)
        lhs = zscale(pt.r0[1], pt.rho)
        bad, mp = zdiff(lhs, rhs, range(0, Z + 1))
        rec(pt, "(12) control %s: KEY FORMULA must FAIL" % tag, bad is not None, "first diff (z^%d, t^%d) = %s" % bad if bad else "UNEXPECTEDLY PASSES")
        g1 = geom(pt.E, 1, Z, N)
        f2bad = None
        for m in range(Z + 1):
            if zcoef(pt.r0[1], m, prec).coef(0) != zcoef(g1, m, N).coef(0) * pt.E % P:
                f2bad = m; break
        rec(pt, "(12) control %s: F2 r_{0,1}(z;0) = E/(E+z) must FAIL" % tag, f2bad is not None, "first diff z^%d" % f2bad if f2bad is not None else "UNEXPECTEDLY PASSES")
        rq = {m + 1: (lhs[m] if (m + 1) % 2 == 0 else -lhs[m]) for m in range(Z + 1)}
        rq[0] = const(0, prec)
        l1bad = []
        for k in range(0, 4):
            coeffs = [rq[a].coef(k) if a in rq else 0 for a in range(Z + 2)]
            c, bad, M = leaf_fit(coeffs, k, pt.E, 1)
            if bad is not None or (k >= 1 and c[0] != 0):
                l1bad.append(k)
        rec(pt, "(12) control %s: L1 (k<=3) must FAIL" % tag, len(l1bad) > 0, "fails at k = %s" % l1bad if l1bad else "UNEXPECTEDLY PASSES")

def cor12(a1, a2, e1, e2, N, Z):
    """Corollary 12.1: a0 = 0 frame-B objects vs the fibrewise point (0, x, 1, x, 1), x = a1/a2, z -> z/a2."""
    a1, a2, e1, e2 = Fr(a1), Fr(a2), Fr(e1), Fr(e2)
    x = a1 / a2
    pt = Point(0, a1, a2, e1, e2, N, Z, tag="cor12-gen")
    tf = Point(0, x, 1, x, 1, N, Z, tag="cor12-TF")
    print("\n=== Corollary 12.1: (0,%s,%s,%s,%s) vs T_F x = %s   CY=%s" % (a1, a2, e1, e2, x, pt.isCY), flush=True)
    A2 = F(a2)
    kappa = F(e1 * e2 * a2 ** 4 / a1)
    prec = min(pt.rowprec, tf.rowprec)
    ok = True; det = ""
    for i in (1, 2):
        for m in range(Z + 1):
            # r_{0,i}(z) = r^TF_{0,i}(z/a2): [z^m] r = a2^-m [z^m] r^TF
            x0 = zcoef(pt.r0[i], m, prec); y0 = zcoef(tf.r0[i], m, prec).scale(pow(inv(A2), m, P))
            x1 = zcoef(pt.r1[i], m, prec); y1 = zcoef(tf.r1[i], m, prec).scale(pow(inv(A2), m, P) * A2 % P)
            k0, _ = s_cmp(x0, y0); k1, _ = s_cmp(x1, y1)
            if k0 is not None or k1 is not None:
                ok = False; det = "colour %d z^%d: r0 diff t^%s, r1 diff t^%s" % (i, m, k0, k1); break
        if not ok: break
    rec(pt, "(10) Cor 12.1: r_{0,i}(z) = r^TF_{0,i}(z/a2), r_{1,i}(z) = a2 r^TF_{1,i}(z/a2), z^<=%d" % Z, ok, det or "t-prec %d" % prec)
    ok = True; det = []
    for i in (1, 2):
        if s_cmp(pt.L[i], tf.L[i].scale(A2))[0] is not None: ok = False; det.append("L_%d" % i)
        if s_cmp(pt.Delta[i], tf.Delta[i].scale(kappa))[0] is not None: ok = False; det.append("Delta_%d" % i)
    if s_cmp(pt.rho, tf.rho)[0] is not None: ok = False; det.append("rho")
    for i in (1, 2):
        for d in range(1, min(6, Z) + 1):
            if s_cmp(pt.dil(i, d), tf.dil(i, d).scale(pow(inv(A2), d - 1, P)))[0] is not None: ok = False; det.append("Dil_%d(%d)" % (i, d))
    for i in (1, 2):
        for j in (1, 2):
            for m in range(min(6, Z - 1) + 1):
                for n in range(min(6, Z - 1 - m) + 1):
                    if s_cmp(pt.edge(i, j, m, n), tf.edge(i, j, m, n).scale(kappa * pow(inv(A2), 1 + m + n, P) % P))[0] is not None:
                        ok = False; det.append("E_%d%d(%d,%d)" % (i, j, m, n))
    rec(pt, "(10) Cor 12.1: L_i = a2 L^TF, Delta_i = kappa Delta^TF, rho = rho^TF, Dil_i(d) = a2^(1-d) Dil^TF, E_ij(m,n) = kappa a2^(-1-m-n) E^TF", ok, "; ".join(det[:5]))
    # h_{(a1,a2,-E)}(z) = h_{(c,-x,-1)}(z/a2)^-1
    hA = h_series([F(a1), F(a2), F(-(a1 + a2))], Z, N)
    hT = h_series([F(1 + x), F(-x), F(-1)], Z, N)
    hTi = zinv_scalar(hT, Z, N)
    ok = all(zcoef(hA, m, N).coef(0) == zcoef(hTi, m, N).coef(0) * pow(inv(A2), m, P) % P for m in range(Z + 1))
    rec(pt, "(10) h_{(a1,a2,-E)}(z) = h_{(c,-x,-1)}(z/a2)^-1", ok)

def remark01(a0, a1, a2, e1, e2, e1b, e2b, N, Z):
    """Remark 0.1: the frame-B rows at (a0,a1,a2,e1,e2) and (a0,a1,a2,e1b,e2b) coincide; Delta_i, E_ij, pairing scale by e1 e2."""
    pa = Point(a0, a1, a2, e1, e2, N, Z, tag="rem01-a")
    pb = Point(a0, a1, a2, e1b, e2b, N, Z, tag="rem01-b")
    print("\n=== Remark 0.1: (%s,%s,%s | e=%s,%s  CY=%s) vs (e=%s,%s  CY=%s)" % (a0, a1, a2, e1, e2, pa.isCY, e1b, e2b, pb.isCY), flush=True)
    prec = min(pa.rowprec, pb.rowprec)
    ok = True; det = ""
    for i in (1, 2):
        for m in range(Z + 1):
            k0, _ = s_cmp(zcoef(pa.r0[i], m, prec), zcoef(pb.r0[i], m, prec))
            k1, _ = s_cmp(zcoef(pa.r1[i], m, prec), zcoef(pb.r1[i], m, prec))
            if k0 is not None or k1 is not None:
                ok = False; det = "colour %d z^%d differs (t^%s / t^%s)" % (i, m, k0, k1); break
        if not ok: break
    rec(pa, "(11) Remark 0.1: frame-B rows r_{0,i}, r_{1,i} independent of (e1, e2)", ok, det or "z^<=%d, t-prec %d" % (Z, prec))
    ratio = F(Fr(e1) * Fr(e2) / (Fr(e1b) * Fr(e2b)))
    ok = True; det = []
    for i in (1, 2):
        if s_cmp(pa.Delta[i], pb.Delta[i].scale(ratio))[0] is not None: ok = False; det.append("Delta_%d" % i)
    if s_cmp(pa.rho, pb.rho)[0] is not None: ok = False; det.append("rho")
    for i in (1, 2):
        for j in (1, 2):
            for m in range(min(6, Z - 1) + 1):
                for n in range(min(6, Z - 1 - m) + 1):
                    if s_cmp(pa.edge(i, j, m, n), pb.edge(i, j, m, n).scale(ratio))[0] is not None:
                        ok = False; det.append("E_%d%d(%d,%d)" % (i, j, m, n))
    rec(pa, "(11) Remark 0.1: Delta_i, E_ij scale by e1 e2; rho, Dil unchanged", ok, "; ".join(det[:5]))

def wick_partition_check(pt, mmax=4):
    """(C7) partition formula wickZ_m = Sum_{n=m}^{3m} ((2n)!(-1)^n/(n! 2^n W_2^n)) Sum_p prod_k W_k^(m_k)/((k!)^(m_k) m_k!),
    p over partitions of 2n into exactly n - m parts >= 3; compared with the exp-route rows rpt[i][m]."""
    def parts(total, nparts, minp):
        if nparts == 0:
            if total == 0:
                yield []
            return
        for first in range(minp, total // nparts + 1):
            for rest in parts(total - first, nparts - 1, first):
                yield [first] + rest
    ok = True; det = ""
    for i in (1, 2):
        Wk = pt.W[i]
        W2inv = Wk[2].inverse()
        for m in range(0, mmax + 1):
            acc = const(0, pt.N) if m > 0 else const(1, pt.N)
            if m > 0:
                acc = None
                for n in range(m, 3 * m + 1):
                    pref = F(Fr(math.factorial(2 * n) * (-1) ** n, math.factorial(n) * 2 ** n))
                    W2p = const(1, pt.N)
                    for _ in range(n):
                        W2p = W2p * W2inv
                    inner = None
                    for p in parts(2 * n, n - m, 3):
                        mult = {}
                        for k in p:
                            mult[k] = mult.get(k, 0) + 1
                        term = const(1, pt.N)
                        cf = Fr(1)
                        for k, mk in mult.items():
                            for _ in range(mk):
                                term = term * Wk[k]
                            cf /= Fr(math.factorial(k)) ** mk * math.factorial(mk)
                        term = term.scale(F(cf))
                        inner = term if inner is None else inner + term
                    if inner is None:
                        continue
                    t = (inner * W2p).scale(pref)
                    acc = t if acc is None else acc + t
            k, pr = s_cmp(acc, pt.rpt[i][m])
            if k is not None:
                ok = False; det = "colour %d m=%d first diff t^%d" % (i, m, k); break
        if not ok: break
    rec(pt, "sanity: (C7) partition formula = exp-route Wick row, m<=%d, both colours" % mmax, ok, det)

# ---------------------------------------------------------------- points
POINTS = {
    # tag: (a0, a1, a2, e1, e2)   e2 = None -> CY (e2 = E - e1)
    "CY1":      (Fr(3, 2), Fr(-5, 7), Fr(11, 4), Fr(2, 3), None),
    "CY2":      (Fr(-7, 3), Fr(9, 5), Fr(4, 7), Fr(-5, 2), None),
    "CY3":      (Fr(7, 5), Fr(13, 9), Fr(-4, 11), Fr(3, 8), None),
    "nonCY1":   (Fr(2, 3), Fr(7, 5), Fr(-3, 4), Fr(5, 2), Fr(1, 7)),
    "nonCY2":   (Fr(-5, 4), Fr(2, 9), Fr(3), Fr(1, 3), Fr(-4, 5)),
    "a1=a2":    (Fr(1, 2), Fr(3, 5), Fr(3, 5), Fr(4, 7), None),
    "u=0":      (Fr(1, 3), Fr(11, 15), Fr(2, 5), Fr(3, 7), None),      # E = 2 a1
    "e1=e2":    (Fr(3, 4), Fr(-2, 5), Fr(7, 3), None, None),           # e1 = e2 = E/2
    "a1=a2-nCY": (Fr(-2, 3), Fr(5, 4), Fr(5, 4), Fr(1, 6), Fr(3, 11)),
}

def mkpoint(tag, N, Z):
    a0, a1, a2, e1, e2 = POINTS[tag]
    E = a0 + a1 + a2
    if e1 is None:
        e1 = E / 2
    if e2 is None:
        e2 = E - e1
    return Point(a0, a1, a2, e1, e2, N, Z, tag=tag)

if __name__ == "__main__":
    import reflib
    if len(sys.argv) > 4:
        reflib.P = int(sys.argv[4])
        P = reflib.P
        print("using prime p = %d" % P)
    t0 = time.time()
    tags = list(POINTS) if which == "all" else which.split(",")
    for tag in tags:
        if tag in POINTS:
            pt = mkpoint(tag, N, Z)
            run_point(pt)
            if tag in ("CY1", "nonCY1"):
                wick_partition_check(pt, 4)
            print("  (elapsed %.1fs)" % (time.time() - t0), flush=True)
    if which in ("all", "controls"):
        control_frames(Fr(3, 2), Fr(-5, 7), Fr(11, 4), Fr(2, 3), Fr(3, 2) + Fr(-5, 7) + Fr(11, 4) - Fr(2, 3), N, Z)
    if which in ("all", "cor12"):
        cor12(Fr(5, 3), Fr(-7, 4), Fr(2, 9), Fr(5, 3) + Fr(-7, 4) - Fr(2, 9), N, Z)   # CY, a0 = 0
        cor12(Fr(-3, 8), Fr(5, 2), Fr(7, 5), Fr(1, 3), N, Z)                        # non-CY, a0 = 0
    if which in ("all", "rem01"):
        a0, a1, a2 = Fr(3, 2), Fr(-5, 7), Fr(11, 4)
        E = a0 + a1 + a2
        remark01(a0, a1, a2, Fr(2, 3), E - Fr(2, 3), Fr(7, 9), Fr(-13, 5), N, Z)   # CY1 vs a non-CY e-pair
    print("\n\nSUMMARY (N=%d, Z=%d): %d checks, %d failed, %.1fs" % (N, Z, len(RESULTS), sum(1 for r in RESULTS if not r[2]), time.time() - t0))
    for r in RESULTS:
        if not r[2]:
            print("  FAIL %s: %s  %s" % r[:2] + (r[3],) if False else "  FAIL %-10s %s   %s" % (r[0], r[1], r[3]))
