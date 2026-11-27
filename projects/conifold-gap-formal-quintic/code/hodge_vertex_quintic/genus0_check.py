# -*- coding: utf-8 -*-
"""genus0_check.py -- PLAN_2026-11-03.md step 1 (S1) driver.  Run:  py -I genus0_check.py   (writes logs/genus0_check.log).
Checks (each logged PASS/FAIL; the script stops at the first FAIL):
  I-function at generic weights (two samples) and at the zeta point; zeta covariance; PF operator (generic product form
  and the explicit polynomial form at the zeta point) on the q-series; tau0 propto e_1(lambda), zero at the zeta point;
  mirror map b = Itilde_1; structure constants kappa_a (flat basis), Yukawa (closed form and instanton numbers);
  Y_0..Y_8; eigenvalues; Delta_i; the Wronskian relation (derived, printed, verified on exact rational q-series to
  q^20); conifold Frobenius basis, gamma, pairing table, t_c, Cf, Cf preserves the relation; two primes.
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, QRing, P1, P2
from fq_series import Ser
from fq_ifun import IFun, vandermonde, components, lin
from fq_genus0 import Genus0, quintic_C0, quintic_It1
from fq_ring import Ring, NEl, NAMES, IDX
from fq_conifold import Conifold

LOG = open(os.path.join(HERE, "logs", "genus0_check.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok:
        log("STOP at first failure"); LOG.close(); sys.exit(1)

def apply_pf_explicit(R, F, lam, z_prec):
    """[(M+lam)^5 - 1 - q(3125 (M+lam)^5 + 9375 z (M+lam)^4 + 10625 z^2 (M+lam)^3 + 5625 z^3 (M+lam)^2 + 1370 z^4 (M+lam)
       + 120 z^5)] on the q-series F (list over d of Laurent Ser in z); M q^d = d z q^d."""
    D = len(F) - 1; out = []
    def Mpow(f, d, k):
        for _ in range(k): f = f.mul(lin(R, lam, d, f.prec + 2))
        return f
    for d in range(D + 1):
        t = Mpow(F[d], d, 5).sub(F[d])
        if d >= 1:
            g = F[d - 1]; e = d - 1
            for k, c in ((5, 3125), (4, 9375), (3, 10625), (2, 5625), (1, 1370), (0, 120)):
                t = t.sub(Mpow(g, e, k).shift(5 - k).scale(c))
        out.append(t.truncate(z_prec))
    return out

def run(p, second=False):
    R = FpRing(p); t0 = time.time()
    log("=" * 78); log("prime p = %d, zeta = %d" % (p, R.zeta))
    D = 12
    # ---- generic weights: two samples
    random.seed(11 + p % 7)
    e1_ratio = None; b_prev = None
    for smp in range(2):
        lams = [random.randrange(2, 10 ** 6) for _ in range(5)]
        G = Genus0.__new__(Genus0)       # only the I-function part is weight-generic: do it by hand
        IF = IFun(R, lams, D)
        for i in range(5):
            res = IF.pf_check(i, 8)
            check("PF (product form) on I_%d, generic sample %d, q^<=%d" % (i, smp, D), all(r.iszero() for r in res))
        V, Vinv = vandermonde(R, lams)
        N = D + 1; C0 = quintic_C0(R, D); C0inv = C0.inv()
        from fq_ifun import SL
        comp1 = None
        J = []
        for i in range(5):
            acc = None
            for d in range(D + 1):
                t = SL.from_sser(R, IF.c_inf(d, i, 3), N).mulq(Ser.mono(R, d, 1, N))
                acc = t if acc is None else acc.add(t)
            J.append(SL.exp_ls(R, lams[i], 2, N).mul(acc, 2).mulq(C0inv))
        c0comp = components(R, Vinv, [J[i].part(0, 0, N) for i in range(5)])
        check("s^0 term of I/(zC0) = 1 (generic %d)" % smp, c0comp[0].eq(Ser.const(R, 1, N)) and all(c.iszero() for c in c0comp[1:]))
        comp1 = components(R, Vinv, [J[i].part(1, 0, N) for i in range(5)])
        comp1l = components(R, Vinv, [J[i].part(1, 1, N) for i in range(5)])
        check("s^1 log q term = H (generic %d)" % smp, comp1l[1].eq(Ser.const(R, 1, N)) and all(comp1l[k].iszero() for k in (0, 2, 3, 4)))
        check("s^1 term has no H^2,H^3,H^4 components (generic %d)" % smp, all(comp1[k].iszero() for k in (2, 3, 4)))
        e1 = sum(lams) % p
        ratio = comp1[0].scale(R.inv(e1))
        if e1_ratio is None: e1_ratio = ratio
        else: check("tau0/e_1(lambda) the same for both weight samples", ratio.eq(e1_ratio))
        log("   tau0/e_1 =", ratio.fmt(3))
        b = comp1[1].mul(C0)
        check("mirror-map coefficient b = Itilde_1 (generic %d)" % smp, b.eq(quintic_It1(R, D)))
    # ---- zeta point
    lams = [R.z(i) for i in range(5)]
    IF = IFun(R, lams, D)
    for i in range(5):
        res = IF.pf_check(i, 8)
        check("PF (product form) on I_%d at the zeta point" % i, all(r.iszero() for r in res))
        F = IF.F_zero(i, 14)
        res2 = apply_pf_explicit(R, [f.shift(1) for f in F], R.z(i), 8)
        check("PF explicit polynomial form [M^5-1-q(3125M^5+9375zM^4+10625z^2M^3+5625z^3M^2+1370z^4M+120z^5)] on I_%d" % i,
              all(r.iszero() for r in res2))
        ok = all(R.eq(IF.c_zero(d, i, 6).coef(k), R.mul(IF.c_zero(d, 0, 6).coef(k), R.z(-i * k))) for d in range(1, 5) for k in range(-d, 6))
        check("I_%d(q,z) = zeta^%d I_0(q, zeta^-%d z) (q^<=4, z^<6)" % (i, i, i), ok)
    G = Genus0(R, lams, D, Ns=7)
    log("C0 =", G.C0.fmt(4)); log("Itilde_1 =", G.It1.fmt(4))
    check("tau0 = 0 at the zeta point", G.tau0.iszero())
    check("b = Itilde_1 at the zeta point (T = log q + Itilde_1/C0 = I_1^Q/I_0^Q)", G.b.eq(G.It1))
    log("T - log q =", G.Tq.fmt(3)); log("C1 = DT =", G.C1.fmt(3))
    check("z^0 terms of z d_T S(H^a) log-free", G.lfree)
    check("H*H^a proportional to H^(a+1 mod 5) (Z/5 charge)", G.offdiag_ok)
    check("closure z d_T S(H^4) = kappa_4 S(1)", G.closure_ok)
    for a in range(5): log("   kappa_%d =" % a, G.kappa[a].fmt(3))
    C1inv = G.C1.inv(); LC = G.L.mul(C1inv)
    check("kappa_0 = kappa_2 = 1", G.kappa[0].eq(Ser.const(R, 1, G.N)) and G.kappa[2].eq(Ser.const(R, 1, G.N)))
    check("kappa_1 = L^5/(C0^2 C1^3) (= Y_3/5)", G.kappa[1].eq(G.L.pow(5).mul(G.C0.pow(2).mul(G.C1.pow(3)).inv())))
    check("kappa_3 = kappa_4 = C0/C1", G.kappa[3].eq(G.C0.mul(C1inv)) and G.kappa[4].eq(G.C0.mul(C1inv)))
    check("Gram eta(v_k,v_l) depends on k+l only", G.gram_ok)
    for m in range(9): log("   Y_%d =" % m, G.Y(m).fmt(3))
    check("Y_0 = Y_1 = Y_2 = Y_4 = Y_5 = Y_6 = Y_7 = 0", all(G.Y(m).iszero() for m in (0, 1, 2, 4, 5, 6, 7)))
    check("Y_3 = 5 L^5/(C0^2 C1^3) to q^%d" % D, G.Y(3).eq(G.L.pow(5).scale(5).mul(G.C0.pow(2).mul(G.C1.pow(3)).inv())))
    check("Y_8 = (L/C1)^5 Y_3", G.Y(8).eq(LC.pow(5).mul(G.Y(3))))
    check("characteristic polynomial x^5 - (L/C1)^5 (prod kappa = (L/C1)^5)", G.kprod.eq(LC.pow(5)))
    Y3Q = G.Qseries(G.Y(3))
    nd = {1: 2875, 2: 609250, 3: 317206375, 4: 242467530000, 5: 229305888887625}
    for d in range(1, 6):
        N0d = sum(Fr(nd[d // k], k ** 3) for k in range(1, d + 1) if d % k == 0)
        check("Yukawa Y_3 = 5 + Sum d^3 N_{0,d} Q^d: d=%d, N_{0,d} = %s" % (d, N0d), R.eq(Y3Q.coef(d), R.of(N0d * d ** 3)))
    for i in range(5):
        check("eigenvalue mu_%d = zeta^%d L/C1" % (i, i), G.mu[i].eq(LC.scale(R.z(i))))
        check("Delta_%d = C0^2 zeta^(3*%d)/L^2, Delta_%d(0) = %s = W_%d^0" % (i, i, i, R.fmt(G.Delta[i].coef(0)), i),
              G.Delta[i].eq(G.C0.pow(2).scale(R.z(3 * i)).mul(G.L.pow(2).inv())) and R.eq(G.Delta[i].coef(0), G.W0[i]))
    for i in range(5):
        for j in range(5):
            if i != j: check("eta(e_%d, e_%d) = 0" % (i, j), G.eta(G.idem[i], G.idem[j]).iszero()) if (i, j) == (0, 1) else None
    # ---- differential ring and the relation
    ring = Ring(R)
    G20 = Genus0(R, lams, 22, Ns=5)
    img = ring.images_from_qseries(G20)
    check("PF reduction D^4 C0 formula on q-series (q^20)", img['c3'].theta().eq(ring.gen('c3').D(ring).subst(img, ring), 21))
    check("PF reduction D^4 Itilde_1 formula on q-series (q^20)", img['i3'].theta().eq(ring.gen('i3').D(ring).subst(img, ring), 21))
    check("DL = (L^6-L)/5 on q-series", G20.L.theta().eq(ring.gen('L').D(ring).subst(img, ring), 21))
    check("C1 = W/c0^2 on q-series", G20.C1.eq(ring.C1().subst(img, ring), 21))
    ring.derive_relation()
    check("relation B(varpi_0, varpi_1) is log-free (Lam-part identically zero)", not ring.rel_lam)
    check("B(varpi_0, varpi_0) = 0 identically (antisymmetry)", ring.concomitant([ring.gen('c%d' % m) for m in range(4)], ring.gen('c0')).iszero())
    terms = []
    relQ = {}
    for k, v in sorted(ring.rel0.num.items()):
        r = R.ratrec(v); assert r is not None
        relQ[k] = r
        mono = " ".join(("%s^%d" % (NAMES[j], k[j]) if k[j] != 1 else NAMES[j]) for j in range(len(k)) if k[j])
        terms.append("(%s) %s" % (r, mono))
    log("RELATION (= 0), %d terms:  " % len(terms) + "  +  ".join(terms))
    log("   i3 elimination: i3 = L^5 * Rest / c0^%d" % ring.i3_c0pow)
    # exact rational verification to q^20
    Q = QRing(); NQ = 21
    C0Q = quintic_C0(Q, NQ); I1Q = quintic_It1(Q, NQ); LQ = Ser.from_list(Q, [1, -3125], 0, NQ + 1).root(5).inv()
    jets = {'L': LQ}
    c = C0Q; i = I1Q
    for j in range(4):
        jets['c%d' % j] = c; jets['i%d' % j] = i; c = c.theta(); i = i.theta()
    tot = Ser.zero(Q, NQ + 1)
    for k, r in relQ.items():
        t = Ser.const(Q, r, NQ + 1)
        for j in range(len(k)):
            if k[j]: t = t.mul(jets[NAMES[j]].pow(k[j]))
        tot = tot.add(t)
    check("relation holds on the EXACT rational q-series to q^20", tot.truncate(21).iszero(), tot.truncate(21).fmt(2))
    # ---- conifold
    C = Conifold(R, Nu=40)
    log("conifold: exponents {0,1,1,2}; indicial P(n) = n(n-1)^2(n-2) checked")
    log("   f0 =", C.f[0].fmt(6)); log("   f1 =", C.f[1].fmt(5)); log("   f2 =", C.f[2].fmt(5))
    log("   gamma = %s;  varpi^cf_1 = g = f1 + gamma f2 =" % C.gamma, C.g.fmt(5))
    log("   log solution g log u + h, h =", C.h.fmt(5))
    check("log-solution obstructions (u^-3..u^-1) vanish for g", all(x == 0 for x in C.log_obstruction))
    check("pairings B(f_a,f_b) are D-constant", C.Bconst_ok)
    log("   B in basis (f0, g, f2):", C.Bmat, ";  B(g, log solution) =", C.Blog.fmt(0))
    check("B(g, analytic) = 0 for all analytic solutions (vanishing period; pairing normalisation vacuous)", C.pairing_vacuous)
    log("   FREEDOM: varpi^cf_0 = f0 + alpha g + beta f2 with (alpha, beta) free and the scale of g free; chosen alpha = beta = 0")
    log("   t_c = g/varpi^cf_0 =", C.tc.fmt(5))
    cf0 = C.cf_images(0); cf1 = C.cf_images(1)
    check("Cf(c1) = D_s Cf(c0) with D = -((1-s^5)/5) s^-4 d/ds", C.Ds(cf0['c0']).eq(cf0['c1'], cf0['c1'].prec - 6))
    check("Cf(i1) = D_s Cf(i0) - Cf(c0)", C.Ds(cf0['i0']).sub(cf0['c0']).eq(cf0['i1'], cf0['i1'].prec - 6))
    v = ring.rel0.subst(cf0, ring)
    check("Cf preserves the relation (Cf(Rel) = 0 as s-series, prec O(s^%d))" % v.prec, v.iszero())
    a0 = ring.C1().subst(cf0, ring); a1 = ring.C1().subst(cf1, ring)
    check("Cf(C1) independent of the Lam value (log-free element)", a0.eq(a1))
    log("   Cf(C1) = Cf(D T) = D t_c =", a0.fmt(25))
    log("time %.1f s" % (time.time() - t0))

if __name__ == "__main__":
    run(P2); run(P1)
    log("ALL S1 CHECKS PASSED"); LOG.close()
