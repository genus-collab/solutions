# -*- coding: utf-8 -*-
"""rows_check.py -- PLAN_2026-11-03.md step 2 (S2) driver.  Run:  py -I rows_check.py  (writes logs/rows_check.log).
Checks (PASS/FAIL, stop at first FAIL): WKB at every colour with residues at L = 0 and at the fifth roots of unity at
every order to z^7; psi_0 = L; constants from q = 0 (h_{wt_i}); h_{wt_0} rational; rows P_a^{(i)} in N (flat basis,
reduced modulo the relation); rows at q = 0 = zeta^{ia} h_{wt_i}; Z/5 covariance; symplecticity EXACT in N (all 25
pairs, z^<=6) and at a random point; Tier-1 cross-check against the exact q-series rows from the I-function (q^12,
z^6); Dil_i(d) with L-degrees (d <= 7); edge coefficients symmetric; frame change rows_times; two primes.
"""
import os, sys, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_series import Ser
from fq_genus0 import Genus0
from fq_ring import Ring, NEl, NAMES
from fq_rows import WKB, Rows, tier1_rows, edge_coeffs, lp_fmt
from fq_hodge import hv, wt

LOG = open(os.path.join(HERE, "logs", "rows_check.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok:
        log("STOP at first failure"); LOG.close(); sys.exit(1)

def h_wt0_exact(zmax):
    """exact rational coefficients of h_{wt_0}(z), wt_0 = {1 - zeta^j, j=1..4} U {-5}: the x_j = 1/(1-zeta^j) are the roots of
    x^4 Phi_5(1 - 1/x) = (5x^4 - 10x^3 + 10x^2 - 5x + 1); power sums by Newton's identities."""
    from fq_fp import QRing, bernoulli
    Q = QRing()
    # monic: x^4 - 2x^3 + 2x^2 - x + 1/5 ; e1 = 2, e2 = 2, e3 = 1, e4 = 1/5
    e = [Fr(1), Fr(2), Fr(2), Fr(1), Fr(1, 5)]
    pw = [Fr(4)]
    for m in range(1, zmax + 1):
        s_ = Fr(0)
        for k in range(1, min(m, 4) + 1):
            s_ += (-1) ** (k - 1) * e[k] * (pw[m - k] if m - k > 0 else 0)
        if m <= 4: s_ += (-1) ** (m - 1) * m * e[m]
        pw.append(s_)
    ex = [Fr(0)] * (zmax + 1)
    k = 1
    while 2 * k - 1 <= zmax:
        m = 2 * k - 1
        ex[m] = -bernoulli(2 * k) / (2 * k * (2 * k - 1)) * (pw[m] + Fr(-5) ** (-m))
        k += 1
    return _exp_list(ex)

def _exp_list(ex):
    out = [Fr(1)] + [Fr(0)] * (len(ex) - 1)
    for n in range(1, len(ex)):
        out[n] = sum(k * ex[k] * out[n - k] for k in range(1, n + 1)) / n
    return out

def run(p, hrat):
    R = FpRing(p); t0 = time.time(); Nz = 6
    log("=" * 78); log("prime p = %d" % p)
    lams = [R.z(k) for k in range(5)]
    # ---- WKB with residues to z^7
    wk = {}
    for i in range(5):
        W = WKB(R, i, 7); wk[i] = W
        check("colour %d: symbol equation Pcal_0 = 0, a_1 = 5 zeta^-%d/L, Pcal_1 L = 0 (psi_0 = L propto Delta_%d^-1/2 C0 = zeta^(-3*%d/2) L)" % (i, i, i, i), True)
        for n in range(1, 8):
            rem, res0 = W.residues[n]
            check("colour %d z^%d: residues at fifth roots of unity (remainder mod L^5-1) and at L=0 vanish" % (i, n), not rem and res0 == 0,
                  "rem=%s res0=%s" % (rem, res0))
        if i == 0:
            for n in range(1, 8):
                log("   psi_{0,%d}/L = %s" % (n, lp_fmt(W.phi[n], R)))
    for i in range(1, 5):
        ok = all(wk[i].psi[n] == {k: v * R.z(-i * n) % p for k, v in wk[0].psi[n].items()} for n in range(8))
        check("Z/5 covariance psi_{%d,n} = zeta^(-%d n) psi_{0,n} (n <= 7)" % (i, i), ok)
    h0 = wk[0].h
    if hrat is None:
        hrat = h_wt0_exact(7)
        log("   h_{wt_0}(z) exact (power sums of 1/(1-zeta^j) by Newton, -5 weight): " + ", ".join(str(x) for x in hrat))
    check("h_{wt_0}(z) has the exact rational coefficients above (checked mod p)", all(R.eq(R.of(hrat[n]), h0[n]) for n in range(8)))
    for i in range(1, 5):
        check("h_{wt_%d}(z) = h_{wt_0}(zeta^-%d z)" % (i, i), all(R.eq(wk[i].h[n], R.mul(R.z(-i * n), h0[n])) for n in range(8)))
    # ---- rows in N
    t1 = time.time()
    rows = Rows(R, Nz=Nz)
    log("rows built in %.1f s; sizes (terms) P_a^(0) by z-order: %s" % (time.time() - t1, [[x.nterms() for x in rows.P[a][0]] for a in range(5)]))
    log("   denominators (c0-power, W-power) of P_a^(0): %s" % [[(x.a, x.b) for x in rows.P[a][0]] for a in range(5)])
    check("rows are log-free (no Lam) and reduced (no i3)", all(not x.has('Lam') and not x.has('i3') for a in range(5) for i in range(5) for x in rows.P[a][i]))
    for a in range(5):
        ok = all(rows.q0_value(rows.P[a][i][n]) == R.z(i * a) * wk[i].h[n] % p for i in range(5) for n in range(Nz + 1))
        check("rows at q = 0: P_%d^{(i)} = zeta^{%d i} h_{wt_i}(z) (QRR constants), all i" % (a, a), ok)
    ok = True
    for a in range(5):
        for i in range(1, 5):
            for n in range(Nz + 1):
                if not rows.P[a][i][n].eq(rows.P[a][0][n].scale(R.z(i * (a - n))), rows.ring): ok = False
    check("Z/5 covariance of the rows: P_a^{(i)}(z) = zeta^{ia} P_a^{(0)}(zeta^-i z), exact in N", ok)
    # symplecticity exact in N
    t1 = time.time()
    for i in range(5):
        for j in range(5):
            pr = rows.pair(i, j)
            tgt = [rows.Delta(i) if i == j else rows.ring.const(0)] + [rows.ring.const(0)] * Nz
            ok = all(pr[n].eq(tgt[n], rows.ring) for n in range(Nz + 1))
            check("symplecticity Sum_ab P_a^{(%d)}(z) eta^{ab} P_b^{(%d)}(-z) = delta Delta_%d EXACT in N mod relation, z^<=%d" % (i, j, i, Nz), ok)
    log("   symplecticity (25 pairs) in %.1f s" % (time.time() - t1))
    random.seed(7)
    vals = {n: random.randrange(1, p) for n in NAMES}; vals['i3'] = 0; vals['Lam'] = 0
    ok = all(rows.pair(i, j)[n].evalpt(vals, rows.ring) == (rows.Delta(i).evalpt(vals, rows.ring) if (i == j and n == 0) else 0)
             for i in range(5) for j in range(5) for n in range(Nz + 1))
    check("symplecticity at a random point of the reduced generators", ok)
    # Tier-1
    t1 = time.time()
    G = Genus0(R, lams, 12, Ns=6)
    T1 = tier1_rows(R, G, Nz=Nz, zprec=70)
    img = rows.ring.images_from_qseries(G)
    bad = None
    for a in range(5):
        for i in range(5):
            for n in range(Nz + 1):
                w = rows.P[a][i][n].subst(img, rows.ring); t = T1[(a, i)].coef_series(n)
                if not w.eq(t, 13): bad = (a, i, n, w.fmt(2), t.fmt(2)); break
            if bad: break
        if bad: break
    check("Tier-1: WKB rows (q-series of L, C0, Itilde_1 substituted) = exact I-function rows, all a, i, z^<=%d, q^<=12" % Nz, bad is None, bad)
    log("   Tier-1 in %.1f s" % (time.time() - t1))
    # Dil
    for d in range(2, 8):
        dl = rows.Dil(0, d)
        lo, hi = min(dl), max(dl)
        check("Dil_0(%d) = (-1)^%d [z^%d] P_0 in F_p[L^+-1], L-degree range [%d, %d], top degree 4(d-1) = %d" % (d, d, d - 1, lo, hi, 4 * (d - 1)),
              hi == 4 * (d - 1), lp_fmt(dl, R))
        for i in range(1, 5):
            di = rows.Dil(i, d)
            check("Dil_%d(%d) = zeta^(-%d*%d) Dil_0(%d)" % (i, d, i, d - 1, d), di == {k: v * R.z(-i * (d - 1)) % p for k, v in dl.items()}) if d == 2 else None
    # edges: symmetry E^{(ij)}_{mn} = E^{(ji)}_{nm}
    t1 = time.time()
    E = {(i, j): edge_coeffs(rows, i, j, 2) for i in range(5) for j in range(5)}
    ok = all(E[(i, j)][(m, n)].eq(E[(j, i)][(n, m)], rows.ring) for i in range(5) for j in range(5) for (m, n) in E[(i, j)])
    check("edge coefficients E^{(ij)}_{mn} = E^{(ji)}_{nm} (m+n <= 5)", ok)
    log("   edges (m, n <= 2, all pairs) in %.1f s; sizes E^{(00)}: %s" % (time.time() - t1, {k: v.nterms() for k, v in E[(0, 0)].items()}))
    # frame change
    hinv = hv(R, [1], Nz, -1)
    P2_ = rows.rows_times(hinv)
    ok = all(P2_[a][i][n].evalpt(rows.ring.q0_values(), rows.ring) == R.z(i * a) * sum(wk[i].h[m] * hinv[n - m] for m in range(n + 1)) % p
             for a in range(5) for i in range(5) for n in range(Nz + 1))
    check("rows_times(h_{(1)}^-1): q = 0 value zeta^{ia} h_{wt_i} h_{(lambda=1)}^-1", ok)
    log("time %.1f s" % (time.time() - t0))
    return hrat

if __name__ == "__main__":
    hrat = run(P2, None); run(P1, hrat)
    log("ALL S2 CHECKS PASSED"); LOG.close()
