# -*- coding: utf-8 -*-
"""genus2_symbol.py -- the CONIFOLD SYMBOL R^(zhat) of the formal-quintic R-matrix and the leading-order system for a q-dependent
factorisation R = R' M_0(zhat), M_0 = exp(mu_1 zhat + mu_2 zhat^2 + mu_3 zhat^3), zhat = z/s^4, at genus 2.
Modules: fq_symbol (the symbol, colour-colour form, logs), fq_leadfit (weight <= 3 fit, 991 monomials), fq_leadsolve (L_Gamma, frozen theory,
classes), fq_leadclass / fq_leadmodes / fq_leadline (exact solving: resultants, line, blow-up), fq_leadcol (per-coloured diagonal
systems), fq_wsearch (Hodge weight search), fq_symreport1 / fq_symreport3 (report blocks).
Run:  py -I genus2_symbol.py   (writes logs/genus2_symbol.log and symbol_result.json).  Two primes (2^89-1 first, then 2^61-1), then a
cross-prime block.  Steps: 1 the symbol (valuations, symplecticity, covariance phase, Teleman form, exp(odd) tests, f, QRR relation);
2 the s-series fit over all weight <= 3 monomials (1010 + 20 random symplectic M); 3 the leading system (frozen-theory check, tautological
point, Jacobians, W, the classes diagonal / covariant diagonal / circulant / per-coloured).
"""
import os, sys, time, random, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_conifold import Conifold
from fq_sser import SRows, SympM
from fq_msum import GraphSummer, outer_graphs
from fq_leadfit import MONOS3, MONO3_W, mono_name3
from fq_scalar_lib import qzeta_ratrec, qzeta_coords
import fq_symreport1 as rep1
import fq_symreport3 as rep3

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_symbol.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

NZ = 6; PREC_ABS = 24


def run(p, results):
    R = FpRing(p); T0 = time.time(); rng = random.Random(p % 1000003 + 2026)
    log("=" * 100); log("prime p = %d" % p)
    rows = Rows(R, Nz=NZ); C = Conifold(R, Nu=10); S = SRows(rows, C, prec_abs=PREC_ABS); graphs = outer_graphs()
    log("rows (Nz = %d), conifold periods, rows truncated at O(s^%d): %.1f s" % (NZ, PREC_ABS, time.time() - T0))
    log("-" * 100); log("(1) the conifold symbol R^(zhat), zhat = z/s^4: R^_{a,n}^{(i)} = [s^{-4n+c_a-1}] Q^{(i)}_{a,n}, through zhat^%d" % NZ)
    sym, d1 = rep1.step1(S, log, check, R)
    log("-" * 100); log("(2) s-series fit over ALL monomials of weight <= 3 in (m_1, m_2, m_3): 1 + 40 + 120 + 150 + 680 = 991 unknowns")
    polys, qzd, ver = rep1.step2(S, graphs, log, check, R, rng, SympM=SympM, GraphSummer=GraphSummer)
    log("-" * 100); log("(3) the leading-order system L_Gamma(mu) for M_0(zhat) = exp(mu_1 zhat + mu_2 zhat^2 + mu_3 zhat^3)")
    d3 = rep3.step3(polys, sym, graphs, log, check, R, rng)
    log("time for prime %d: %.1f s" % (p, time.time() - T0))
    results[p] = dict(sym=sym, d1=d1, polys=polys, qzd=qzd, d3=d3, R=R)


def cross(results):
    A = results[P2]; B = results[P1]; RA = A['R']; RB = B['R']; pB = RB.p
    log("=" * 100); log("(4) cross-prime comparison (2^89-1 reconstructed Q(zeta)/Q values vs 2^61-1 residues)")
    bad = 0; n = 0
    for key, poly in B['polys'].items():
        if set(poly) != set(A['polys'][key]): bad += 1
        for k, c in poly.items():
            cs = qzeta_ratrec(RA, A['qzd'][key][k]); n += 1
            if cs is None or sum(RB.of(cs[j]) * RB.z(j) for j in range(4)) % pB != c: bad += 1
    check("all %d fitted polynomial coefficients (weight <= 3, 13 orders, 7 graphs) agree at both primes" % n, bad == 0, bad)
    CA = A['d1']['C']; CB = B['d1']['C']; bad = 0; n = 0
    for nn in range(6):
        for j in range(5):
            for i in range(5):
                cs = qzeta_ratrec(RA, qzeta_coords(RA, [CA[nn][(m * j) % 5][(m * i) % 5] for m in range(1, 5)])); n += 1
                if cs is None or sum(RB.of(cs[k]) * RB.z(k) for k in range(4)) % pB != CB[nn][j][i]: bad += 1
    check("symbol colour-colour matrices C_0..C_5 (%d entries) agree at both primes" % n, bad == 0, bad)
    fA = A['d1']['f']; fB = B['d1']['f']; bad = 0
    for k in range(len(fA)):
        r = RA.ratrec(fA[k], 10 ** 13)
        if r is None or RB.of(r) != fB[k]: bad += 1
    check("identity-row symbol f (zhat^0..%d, large-bound reconstruction at 2^89-1) agrees at 2^61-1" % (len(fA) - 1), bad == 0, bad)
    same = (A['d3'].get('gen') == B['d3'].get('gen') and A['d3']['rk_taut'] == B['d3']['rk_taut'] and A['d3']['rBT'] == B['d3']['rBT']
            and A['d3']['cov_gcd_deg'] == B['d3']['cov_gcd_deg'] and A['d3'].get('fib') == B['d3'].get('fib') and A['d3'].get('percolour') == B['d3'].get('percolour')
            and A['d3'].get('ell_dir') == B['d3'].get('ell_dir') and A['d3'].get('X') == B['d3'].get('X') and A['d3'].get('X_fibre') == B['d3'].get('X_fibre'))
    check("step-3 ranks, dimensions, gcd degrees, the line direction and the point X (Q(zeta)) identical at both primes", same)
    out = {"prime": P2, "symbol_colour_matrices_qzeta": [[[rep1.qz(RA, [CA[nn][(m * j) % 5][(m * i) % 5] for m in range(1, 5)]) for i in range(5)] for j in range(5)] for nn in range(6)],
           "f": [rep1.ratbig(RA, x) for x in fA], "log_f_odd_even": {k: rep1.ratbig(RA, v) for k, v in rep1.exp_odd_test(fA, RA)[0].items()},
           "lead_polys_qzeta": {str(G): {mono_name3(MONOS3[k]): [str(x) for x in (qzeta_ratrec(RA, A['qzd'][(G, {0: -10, 1: -6, 2: -2, 3: 2}[MONO3_W[k]])][k]) or ())] for k in A['d3']['L'][G]} for G in range(1, 8)},
           "ell_direction": A['d3'].get('ell_dir'), "X": A['d3'].get('X'), "X_fibre": A['d3'].get('X_fibre'), "fibres_over_ell": A['d3'].get('fib'),
           "diag_conditions": [[RA.fmt(x) for x in c] for c in A['d3']['diag_conditions']], "generic_W": A['d3'].get('gen')}
    json.dump(out, open(os.path.join(HERE, "symbol_result.json"), "w"), indent=1)


if __name__ == "__main__":
    results = {}; T0 = time.time()
    run(P2, results); run(P1, results); cross(results)
    log("total time %.1f s" % (time.time() - T0))
    log("ALL SYMBOL CHECKS PASSED" if not FAILED else "SYMBOL FAILURES: %s" % FAILED)
    LOG.close()
