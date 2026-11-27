# -*- coding: utf-8 -*-
"""genus2_ou.py -- the O(u) correction of the tautological factorisation of the formal-quintic R-matrix at genus 2:
M(zhat; u) = M_0 exp(u N(zhat)) with M_0 = R~ mod zhat^4 (the conifold symbol in Teleman's colour-colour form), N in sp (40 coordinates).
Modules: fq_hat (the hatted theory: HatRows, RawSympM, products, t_c conversion, interpolation), fq_ou (directions, the 2 x 40 system J,
classes right/left, solving, verification by eps -> u), fq_ourec (Q / Q(zeta) reconstruction, cross-prime record), fq_oureport /
fq_oureport2 (report blocks).  Run:  py -I genus2_ou.py  (writes logs/genus2_ou.log and ou_result.json).  Two primes (2^89-1, then
2^61-1), then a cross-prime block.  Steps per prime: 1 HatRows validation and the point frame in t_c; 2 M_0 and g1_Gamma = [u^1] G_Gamma[M_0];
3 the directional derivatives D_N L_Gamma (40 directions, eps = 0..9, frozen theory) -> J, cross-checked against the stored leading
polynomials via the dexp map; 4 the linear system J n = -(g1_2, g1_4) on the classes (i)-(iv), right and left placement; 5 verification
of one solution per solvable class by the eps = u substitution (all six edge graphs O(u^2), vertex graph (1/1200) t_c^-2 + 0 t_c^-1);
6 the same at the conifold-frame ambiguity (alpha, beta) = (1/3, -2/7).  Genus 2 only; exact finite computation; AI-checked only.
"""
import os, sys, time, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_conifold import Conifold
from fq_sser import SRows
from fq_msum import outer_graphs
from fq_symbol import Symbol
from fq_hat import HatRows, tc_series
from fq_ou import GRAPHS, frozen_eps_poly, direction_N
from fq_leadsolve import FrozenRows
from fq_ourec import Rec, cross, qzfmt
import fq_oureport as rep
import fq_oureport2 as rep2

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_ou.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

NZ = 6; PREC_ABS = 24; AB2 = (Fr(1, 3), Fr(-2, 7))
JSON_IN = os.path.join(HERE, "symbol_result.json")


def run(p, results):
    R = FpRing(p); T0 = time.time(); rec = Rec(R)
    log("=" * 100); log("prime p = %d" % p)
    rows = Rows(R, Nz=NZ); C = Conifold(R, Nu=10); S = SRows(rows, C, prec_abs=PREC_ABS); graphs = outer_graphs()
    sym = Symbol(S); H = HatRows(S); tc = tc_series(C, R)
    log("rows (Nz = %d), conifold periods (alpha = beta = 0), rows truncated at O(s^%d), symbol, HatRows: %.1f s" % (NZ, PREC_ABS, time.time() - T0))
    log("-" * 100); log("(1) the hatted theory at (alpha, beta) = (0, 0): Qhat = Q s^{4n+1-c_a}, xhat = x/s; checks (a), (b), (c); point frame in t_c")
    d1 = rep.step1(S, H, sym, tc, graphs, log, check, R, rec, "ab00")
    log("-" * 100); log("(2) M_0 = R~ mod zhat^4 and the hatted sums G_Gamma[M_0](u): g1_Gamma = [u^1]")
    d2 = rep.step2(H, sym, tc, graphs, log, check, R, rec, "ab00", d1['tot'])
    log("-" * 100); log("(3) directional derivatives D_N L_Gamma = d/d eps|_0 G^0_Gamma[M_0 exp(eps N)]: the 2 x 40 matrix J (graphs 2, 4)")
    d3 = rep.step3(sym, d2, graphs, log, check, R, rec, JSON_IN)
    log("-" * 100); log("(4) solving J n = -(g1_2, g1_4): full class and the classes (i)-(iii), RIGHT (M_0 exp(uN)) and LEFT (exp(uN') M_0) placement")
    sols = rep2.step4(d3['J'], d2['g1'], d2['M0'], log, check, R, rec, "ab00")
    log("-" * 100); log("(5) verification by eps -> u: hatted sums with M_0 exp(eps N) at eps = 0..9, P_Gamma(eps; u), eps = u")
    v5 = rep2.step5(H, d2['M0'], sols, graphs, d3['Vinv'], log, check, R, rec, "ab00", d1['tot'], tc.coef(2))
    log("   full s-series picture: every hatted factor is a u-series, so after eps = u every colour-summed graph term is s^-10 x (u-series):")
    log("   the s^k coefficients with 5 not dividing k vanish identically for every graph and every N (checked: is_u_series on all series above).")
    log("-" * 100); log("(6) the conifold-frame ambiguity: varpi^cf_0 = f_0 + alpha g + beta f_2 at (alpha, beta) = (%s, %s)" % AB2)
    C2 = Conifold(R, Nu=10, alpha=AB2[0], beta=AB2[1]); S2 = SRows(rows, C2, prec_abs=PREC_ABS); sym2 = Symbol(S2); H2 = HatRows(S2); tc2 = tc_series(C2, R)
    log("   varpi^cf_0 = %s;  t_c = %s" % (C2.varpi0.fmt(3), tc2.fmt(3)))
    check("ab13 the conifold SYMBOL is independent of (alpha, beta): R^_{a,n}^{(i)} identical for all (a, i, n), n <= %d (so J is unchanged)" % NZ, sym2.Rh == sym.Rh)
    check("ab13 xhat = varpi^cf_0(u) changes at O(u) (xhat_1 = alpha = %s) while xhat(0) = 1" % R.fmt(H2.x.coef(5)), H2.x.coef(0) == 1 and H2.x.coef(5) == R.of(AB2[0]) and H.x.coef(5) == 0)
    d1b = rep.step1(S2, H2, sym2, tc2, graphs, log, check, R, rec, "ab13", brief=True)
    d2b = rep.step2(H2, sym2, tc2, graphs, log, check, R, rec, "ab13", d1b['tot'])
    check("ab13 M_0 (colour matrices C_0..C_3) identical to that of (0, 0)", d2b['C'][:4] == d2['C'][:4])
    pol = frozen_eps_poly(FrozenRows(sym2), d2b['M0'].M, d2b['M0'].Minv, direction_N(7, p), graphs, R, d3['Vinv'])
    check("ab13 J column 7 re-derived from FrozenRows(sym2) equals the (0, 0) one (J depends on the symbol only)", all(pol[G][1] == d3['J'][G][7] for G in GRAPHS))
    sols_b = rep2.step4(d3['J'], d2b['g1'], d2b['M0'], log, check, R, rec, "ab13")
    v5b = rep2.step5(H2, d2b['M0'], sols_b, graphs, d3['Vinv'], log, check, R, rec, "ab13", d1b['tot'], tc2.coef(2), keys=[('cov', 'right'), ('cov', 'left')], random_full=False)
    log("   summary of the (alpha, beta) dependence:")
    rep2.frame_summary(dict(g1=d2['g1'], A=d2['A'], A1=d1['A'], sols=sols, c1=tc.coef(2)), dict(g1=d2b['g1'], A=d2b['A'], A1=d1b['A'], sols=sols_b, c1=tc2.coef(2), alpha=R.of(AB2[0])), log, check, R, rec)
    log("time for prime %d: %.1f s" % (p, time.time() - T0))
    struct = dict(J_rank=d3['rank'], J2_vars=d3.get('J2_vars'), J4_vars=d3.get('J4_vars'), dexp_ok=d3.get('dexp_ok'),
                  classes={"%s:%s:%s" % (tag, k[0], k[1]): dict(rank=s['rank'], dim=s['dim'], consistent=s['N'] is not None) for tag, ss in (("ab00", sols), ("ab13", sols_b)) for k, s in ss.items()})
    results[p] = dict(rec=rec, struct=struct, R=R)


def cross_block(results):
    A = results[P2]; B = results[P1]
    log("=" * 100); log("(7) cross-prime comparison: every Q / Q(zeta) value reconstructed at 2^89-1 reduced mod 2^61-1 against the residue there")
    n, bad, norec, missing = cross(A['rec'], B['rec'])
    check("all %d reconstructed values (g1, u- and t_c-coefficients, J rows in Q(zeta), restricted matrices, solutions (b_1, b_3), verification coefficients) agree at both primes" % n, not bad and not missing, (bad[:5], missing[:5]))
    log("   values not reconstructible at 2^89-1 (bound 1e10; compared by nothing, residues logged only): %d %s" % (len(norec), norec[:8]))
    check("step-4 ranks, dimensions and consistency identical at both primes; J rank 2 at both", A['struct']['classes'] == B['struct']['classes'] and A['struct']['J_rank'] == B['struct']['J_rank'] == 2)
    rec = A['rec']
    nbig = sum(1 for name in rec.names() if rec.kind[name] == 'Qbig')
    log("   %d values used the 1e13 fallback bound at 2^89-1; all of them are confirmed by their 2^61-1 residues above" % nbig)
    out = dict(prime=P2, values={name: (str(rec.rec[name]) if rec.kind[name] in ('Q', 'Qbig') else qzfmt(rec.rec[name])) if rec.rec.get(name) is not None else "[%d mod p]" % rec.res[name] for name in rec.names()},
               structure=A['struct'], cross_prime=dict(compared=n, mismatches=bad, not_reconstructed=norec, missing=missing), failed=FAILED)
    json.dump(out, open(os.path.join(HERE, "ou_result.json"), "w"), indent=1)


if __name__ == "__main__":
    results = {}; T0 = time.time()
    run(P2, results); run(P1, results); cross_block(results)
    log("total time %.1f s" % (time.time() - T0))
    log("ALL OU CHECKS PASSED" if not FAILED else "OU FAILURES: %s" % FAILED)
    LOG.close()
