# -*- coding: utf-8 -*-
"""genus2_symplectic.py -- constant symplectic (colour-mixing) factorisation R = R' M of the formal-quintic R-matrix at genus 2:
can the polar part of each stable graph be reduced to c_Gamma t_c^-2?  Modules fq_sser (rows in s, M), fq_mvertex (Omega^M = M.omega),
fq_msum (outer sum R'.Omega^M), fq_mfit (interpolation), fq_msolve / fq_mfourier (solve, Z/5 structure), fq_mreport (validation).
Run:  py -I genus2_symplectic.py   (writes logs/genus2_symplectic.log and symplectic_result.json).  Two primes.
Steps: 1 s-series pipeline validated against the point frame; 2 vertex tensors; 3 outer sum (F4 point, random diagonal frames via the
N-ring route, random full M: M-independent total); 4 polar coefficients s^-10..s^-1 of the 7 colour-summed graphs as polynomials in the
25 (+15 m3) entries of M (161 monomials, 170 + 10 points); 5 solves (full, edge graphs, s^-5 only, subleading only) and the Fourier structure.
"""
import os, sys, time, random, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_rows import Rows
from fq_graphs import Decorated
from fq_conifold import Conifold
from gsum_q import ser_fmt
from fq_sser import SRows, SympM, rows_prime, EdgeLeaf
from fq_msum import GraphSummer, outer_graphs, sum_series
from fq_mfit import MONOS, MONO_W, VNAME, PAIRS1, PAIRS2, params_of, fit, poly_eval, mono_name, conj_coeffs
from fq_msolve import solve, basis_qzeta, fmt_basis_qzeta, poly_qzeta, M1, M2
from fq_mfourier import structure, fmt_mat
from fq_scalar_lib import qzeta_coords, qzeta_ratrec, qzeta_fmt
import fq_mreport as rep

os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
LOG = open(os.path.join(HERE, "logs", "genus2_symplectic.log"), "w", encoding="utf-8")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
FAILED = []
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)

ORDERS = list(range(-10, 0)); NSAMP = 170; NVER = 10; PREC_ABS = 24
SYSTEMS = [("full (all 7 graphs, s^-9..s^-1, target_-5 = -(7/5) c_Gamma)", range(1, 8), range(-9, 0)),
           ("(a) edge graphs 2..7 only", range(2, 8), range(-9, 0)),
           ("(b) s^-5 targets only, all graphs", range(1, 8), [-5]),
           ("subleading only: s^-6, s^-2, s^-1 = 0 for all 7 graphs (no s^-5)", range(1, 8), [k for k in range(-9, 0) if k != -5]),
           ("subleading only, edge graphs 2..7", range(2, 8), [k for k in range(-9, 0) if k != -5])]


def shift_orbit_reps():
    """representatives of the Z/5 colour-shift orbits of the monomials (weight <= 2)"""
    def shift(v):
        if v < 15: i, j = PAIRS1[v]; return PAIRS1.index(tuple(sorted(((i + 1) % 5, (j + 1) % 5))))
        if v < 25: i, j = PAIRS2[v - 15]; a, b = (i + 1) % 5, (j + 1) % 5; return 15 + PAIRS2.index((min(a, b), max(a, b)))
        i, j = PAIRS1[v - 25]; return 25 + PAIRS1.index(tuple(sorted(((i + 1) % 5, (j + 1) % 5))))
    reps = []; seen = set()
    for k, mu in enumerate(MONOS):
        if k in seen: continue
        orb = []; cur = mu
        for _ in range(5):
            kk = MONOS.index(tuple(sorted(shift(v) for v in cur))) if cur else 0
            orb.append(kk); cur = MONOS[kk]
        seen.update(orb); reps.append(k)
    return reps


def run(p, results):
    R = FpRing(p); T0 = time.time(); rng = random.Random(p % 1000003 + 7)
    log("=" * 100); log("prime p = %d" % p)
    rows = Rows(R, Nz=6); ring = rows.ring; C = Conifold(R, Nu=10); S = SRows(rows, C, prec_abs=PREC_ABS); graphs = outer_graphs()
    log("rows, conifold periods (u-precision %d), rows truncated at O(s^%d): %.1f s" % (min(x.prec for x in S.imu.values()), PREC_ABS, time.time() - T0))
    log("-" * 100); log("(1) s-series pipeline (normalised basis), point frame M = 1")
    res1, Tref = rep.step1(S, rows, graphs, log, check)
    log("-" * 100); log("(2) vertex tensors Omega^M = M.omega (fq_mvertex)")
    rep.step2(R, log, check, rng)
    log("-" * 100); log("(3) outer sum R'.Omega^M (fq_msum)")
    deco = Decorated(2, hodge=True)
    rep.step3(S, rows, deco, graphs, log, check, rng)
    M = SympM.random(R, rng); Qp = rows_prime(S, M); EL = EdgeLeaf(S, Qp, M)
    check("random full M: T'_d = (-1)^d Q'_{0,d-1} equals -[z^{d-1}] M(z) Q_0(-z) (R'^-1 = M R^-1), d = 2..4",
          all(EL.T(d)[i].eq(EL.T_alt(d)[i], min(EL.T(d)[i].prec, EL.T_alt(d)[i].prec) - 1) for d in (2, 3, 4) for i in range(5)))
    S2 = SRows.__new__(SRows); S2.__dict__.update(S.__dict__); S2.Q = Qp
    check("random full M: the rows Q' = Q M^-1 are symplectic to z^3", not S2.symplectic_check(3))
    t0 = time.time(); resM = GraphSummer(S, M).all_graphs(graphs); TM = sum_series(list(resM.values()))
    log("   random full M: graph sums in %.1f s; total prec O(s^%d)" % (time.time() - t0, TM.prec))
    for idx in sorted(resM): log("   graph %d: %s" % (idx, ser_fmt(resM[idx], R, kmax=-1)))
    check("random full M: TOTAL F_2 series M-independent through s^10 (= point-frame total)", TM.eq(Tref, 11))
    check("random full M: leading s^-10 coefficient of every graph M-independent", all(resM[i].coef(-10) == res1[i].coef(-10) for i in resM))
    check("random full M: s^-9, s^-8, s^-7 coefficients of every graph zero", all(resM[i].coef(k) == 0 for i in resM for k in (-9, -8, -7)))
    check("random full M: s^-5 coefficient of every graph M-independent (weight classes), s^-4, s^-3 zero",
          all(resM[i].coef(-5) == res1[i].coef(-5) and resM[i].coef(-4) == 0 and resM[i].coef(-3) == 0 for i in resM))
    # ---- (4) interpolation
    log("-" * 100); log("(4) polar coefficients s^-10..s^-1 as polynomials in M: %d monomials (1, m1 (15), m2 (10), m3 (15), m1 m1 (120)); %d + %d random symplectic M" % (len(MONOS), NSAMP, NVER))
    keys = [(G['idx'], k) for G in graphs for k in ORDERS]; t0 = time.time(); samples = []
    for t in range(NSAMP + NVER):
        Mt = SympM.random(R, rng); r = GraphSummer(S, Mt).all_graphs(graphs)
        samples.append((params_of(Mt), {(idx, k): r[idx].coef(k) for idx in r for k in ORDERS}))
    log("   %d evaluations in %.1f s" % (len(samples), time.time() - t0))
    polys, rank, nbad = fit(samples[:NSAMP], keys, p)
    check("fit: Vandermonde rank %d = number of monomials, %d inconsistent rows" % (rank, nbad), rank == len(MONOS) and nbad == 0)
    check("fit verified on %d extra random M for all 70 (graph, order) coefficients" % NVER, all(poly_eval(polys[k], par, p) == vals[k] for (par, vals) in samples[NSAMP:] for k in keys))
    kinds = {}
    for (G, k), poly in polys.items():
        ks = sorted(set(("1" if not MONOS[m] else "m3" if any(v >= 25 for v in MONOS[m]) else "m2" if any(15 <= v < 25 for v in MONOS[m]) else "m1" if len(MONOS[m]) == 1 else "m1m1") for m in poly))
        kinds[(G, k)] = "+".join(ks) if ks else "0"
    for G in range(1, 8): log("   graph %d: " % G + "  ".join("s^%d: %s" % (k, kinds[(G, k)]) for k in ORDERS if kinds[(G, k)] != "0") + "  (other orders identically 0)")
    check("monomial content: s^-10, s^-5 constant; s^-6, s^-1 linear in m1; s^-2 = m1 m1 + m2; s^-9,-8,-7,-4,-3 zero; m3 never",
          all(kinds[(G, -10)] == "1" and kinds[(G, -5)] == "1" and kinds[(G, -6)] == "m1" and kinds[(G, -1)] in ("m1", "0") and kinds[(G, -2)] in ("m1m1+m2", "m1m1")
              and all(kinds[(G, k)] == "0" for k in (-9, -8, -7, -4, -3)) for G in range(1, 8)))
    # Galois reconstruction of every coefficient in Q(zeta); orbit structure c_{shift mu} = zeta^w c_mu
    qz = {}; nrec = 0; nshift = 0
    for key, poly in polys.items():
        qz[key] = poly_qzeta(R, polys, key)
        for k, cs in qz[key].items():
            if qzeta_ratrec(R, cs) is None: nrec += 1
    reps_ = shift_orbit_reps()
    if p == P2: check("all %d polynomial coefficients reconstruct in Q(zeta) (Galois conjugates sigma_m(c_mu) = c_{pi_m mu})" % sum(len(v) for v in qz.values()), nrec == 0, nrec)
    else: log("   (p = 2^61-1: %d of %d coefficient coordinates beyond the reconstruction bound ~1.07e7 of fq_fp.ratrec; compared by residues in (6))" % (nrec, sum(len(v) for v in qz.values())))
    # ---- (5) solves and structure
    log("-" * 100); log("(5) the solve.  c_Gamma (s^-10), s^-5 coefficient, and the t_c^-1 coefficient c' = c_-5 + (7/5) c_Gamma (M-INDEPENDENT):")
    c10 = {G: polys[(G, -10)].get(0, 0) for G in range(1, 8)}
    for G in range(1, 8):
        c5 = polys[(G, -5)].get(0, 0)
        log("   graph %d: c = %s, c_-5 = %s, t_c^-1 coefficient %s, target c_-5 = -(7/5)c = %s" % (G, R.fmt(c10[G]), R.fmt(c5), R.fmt((c5 + R.of(Fr(7, 5)) * c10[G]) % p), R.fmt((-R.of(Fr(7, 5)) * c10[G]) % p)))
    st = structure(polys, R)
    check("Z/5 structure: weight-1 coefficients in the modes L_{d,1} only, m2 part in B_{d,2} only, m1 m1 part only in L_{d,a} L_{d',a'} with a+a' = 2 mod 5", st['classes_ok'])
    log("   rank of the 14 weight-1 forms (s^-6, s^-1 of 7 graphs) in the 3 modes L_{d,1} = Sum_i zeta^i m1_{i,i+d}: %d;  m2 forms beta^Gamma = (coeff of B_{1,2}, B_{2,2}): rank %d" % (st['lin_rank'], st['beta_rank']))
    log("   coefficients (Q(zeta)) of m1_00, m1_01, m1_02 in the s^-6 and s^-1 forms  [L^Gamma = Sum_d c_{0d} L_{d,1}]; of m2_01, m2_02 in s^-2:")
    for G in range(1, 8):
        parts = []
        for k in (-6, -1):
            cs = [qzeta_fmt(R, qz[(G, k)].get(MONOS.index((PAIRS1.index((0, d)),)), (0, 0, 0, 0))) for d in range(3)]
            parts.append("s^%d: (%s)" % (k, ", ".join(cs)))
        cs2 = [qzeta_fmt(R, qz[(G, -2)].get(MONOS.index((15 + PAIRS2.index((0, d)),)), (0, 0, 0, 0))) for d in (1, 2)]
        log("   graph %d: %s;  s^-2 m2: (%s)" % (G, "; ".join(parts), ", ".join(cs2)))
    if st['beta_rank'] == 1:
        g0 = st['pivot']; b = st['beta'][g0]; kap = b[1] * R.inv(b[0]) % p if b[0] else None
        log("   all m2 forms proportional to B_{1,2} + kappa B_{2,2}, kappa = %s (conjugates: %s)" % (R.fmt(kap), [R.fmt(st['beta'][G][1] * R.inv(st['beta'][G][0]) % p) for G in range(1, 8) if st['beta'][G][0]]))
    log("   quadratic part of s^-2 in the modes x_d = L_{d,0}, y_d = L_{d,2}, z_d = L_{d,3}, w_d = L_{d,4}: Q^Gamma = x^T C y + z^T E w (F_p residues)")
    for G in range(1, 8): log("   graph %d: C = %s  E = %s" % (G, fmt_mat(R, st['C'][G]), fmt_mat(R, st['E'][G])))
    log("   Q(zeta) coefficients of the orbit representatives of the m1 m1 monomials (c_{shift^s mu} = zeta^{2s} c_mu):")
    for G in range(1, 8):
        ent = ["%s: %s" % (mono_name(MONOS[k]), qzeta_fmt(R, qz[(G, -2)][k])) for k in reps_ if k in qz[(G, -2)] and len(MONOS[k]) == 2]
        log("   graph %d (%d): " % (G, len(ent)) + "; ".join(ent))
    log("   reduced quadrics on K (pivot graph %s): %s independent; local dimension of their common zero locus in K (dim 12) at random points of the linear subspaces: %s"
        % (st.get('pivot'), st.get('n_quadrics'), st.get('jacobian_dims')))
    sols = {}
    for label, gs, ords in SYSTEMS:
        r = solve(polys, list(gs), list(ords), c10, R); sols[label] = r
        log("   == system %s" % label)
        if r['obstructed']: log("      INCONSISTENT at the constant (M-independent) equations: %s" % ", ".join("(graph %d, s^%d, residual %s)" % (G, k, R.fmt(x)) for (G, k, x) in r['obstructed']))
        else: log("      no constant obstruction")
        log("      linear part (s^-6, s^-1): rank %d, consistent %s, kernel K dim %d;  s^-2: rank of the m2 coefficient matrix %s, quadratic obstruction on K: %s" %
            (r['lin_rank'], r['lin_consistent'], r['dimK'], r['rankB'], "none" if r['quad_ok'] else "%d non-zero restricted quadratic forms" % len(r.get('quad_obstruction', []))))
        if r['dimK'] < 15:
            ok, rws = basis_qzeta(R, r['K'], M1); log("      RREF basis of K in Q(zeta) (Galois-consistent pivots: %s):" % ok)
            if rws:
                for line in fmt_basis_qzeta(R, rws, M1): log(line)
        if r.get('N2') is not None and r['rankB']:
            ok2, rws2 = basis_qzeta(R, r['N2'], M2); log("      m2 kernel (b . m2 = 0), dim %d, RREF basis in Q(zeta) (pivots ok %s):" % (len(r['N2']), ok2))
            if rws2:
                for line in fmt_basis_qzeta(R, rws2, M2): log(line)
    log("time for prime %d: %.1f s" % (p, time.time() - T0))
    results[p] = dict(polys=polys, qz=qz, st=st, sols=sols, c10=c10)


if __name__ == "__main__":
    results = {}; T0 = time.time()
    run(P2, results); run(P1, results)
    log("=" * 100); log("(6) cross-prime comparison (2^89-1 reconstructed Q(zeta) coordinates vs 2^61-1 residues)")
    A, B = results[P2], results[P1]; R1 = FpRing(P1); bad = 0; n = 0
    for key, poly in B['polys'].items():
        for k, c in poly.items():
            cs = qzeta_ratrec(FpRing(P2), A['qz'][key][k]); n += 1
            if cs is None or sum(R1.of(cs[j]) * R1.z(j) for j in range(4)) % P1 != c: bad += 1
        if set(poly) != set(A['polys'][key]): bad += 1
    check("all %d polynomial coefficients agree at both primes (by residues of the reconstructed Q(zeta) values)" % n, bad == 0, bad)
    check("ranks, dimensions and obstruction sets identical at both primes",
          all(A['sols'][l]['dimK'] == B['sols'][l]['dimK'] and A['sols'][l]['rankB'] == B['sols'][l]['rankB'] and A['sols'][l]['lin_rank'] == B['sols'][l]['lin_rank']
              and [(G, k) for (G, k, x) in A['sols'][l]['obstructed']] == [(G, k) for (G, k, x) in B['sols'][l]['obstructed']] and A['sols'][l]['quad_ok'] == B['sols'][l]['quad_ok'] for l in A['sols'])
          and A['st']['lin_rank'] == B['st']['lin_rank'] and A['st']['beta_rank'] == B['st']['beta_rank'] and A['st']['n_quadrics'] == B['st']['n_quadrics'] and A['st']['jacobian_dims'] == B['st']['jacobian_dims'])
    out = {"prime": P2, "polys_qzeta": {"%d,%d" % key: {mono_name(MONOS[k]): [str(x) for x in (qzeta_ratrec(FpRing(P2), cs) or cs)] for k, cs in v.items()} for key, v in A['qz'].items()},
           "structure": {"lin_rank": A['st']['lin_rank'], "beta_rank": A['st']['beta_rank'], "n_quadrics": A['st']['n_quadrics'], "jacobian_dims": {str(k): v for k, v in A['st']['jacobian_dims'].items()}},
           "systems": {l: {"obstructed": [(G, k) for (G, k, x) in s['obstructed']], "dimK": s['dimK'], "rankB": s['rankB'], "quad_ok": s['quad_ok'], "dim_solution": s['dim_solution']} for l, s in A['sols'].items()}}
    json.dump(out, open(os.path.join(HERE, "symplectic_result.json"), "w"), indent=1)
    log("total time %.1f s" % (time.time() - T0))
    log("ALL SYMPLECTIC CHECKS PASSED" if not FAILED else "SYMPLECTIC FAILURES: %s" % FAILED)
    LOG.close()
