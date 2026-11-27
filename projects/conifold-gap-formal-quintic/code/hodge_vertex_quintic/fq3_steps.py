# -*- coding: utf-8 -*-
"""fq3_steps.py -- the steps of the genus-3 O(u^k) programme (genus3_ou.py), each reading/writing the JSON state of one prime.

State file  g3_state_<p>[_<frame>].json : residues mod p of every computed quantity (ints), plus structural data; the driver prints
Q / Q(zeta) reconstructions at the reporting prime and compares residues across primes (genus3_ou.py cross block).
Graph indices: Mathematica order of StableGraph-3.bin (fq3_graphs.outer_graphs3): 1 = genus-3 vertex (E = 0); E = 1: 2 (loop on a
genus-2 vertex), 5 (bridge g1 - g2); E = 2: 3, 6, 7, 8, 17; E = 3: 4, 9-12, 18-20, 32; E >= 4: the remaining 25.
Steps (per prime p and frame (alpha, beta)):
  step1  setup checks; point frame M = 1: hatted sums of the 42 graphs (pool), frozen sums = [u^0], total and its t_c expansion;
         M-independence of the total at a random constant M (and [u^0] = frozen sums there); genus-2 graphs through the new machinery
         against the genus-2 values of record; s-series check for the graphs with E <= 1.
  step2  pass A: hatted sums at M_0 (constant): g_k[Gamma] = [u^k] G_Gamma[M_0]; [u^0] = 0 on edge graphs; O(u^E); regular set; vertex in t_c.
  step3  frozen Jacobian J[Gamma][v] = d/d eps|_0 G^0_Gamma[M_0 exp(eps N_v)], v = 0..74 (RIGHT placement), graphs with E <= 1 and two
         E = 2 controls; genus-2 J recomputed with the new machinery against fq_ou.jacobian_N at the same prime.
  step4  the [u^1] systems: genus 3 alone (equations for graphs 2 and 5) and genus 2 + 3 combined (4 equations) on the classes
         cov / diag / circ / full, right and left placement (left: columns conjugated by M_0); explicit cov solutions.
  step5  u-dependent passes M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3): tables [u^k], the next-order equations for N_{k+1} from the 1-edge
         graphs (same J), the tests on the >= 2-edge graphs, the vertex graph in t_c.
"""
import os, sys, json, time, random
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_hat import vander_inv, interp_fp, tc_series
from fq_mfit import fp_solve, fp_rref
from fq3_graphs import outer_graphs3
from fq3_m import NSP, SP_NAMES, sp_direction, sp_mats, sp_params, class_columns, expand, conj, CLASSES, MU
from fq3_msum import tc_expansion3
from fq3_par import run_pool, job_hat, job_frozen, job_jac, setup

GRAPHS3 = outer_graphs3()
GD = {G['idx']: G for G in GRAPHS3}
IDX_BY_E = {E: [G['idx'] for G in GRAPHS3 if G['E'] == E] for E in range(7)}
IDX_ALL = [G['idx'] for G in GRAPHS3]
IDX_VERTEX = 1
IDX_1EDGE = IDX_BY_E[1]                 # [2, 5]
assert IDX_1EDGE == [2, 5] and GD[2]['gv'] == (2,) and GD[5]['gv'] == (1, 2) and GD[1]['E'] == 0


def frame_of(p, ab, prec_abs, K): return (p, str(Fr(ab[0])), str(Fr(ab[1])), prec_abs, K)


class State(object):
    def __init__(self, path):
        self.path = path
        self.d = json.load(open(path)) if os.path.exists(path) else {}
    def save(self): json.dump(self.d, open(self.path, "w"), indent=1)
    def put(self, key, val): self.d[key] = val; self.save()
    def get(self, key, default=None): return self.d.get(key, default)


def hat_pass(frame, mspec, idxs, nproc, log, NN=21, label=""):
    """({idx: [u^0..u^{K-1}]}, {idx: info})"""
    res = run_pool(job_hat, [(frame, mspec, i, NN) for i in idxs], nproc, log, label)
    out = {}; info = {}
    for idx, coefs, inf in res: out[idx] = coefs; info[idx] = inf
    return out, info


def fmt_u(R, coefs, kmax=None, rec=None):
    terms = []
    for k, c in enumerate(coefs):
        if c is None or (kmax is not None and k > kmax): continue
        terms.append("(%s) u^%d" % (rec(c) if rec else R.fmt(c), k))
    return " + ".join(terms) if terms else "0"


def total_of(R, table, idxs=None):
    p = R.p; K = len(next(iter(table.values())))
    idxs = list(table) if idxs is None else idxs
    return [sum(table[i][k] for i in idxs) % p for k in range(K)]


def tc_fmt(rec, T, nmax=5):
    return " + ".join("(%s) t_c^%d" % (rec.fmt(T[k]), k - 4) for k in range(min(nmax, len(T))))


def log_table(log, R, rec, A, info, idxs, kmax=None):
    for i in idxs: log("      graph %2d (E=%d, gv=%s): %s" % (i, GD[i]['E'], GD[i]['gv'], fmt_u(R, A[i], kmax, rec.fmt)))


# ------------------------------------------------------------------ step 1
def step1(frame, st, state, nproc, log, check, rec):
    R = st['R']; p = R.p; K = st['K']; H = st['H']; sym = st['sym']
    check("HatRowsU: every Qhat has valuation 0 with non-zero constant term, xhat = 1 + O(u); u-precision %d >= K = %d" % (H.min_prec, K), H.val_ok and H.min_prec >= K)
    check("HatRowsU: Qhat|_{u=0} = the conifold symbol R^ for all (a, i, n), n <= 6", H.symbol_check(sym))
    check("M_0 = R~ mod zhat^7: symplectic to order 6, M_0 M_0^-1 = 1, exp(log M_0) = M_0, log in sp, symbol symplectic to order 6", all(st['M0checks'].values()), st['M0checks'])
    log("   1-edge graphs: " + "; ".join("idx %d gv=%s ed=%s" % (G['idx'], G['gv'], G['ed']) for G in GRAPHS3 if G['E'] == 1))
    tc = tc_series(st['C'], R, 10); state.put('tc', [tc.coef(k) for k in range(1, 8)])
    log("   t_c(u) = %s" % tc.fmt(4))
    t = time.time(); A, info = hat_pass(frame, ('id',), IDX_ALL, nproc, log, label="point frame")
    log("   point-frame pass: %.0f s (slowest graph %.0f s)" % (time.time() - t, max(i['time'] for i in info.values())))
    state.put('A_id', {str(i): A[i] for i in A}); state.put('time_id', {str(i): round(info[i]['time'], 1) for i in info})
    log_table(log, R, rec, A, info, IDX_ALL)
    tot = total_of(R, A); state.put('total_id', tot)
    log("   total F_3^cf u^4 = %s" % fmt_u(R, tot, rec=rec.fmt))
    T = tc_expansion3(tot, tc, R); state.put('total_tc', T[:6])
    log("   total in t_c: " + tc_fmt(rec, T, 6))
    check("point frame total: polar part in t_c is (1/25200) t_c^-4 only (gap: t_c^-3, t_c^-2, t_c^-1 coefficients vanish)", T[0] == R.of(Fr(1, 25200)) and T[1] == 0 and T[2] == 0 and T[3] == 0, tc_fmt(rec, T, 4))
    t = time.time(); fr = {}
    for d in run_pool(job_frozen, [(frame, ('id',), [i]) for i in IDX_ALL], nproc, log, label="frozen id"): fr.update(d)
    check("point frame: frozen sums (u = 0 theory, rows R^, x = 1) = [u^0] of the hatted sums, all 42 graphs (%.0f s)" % (time.time() - t), all(fr[i] == A[i][0] for i in IDX_ALL))
    rng = random.Random(2026); par = [rng.randrange(p) for _ in range(NSP)]; mu = sp_mats(par, p)
    t = time.time(); B, infoB = hat_pass(frame, ('const', mu), IDX_ALL, nproc, log, label="random M")
    totB = total_of(R, B); state.put('A_rand', {str(i): B[i] for i in B}); state.put('total_rand', totB)
    check("random constant symplectic M (order 6): total Sum_Gamma G_Gamma = point-frame total through u^%d (M-independence; %.0f s)" % (K - 1, time.time() - t), totB == tot)
    frB = {}
    for d in run_pool(job_frozen, [(frame, ('const', mu), [i]) for i in IDX_ALL], nproc, log, label="frozen rand"): frB.update(d)
    check("random M: frozen sums = [u^0] of the hatted sums, all 42 graphs", all(frB[i] == B[i][0] for i in IDX_ALL))
    return dict(A=A, tot=tot, tc=tc)


def step1_genus2(frame, st, state, log, check, rec):
    """the genus-2 graphs through the new machinery (point frame and M_0) against the genus-2 values of record"""
    from fq_msum import outer_graphs
    from fq3_mvertex import MVertex6, MConst
    from fq3_msum import GraphSummer6
    R = st['R']; p = R.p; H = st['H']; G2 = outer_graphs()
    ref_id0 = [Fr(1807, 3072000), Fr(-811, 256000), Fr(-57, 64000), Fr(843, 1024000), Fr(3, 51200), Fr(87, 32000), Fr(9, 12800)]
    ref_M0 = {1: [Fr(1, 1200), Fr(-7, 8000), Fr(67, 288000), Fr(1, 187500)], 2: [0, Fr(1, 750), Fr(-1, 480), Fr(-161, 150000)], 3: [0, 0, Fr(1, 1000), Fr(2, 625)],
              4: [0, Fr(-13, 8000), Fr(19, 96000), Fr(203, 3600000)], 5: [0, 0, Fr(1, 500), Fr(1, 6000)], 6: [0, 0, 0, Fr(-7, 7500)], 7: [0, 0, 0, Fr(-9, 5000)]}
    t = time.time(); I = MConst.identity(R); A = GraphSummer6(H, I, MVertex6(R, I)).all_graphs(G2)
    check("genus-2 graphs with the genus-3 code (GraphSummer6, order-6 M = 1): [u^0] = the 7 point-frame values of genus2_ou (%.0f s)" % (time.time() - t), all(A[i].coef(0) == R.of(ref_id0[i - 1]) for i in range(1, 8)))
    A0 = GraphSummer6(H, st['M0'], MVertex6(R, st['M0'])).all_graphs(G2)
    check("genus-2 graphs with M_0 (order 6; orders > 3 cannot enter at genus 2): [u^0..u^3] = the table of OU_EXTENSION_genus2 s.3", all(A0[i].coef(k) == R.of(ref_M0[i][k]) for i in range(1, 8) for k in range(4)))
    tot2 = None
    for i in A0: tot2 = A0[i] if tot2 is None else tot2.add(A0[i])
    check("genus-2 total through u^4 = 1/1200 - 7/6000 u + 97/72000 u^2 - 6809/18000000 u^3 - 2551/15000000 u^4",
          [tot2.coef(k) for k in range(5)] == [R.of(x) for x in (Fr(1, 1200), Fr(-7, 6000), Fr(97, 72000), Fr(-6809, 18000000), Fr(-2551, 15000000))])


def step1_sseries(frame, st, state, log, check, idxs=None, kmax=6):
    """s-series check: GraphSummer6 on the s-series rows (SRows) at M = 1 equals s^{-20} x hatted sums (u = s^5) for the listed graphs"""
    from fq3_mvertex import MVertex6, MConst
    from fq3_msum import GraphSummer6
    from gsum_q import u_to_s
    R = st['R']; S = st['S']; H = st['H']; idxs = idxs or (IDX_BY_E[0] + IDX_BY_E[1])
    I = MConst.identity(R); mv = MVertex6(R, I); GSs = GraphSummer6(S, I, mv); GSu = GraphSummer6(H, I, mv)
    t = time.time(); bad = []
    for i in idxs:
        G = GD[i]; ss = GSs.graph(G); su = u_to_s(GSu.graph(G)).shift(-20)
        N = min(ss.prec, su.prec, -20 + kmax + 1)
        if ss.v != -20 or not ss.truncate(N).eq(su.truncate(N)): bad.append((i, ss.v, ss.fmt(-16), su.fmt(-16)))
    check("s-series check: direct s-series sums (SRows, M = 1) have valuation -20 and equal s^-20 G_Gamma(s^5) through s^%d for graphs %s (%.0f s)" % (-20 + kmax, idxs, time.time() - t), not bad, bad[:2])


# ------------------------------------------------------------------ step 2
def step2(frame, st, state, nproc, log, check, rec, tot_id, tc):
    R = st['R']; p = R.p; K = st['K']
    t = time.time(); A, info = hat_pass(frame, ('M0',), IDX_ALL, nproc, log, label="pass A (M0)")
    log("   pass A (M_0): %.0f s" % (time.time() - t)); state.put('A_M0', {str(i): A[i] for i in A})
    log_table(log, R, rec, A, info, IDX_ALL)
    check("pass A: [u^0] G_Gamma[M_0] = 0 for every edge graph (M_0 kills the leading t_c^-4 pole of every edge graph)", all(A[i][0] == 0 for i in IDX_ALL if i != IDX_VERTEX))
    check("pass A: [u^k] G_Gamma[M_0] = 0 for k < E (an E-edge graph is O(u^E) since R' = 1 + O(u))", all(A[i][k] == 0 for i in IDX_ALL for k in range(GD[i]['E'])))
    reg = [i for i in IDX_ALL if i != IDX_VERTEX and all(A[i][k] == 0 for k in range(4))]
    nonreg = [i for i in IDX_ALL if i != IDX_VERTEX and i not in reg]
    log("   with M_0 alone: REGULAR edge graphs (O(u^4)): %s; NOT regular: %s" % (reg, nonreg)); state.put('M0_regular', reg); state.put('M0_nonregular', nonreg)
    check("pass A: every graph with >= 4 edges is regular with M_0 alone", all(i in reg for E in (4, 5, 6) for i in IDX_BY_E[E]))
    log("   non-regular graphs with M_0 (u^0..u^3): " + "; ".join("%d (E=%d): %s" % (i, GD[i]['E'], fmt_u(R, A[i], 3, rec.fmt)) for i in nonreg))
    check("pass A: total = point-frame total through u^%d (M-independence)" % (K - 1), total_of(R, A) == tot_id)
    check("pass A: vertex graph [u^0] = 1/25200 = (B_6/(6 . 4))/5 (the genus-3 gap coefficient)", A[IDX_VERTEX][0] == R.of(Fr(1, 25200)), rec.fmt(A[IDX_VERTEX][0]))
    V = tc_expansion3(A[IDX_VERTEX], tc, R); state.put('vertex_M0_tc', V[:5])
    log("   vertex graph with M_0 in t_c: %s" % tc_fmt(rec, V, 5))
    return dict(A=A, reg=reg, nonreg=nonreg)


# ------------------------------------------------------------------ step 3
def step3(frame, st, state, nproc, log, check, rec, nodes=7, control_idxs=(3, 6)):
    R = st['R']; p = R.p
    idxs = IDX_BY_E[0] + IDX_BY_E[1] + list(control_idxs)
    Vinv = vander_inv(list(range(nodes)), p)
    t = time.time()
    jobs = [(frame, sp_direction(v, p), idxs, list(range(nodes)), 'right', v) for v in range(NSP)]
    res = dict(run_pool(job_jac, jobs, nproc, log, label="J directions"))
    J = {i: [0] * NSP for i in idxs}; maxdeg = 0; consts = {i: set() for i in idxs}
    for v in range(NSP):
        for i in idxs:
            c = interp_fp(res[v][i], Vinv, p)
            J[i][v] = c[1]; consts[i].add(c[0]); maxdeg = max(maxdeg, max([k for k, x in enumerate(c) if x] or [0]))
    log("   J: %d directions x %d nodes x %d graphs in %.0f s; maximal eps-degree seen %d" % (NSP, nodes, len(idxs), time.time() - t, maxdeg))
    state.put('J', {str(i): J[i] for i in idxs}); state.put('J_maxdeg', maxdeg)
    check("J: eps-polynomials of degree <= 6 (= 3g-3; eps^6 is the last node-resolved power with 7 nodes: degree bound consistent)", maxdeg <= 6)
    A0 = state.get('A_M0')
    check("J: constant terms = [u^0] G_Gamma[M_0] (pass A) for every direction and graph", all(consts[i] == {A0[str(i)][0] % p} for i in idxs))
    check("J: D_N L_Gamma = 0 identically for the E = 2 control graphs %s" % (list(control_idxs),), all(not any(J[i]) for i in control_idxs))
    check("J: Sum over graphs (vertex + the two 1-edge graphs) = 0 (the frozen total is M-independent)", all((J[1][v] + J[2][v] + J[5][v]) % p == 0 for v in range(NSP)))
    rank = len(fp_rref([J[2], J[5]], NSP, p)[1])
    check("J: rank of the 2 x 75 matrix (rows graphs 2, 5) = 2", rank == 2); state.put('J_rank', rank)
    for i in (2, 5): log("      J[%d]: %d non-zero entries; first-block (k = 1) entries: %s" % (i, sum(1 for x in J[i] if x), [rec.fmt(x) for x in J[i][:15]]))
    return dict(J=J, Vinv=Vinv)


def step3_genus2(frame, st, state, log, check, rec, nodes=7):
    """the genus-2 J (40 coordinates, graphs 2 and 4 of genus 2) with the NEW machinery against fq_ou.jacobian_N (old code), and extended to 75"""
    from fq_msum import outer_graphs as outer2
    from fq_leadsolve import FrozenRows
    from fq_sser import SympM
    from fq_hat import tautological_M0
    from fq_ou import jacobian_N
    from fq3_msum import frozen_sums6
    from fq3_par import make_M
    R = st['R']; p = R.p; G2 = outer2(); fr = st['frozen']
    t = time.time(); Vinv = vander_inv(list(range(nodes)), p); J2 = {i: [0] * NSP for i in (1, 2, 4)}
    for v in range(NSP):
        N = sp_direction(v, p); vals = {i: [] for i in (1, 2, 4)}
        for eps in range(nodes):
            M, _ = make_M(st, ('M0exp_eps', N, eps)); res = frozen_sums6(fr, M, [G2[i - 1] for i in (1, 2, 4)])
            for i in (1, 2, 4): vals[i].append(res[i])
        for i in (1, 2, 4): J2[i][v] = interp_fp(vals[i], Vinv, p)[1]
    log("   genus-2 J with the new machinery: %.0f s" % (time.time() - t))
    M0s, mu, C, ok = tautological_M0(st['sym'], R); Vinv10 = vander_inv(list(range(10)), p)
    t = time.time(); Jold, maxdeg, consts = jacobian_N(FrozenRows(st['sym']), M0s.M, M0s.Minv, G2, R, Vinv10)
    log("   genus-2 J with the old code (fq_ou.jacobian_N): %.0f s" % (time.time() - t))
    check("genus-2 J: new machinery (75 coordinates) = old fq_ou.jacobian_N on the 40 genus-2 coordinates, graphs 1, 2, 4", all(J2[i][v] == Jold[i][v] % p for i in (1, 2, 4) for v in range(40)))
    check("genus-2 J: the 35 coordinates of zhat-orders 4, 5, 6 have zero entries (orders > 3 cannot enter at genus 2)", all(J2[i][v] == 0 for i in (1, 2, 4) for v in range(40, NSP)))
    state.put('J2', {str(i): J2[i] for i in (2, 4)})
    return J2


# ------------------------------------------------------------------ step 4
def solve_system(R, rows, rhs, cols):
    """restricted system (rows . col_j) t = rhs on the class spanned by cols: (matrix, particular t, nullspace, rank, consistent)"""
    p = R.p
    Jr = [[sum(row[v] * col[v] for v in range(NSP)) % p for col in cols] for row in rows]
    sol = fp_solve(Jr, [r % p for r in rhs], p); x, ns = sol[0], sol[1]
    rank = len(fp_rref(Jr, len(cols), p)[1])
    return Jr, x, ns, rank, x is not None


def step4(frame, st, state, log, check, rec, J3, g1_3, J2=None, g1_2=None):
    """[u^1] systems; returns {(cls, placement, which): dict}"""
    R = st['R']; p = R.p; M0 = st['M0']; out = {}
    systems = {'g3': ([J3[2], J3[5]], [(-g1_3[2]) % p, (-g1_3[5]) % p], ["graph 2 (loop on g2)", "graph 5 (bridge g1-g2)"])}
    if J2 is not None:
        systems['g2+g3'] = ([J2[2], J2[4], J3[2], J3[5]], [(-g1_2[2]) % p, (-g1_2[4]) % p, (-g1_3[2]) % p, (-g1_3[5]) % p],
                            ["genus-2 graph 2 (g1 loop)", "genus-2 graph 4 (g1-g1)", "graph 2 (loop on g2)", "graph 5 (bridge g1-g2)"])
    for which, (rows, rhs, names) in systems.items():
        log("   --- system %s: %d equations %s" % (which, len(rows), names))
        for cls in ('cov', 'diag', 'circ', 'full'):
            for placement in ('right', 'left'):
                cols, cnames = class_columns(R, cls, M0, placement)
                Jr, x, ns, rank, cons = solve_system(R, rows, rhs, cols)
                dim = (len(cols) - rank) if cons else None
                key = "%s:%s:%s" % (which, cls, placement)
                out[(which, cls, placement)] = dict(Jr=Jr, x=x, ns=ns, rank=rank, consistent=cons, dim=dim, cols=cols, names=cnames)
                state.put('sys4_' + key, dict(Jr=Jr, rhs=rhs, x=x, rank=rank, consistent=cons, dim=dim, names=cnames))
                msg = "rank %d, %s" % (rank, ("consistent, solution set of dimension %d" % dim) if cons else "INCONSISTENT")
                if cls == 'cov':
                    log("      %-22s matrix rows (b1, b3, b5): %s; rhs %s" % (key, [[rec.fmt(c) for c in row] for row in Jr], [rec.fmt(r) for r in rhs]))
                    if cons:
                        log("      %-22s solution (b1, b3, b5) = %s%s" % (key, [rec.fmt(c) for c in x], ("; free directions %s" % [[rec.fmt(c) for c in nvec] for nvec in ns]) if ns else ""))
                check("%s: %s" % (key, msg), True)
                log("      %-22s %s" % (key, msg))
    return out


# ------------------------------------------------------------------ step 5
def pass_u(frame, st, state, nproc, log, check, rec, Ns, placement, tag, tot_id, NN=21):
    """one u-dependent pass M = M_0 exp(Sum u^k N_k); returns the table"""
    R = st['R']; p = R.p; K = st['K']
    t = time.time(); A, info = hat_pass(frame, ('M0exp', Ns, placement), IDX_ALL, nproc, log, NN, label="pass %s" % tag)
    maxdeg = max(i['maxdeg'] or 0 for i in info.values())
    log("   pass %s (%s, %d N's): %.0f s; maximal eps-degree of the vertex tensors %d (bound 18; NN = %d nodes)" % (tag, placement, len(Ns), time.time() - t, maxdeg, NN))
    state.put('A_' + tag, {str(i): A[i] for i in A}); state.put('Ns_' + tag, dict(placement=placement, Ns=[sp_params(N, p) for N in Ns]))
    log_table(log, R, rec, A, info, IDX_ALL)
    check("pass %s: eps-degree <= 18 for every vertex tensor (eps^19, eps^20 coefficients vanish)" % tag, maxdeg <= 18)
    check("pass %s: [u^k] = 0 for k < E for every graph" % tag, all(A[i][k] == 0 for i in IDX_ALL for k in range(GD[i]['E'])))
    check("pass %s: total = point-frame total through u^%d (M-independence for u-dependent M)" % (tag, K - 1), total_of(R, A) == tot_id)
    return A


def status(R, A, kmax=3):
    """per graph: list of the u-orders k <= kmax with non-zero coefficient (edge graphs only)"""
    return {i: [k for k in range(kmax + 1) if A[i][k]] for i in IDX_ALL if i != IDX_VERTEX}


def next_N(R, J3, A, order, cls, placement, M0, log, check, rec, tag):
    """solve J . n = -[u^order] A for the 1-edge graphs on the class; returns (N as sp series or None, info)"""
    p = R.p; cols, names = class_columns(R, cls, M0, placement)
    rows = [J3[2], J3[5]]; rhs = [(-A[2][order]) % p, (-A[5][order]) % p]
    Jr, x, ns, rank, cons = solve_system(R, rows, rhs, cols)
    log("   %s: equations for N_%d on class %s/%s: rhs (graphs 2, 5) = %s; rank %d; %s" % (tag, order, cls, placement, [rec.fmt(r) for r in rhs], rank, "consistent (dim %d)" % (len(cols) - rank) if cons else "INCONSISTENT"))
    if not cons: return None, dict(Jr=Jr, rhs=rhs, rank=rank, consistent=False)
    # canonical choice: the particular solution of fp_solve (free parameters = 0); report in class coordinates
    log("      N_%d class coordinates %s = %s" % (order, names, [rec.fmt(c) for c in x]))
    N = sp_mats(expand(cols, x, p), p)
    return N, dict(Jr=Jr, rhs=rhs, rank=rank, consistent=True, t=x, ns=ns, names=names)
