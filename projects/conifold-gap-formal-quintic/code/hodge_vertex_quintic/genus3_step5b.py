# -*- coding: utf-8 -*-
"""genus3_step5b.py -- the TWO-SIDED covariant-scalar dressing  M = exp(u S_L) M_0 exp(u S_R),  S_L = Sum_{k odd} b'_k diag(zeta^{-ki}) zhat^k,
S_R = Sum_{k odd} b_k diag(zeta^{-ki}) zhat^k  (6 parameters; pointed out by the results referee): at O(u) its N_1 is M_0^{-1} S_L M_0 + S_R,
and the four O(u) equations (genus 2: graphs 2, 4; genus 3: graphs 2, 5) have rank 4 on the 6 parameters, leaving a 2-dimensional affine
family (t_1, t_2).  Along it the [u^2] coefficient of every two-edge graph is a quadratic form Q_Gamma(t_1, t_2) (two O(u) edge factors, each
affine in (t_1, t_2)); it is interpolated exactly from the 6 passes at (t_1, t_2) = (0,0), (1,0), (2,0), (0,1), (0,2), (1,1) and checked at
(2,1).  Common zeros of the five Q_Gamma over C are decided by resultants: Res_{t_2}(Q_3, Q_j) are polynomials of degree <= 4 in t_1 whose
gcd over F_p (two primes) must have a root for a common zero to exist.  Graphs: E <= 2 only (1, 2, 5, 3, 6, 7, 8, 17; the E >= 3 graphs are
O(u^3) for every such M).  Run:  py -I genus3_step5b.py --prime 89|61 [--nproc 24]  -> appends to logs/genus3_ou_<p>.log, state g3_state_<p>.json
"""
import os, sys, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import P1, P2
from fq_mfit import fp_solve, fp_rref
from fq3_steps import frame_of, State, hat_pass, IDX_BY_E, GD, total_of, fmt_u, log_table
from fq3_m import sp_mats, sp_params, expand, class_columns, conj, NSP, sp_add, sp_scale
from fq3_par import setup
from fq3_rec import Rec

PTS = [(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1)]; CHK = (2, 1)
MONOS = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]


def mono_eval(t, m, p): return pow(t[0], m[0], p) * pow(t[1], m[1], p) % p


def fit_quadric(vals, p):
    """coefficients c_m (m in MONOS) with Sum_m c_m t^m = vals[t] at the 6 points PTS"""
    A = [[mono_eval(t, m, p) for m in MONOS] + [v] for t, v in zip(PTS, vals)]
    n = 6
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c]); A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p); A[c] = [x * inv % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]: f = A[r][c]; A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [A[m][n] for m in range(n)]


def quadric_eval(c, t, p): return sum(ci * mono_eval(t, m, p) for ci, m in zip(c, MONOS)) % p


# ---- polynomials in t_1 over F_p as coefficient lists (low to high)
def padd(a, b, p):
    n = max(len(a), len(b)); return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)]
def psub(a, b, p): return padd(a, [(-x) % p for x in b], p)
def pmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] = (out[i + j] + x * y) % p
    return out
def ptrim(a):
    while len(a) > 1 and a[-1] == 0: a = a[:-1]
    return a
def pgcd(a, b, p):
    a = ptrim(a); b = ptrim(b)
    while b != [0]:
        # a mod b
        a = list(a)
        while len(a) >= len(b) and a != [0]:
            f = a[-1] * pow(b[-1], p - 2, p) % p; sh = len(a) - len(b)
            for i, y in enumerate(b): a[sh + i] = (a[sh + i] - f * y) % p
            a = ptrim(a)
            if len(a) < len(b): break
        a, b = b, a
    return ptrim(a)


def as_t2_quadratic(c, p):
    """Q(t_1, t_2) = A(t_1) t_2^2 + B(t_1) t_2 + C(t_1), A, B, C polynomials in t_1"""
    d = dict(zip(MONOS, c))
    return [d[(0, 2)]], [d[(0, 1)], d[(1, 1)]], [d[(0, 0)], d[(1, 0)], d[(2, 0)]]


def resultant_t2(c1, c2, p):
    """Res_{t_2} of two quadratics in t_2: (a2 c0 - a0 c2)^2 - (a2 b1 - a1 b2)(a1 c0 - a0 c1) with (a, b, c) = (A, B, C) of each"""
    A1, B1, C1 = as_t2_quadratic(c1, p); A2, B2, C2 = as_t2_quadratic(c2, p)
    X = psub(pmul(A1, C2, p), pmul(C1, A2, p), p); Y = psub(pmul(A1, B2, p), pmul(B1, A2, p), p); Z = psub(pmul(B1, C2, p), pmul(C1, B2, p), p)
    return ptrim(psub(pmul(X, X, p), pmul(Y, Z, p), p))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--nproc", type=int, default=24); ap.add_argument("--NN", type=int, default=21)
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1
    frame = frame_of(p, (0, 0), 34, 6); state = State(os.path.join(HERE, "g3_state_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_ou_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    st = setup(frame); R = st['R']; rec = Rec(R); M0 = st['M0']; T0 = time.time()
    tot_id = state.get('total_id'); J3 = {int(i): v for i, v in state.get('J').items()}; J2 = {int(i): v for i, v in state.get('J2').items()}
    A0 = {int(i): v for i, v in state.get('A_M0').items()}
    log("-" * 110); log("(5b) two-sided covariant-scalar dressing exp(u S_L) M_0 exp(u S_R), genus 2 + 3 O(u) system  (%s)" % time.strftime("%Y-%m-%d %H:%M"))
    colsR, namesR = class_columns(R, 'cov', M0, 'right'); colsL, namesL = class_columns(R, 'cov', M0, 'left')   # left columns = conj by M_0
    cols = colsR + colsL; names = ["b1", "b3", "b5", "b1'", "b3'", "b5'"]
    rows = [J2[2], J2[4], J3[2], J3[5]]; rhs = [(-R.of(Fr(1, 750))) % p, (-R.of(Fr(-13, 8000))) % p, (-A0[2][1]) % p, (-A0[5][1]) % p]
    Jr = [[sum(row[v] * col[v] for v in range(NSP)) % p for col in cols] for row in rows]
    x, ns, rank = fp_solve(Jr, rhs, p)
    log("   4 x 6 matrix (rows: genus-2 graphs 2, 4; genus-3 graphs 2, 5; columns %s):" % names)
    for row in Jr: log("      " + str([rec.fmt(c) for c in row]))
    check("two-sided scalar class: rank 4, consistent, 2-dimensional solution set", x is not None and rank == 4 and len(ns) == 2, (rank, x is not None, len(ns)))
    state.put('sys5b_Jr', Jr); state.put('sys5b_rank', rank); state.put('sys5b_x', x); state.put('sys5b_ns', ns)
    log("   particular solution %s; directions %s, %s" % ([rec.fmt(c) for c in x], [rec.fmt(c) for c in ns[0]], [rec.fmt(c) for c in ns[1]]))
    def SLR(t):
        par = [(x[j] + t[0] * ns[0][j] + t[1] * ns[1][j]) % p for j in range(6)]
        SR = sp_mats(expand(colsR, par[:3], p), p); SL = sp_mats(expand(colsR, par[3:], p), p)      # S_L in its own (scalar) coordinates
        return SL, SR
    idxs = IDX_BY_E[0] + IDX_BY_E[1] + IDX_BY_E[2]
    vals = {i: [] for i in idxs}; tables = {}
    for t in PTS + [CHK]:
        SL, SR = SLR(t)
        t0 = time.time(); A, info = hat_pass(frame, ('M0exp2', [SL], [SR]), idxs, a.nproc, log, a.NN, label="t=%s" % (t,))
        log("   pass (t_1, t_2) = %s: %.0f s; max eps-degree %d" % (t, time.time() - t0, max(i['maxdeg'] or 0 for i in info.values())))
        tables[t] = A; state.put('A_cov5b_t%d%d' % t, {str(i): A[i] for i in A})
        check("(%d,%d): 1-edge graphs vanish through u^1, every graph O(u^E)" % t, all(A[i][k] == 0 for i in (2, 5) for k in range(2)) and all(A[i][k] == 0 for i in idxs for k in range(GD[i]['E'])))
        if t in PTS:
            for i in idxs: vals[i].append(A[i][2])
    log("   quadratic forms Q_Gamma(t_1, t_2) = [u^2] of the two-edge graphs (coefficients of %s):" % (MONOS,))
    Q = {}
    for i in IDX_BY_E[2]:
        c = fit_quadric(vals[i], p); Q[i] = c; state.put('obs5b_u2_%d' % i, c)
        check("graph %d: the fitted quadric reproduces the check point (2,1) (degree <= 2 confirmed)" % i, quadric_eval(c, CHK, p) == tables[CHK][i][2])
        log("      graph %2d: %s" % (i, [rec.fmt(z) for z in c]))
    # common zeros: resultants in t_2 against graph 3, gcd in t_1
    g = None
    for i in IDX_BY_E[2]:
        if i == 3: continue
        res = resultant_t2(Q[3], Q[i], p)
        log("      Res_t2(Q_3, Q_%d): degree %d in t_1" % (i, len(res) - 1))
        g = res if g is None else pgcd(g, res, p)
    log("   gcd of the four resultants Res_t2(Q_3, Q_j), j = 6, 7, 8, 17: degree %d (%s)" % (len(g) - 1, [rec.fmt(z) for z in g]))
    state.put('obs5b_gcd_degree', len(g) - 1)
    check("KEY (two-sided scalar class): the five [u^2] quadrics have NO common zero over the algebraic closure (gcd of resultants = constant)", len(g) == 1 and g[0] != 0)
    # the same with the roles of graphs exchanged (robustness: leading coefficients)
    g2 = None
    for i in (3, 7, 8, 17):
        res = resultant_t2(Q[6], Q[i], p); g2 = res if g2 is None else pgcd(g2, res, p)
    log("   gcd of Res_t2(Q_6, Q_j), j = 3, 7, 8, 17: degree %d" % (len(g2) - 1))
    log("step 5b time %.0f s" % (time.time() - T0)); log("STEP 5b ALL CHECKS PASSED" if not FAILED else "STEP 5b FAILURES: %s" % FAILED); LOG.close()


if __name__ == "__main__":
    main()
