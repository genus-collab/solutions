# -*- coding: utf-8 -*-
"""genus3_natural.py -- the NATURAL-VARIABLE form of the conifold factorisation of the formal quintic (genus 2 and 3, O(u^3)).

Statement tested (finite orders, modulo the prime):  with tau(u) defined by
        1/tau(u) = 2F1(1, 4/5; 9/5; u) = Sum_n 4 u^n/(5n + 4) = (4/5) u^{-4/5} Int_0^u t^{-1/5} (1 - t)^{-1} dt ,
i.e. tau(u) zhat = (5/4) z / Int_{q_c}^{q} L dq/q = z over the canonical coordinate measured from the conifold point (up to the
constant zeta^i 5/4), the Teleman form of the hatted R-matrix factorises as
        R~(zhat; u) = G(u, zhat) . M_0(tau(u) zhat),      G = 1 + Sum_{n >= 1} u^n G_n(zhat),   G_n a POLYNOMIAL in zhat,
with M_0 = R~(zhat; 0) the conifold symbol (u-independent).  Checked here through u^4 on the zhat-orders <= 6 (the degrees of G_1..G_4 are
1, 3, 4, 5, so the orders above them are genuine tests: 5 + 3 + 2 + 1 matrix identities) and G_1 = zhat X_1 with X_1 a combination of the
two propagators 1 (x) H^2 + H^2 (x) 1 and H^3 (x) H^4 + H^4 (x) H^3 only (no H (x) H): a point of the kinematic plane Pi of
OU_EXTENSION_genus3.md s.6d.

M = A(u, zhat) M_0(tau(u) zhat) with A = exp(Sum_n u^n a_n(zhat)) a left dressing (--variant nat: A = 1; other variants read the a_n from
g3_nat_dress_<variant>_<p>.json, key 'a': list over n = 1.. of 7 matrices).  M = M_0 exp(Sum u^k N_k) defines N_k (right placement); the
pass over all 42 genus-3 graphs uses (N_1, N_2, N_3) (M mod u^4) and reports [u^0..u^3] of every graph, the vertex graph in t_c, and the
genus-2 graphs with the same M.
Run: py -I genus3_natural.py --prime 89|61 [--variant nat] [--nproc 12] [--nopass] [--ab 1/3,-2/7]
(--ab: conifold frame (alpha, beta) of OU_EXTENSION_genus3.md s.7; default (0, 0); the stored N_1-conditions and the value of the string
component of G_1 are frame-(0, 0) data and are only checked there).
"""
import os, sys, json, argparse, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq_sser import mat_id
from fq3_steps import State, frame_of, hat_pass, IDX_ALL, GD, total_of
from fq3_par import setup
from fq3_m import NSP, sp_params, in_sp, MU, _biv_mul, _biv_exp, Z5, KMAX
from fq3_taut import colour_series_u
from fq3_rec import Rec
NU = 6


def tau_pred(n):
    """coefficients of tau(u) = 1 / Sum_k 4 u^k/(5k+4) as Fractions"""
    inv = [Fr(4, 5 * k + 4) for k in range(n)]; tau = [Fr(1)]
    for k in range(1, n): tau.append(-sum(inv[j] * tau[k - j] for j in range(1, k + 1)))
    return tau


def biv_clean(X, p): return {key: A for key, A in X.items() if any(x % p for r in A for x in r)}


def biv_log(X, p, K, R):
    """log of a bivariate series X = 1 + Y, every term of Y of zhat-order >= 1 and u-order >= 1"""
    Y = biv_clean({key: A for key, A in X.items() if key != (0, 0)}, p)
    assert X[(0, 0)] == mat_id(5, p) and all(k >= 1 and j >= 1 for (k, j) in Y), sorted(Y)[:5]
    out = {}; power = {(0, 0): mat_id(5, p)}
    for m in range(1, K):
        power = _biv_mul(power, Y, p, K); c = R.inv(m) * (1 if m % 2 == 1 else p - 1) % p
        for key, A in power.items():
            B = out.get(key, Z5()); out[key] = [[(B[i][j] + c * A[i][j]) % p for j in range(5)] for i in range(5)]
    return out


def fourier(X, z, p):
    inv25 = pow(25, p - 2, p); out = {}
    for m in range(5):
        for mp in range(5):
            v = sum(pow(z, m * i, p) * X[i][j] * pow(z, mp * j, p) for i in range(5) for j in range(5)) % p
            if v: out[(m, mp)] = v * inv25 % p
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--variant", default="nat")
    ap.add_argument("--nproc", type=int, default=12); ap.add_argument("--nopass", action="store_true"); ap.add_argument("--ab", default="0,0")
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R); var = a.variant
    ab = tuple(Fr(x) for x in a.ab.split(",")); f00 = ab == (0, 0)
    ftag = "" if f00 else "_ab" + "_".join(str(x).replace("/", "o").replace("-", "m") for x in ab)
    us = State(os.path.join(HERE, "g3_nat_%s_%s%s.json" % (a.prime, var, ftag)))
    LOG = open(os.path.join(HERE, "logs", "genus3_natural_%s_%s%s.log" % (a.prime, var, ftag)), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("=" * 110); log("genus3_natural: variant %s, p = 2^%s-1, frame (alpha, beta) = %s" % (var, a.prime, a.ab))
    frame6 = frame_of(p, ab, 34, NU); st = setup(frame6); H = st['H']; M0 = st['M0']; z = R.zeta
    # ------------------------------------------------------------ tau, M_nat, G
    tau = [R.of(x) for x in tau_pred(NU)]
    tp = [[1] + [0] * (NU - 1)]
    for k in range(1, KMAX + 1): tp.append([sum(tp[-1][i] * tau[n - i] for i in range(n + 1)) % p for n in range(NU)])
    sc = lambda A, c: [[x * c % p for x in r] for r in A]
    Mnat = biv_clean({(k, j): sc(M0.M[k], tp[k][j]) for k in range(KMAX + 1) for j in range(NU)}, p)
    Mnat_inv = biv_clean({(k, j): sc(M0.Minv[k], tp[k][j]) for k in range(KMAX + 1) for j in range(NU)}, p)
    Pchk = _biv_mul(Mnat, Mnat_inv, p, NU)
    check("M_0(tau zhat) M_0(tau zhat)^-1 = 1 mod (zhat^7, u^6)", all((A == mat_id(5, p)) if key == (0, 0) else not any(x for r in A for x in r) for key, A in Pchk.items()))
    C = [colour_series_u(H, n) for n in range(NU)]
    Rt = biv_clean({(k, n): C[n][k] for n in range(NU) for k in range(KMAX + 1)}, p)
    G = biv_clean(_biv_mul(Rt, Mnat_inv, p, NU), p)
    deg = {n: max([k for (k, j) in G if j == n] or [0]) for n in range(NU)}
    log("   G = R~(zhat; u) M_0(tau(u) zhat)^-1: zhat-degree of [u^n] G for n = 0..5: %s  (truncation: zhat^6)" % [deg[n] for n in range(NU)])
    check("FACTORISATION: [u^0] G = 1 and deg_zhat [u^n] G = 1, 3, 4, 5 for n = 1, 2, 3, 4 with tau = 1/2F1(1, 4/5; 9/5; u) (orders up to 6 tested)",
          deg[0] == 0 and [deg[n] for n in (1, 2, 3, 4)] == [1, 3, 4, 5])
    for n in (1, 2, 3):
        for k in range(1, deg[n] + 1):
            if (k, n) in G: log("      [u^%d zhat^%d] G, Fourier components Xhat(m, m')/25: %s" % (n, k, {key: rec.fmt(v) for key, v in fourier(G[(k, n)], z, p).items()}))
    F1 = fourier(G[(1, 1)], z, p)
    check("G_1 = zhat X_1, X_1 = (%s) [1 (x) H^2 + H^2 (x) 1 component] + (-1/75) [H^3 (x) H^4 + H^4 (x) H^3 component], NO H (x) H component" % rec.fmt(F1.get((2, 4), 0)),
          set(F1) == {(0, 1), (1, 0), (2, 4), (4, 2)} and (not f00 or F1[(2, 4)] == R.of(Fr(-8, 225))) and F1[(0, 1)] == R.of(Fr(-1, 75)))
    us.put('G', {"%d,%d" % key: A for key, A in G.items() if key[1] <= 4}); us.put('tau', tau)
    # ------------------------------------------------------------ the dressing and the N_k
    M = Mnat
    if var != 'nat':
        dr = json.load(open(os.path.join(HERE, "g3_nat_dress_%s_%s.json" % (var, a.prime))))
        Xd = biv_clean({(k, n + 1): dr['a'][n][k] for n in range(len(dr['a'])) for k in range(1, KMAX + 1)}, p)
        M = _biv_mul(_biv_exp(Xd, p, NU, R), Mnat, p, NU)
        log("   left dressing A = exp(Sum u^n a_n): %s" % dr.get('desc', ''))
    X = _biv_mul({(k, 0): M0.Minv[k] for k in range(KMAX + 1)}, M, p, NU)
    Lg = biv_log(X, p, NU, R)
    Ns = [[Lg.get((k, j), Z5()) for k in range(KMAX + 1)] for j in range(1, NU)]
    check("the N_k = [u^k] log(M_0^-1 M) are in sp (k = 1..5)", all(in_sp(N, p) for N in Ns))
    MUo = MU(R, M0, Ns, NU, 'right')
    check("M_0 exp(Sum_k u^k N_k) = M exactly mod (zhat^7, u^6)", all(MUo.M[k][i][j][n] == M.get((k, n), Z5())[i][j] for k in range(KMAX + 1) for n in range(NU) for i in range(5) for j in range(5)))
    n123 = [sp_params(N, p) for N in Ns[:3]]; us.put('n1', n123[0]); us.put('n2', n123[1]); us.put('n3', n123[2])
    # N_1-conditions
    from fq3_circsys import CircSys
    S = CircSys(a.prime)        # frame-(0, 0) data: N_1^taut and the fitted quadrics/cubics; the O(u) rows are frame independent
    d1 = [(n123[0][v] - S.nt[v]) % p for v in range(NSP)]
    c15 = fp_solve([[S.cols[j][v] for j in range(15)] for v in range(NSP)], d1, p)[0]
    check("N_1 - N_1^taut is circulant", c15 is not None)
    res = S.residuals(c15); bad = [k for k, v in res.items() if v]
    if f00: check("N_1 satisfies EVERY N_1-condition of genus 2 and 3 (4 linear rows, 5 two-edge quadrics, 9 three-edge cubics)", not bad, bad)
    log("   circulant coordinates of N_1 - N_1^taut: %s" % ["%s: %s" % (S.names[j], rec.fmt(c15[j])) for j in range(15) if c15[j]])
    if a.nopass:
        us.put('failed', FAILED); log("GENUS3_NATURAL (%s, no pass) %s" % (var, "ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); LOG.close(); return
    # ------------------------------------------------------------ the pass over all 42 graphs
    A3 = us.get('A3')
    if A3 is None:
        t = time.time(); A3, info = hat_pass(frame6, ('M0exp', Ns[:3], 'right'), IDX_ALL, a.nproc, log, 21, label="pass (N1, N2, N3)")
        log("   pass (42 graphs): %.0f s; maximal eps-degree %s" % (time.time() - t, max(i['maxdeg'] or 0 for i in info.values())))
        us.put('A3', {str(i): A3[i] for i in A3})
    else: A3 = {int(k): v for k, v in A3.items()}
    edge = [i for i in IDX_ALL if i != 1]
    for i in IDX_ALL:
        if i == 1 or any(A3[i][k] for k in range(4)): log("      graph %2d (E=%d, gv=%s): %s" % (i, GD[i]['E'], GD[i]['gv'], [rec.fmt(x) for x in A3[i][:5]]))
    nonreg = {i: [k for k in range(4) if A3[i][k]] for i in edge if any(A3[i][k] for k in range(4))}
    us.put('nonregular', {str(i): v for i, v in nonreg.items()})
    check("FINAL: all 41 edge graphs of genus 3 are O(u^4), i.e. REGULAR at the conifold point", not nonreg, nonreg)
    tot = total_of(R, A3); g3 = json.load(open(os.path.join(HERE, "g3_state_%s%s.json" % (a.prime, "" if f00 else "_ab1o3_m2o7"))))
    check("FINAL: total over the 42 graphs = point-frame total through u^3 (M-independence)", tot[:4] == [x % p for x in g3['total_id']][:4])
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3, GraphSummer6, UVertex
    tc = tc_series(st['C'], R, 10); V = tc_expansion3(A3[1], tc, R)
    log("   vertex graph in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(V[k]), k - 4) for k in range(4)))
    check("FINAL: the vertex graph alone carries the polar part: (1/25200) t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1", V[0] == R.of(Fr(1, 25200)) and V[1] == 0 and V[2] == 0 and V[3] == 0)
    us.put('vertex_tc', V[:4])
    from fq_msum import outer_graphs
    MU3 = MU(R, M0, Ns[:3], NU, 'right'); t = time.time()
    G2 = GraphSummer6(H, MU3, UVertex(R, MU3, NU, 21)).all_graphs(outer_graphs())
    log("   genus-2 graphs with the same M (%.0f s): %s" % (time.time() - t, {i: [rec.fmt(G2[i].coef(k)) for k in range(3)] for i in sorted(G2)}))
    check("GENUS 2: all six edge graphs are O(u^2) (regular) with the same M; vertex graph [u^0] = 1/1200", all(G2[i].coef(k) == 0 for i in range(2, 8) for k in range(2)) and G2[1].coef(0) == R.of(Fr(1, 1200)))
    us.put('failed', FAILED)
    log("GENUS3_NATURAL (%s) ALL CHECKS PASSED" % var if not FAILED else "GENUS3_NATURAL (%s) FAILURES: %s" % (var, FAILED)); LOG.close()


if __name__ == "__main__":
    main()
