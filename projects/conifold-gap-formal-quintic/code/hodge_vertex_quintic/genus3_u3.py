# -*- coding: utf-8 -*-
"""genus3_u3.py -- explicit factorisations M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3) regularising EVERY edge graph of genus 2 and 3 (O(u^3)).

Conditions (right placement, frame (0, 0)); a genus-3 graph term is u^-4 G_Gamma(u), regular iff [u^0..u^3] G_Gamma = 0:
  N_1:  O(u) rows of genus 2 (2) and genus 3 (2); the five two-edge quadrics [u^2]; the nine three-edge cubics [u^3]  (fq3_quad, fq3_cubic);
  N_2:  one-edge graphs at u^2:  J . n_2 = -[u^2] G_Gamma(N_1, 0)   (Gamma = 2, 5; J = frozen Jacobian);
        two-edge graphs at u^3:  grad Q_Gamma(N_1) . n_2 = -[u^3] G_Gamma(N_1, 0)   (the second edge factor at order u^2 is linear in N_2
        through the same map as E'_1 in N_1, so the derivative is the polar form of the quadric);
  N_3:  one-edge graphs at u^3:  J . n_3 = -[u^3] G_Gamma(N_1, N_2, 0).
Classes:
  --cls full : no covariance; Delta_L = N_1 - N_1^taut non-zero in the orders 1, 2, 3 (fq3_fullsol.full_n1: orders 1, 2 on the linear locus of
               the star-graph cubic, order 3 and then the orders 4..6 from linear equations);
  --cls circpi: Z/5-covariant; the point A + 2 B of the plane Pi of solutions (genus3_pi.py);
  --cls circb: Z/5-covariant; the isolated solution with non-zero low-order part on {c_1(0) = 0, m_1 = 0} (genus3_cubicb);
  --cls circ : Z/5-covariant; Delta pure-high (orders 4..6) in the kernel of the genus-3 rows;
  --cls circstar: Z/5-covariant; the isolated solution Delta* with c_1(0) != 0 (key 'circ_star').
The full classification of the circulant N_1-solutions is in OU_EXTENSION_genus3.md s.6d.  --cls full is an abandoned attempt (fq3_fullsol.py,
inconsistent at its last step); no full-class pass exists.
Three successive passes per case (hatted theory, 21-node vertex interpolation): (N_1) and (N_1, N_2) on the 17 graphs with at most three
edges, used to solve N_2 and then N_3; (N_1, N_2, N_3) on all 42 graphs.
Final checks: all 41 edge graphs O(u^4); total = point-frame total (hence the vertex graph carries the polar part (1/25200) t_c^-4);
genus-2 edge graphs O(u^2) with the same M (a condition on N_1 only).
Run: py -I genus3_u3.py --prime 89|61 --cls circ|circpi|circb|circstar [--nproc 12]   (resumable: n1, n2, n3 and the pass tables are stored
in g3_u3_<p>_<cls>.json).  The tables of record for circ, circb, circpi were produced before the last docstring edits of this file.
"""
import os, sys, json, argparse, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq_mfit import fp_solve
from fq3_steps import State, frame_of, hat_pass, IDX_BY_E, IDX_ALL, GD, solve_system, total_of
from fq3_quad import QuadForm, rank_fp, IDX_2EDGE
from fq3_cubic import job_cubic, IDX_3EDGE
from fq3_par import run_pool, setup
from fq3_m import NSP, class_columns, sp_mats, expand, BLOCKS
from fq3_rec import Rec
import fq3_quadlib as L


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--cls", default="circ"); ap.add_argument("--nproc", type=int, default=12)
    a = ap.parse_args(); p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R); cls = a.cls
    qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % a.prime))); g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % a.prime)))
    us = State(os.path.join(HERE, "g3_u3_%s_%s.json" % (a.prime, cls)))
    LOG = open(os.path.join(HERE, "logs", "genus3_u3_%s_%s.log" % (a.prime, cls)), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    FAILED = []
    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)
    log("=" * 110); log("genus3_u3: class %s, p = 2^%s-1" % (cls, a.prime))
    Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}; nt = qs['n_taut']
    J2 = [g3['J2']['2'], g3['J2']['4']]; J3 = [g3['J']['2'], g3['J']['5']]
    rhs4 = [(-R.of(Fr(1, 750))) % p, (-R.of(Fr(-13, 8000))) % p, (-g3['A_M0']['2'][1]) % p, (-g3['A_M0']['5'][1]) % p]
    dot = lambda r, v: sum(r[k] * v[k] for k in range(NSP)) % p
    frame6 = frame_of(p, (0, 0), 34, 6); frame4 = frame_of(p, (0, 0), 34, 4)
    idx17 = IDX_BY_E[0] + IDX_BY_E[1] + IDX_BY_E[2] + IDX_BY_E[3]
    def unit(v): e = [0] * NSP; e[v] = 1; return e
    if cls == 'full': cols = [unit(v) for v in range(NSP)]
    else: cols, cnames = class_columns(R, 'circ')
    # ------------------------------------------------------------ N_1
    n1 = us.get('n1')
    if n1 is None:
        if cls == 'full':
            from fq3_fullsol import full_n1
            n1 = full_n1(R, Q, nt, J2, J3, frame4, a.nproc, log, check, us)
        elif cls == 'circstar':
            # the isolated solution Delta* with c_1(0) != 0, over the point where the c1-line of Pi touches the conic a_4 = 0 (found by the
            # elimination referee, g3_referee_u3.json; stored under 'circ_star')
            from fq3_circsys import CircSys
            S = CircSys(a.prime); v = json.load(open(os.path.join(HERE, "g3_cubic_%s.json" % a.prime)))['circ_star']
            check("circstar: the stored point satisfies every N_1-condition (4 linear rows, 5 quadrics, 9 cubics) and has c1_0 = 1", not any(S.residuals(v).values()) and v[0] == 1)
            n1 = S.to_sp(v)
        elif cls == 'circpi':
            # a point of the plane Pi of circulant solutions (genus3_pi.py): Delta = A + 2 B
            from fq3_circsys import CircSys
            S = CircSys(a.prime); cj = json.load(open(os.path.join(HERE, "g3_cubic_%s.json" % a.prime))); v = [(x + 2 * y) % p for x, y in zip(cj['Pi_A'], cj['Pi_B'])]
            check("circpi: Delta = A + 2 B satisfies every N_1-condition (4 linear rows, 5 quadrics, 9 cubics); c1 = (1, 2, -zeta^4/2 - 2 zeta^2)", not any(S.residuals(v).values()) and v[0] == 1 and v[1] == 2)
            n1 = S.to_sp(v)
        elif cls == 'circb':
            # the isolated Z/5-covariant solution with c_1(0) = 0, m_1 = 0 and non-zero low-order part (genus3_cubicb.py)
            from fq3_circsys import CircSys
            S = CircSys(a.prime); v = json.load(open(os.path.join(HERE, "g3_cubic_%s.json" % a.prime)))['circ_b_solutions'][0]
            check("circb: the stored circulant point satisfies every N_1-condition (4 linear rows, 5 quadrics, 9 cubics) and has c1_0 = 0, c1_1 = 1", not any(S.residuals(v).values()) and v[0] == 0 and v[1] == 1)
            log("   circulant coordinates of Delta: %s" % ["%s: %s" % (cnames[j], rec.fmt(v[j])) for j in range(15)])
            n1 = S.to_sp(v)
        else:
            rows15 = [[dot(r, cols[j]) for j in range(15)] for r in J2 + J3]; Kc = L.kernel_basis(rows15, 15, p); assert len(Kc) == 11
            lo = [j for j, nm in enumerate(cnames) if int(nm[1]) <= 3]
            Tp = L.kernel_basis([[Kc[t][j] for t in range(11)] for j in lo], 11, p); assert len(Tp) == 5
            s = [1, 2, -1, 3, 1]; tt = [sum(s[k] * Tp[k][t] for k in range(5)) % p for t in range(11)]
            c15 = [sum(Kc[t][j] * tt[t] for t in range(11)) % p for j in range(15)]
            check("circulant Delta is pure-high (orders 4, 5, 6) and non-zero", all(c15[j] == 0 for j in lo) and any(c15))
            log("   circulant coordinates of Delta: %s" % ["%s: %s" % (cnames[j], rec.fmt(c15[j])) for j in range(15) if c15[j]])
            n1 = [(nt[v] + sum(cols[j][v] * c15[j] for j in range(15))) % p for v in range(NSP)]
        check("N_1: the four O(u) equations (genus 2 and 3) and the five two-edge quadrics hold", [dot(r, n1) for r in J2 + J3] == rhs4 and all(Q[i].eval(n1) == 0 for i in IDX_2EDGE))
        us.put('n1', n1)
    N1 = sp_mats(n1, p)
    def do_pass(key, Ns, idxs, label):
        A = us.get(key)
        if A is not None: return {int(k): v for k, v in A.items()}
        t = time.time(); A, info = hat_pass(frame6, ('M0exp', Ns, 'right'), idxs, a.nproc, log, 21, label=label)
        log("   pass %s (%d graphs): %.0f s; maximal eps-degree %s" % (label, len(idxs), time.time() - t, max(i['maxdeg'] or 0 for i in info.values())))
        us.put(key, {str(i): A[i] for i in A}); return A
    def table(A, idxs, kmax=4):
        for i in idxs: log("      graph %2d (E=%d): %s" % (i, GD[i]['E'], [rec.fmt(x) for x in A[i][:kmax + 1]]))
    def generic(x, ns, nmax=8):
        t = list(x)
        for k, v in enumerate(ns[:nmax]):
            w = (k + 1) * (-1) ** k
            t = [(t[j] + w * v[j]) % p for j in range(len(t))]
        return t
    # ------------------------------------------------------------ pass 1 and N_2
    A1 = do_pass('A1', [N1], idx17, "P1 (N1)"); table(A1, idx17)
    check("P1: one-edge graphs O(u^2), two-edge graphs O(u^3), THREE-EDGE GRAPHS O(u^4) (the nine cubic conditions hold)",
          all(A1[i][k] == 0 for i in (2, 5) for k in range(2)) and all(A1[i][k] == 0 for i in IDX_2EDGE for k in range(3)) and all(A1[i][k] == 0 for i in IDX_3EDGE for k in range(4)))
    n2 = us.get('n2')
    if n2 is None:
        rows = [J3[0], J3[1]] + [Q[i].gradient(n1) for i in IDX_2EDGE]; rhs = [(-A1[2][2]) % p, (-A1[5][2]) % p] + [(-A1[i][3]) % p for i in IDX_2EDGE]
        Jr, x, ns, rank, cons = solve_system(R, rows, rhs, cols)
        log("   N_2: 7 linear equations on the class (%d parameters): rank %d, %s" % (len(cols), rank, "consistent, %d free parameters" % (len(cols) - rank) if cons else "INCONSISTENT"))
        check("N_2 system consistent", cons)
        n2 = expand(cols, generic(x, ns), p); us.put('n2', n2); us.put('n2_rank', rank)
    N2 = sp_mats(n2, p)
    A2 = do_pass('A2', [N1, N2], idx17, "P2 (N1, N2)"); table(A2, idx17)
    check("P2: one-edge graphs O(u^3), two-edge graphs O(u^4), three-edge graphs O(u^4)",
          all(A2[i][k] == 0 for i in (2, 5) for k in range(3)) and all(A2[i][k] == 0 for i in IDX_2EDGE for k in range(4)) and all(A2[i][k] == 0 for i in IDX_3EDGE for k in range(4)))
    n3 = us.get('n3')
    if n3 is None:
        Jr, x, ns, rank, cons = solve_system(R, [J3[0], J3[1]], [(-A2[2][3]) % p, (-A2[5][3]) % p], cols)
        check("N_3 system (2 equations, rank %d) consistent" % rank, cons)
        n3 = expand(cols, generic(x, ns), p); us.put('n3', n3)
    N3 = sp_mats(n3, p)
    # ------------------------------------------------------------ final pass on all 42 graphs
    A3 = do_pass('A3', [N1, N2, N3], IDX_ALL, "P3 (N1, N2, N3)"); table(A3, IDX_ALL, 5)
    edge = [i for i in IDX_ALL if i != 1]
    check("FINAL: all 41 edge graphs of genus 3 are O(u^4), i.e. REGULAR at the conifold point, for M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3)", all(A3[i][k] == 0 for i in edge for k in range(4)))
    tot = total_of(R, A3); check("FINAL: total over the 42 graphs = point-frame total through u^5 (M-independence)", tot == [x % p for x in g3['total_id']][:len(tot)])
    st = setup(frame6)
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3, GraphSummer6, UVertex
    from fq3_m import MU
    tc = tc_series(st['C'], R, 10); V = tc_expansion3(A3[1], tc, R)
    log("   vertex graph in t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(V[k]), k - 4) for k in range(6)))
    check("FINAL: the vertex graph alone carries the polar part: (1/25200) t_c^-4 + 0 t_c^-3 + 0 t_c^-2 + 0 t_c^-1", V[0] == R.of(Fr(1, 25200)) and V[1] == 0 and V[2] == 0 and V[3] == 0)
    us.put('vertex_tc', V[:6])
    from fq_msum import outer_graphs
    MUo = MU(R, st['M0'], [N1, N2, N3], 6, 'right'); t = time.time()
    G2 = GraphSummer6(st['H'], MUo, UVertex(R, MUo, 6, 21)).all_graphs(outer_graphs())
    log("   genus-2 graphs with the same M (%.0f s): %s" % (time.time() - t, {i: [rec.fmt(G2[i].coef(k)) for k in range(3)] for i in sorted(G2)}))
    check("GENUS 2: all six edge graphs are O(u^2) (regular) with the same M; vertex graph [u^0] = 1/1200", all(G2[i].coef(k) == 0 for i in range(2, 8) for k in range(2)) and G2[1].coef(0) == R.of(Fr(1, 1200)))
    tot2 = None
    for i in G2: tot2 = G2[i] if tot2 is None else tot2.add(G2[i])
    check("GENUS 2: total = 1/1200 - 7/6000 u + 97/72000 u^2 - 6809/18000000 u^3", [tot2.coef(k) for k in range(4)] == [R.of(x) for x in (Fr(1, 1200), Fr(-7, 6000), Fr(97, 72000), Fr(-6809, 18000000))])
    us.put('failed', FAILED)
    log("GENUS3_U3 (%s) ALL CHECKS PASSED" % cls if not FAILED else "GENUS3_U3 (%s) FAILURES: %s" % (cls, FAILED)); LOG.close()


if __name__ == "__main__":
    main()
