# -*- coding: utf-8 -*-
"""p1_audit.py -- referee PASSES: audit of the stored tables g3_u3_<p>_<cls>.json against g3_state_<p>.json, own code (read-only on the project).
Run: py -I p1_audit.py   (writes p1_audit.json next to this script)"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq3_graphs import outer_graphs3
from fq_fp import FpRing, P1, P2
from fq3_circsys import CircSys

PR = {"61": P1, "89": P2}
G3 = outer_graphs3(); EDG = {G['idx']: G['E'] for G in G3}
BYE = {E: sorted(i for i in EDG if EDG[i] == E) for E in range(7)}
OUT = {}


def rat(a, p, B=10 ** 10):
    """own rational reconstruction (half-gcd), None if not within the bound"""
    a %= p
    if a == 0: return Fr(0)
    r0, r1, s0, s1 = p, a, 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1)


def res(fr, p): return fr.numerator * pow(fr.denominator, p - 2, p) % p


def smul(a, b, p, n):
    c = [0] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]): c[i + j] = (c[i + j] + x * y) % p
    return c


def sinv(a, p, n):
    i0 = pow(a[0], p - 2, p); b = [i0] + [0] * (n - 1)
    for k in range(1, n):
        b[k] = (-i0 * sum(a[j] * b[k - j] for j in range(1, min(k, len(a) - 1) + 1))) % p
    return b


def tc_expand(acoefs, tcc, p, n=6):
    """u^-4 Sum a_k u^k in t_c, t_c = Sum tcc[j] u^(j+1): coefficients of t_c^-4 .. t_c^(n-5).  Own reversion: u = t h(t)."""
    # h with u = t h: fixed point of  h = 1 / (Sum_j tcc[j] (t h)^j)
    h = [1] + [0] * (n - 1)
    for _ in range(n + 1):
        f = [0] * n; pw = [1] + [0] * (n - 1)
        for j in range(min(len(tcc), n)):
            f = [(f[k] + tcc[j] * pw[k]) % p for k in range(n)]
            pw = smul(pw, [0] + h, p, n)                    # (t h)^(j+1)
        h = sinv(f, p, n)
    hi = sinv(h, p, n); hm4 = smul(smul(hi, hi, p, n), smul(hi, hi, p, n), p, n)
    T = [0] * n; pw = hm4                                    # pw = h^(k-4)
    for k in range(min(len(acoefs), n)):
        for m in range(n - k): T[k + m] = (T[k + m] + acoefs[k] * pw[m]) % p
        pw = smul(pw, h, p, n)
    return T, h


def rank_aug(rows, rhs, p):
    """(rank of rows, rank of the augmented matrix)"""
    def rk(M):
        M = [list(r) for r in M]; r = 0; nc = len(M[0]) if M else 0
        for c in range(nc):
            piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c] % p:
                    f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
            r += 1
        return r
    return rk(rows), rk([list(r) + [b] for r, b in zip(rows, rhs)])


def circ_coords(S, vec75):
    """solve cols . t = vec75 (75 equations, 15 unknowns): (t, exact?)"""
    p = S.p; M = [[S.cols[j][k] for j in range(15)] + [vec75[k] % p] for k in range(75)]
    r = 0; pivs = []
    for c in range(15):
        piv = next((i for i in range(r, 75) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(75):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        pivs.append(c); r += 1
    ok = all(M[i][15] % p == 0 for i in range(r, 75)); t = [0] * 15
    for i, c in enumerate(pivs): t[c] = M[i][15] % p
    return t, ok and len(pivs) == 15


def lowest(v):
    return next((k for k, x in enumerate(v) if x), None)


def audit(tag, cls):
    p = PR[tag]; out = {}
    f = os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls)); d = json.load(open(f))
    if 'A3' not in d: return {'status': 'no A3 (not finished)', 'keys': sorted(d)}
    g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    A = {k: {int(i): [x % p for x in v] for i, v in d[k].items()} for k in ('A1', 'A2', 'A3')}
    A1, A2, A3 = A['A1'], A['A2'], A['A3']
    out['shape'] = {k: (len(A[k]), sorted(set(len(v) for v in A[k].values()))) for k in A}
    # (1) A3: the 41 edge graphs vanish at u^0..u^3
    edge = [i for i in sorted(A3) if i != 1]
    bad = [(i, k) for i in edge for k in range(4) if A3[i][k]]
    out['c1_n_edge'] = len(edge); out['c1_bad'] = bad; out['c1_ok'] = (len(edge) == 41 and not bad and sorted(A3) == list(range(1, 43)) and EDG[1] == 0)
    low = {i: lowest(A3[i]) for i in edge}
    out['c1_lowest_order_by_E'] = {E: sorted(set(str(low[i]) for i in BYE[E])) for E in range(1, 7)}
    out['c1_u4_nonzero'] = [i for i in edge if A3[i][4]]; out['c1_u5_nonzero'] = [i for i in edge if A3[i][5]]
    out['c1_below_E_nonzero'] = [(i, k) for i in edge for k in range(min(EDG[i], 6)) if A3[i][k]]
    # (2) total
    tot = [sum(A3[i][k] for i in A3) % p for k in range(6)]; tid = [x % p for x in g3['total_id']]
    out['c2_ok'] = (tot == tid and len(tid) == 6)
    out['c2_total_rat'] = [str(rat(x, p)) for x in tot]
    out['c2_vertex_eq_total_u0_3'] = [A3[1][k] == tid[k] for k in range(4)]
    # (3) vertex graph in t_c
    tcc = [x % p for x in g3['tc']]
    out['c3_tc_rat'] = [str(rat(x, p)) for x in tcc]
    V, h = tc_expand(A3[1], tcc, p)
    out['c3_V_rat'] = [str(rat(x, p)) for x in V]
    out['c3_ok'] = (V[0] == res(Fr(1, 25200), p) and V[1] == 0 and V[2] == 0 and V[3] == 0)
    out['c3_matches_stored_vertex_tc'] = (V == [x % p for x in d['vertex_tc']])
    Tt, _ = tc_expand(tid, tcc, p); out['c3_total_tc_rat'] = [str(rat(x, p)) for x in Tt]
    out['c3_total_tc_matches_state'] = (Tt == [x % p for x in g3['total_tc']])
    # exact check of the reversion in Q for the total (if the total is rational)
    # (4) A1, A2 patterns
    def pattern(T):
        r = {}
        for E in (1, 2, 3):
            r["E%d" % E] = {i: [k for k in range(6) if T[i][k]] for i in BYE[E]}
        return r
    out['c4_A1_nonzero_orders'] = pattern(A1); out['c4_A2_nonzero_orders'] = pattern(A2); out['c4_A3_nonzero_orders'] = pattern(A3)
    out['c4_A1_ok'] = (all(A1[i][k] == 0 for i in BYE[1] for k in range(2)) and all(A1[i][k] == 0 for i in BYE[2] for k in range(3))
                       and all(A1[i][k] == 0 for i in BYE[3] for k in range(4)))
    out['c4_A1_not_regular'] = [(i, k) for i in BYE[1] + BYE[2] + BYE[3] for k in range(4) if A1[i][k]]
    out['c4_A2_ok'] = (all(A2[i][k] == 0 for i in BYE[1] for k in range(3)) and all(A2[i][k] == 0 for i in BYE[2] + BYE[3] for k in range(4)))
    out['c4_A2_not_regular'] = [(i, k) for i in BYE[1] + BYE[2] + BYE[3] for k in range(4) if A2[i][k]]
    # cross-pass consistency (independent passes must agree where the later N_k cannot enter)
    out['x_vertex_A1_A2_A3_u01'] = [A1[1][k] == A2[1][k] == A3[1][k] for k in range(2)]
    out['x_vertex_A2_A3_u2'] = (A2[1][2] == A3[1][2])
    out['x_3edge_u4_A2_eq_A3'] = all(A2[i][4] == A3[i][4] for i in BYE[3]); out['x_3edge_u4_nonzero'] = [i for i in BYE[3] if A3[i][4]]
    out['x_2edge_u4_A2_ne_A3'] = [i for i in BYE[2] if A2[i][4] != A3[i][4]]
    out['x_1edge_u3_A2_ne_A3'] = [i for i in BYE[1] if A2[i][3] != A3[i][3]]
    # (5) n1, n2, n3
    S = CircSys(tag); nt = [x % p for x in qs['n_taut']]; n1, n2, n3 = d['n1'], d['n2'], d['n3']
    dl, ok1 = circ_coords(S, [(n1[k] - nt[k]) % p for k in range(75)])
    out['c5_delta_circulant'] = ok1; out['c5_delta_nonzero_coords'] = [S.names[j] for j in range(15) if dl[j]]
    out['c5_delta_orders'] = sorted(set(int(S.names[j][1]) for j in range(15) if dl[j]))
    out['c5_delta_rat'] = {S.names[j]: str(rat(dl[j], p)) for j in range(15) if dl[j]}
    rs = S.residuals(dl); out['c5_N1_conditions_nonzero'] = [k for k, v in rs.items() if v]; out['c5_n_conditions'] = len(rs)
    for nm, v in (('n2', n2), ('n3', n3)):
        t, ok = circ_coords(S, v)
        out['c5_%s_circulant' % nm] = ok; out['c5_%s_nonzero_coords' % nm] = [S.names[j] for j in range(15) if t[j]]
        out['c5_%s_nonzero75' % nm] = sum(1 for x in v if x % p)
    # the N_2 and N_3 systems rebuilt here: rows on the 15 circulant columns
    dot = lambda r, v: sum(x * y for x, y in zip(r, v)) % p
    J3 = [g3['J']['2'], g3['J']['5']]
    Qf = {i: qs['Q'][str(i)] for i in (3, 6, 7, 8, 17)}
    def gradQ(i):
        b, Am = Qf[i]['b'], Qf[i]['A']
        return [(b[a] + 2 * sum(Am[a][c] * n1[c] for c in range(75))) % p for a in range(75)]
    def Qval(i, n):
        b, Am = Qf[i]['b'], Qf[i]['A']
        return (Qf[i]['c0'] + dot(b, n) + sum(n[a] * sum(Am[a][c] * n[c] for c in range(75)) for a in range(75))) % p
    out['c5_Q_at_n1'] = {i: Qval(i, n1) for i in Qf}
    rows = J3 + [gradQ(i) for i in (3, 6, 7, 8, 17)]
    rhs = [(-A1[2][2]) % p, (-A1[5][2]) % p] + [(-A1[i][3]) % p for i in (3, 6, 7, 8, 17)]
    R15 = [[dot(r, S.cols[j]) for j in range(15)] for r in rows]
    out['n2_rank_aug'] = rank_aug(R15, rhs, p); out['n2_rank_stored'] = d.get('n2_rank')
    out['n2_rank_gradQ_only'] = rank_aug(R15[2:], rhs[2:], p)[0]
    out['n2_solves'] = [dot(r, n2) == b for r, b in zip(rows, rhs)]; out['n2_rhs_nonzero'] = [bool(b) for b in rhs]
    rhs3 = [(-A2[2][3]) % p, (-A2[5][3]) % p]
    out['n3_rank_aug'] = rank_aug([R15[0], R15[1]], rhs3, p); out['n3_solves'] = [dot(r, n3) == b for r, b in zip(J3, rhs3)]
    out['n3_rhs_nonzero'] = [bool(b) for b in rhs3]
    # linear responses predicted by the docstring: differences between passes
    out['resp_1edge_u2_A2-A1_eq_J.n2'] = [(A2[i][2] - A1[i][2]) % p == dot(J3[k], n2) for k, i in enumerate((2, 5))]
    out['resp_2edge_u3_A2-A1_eq_gradQ.n2'] = [(A2[i][3] - A1[i][3]) % p == dot(gradQ(i), n2) for i in (3, 6, 7, 8, 17)]
    out['resp_1edge_u3_A3-A2_eq_J.n3'] = [(A3[i][3] - A2[i][3]) % p == dot(J3[k], n3) for k, i in enumerate((2, 5))]
    out['resp_1edge_u1_J.n1_eq_-A_M0'] = [dot(J3[k], n1) == (-g3['A_M0'][str(i)][1]) % p for k, i in enumerate((2, 5))]
    out['failed_stored'] = d.get('failed')
    out['_tables'] = A; out['_n'] = (n1, n2, n3); out['_delta'] = dl; out['_V'] = V
    return out


def main():
    res_all = {}
    for cls in ("circ", "circpi", "circb", "circstar"):
        for tag in ("61", "89"):
            res_all[(cls, tag)] = audit(tag, cls)
    # (6) both primes: structure and rational values
    for cls in ("circ", "circpi", "circb", "circstar"):
        a, b = res_all[(cls, "61")], res_all[(cls, "89")]
        if '_tables' not in a or '_tables' not in b: continue
        same = [k for k in a if not k.startswith('_') and k not in ('c2_total_rat', 'c3_V_rat', 'c3_tc_rat', 'c3_total_tc_rat', 'c5_delta_rat', 'c5_Q_at_n1') and a[k] != b[k]]
        st = dict(agree=0, mismatch=0, unrec=0, zero_both=0, zero_one=0, n=0); mism = []; recs = {}
        for key in ('A1', 'A2', 'A3'):
            for i in a['_tables'][key]:
                for k in range(6):
                    x61, x89 = a['_tables'][key][i][k], b['_tables'][key][i][k]; st['n'] += 1
                    if x61 == 0 and x89 == 0: st['zero_both'] += 1; continue
                    if (x61 == 0) != (x89 == 0): st['zero_one'] += 1; mism.append((key, i, k)); continue
                    r = rat(x89, P2)
                    if r is None or r.denominator % P1 == 0: st['unrec'] += 1; continue
                    if res(r, P1) == x61: st['agree'] += 1; recs["%s[%d][u^%d]" % (key, i, k)] = str(r)
                    else: st['mismatch'] += 1; mism.append((key, i, k, str(r)))
        dd = dict(agree=0, mismatch=0, unrec=0)
        for j in range(15):
            x61, x89 = a['_delta'][j], b['_delta'][j]
            if x61 == 0 and x89 == 0: continue
            r = rat(x89, P2)
            if r is None: dd['unrec'] += 1
            elif res(r, P1) == x61: dd['agree'] += 1
            else: dd['mismatch'] += 1
        vv = [(str(rat(x, P2)), rat(x, P2) is not None and res(rat(x, P2), P1) == y) for x, y in zip(b['_V'], a['_V'])]
        res_all[(cls, "x")] = dict(structural_keys_differing=same, table_entries=st, mismatches=mism[:20], rational_entries=recs,
                                   delta_coords=dd, vertex_tc_cross=vv)
    # distinctness of the tables between classes (no shared cache)
    for tag in ("61", "89"):
        fin = [c for c in ("circ", "circpi", "circb", "circstar") if '_tables' in res_all[(c, tag)]]; dist = {}
        for x in range(len(fin)):
            for y in range(x + 1, len(fin)):
                ta, tb = res_all[(fin[x], tag)]['_tables'], res_all[(fin[y], tag)]['_tables']
                dist["%s|%s" % (fin[x], fin[y])] = {k: sum(1 for i in ta[k] for j in range(6) if ta[k][i][j] != tb[k][i][j]) for k in ('A1', 'A2', 'A3')}
                dist["%s|%s" % (fin[x], fin[y])]['n_equal'] = [res_all[(fin[x], tag)]['_n'][m] == res_all[(fin[y], tag)]['_n'][m] for m in range(3)]
        res_all[("distinct", tag)] = dist
    dump = {"%s_%s" % k: {kk: vv for kk, vv in v.items() if not kk.startswith('_')} for k, v in res_all.items()}
    json.dump(dump, open(os.path.join(HERE, "p1_audit.json"), "w"), indent=1, default=str)
    for k in dump:
        print("=" * 30, k)
        for kk, vv in dump[k].items(): print("  %s: %s" % (kk, str(vv)[:1400]))


if __name__ == "__main__":
    main()
