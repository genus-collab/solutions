# -*- coding: utf-8 -*-
"""p7_star.py -- referee PASSES: partial audit of the unfinished class circstar (Delta*), from whatever g3_u3_<p>_circstar.json holds:
A1 pattern, the N_2 system (rank, consistency, stored n2 solves it), A2 pattern and the N_3 system if present, A3 checks if present,
and cross-prime agreement of every entry that reconstructs rationally at 2^89-1 (bound 1e10).   Run: py -I p7_star.py"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import P1, P2
from fq3_circsys import CircSys
from p1_audit import rat, res, rank_aug, circ_coords, BYE, tc_expand
PR = {"61": P1, "89": P2}; D = {}
for tag in ("61", "89"):
    p = PR[tag]; d = json.load(open(os.path.join(PROJ, "g3_u3_%s_circstar.json" % tag))); D[tag] = d
    g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    S = CircSys(tag); n1 = d['n1']; nt = qs['n_taut']; dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
    dl, ok = circ_coords(S, [(n1[k] - nt[k]) % p for k in range(75)])
    print("p=2^%s-1 keys %s" % (tag, sorted(d)))
    print("   Delta* circulant: %s; non-zero coords %s; N_1-conditions violated: %s" % (ok, [S.names[j] for j in range(15) if dl[j]], [k for k, v in S.residuals(dl).items() if v]))
    if tag == "89": print("   Delta* coordinates (rational reconstruction, bound 1e10): %s" % {S.names[j]: str(rat(dl[j], p)) for j in range(15) if dl[j]})
    D[tag + "dl"] = dl
    if 'A1' in d:
        A1 = {int(i): [x % p for x in v] for i, v in d['A1'].items()}
        okA1 = all(A1[i][k] == 0 for i in BYE[1] for k in range(2)) and all(A1[i][k] == 0 for i in BYE[2] for k in range(3)) and all(A1[i][k] == 0 for i in BYE[3] for k in range(4))
        print("   A1 (%d graphs): pattern ok %s; non-regular coefficients %s" % (len(A1), okA1, [(i, k) for i in BYE[1] + BYE[2] + BYE[3] for k in range(4) if A1[i][k]]))
        def gradQ(i):
            b, Am = qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']
            return [(b[a] + 2 * sum(Am[a][c] * n1[c] for c in range(75))) % p for a in range(75)]
        J3 = [g3['J']['2'], g3['J']['5']]; rows = J3 + [gradQ(i) for i in (3, 6, 7, 8, 17)]
        rhs = [(-A1[2][2]) % p, (-A1[5][2]) % p] + [(-A1[i][3]) % p for i in (3, 6, 7, 8, 17)]
        R15 = [[dot(r, S.cols[j]) for j in range(15)] for r in rows]
        print("   N_2 system: (rank, augmented rank) = %s; rhs non-zero %s; stored n2_rank %s" % (rank_aug(R15, rhs, p), [bool(b) for b in rhs], d.get('n2_rank')))
        if 'n2' in d:
            t, okc = circ_coords(S, d['n2'])
            print("   stored n2 solves the 7 equations: %s; circulant %s; non-zero coords %d" % ([dot(r, d['n2']) == b for r, b in zip(rows, rhs)], okc, sum(1 for x in t if x)))
    if 'A2' in d:
        A2 = {int(i): [x % p for x in v] for i, v in d['A2'].items()}
        print("   A2: pattern ok %s; non-regular %s" % (all(A2[i][k] == 0 for i in BYE[1] for k in range(3)) and all(A2[i][k] == 0 for i in BYE[2] + BYE[3] for k in range(4)),
              [(i, k) for i in BYE[1] + BYE[2] + BYE[3] for k in range(4) if A2[i][k]]))
        if 'n3' in d: print("   stored n3 solves J3 . n3 = -[u^3]: %s" % [dot(J3[k], d['n3']) == (-A2[i][3]) % p for k, i in enumerate((2, 5))])
    if 'A3' in d:
        A3 = {int(i): [x % p for x in v] for i, v in d['A3'].items()}
        bad = [(i, k) for i in A3 if i != 1 for k in range(4) if A3[i][k]]; tot = [sum(A3[i][k] for i in A3) % p for k in range(6)]
        V, _ = tc_expand(A3[1], [x % p for x in g3['tc']], p)
        print("   A3 (%d graphs): edge-graph coefficients u^0..u^3 non-zero: %s; total == total_id: %s; vertex t_c: %s, zeros %s" % (
            len(A3), bad, tot == [x % p for x in g3['total_id']], V[0] == res(Fr(1, 25200), p), [V[k] == 0 for k in (1, 2, 3)]))
st = dict(agree=0, mismatch=0, unrec=0, zero_both=0, zero_one=0); ex = []
for key in ('A1', 'A2', 'A3'):
    if key not in D["61"] or key not in D["89"]: continue
    for i in D["89"][key]:
        for k in range(6):
            x89, x61 = D["89"][key][i][k] % P2, D["61"][key][i][k] % P1
            if x89 == 0 and x61 == 0: st['zero_both'] += 1; continue
            if (x89 == 0) != (x61 == 0): st['zero_one'] += 1; continue
            r = rat(x89, P2)
            if r is None: st['unrec'] += 1
            elif res(r, P1) == x61: st['agree'] += 1; ex.append("%s[%s][u^%d]=%s" % (key, i, k, r))
            else: st['mismatch'] += 1; ex.append("MISMATCH %s[%s][u^%d]" % (key, i, k))
print("cross-prime (tables present at both primes):", st); print("   examples:", ex[:14])
dd = [(rat(a, P2), b) for a, b in zip(D["89dl"], D["61dl"])]
print("Delta* coordinates: rational at 2^89-1 and equal mod 2^61-1: %s" % [r is not None and res(r, P1) == b for r, b in dd])
