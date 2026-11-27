# -*- coding: utf-8 -*-
"""p4_extra.py -- referee PASSES, small follow-ups on the stored data (no graph evaluation):
 (a) point-frame total: the rationals quoted in the note against the residues at both primes (u^0..u^5), and the t_c^0 coefficient;
 (b) pure-high class: support of the five polar forms grad Q_Gamma(N_1) in the circulant coordinates of n_2 (which orders of N_2 they see);
 (c) support of the two genus-3 rows J3 in the circulant coordinates; which coordinates of n2, n3 are untouched by the equations;
 (d) Z/5-covariant conjugate test: values of a table entry times zeta^-k rational? (Q(zeta)-monomial entries), cross-prime;
 (e) vertex graph t_c^0, t_c^1 of each class against the total (they must differ by the regular edge-graph contributions).
Run: py -I p4_extra.py"""
import os, sys, json
from fractions import Fraction as Fr
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from fq3_circsys import CircSys

PR = {"61": P1, "89": P2}
res = lambda fr, p: fr.numerator * pow(fr.denominator, p - 2, p) % p


def rat(a, p, B=10 ** 10):
    a %= p
    if a == 0: return Fr(0)
    r0, r1, s0, s1 = p, a, 0, 1
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1)


def main():
    TOT = [Fr(1, 25200), Fr(-1, 9000), Fr(407, 3780000), Fr(-151, 3780000), Fr(-29663, 9072000000), Fr(364981, 75600000000)]
    for tag in ("61", "89"):
        p = PR[tag]; g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
        print("p=2^%s-1 (a) total_id == note rationals: %s; total_tc[4] == -57737/8400000000: %s; total_tc[1..3] zero: %s" % (
            tag, [res(TOT[k], p) == g3['total_id'][k] % p for k in range(6)], res(Fr(-57737, 8400000000), p) == g3['total_tc'][4] % p,
            [g3['total_tc'][k] % p == 0 for k in (1, 2, 3)]))
        S = CircSys(tag); dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
        J3 = [g3['J']['2'], g3['J']['5']]
        J15 = [[dot(r, S.cols[j]) for j in range(15)] for r in J3]
        print("   (c) J3 rows: non-zero circulant coordinates: %s" % [[S.names[j] for j in range(15) if r[j]] for r in J15])
        for cls in ("circ", "circpi", "circb"):
            d = json.load(open(os.path.join(PROJ, "g3_u3_%s_%s.json" % (tag, cls))))
            if 'A3' not in d: continue
            n1 = d['n1']
            def gradQ(i):
                b, Am = qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']
                return [(b[a] + 2 * sum(Am[a][c] * n1[c] for c in range(75))) % p for a in range(75)]
            G15 = {i: [dot(gradQ(i), S.cols[j]) for j in range(15)] for i in (3, 6, 7, 8, 17)}
            sup = sorted(set(j for i in G15 for j in range(15) if G15[i][j]))
            print("   (b) %s: grad Q_Gamma(N_1) rows see the circulant coordinates %s; zero rows: %s" % (cls, [S.names[j] for j in sup], [i for i in G15 if not any(G15[i])]))
            # (e) vertex t_c^0 against the total
            V = [x % p for x in d['vertex_tc']]; T = [x % p for x in g3['total_tc']]
            print("   (e) %s: vertex_tc == total_tc at t_c^-4..t_c^-1: %s; at t_c^0, t_c^1: %s (edge graphs contribute there)" % (cls, V[:4] == T[:4], [V[k] == T[k] for k in (4, 5)]))
            A3 = {int(i): [x % p for x in v] for i, v in d['A3'].items()}
            e4 = sum(A3[i][4] for i in A3 if i != 1) % p
            print("       sum of the 41 edge graphs at u^4 (= their t_c^0 contribution): non-zero %s; vertex t_c^0 + it == total t_c^0: %s" % (bool(e4), (V[4] + e4) % p == T[4]))
    # (d) cross-prime test of zeta-monomial entries r * zeta^k (r rational at 2^89-1 with bound 1e10, same r with some zeta^k' at 2^61-1)
    R1, R2 = FpRing(P1), FpRing(P2)
    for cls in ("circ", "circpi", "circb"):
        d1 = json.load(open(os.path.join(PROJ, "g3_u3_61_%s.json" % cls))); d2 = json.load(open(os.path.join(PROJ, "g3_u3_89_%s.json" % cls)))
        if 'A3' not in d1 or 'A3' not in d2: continue
        n = dict(mono=0, mono_ok=0, mono_bad=0, other=0, rational=0)
        for key in ('A1', 'A2', 'A3'):
            for i in d2[key]:
                for k in range(6):
                    x2, x1 = d2[key][i][k] % P2, d1[key][i][k] % P1
                    if x2 == 0: continue
                    if rat(x2, P2) is not None: n['rational'] += 1; continue
                    hit = None
                    for e in range(1, 5):
                        r = rat(x2 * R2.z(-e) % P2, P2)
                        if r is not None: hit = r; break
                    if hit is None: n['other'] += 1; continue
                    n['mono'] += 1
                    if any(res(hit, P1) * R1.z(e) % P1 == x1 for e in range(5)): n['mono_ok'] += 1
                    else: n['mono_bad'] += 1
        print("(d) %s: non-zero table entries: rational %d, r*zeta^k %d (cross-prime consistent %d, inconsistent %d), neither %d" % (
            cls, n['rational'], n['mono'], n['mono_ok'], n['mono_bad'], n['other']))


if __name__ == "__main__":
    main()
