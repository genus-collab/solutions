"""Lemma-K linear test: response of the [u^2] / [u^3] regularity conditions to left dressings u^n zhat^k E(m, m') of M_nat (K-regular vs not)"""
import sys, json, os
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, HERE)
from fq_fp import P2, P1
from fq3_par import setup
from fq3_rec import Rec
from fq_symbol import mat_series_mul
from fq3_m import KMAX, sp_params, in_sp, NSP, Z5
from fq3_quad import QuadForm, IDX_2EDGE
W = (0, 4, 3, 2, 1)

def kreg(n, k, m, mp):
    a, b = (-m - 1) % 5, (-mp - 1) % 5
    return 5 * n - 4 * k + 2 - W[a] - W[b] >= 0

for prime, tag in ((P1, "61"), (P2, "89")):
    st = setup((prime, "0", "0", 34, 6)); R = st['R']; p = R.p; rec = Rec(R); M0 = st['M0']; z = R.zeta; inv25 = pow(25, p - 2, p)
    g3 = json.load(open(os.path.join(HERE, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(HERE, "g3_quad_%s.json" % tag)))
    n1 = json.load(open(os.path.join(HERE, "g3_nat_%s_nat.json" % tag)))['n1']
    Q = {i: QuadForm(p, qs['Q'][str(i)]['c0'], qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']) for i in IDX_2EDGE}
    rows2 = [g3['J']['2'], g3['J']['5']] + [Q[i].gradient(n1) for i in IDX_2EDGE]; rows3 = [g3['J']['2'], g3['J']['5']]
    dot = lambda r, v: sum(r[j] * v[j] for j in range(NSP)) % p
    def direction(k, m, mp):
        sgn = 1 if k % 2 == 1 else -1
        if m == mp and sgn == -1: return None
        Xh = {(m, mp): 25, (mp, m): 25 * sgn % p} if m != mp else {(m, m): 25}
        X = [[inv25 * sum(pow(z, (-a * i) % 5, p) * v * pow(z, (-b * j) % 5, p) for (a, b), v in Xh.items()) % p for j in range(5)] for i in range(5)]
        S = [Z5() for _ in range(KMAX + 1)]; S[k] = X
        N = mat_series_mul(mat_series_mul(M0.Minv, S, p, KMAX), M0.M, p, KMAX)
        return sp_params(N, p) if in_sp(N, p) else 'not sp'
    print("p = 2^%s-1" % tag)
    for n, rows, label in ((2, rows2, "[u^2] one-edge (2), [u^3] two-edge (5)"), (3, rows3, "[u^3] one-edge (2)")):
        print("  order u^%d: conditions %s" % (n, label))
        for k in range(1, 7):
            for m in range(5):
                mp = (k - m) % 5
                if mp < m: continue
                d = direction(k, m, mp)
                if d is None: continue
                if d == 'not sp': print("     zhat^%d (%d,%d): not in sp" % (k, m, mp)); continue
                resp = [dot(r, d) for r in rows]
                print("     zhat^%d (%d,%d) %-13s response: %s" % (k, m, mp, "K-regular" if kreg(n, k, m, mp) else "NOT K-regular", "ALL ZERO" if not any(resp) else [rec.fmt(x) for x in resp]))
