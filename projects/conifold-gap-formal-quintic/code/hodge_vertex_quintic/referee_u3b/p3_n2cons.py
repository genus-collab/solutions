# -*- coding: utf-8 -*-
"""p3_n2cons.py -- referee PASSES: is the N_2 system (7 equations: J3 . n2 = -[u^2] G_Gamma(N_1, 0) for Gamma = 2, 5 and
grad Q_Gamma(N_1) . n2 = -[u^3] G_Gamma(N_1, 0) for the five two-edge graphs) CONSISTENT at N_1-solutions other than the three points for
which genus3_u3.py was run?  Only the seven graphs with at most two edges are evaluated (21-node interpolation, N_1 only), so this is not a
full pass.  Points (circulant 15-vectors of Delta, all satisfying the 18 N_1-conditions by CircSys.residuals):
   pi_rand  : x A + y B, a random point of the plane Pi (seeded)
   pi_A, pi_B : the stored generators of Pi
   extra0..3: the four stored special points 'circ_b_extra' on c_1(0) = 0 (b = +-zeta^2)
   star     : the stored Delta* ('circ_star')
   ctrl     : control = the verified point A + 2 B of the driver (must reproduce rank 6, consistent, and the stored A1 entries)
Read-only on the project.  Run: py -I p3_n2cons.py <61|89> name[,name..]   -> p3_n2cons_<p>_<names>.json"""
import os, sys, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from fq3_steps import frame_of
from fq3_par import job_hat
from fq3_m import sp_mats, NSP
from fq3_circsys import CircSys

TWO = (3, 6, 7, 8, 17)


def rank(M, p):
    M = [list(r) for r in M]; r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
    return r


def main():
    tag = sys.argv[1]; names = sys.argv[2].split(","); p = P2 if tag == "89" else P1
    S = CircSys(tag); cj = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % tag))); qs = json.load(open(os.path.join(PROJ, "g3_quad_%s.json" % tag)))
    A = [x % p for x in cj['Pi_A']]; B = [x % p for x in cj['Pi_B']]
    rng = random.Random(20261110 + int(tag)); x, y = rng.randrange(2, p), rng.randrange(2, p)
    pts = {'pi_rand': [(x * a + y * b) % p for a, b in zip(A, B)], 'pi_A': A, 'pi_B': B, 'star': [v % p for v in cj['circ_star']],
           'ctrl': [(a + 2 * b) % p for a, b in zip(A, B)]}
    for k, v in enumerate(cj['circ_b_extra']): pts['extra%d' % k] = [t % p for t in v]
    frame6 = frame_of(p, (0, 0), 34, 6); J3 = [g3['J']['2'], g3['J']['5']]; out = {}
    dot = lambda r, v: sum(a * b for a, b in zip(r, v)) % p
    for nm in names:
        t0 = time.time(); v = pts[nm]; res = S.residuals(v); bad = [k for k, r in res.items() if r]
        n1 = S.to_sp(v); N1 = sp_mats(n1, p); G = {}
        for i in (2, 5) + TWO:
            idx, coefs, info = job_hat((frame6, ('M0exp', [N1], 'right'), i, 21)); G[i] = [c % p for c in coefs]
        def gradQ(i):
            b, Am = qs['Q'][str(i)]['b'], qs['Q'][str(i)]['A']
            return [(b[a] + 2 * sum(Am[a][c] * n1[c] for c in range(NSP))) % p for a in range(NSP)]
        rows = J3 + [gradQ(i) for i in TWO]; rhs = [(-G[2][2]) % p, (-G[5][2]) % p] + [(-G[i][3]) % p for i in TWO]
        R15 = [[dot(r, S.cols[j]) for j in range(15)] for r in rows]
        rk = rank(R15, p); rka = rank([r + [b] for r, b in zip(R15, rhs)], p)
        o = dict(c1=[v[0], v[1], v[2]] if max(v[:3]) < 10 ** 6 else "c1_0=%d, c1_1=%d" % (v[0] if v[0] < 10 ** 6 else -1, v[1] if v[1] < 10 ** 6 else -1),
                 N1_conditions_violated=bad,
                 low_orders_vanish=(all(G[i][k] == 0 for i in (2, 5) for k in range(2)) and all(G[i][k] == 0 for i in TWO for k in range(3))),
                 rhs_nonzero=[bool(b) for b in rhs], rank=rk, rank_augmented=rka, consistent=(rk == rka), seconds=round(time.time() - t0))
        if nm == 'ctrl':
            d = json.load(open(os.path.join(PROJ, "g3_u3_%s_circpi.json" % tag)))
            o['ctrl_n1_eq_stored'] = (n1 == [t % p for t in d['n1']])
            o['ctrl_A1_eq_stored'] = all(G[i] == [t % p for t in d['A1'][str(i)]] for i in G)
        out[nm] = o
        print("p=2^%s-1 %-8s: %s" % (tag, nm, o), flush=True)
        json.dump(out, open(os.path.join(HERE, "p3_n2cons_%s_%s.json" % (tag, "_".join(names))), "w"), indent=1)


if __name__ == "__main__":
    main()
