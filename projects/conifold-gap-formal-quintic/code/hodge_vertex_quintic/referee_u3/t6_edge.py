# -*- coding: utf-8 -*-
"""the structural identity behind C1, C2 and Pi, checked on the code's own edge / leaf tables at a random FULL-class (N_1, N_2):
   [u^1] E'(m, n) = (-1)^(n+1) D_{m+n+1},  D = M0 (N_1 - N_1^taut) M0^-1   (left-placed generator: exact grading by m + n + 1);
   [u^1] T'_d    = (-1)^(d+1) (row 0 of Psi_0) D_{d-1}  (d >= 2);   [u^0] of all edges and leaves = 0;
   d/dN_2 of [u^2] E'(m, n) = (-1)^(n+1) (M0 n_2 M0^-1)_{m+n+1}   (same linear map: the origin of the polar form in C2)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import random
from evlib import *
from rlib import Ad, load
from fq3_steps import frame_of
for tag in ("61", "89"):
    S, d, cs = load(tag); p = S.p; R = S.R; rng = random.Random(31); M = d['M']; Mi = d['Minv']
    st = setup(frame_of(p, (0, 0), 34, 4)); H = st['H']
    n1 = [rng.randrange(p) for _ in range(NSP)]; n2 = [rng.randrange(p) for _ in range(NSP)]; n2b = [rng.randrange(p) for _ in range(NSP)]
    Dl = Ad(M, Mi, sp_mats([(n1[k] - S.nt[k]) % p for k in range(NSP)], p), p)
    def EL(Ns): return GraphSummer6(H, MU(R, st['M0'], [sp_mats(n, p) for n in Ns], 4, 'right'), None).EL
    e1 = EL([n1]); e2 = EL([n1, n2]); e2b = EL([n1, n2b])
    mn = [(m, n) for m in range(6) for n in range(6 - m)]
    ok0 = all(e1.E(m, n)[i][j].coef(0) == 0 for (m, n) in mn for i in range(5) for j in range(5)) and all(e1.T(dd)[i].coef(0) == 0 for dd in range(2, 8) for i in range(5))
    ok1 = all((e1.E(m, n)[i][j].coef(1) - (-1) ** (n + 1) * Dl[m + n + 1][i][j]) % p == 0 for (m, n) in mn for i in range(5) for j in range(5))
    psi0 = [H.Q[0][i][0].coef(0) for i in range(5)]
    okT = all((e1.T(dd)[i].coef(1) - (-1) ** (dd + 1) * sum(psi0[j] * Dl[dd - 1][j][i] for j in range(5))) % p == 0 for dd in range(2, 8) for i in range(5))
    def lin2(nn, e):
        D2 = Ad(M, Mi, sp_mats(nn, p), p)
        return all((e.E(m, n)[i][j].coef(2) - e1.E(m, n)[i][j].coef(2) - (-1) ** (n + 1) * D2[m + n + 1][i][j]) % p == 0 for (m, n) in mn for i in range(5) for j in range(5))
    okN1 = all(e2.E(m, n)[i][j].coef(1) == e1.E(m, n)[i][j].coef(1) for (m, n) in mn for i in range(5) for j in range(5))
    xs = [H.x.coef(k) for k in range(4)]
    print("p = 2^%s-1: [u^0] edges, leaves = 0: %s | [u^1] E'(m,n) = (-1)^(n+1) D_(m+n+1): %s | [u^1] T'_d = (-1)^(d+1) Psi_0[0,.] D_(d-1): %s | [u^1] edges independent of N_2: %s | [u^2] E'(N_1,N_2) - [u^2] E'(N_1,0) = (-1)^(n+1) (Ad_M0 n_2)_(m+n+1): %s, %s | xhat = %s"
          % (tag, ok0, ok1, okT, okN1, lin2(n2, e2), lin2(n2b, e2b), xs))
