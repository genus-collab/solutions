# -*- coding: utf-8 -*-
"""p10_exact.py -- referee PASSES: exact (Fraction) t_c expansion of the vertex graph.  In every finished class [u^0..u^3] of the vertex
graph are the rationals 1/25200, -1/9000, 407/3780000, -151/3780000 (= point-frame total, checked mod both primes in p1_audit.py), and
t_c = u + 7/10 u^2 + 41/75 u^3 + 9/20 u^4 + O(u^5).  Polar part of u^-4 (a0 + a1 u + a2 u^2 + a3 u^3) in t_c, by direct substitution
t_c^-k = u^-k (1 + 7/10 u + 41/75 u^2 + 9/20 u^3)^-k (no series reversion).   Run: py -I p10_exact.py"""
from fractions import Fraction as Fr
a = [Fr(1, 25200), Fr(-1, 9000), Fr(407, 3780000), Fr(-151, 3780000)]
c = [Fr(1), Fr(7, 10), Fr(41, 75), Fr(9, 20)]
def mul(x, y, n=4): return [sum(x[i] * y[k - i] for i in range(k + 1) if i < len(x) and k - i < len(y)) for k in range(n)]
def inv(x, n=4):
    b = [1 / x[0]] + [Fr(0)] * (n - 1)
    for k in range(1, n): b[k] = -b[0] * sum(x[j] * b[k - j] for j in range(1, k + 1) if j < len(x))
    return b
ci = inv(c); pw = {0: [Fr(1), 0, 0, 0]}
for k in range(1, 5): pw[k] = mul(pw[k - 1], ci)            # (t_c/u)^-k
# solve u^-4 Sum a_k u^k = Sum_j V_j t_c^(j-4) mod u^0:  V_j t_c^(j-4) = V_j u^(j-4) (t_c/u)^-(4-j)
V = []; rem = list(a)
for j in range(4):
    V.append(rem[j]); s = pw[4 - j]
    rem = [rem[k] - (V[j] * s[k - j] if k >= j else 0) for k in range(4)]
print("exact polar part: t_c^-4: %s, t_c^-3: %s, t_c^-2: %s, t_c^-1: %s" % tuple(V))
