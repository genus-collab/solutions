# Rank-1 sanity check of the sign of the PPZ translation T_b = z(1 - R^{-1}(z) 1).
# V = Q, eta = 1, Omega = trivial TFT (all classes 1), R(z) = exp(a z + ...).
# Unit axiom at (0,4): the degree-1 part of (R T Omega)_{0,4}(1,1,1,1) must integrate to 0.
# Contributions to the integral over M_{0,4}bar:
#  legs: R^{-1}(psi_l) = 1 - a psi_l + ..., sum_l int psi_l = 4  -> -4a
#  edges: 3 one-edge graphs, edge = (1 - R^{-1}(x)R^{-1}(y))/(x+y) = a + O(x,y)  -> +3a
#  translation: T(z) = s*a z^2 + ..., p_*(psi_5^2) = kappa_1, int kappa_1 = int_{M05} psi_5^2 = 1 -> +s*a
from fractions import Fraction as F
a = F(7, 3)
for s, name in [(1, "PPZ sign T_b"), (-1, "opposite sign")]:
    total = -4*a + 3*a + s*a
    print(name, "degree-1 integral =", total)
