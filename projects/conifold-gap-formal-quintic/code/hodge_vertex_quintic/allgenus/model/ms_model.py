# -*- coding: utf-8 -*-
"""ms_model.py -- genus-0 structure of the LG model lambda = -(x^4/4 + v/x), primitive form dx, on its Z/5 line, in the flat
basis H^a (a in Z/5; H = d/dv), as needed to organise the correlators of ms_engine (which are the model at v = 1).

Homogeneity: weight(v) = 5/4, weight(z) = 1, q_a = c_a/4, c = (0,-1,3,2,1); a correlator <prod tau_k(H^a)>_g is
  const * v^e,  e = 2 - 2g - (1/5) Sum (4 - c_a - 4k)   (ms_engine.vexp), up to log v at genus 0 (n <= 2) and for F_1.
Product (PROVED from <h^a h^b h^c>_0 = 5 [a+b+c = 3] at v = 1 and homogeneity; h^a = N_a H^a):
  H^a * H^b = N_a N_b N_(a+b) v^((c_a + c_b - c_(a+b))/5) H^(a+b),      eta(H^a, H^b) = 5 N_a N_b [a + b = 3 mod 5];
  in particular H * H^b = beta_b H^(b+1), beta = (1, -1/v, 1, -1, -1): the Gauss-Manin system of UNIFORM_REDUCTION s.1.
Calibration S(z) = 1 + S_1/z + .., z d/dv S = (H *) S, homogeneous:  S_j H^b = s(j,b) v^x(j,b) H^(b+j),
  x(j+1,b) = x(j,b) + 1 - [b+j = 1 mod 5],   s(j+1,b) = s(j,b) betac(b+j)/x(j+1,b)      (b != 1: no logarithm occurs;
  the column b = 1 has S_1 H = -log(v) H^2 and is not used here).
Vanishing-period functionals (the H-components of the DESCENDANT times t_i = [S^-1(z) tbar(z)]_+ at z^i):
  Pi_i(tau_k(H^b)) = PiC(tau_(k-i)(H^b)),    PiC(tau_k(H^b)) = (1/5) eta(S^-1_k H^b, H^2) = (-1)^k (1/5) eta(H^b, S_k H^2)
                   = (-1)^k s(k,2) v^x(k,2) N_b N_(2+k) [b + k = 1 mod 5]:
  PiC = 1 on H, -v on tau_1(1), -v^2/2 on tau_2(H^4), -v^3/6 on tau_3(H^3), v^4/24 on tau_4(H^2), v^4/96 on tau_5(H), ..
  (minus the pure-1/z solution Y^C of the common brief, read at -z).
Subspace K = {tbar : Pi_i(tbar) = 0 for all i >= 1} (descendant times t^H_k = 0, k >= 1); it is spanned by tau_0(H) and
  kappa(X, k) = [S(z) X z^k]_+ = Sum_{j <= k} S_j X z^(k-j),  X in {1, H^2, H^3, H^4}   (descendant insertion tau_k(X)).
"""
import os, sys
from fractions import Fraction as Fr
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ms_rows as mr

C = mr.C_EXP; N = mr.N_SIGN
BETAC = (1, -1, 1, -1, -1)            # beta_b = BETAC[b] v^(-[b = 1])


def prod_coef(a, b):
    """H^a * H^b = coef * v^pw * H^(a+b): returns (coef, pw)"""
    a %= 5; b %= 5; c = (a + b) % 5
    pw = Fr(C[a] + C[b] - C[c], 5); assert pw.denominator == 1
    return N[a] * N[b] * N[c], int(pw)


@lru_cache(maxsize=None)
def s_coef(j, b):
    """S_j H^b = s v^x H^(b+j): returns (s, x); b != 1 (the column H has a logarithm)"""
    b %= 5
    assert b != 1, "column H of S has a logarithm"
    s = Fr(1); x = 0
    for i in range(j):
        idx = (b + i) % 5
        x = x + 1 - (1 if idx == 1 else 0)
        assert x != 0
        s = s * BETAC[idx] / x
    return s, x


def pic(b, k):
    """PiC(tau_k(H^b)) = coefficient * v^power: returns (coef, power); coef = 0 unless b + k = 1 mod 5"""
    b %= 5
    if k < 0 or (b + k) % 5 != 1: return Fr(0), 0
    s, x = s_coef(k, 2)
    return (-1) ** k * s * N[b] * N[(2 + k) % 5], x


def pi_i(i, b, k):
    return pic(b, k - i)


def kappa(b, k):
    """kappa(H^b, k) at v = 1 as a list of (coef, (index, psi-power)) : descendant insertion tau_k(H^b), b != 1"""
    out = []
    for j in range(k + 1):
        s, x = s_coef(j, b)
        out.append((s, ((b + j) % 5, k - j)))
    return out


def weight(b, k):
    """w(b, k) = (4 - c_b - 4k)/5: an insertion tau_k(H^b) lowers the exponent of v by w"""
    return Fr(4 - C[b % 5] - 4 * k, 5)


def charge(b, k):
    return (b + k - 1) % 5


if __name__ == "__main__":
    print("H * H^b:", [prod_coef(1, b) for b in range(5)])
    for b in (0, 2, 3, 4):
        print("S_j H^%d, j = 0..6:" % b, ["%s v^%d H^%d" % (s_coef(j, b)[0], s_coef(j, b)[1], (b + j) % 5) for j in range(7)])
    print("PiC nonzero components (b, k): coef, power of v")
    for k in range(0, 12):
        for b in range(5):
            c, x = pic(b, k)
            if c: print("   tau_%d(H^%d): %s v^%d" % (k, b, c, x))
    # checks: kappa(X, k) is annihilated by Pi_i for i >= 1 (X != H) and Pi_0(kappa(X, k)) = 0 unless ..
    ok = True
    for b in (0, 2, 3, 4):
        for k in range(0, 12):
            for i in range(0, k + 1):
                tot = sum(c * pi_i(i, bb, kk)[0] for c, (bb, kk) in kappa(b, k))
                if tot != 0: ok = False; print("   Pi_%d(kappa(H^%d, %d)) = %s" % (i, b, k, tot))
    print("Pi_i(kappa(X, k)) = 0 for all i >= 0, X in {1, H^2, H^3, H^4}, k < 12:", ok)
