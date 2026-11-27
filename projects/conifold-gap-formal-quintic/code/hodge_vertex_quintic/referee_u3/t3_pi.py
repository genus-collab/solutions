# -*- coding: utf-8 -*-
"""tests of the characterisation Pi = Ad_{M_0^-1}( zhat D_1 ), D_1 symmetric covariant with D_1 phi_3 = 0, and of rival hypotheses"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import random
from rlib import *
from fq3_rec import Rec
from fq3_quad import rank_fp
tag = sys.argv[1] if len(sys.argv) > 1 else "61"
S, d, cs = load(tag); p = S.p; rec = Rec(S.R); rng = random.Random(5)
M = d['M']; Mi = d['Minv']; mu0 = d['mu0']; names = S.names
A = cs['Pi_A']; B = cs['Pi_B']
def order1(al, be, ga):
    """zhat D_1 with hat c_1(3) = al, hat c_1(2) = hat c_1(4) = be, hat c_1(0) = hat c_1(1) = ga"""
    return from_fourier(S, {(1, 3): al, (1, 2): be, (1, 4): be, (1, 0): ga, (1, 1): ga})
def res_of(Delta):
    c = circ_of(S, Delta)
    if c is None: return None, None
    r = S.residuals(c); return c, [k for k, v in r.items() if v]
ok = lambda s, v, det="": print(("PASS " if v else "FAIL ") + s + ((": " + str(det)) if det else ""))
# (1) pure order one, no (H, H) component
bad = []
for t in range(6):
    be, ga = rng.randrange(p), rng.randrange(p)
    if t == 0: be, ga = 1, 0
    if t == 1: be, ga = 0, 1
    D = order1(0, be, ga); assert in_sp(D, p)
    Delta = Ad(Mi, M, D, p); c, nz = res_of(Delta)
    if nz: bad.append((t, nz))
ok("(1) Delta = M0^-1 (zhat D_1) M0 with D_1 covariant symmetric, hat c_1(3) = 0 (2 parameters): ALL 18 conditions vanish, 6 points (no order-6 correction needed)", not bad, bad)
# (2) the stored A, B: Ad_{M0} is pure order 1 up to order 6; the order-6 remainder
for nm, v in (("A", A), ("B", B)):
    DA = Ad(M, Mi, ser_from_circ(S, v), p); F, cov = fourier(S, DA)
    mid = all(F[(k, m)] == 0 for k in range(2, 6) for m in range(5))
    D1 = zeros(); D1[1] = DA[1]; rem = add(ser_from_circ(S, v), Ad(Mi, M, D1, p), p, p - 1)
    cr = circ_of(S, rem); r = S.residuals(cr)
    ok("(2) %s: Ad_{M0}(%s) has no components of orders 2..5 and hat c_1(3) = 0; %s - M0^-1 (zhat D_1) M0 is a pure order-6 element (coords %s) killed by the genus-3 rows" % (nm, nm, nm, [names[j] for j in range(15) if cr[j]]),
       mid and F[(1, 3)] == 0 and all(cr[j] == 0 for j in range(13)) and r['J3_0'] == 0 and r['J3_1'] == 0)
# span of Pi (mod order 6) = image of the 2-parameter family
V = [circ_of(S, Ad(Mi, M, order1(0, 1, 0), p))[:13], circ_of(S, Ad(Mi, M, order1(0, 0, 1), p))[:13], A[:13], B[:13]]
ok("(2) span(A, B) = { M0^-1 zhat D_1 M0 } modulo order 6: rank of the four 13-vectors", rank_fp(V, p) == 2, rank_fp(V, p))
# (3) with an (H, H) component
c, nz = res_of(Ad(Mi, M, order1(1, 0, 0), p)); print("   (3) D_1 = pure (H, H) propagator (hat c_1(3) = 1): non-zero conditions: %s" % nz)
c, nz = res_of(Ad(Mi, M, order1(1, rng.randrange(p), rng.randrange(p)), p)); print("   (3) generic covariant D_1 (hat c_1(3) = 1): non-zero conditions: %s" % nz)
# (4) pure single orders k >= 2 in the left placement: do any vanish on everything?
for k in range(2, 7):
    F0 = {}
    ms = {2: [(0, 2, -1), (3, 4, -1)], 3: [(0, 3, 1), (1, 2, 1), (4, 4, 1)], 4: [(0, 4, -1), (1, 3, -1)], 5: [(0, 0, 1), (1, 4, 1), (2, 3, 1)], 6: [(0, 1, -1), (2, 4, -1)]}[k]
    for (m1, m2, sg) in ms:
        Fk = {(k, m1): 1, (k, m2): sg % p} if m1 != m2 else {(k, m1): 1}
        D = from_fourier(S, Fk)
        if not in_sp(D, p): print("   not in sp", k, m1, m2); continue
        c, nz = res_of(Ad(Mi, M, D, p))
        print("   (4) left-placed pure order %d, Fourier pair (%d, %d) [H^%d -> H^%d]: non-zero conditions: %s" % (k, m1, m2, (m1 - 1) % 5, (m1 - 1 - k) % 5, nz))
# (5) rival hypotheses: membership in Pi (mod order 6)  <=>  Ad_{M0} X = zhat D_1 + O(zhat^6), hat c_1(3) = 0
def in_pi(X):
    F, cov = fourier(S, Ad(M, Mi, X, p))
    return cov and F[(1, 3)] == 0 and all(F[(k, m)] == 0 for k in range(2, 6) for m in range(5))
def resid(X):
    c = circ_of(S, X)
    if c is None: return "not circulant"
    return [k for k, v in S.residuals(c).items() if v]
zM = [[[k * x % p for x in r] for r in M[k]] for k in range(KMAX + 1)]          # zhat dM0/dzhat
XR = mul(Mi, zM, p); XL = mul(zM, Mi, p)
cands = [("zhat-rescaling, right: M0^-1 zhat M0'", XR), ("zhat-rescaling, left-placed generator zhat M0' M0^-1 used as Delta", XL), ("mu0 = log M0", mu0), ("N1taut (u-rescaling u -> lambda u)", sp_mats(S.nt, p))]
for k in range(1, 4):
    g = zeros(); g[k] = mu0[k]; cands.append(("graded piece (mu0)_%d zhat^%d" % (k, k), g))
g = zeros(); g[1] = M[1]; cands.append(("(M0)_1 zhat (order-1 part of M0, = (mu0)_1 zhat)", g))
cands.append(("Ad_{M0^-1} of (mu0)_1 zhat", Ad(Mi, M, g, p)))
g = zeros(); g[1] = sp_mats(S.nt, p)[1]; cands.append(("Ad_{M0^-1} of (N1taut)_1 zhat", Ad(Mi, M, g, p)))
for nm, X in cands:
    F, cov = fourier(S, X)
    print("   (5) %-70s in sp %s, covariant %s, hat c_1(3) = %s, in Pi: %s; failing conditions as Delta: %s" % (nm, in_sp(X, p), cov, rec.fmt(F[(1, 3)]), in_pi(X), resid(X)))
# (6) vectors annihilated: the unit phi_1, phi_3 = H^2, in both placements
for nm, v in (("A", A), ("B", B)):
    X = ser_from_circ(S, v); F, _ = fourier(S, X); FD, _ = fourier(S, Ad(M, Mi, X, p))
    print("   (6) %s right-placed: orders k with hat c_k(1) = 0 (unit killed): %s; with hat c_k(3) = 0 (H^2 killed): %s | left-placed (orders 1..5): unit %s, H^2 %s" % (nm,
          [k for k in range(1, 7) if F[(k, 1)] == 0], [k for k in range(1, 7) if F[(k, 3)] == 0], [k for k in range(1, 6) if FD[(k, 1)] == 0], [k for k in range(1, 6) if FD[(k, 3)] == 0]))
