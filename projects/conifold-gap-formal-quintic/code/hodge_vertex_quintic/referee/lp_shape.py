# Referee: inspect the ring structure of the genus-2 total (point frame F5) computed with the PROJECT code (read-only use):
# joint degree in (c0..c3, i0..i2) of every numerator monomial (C0^{2g-2} homogeneity), L-range, Lam/i3 absence.
import os
import sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fq_fp import FpRing, P2
from fq_rows import Rows
from fq_graphs import Decorated
from gsum_q import make_frames, genus2_sums, total
from fq_ring import NAMES
R = FpRing(P2); rows = Rows(R, Nz=6); ring = rows.ring; deco = Decorated(2)
t = time.time()
F = make_frames(rows, ['F5'])['F5']
S = genus2_sums(F, deco); T = total(S, ring)
print("F5 total in %.0f s: terms %d, denominators c0^%d W^%d (= c0^%d C1^%d), L-range %s" % (time.time() - t, T.nterms(), T.a, T.b, T.a + 2 * T.b, T.b, T.Ldeg()))
degs = {}
for k, v in T.num.items():
    jd = sum(k[1:9])            # c0..c3, i0..i3
    degs[jd] = degs.get(jd, 0) + 1
print("joint (c,i)-degree of the numerator monomials: %s ; expected 2g-2 + a + 2b = %d" % (degs, 2 + T.a + 2 * T.b))
print("i3 present:", any(k[8] for k in T.num), " Lam present:", any(k[9] for k in T.num))
# which generators appear
used = {n: any(k[j] for k in T.num) for j, n in enumerate(NAMES)}
print("generators appearing:", used)
# degree in i's and in L per monomial: table (Lpow -> max i-degree)
by_i = {}
for k in T.num:
    idg = sum(k[5:9]); by_i[idg] = by_i.get(idg, 0) + 1
print("i-degree distribution:", by_i)
# per coloured graph: denominators and L-range of a few edge graphs
for key in [(4, (0, 0)), (4, (0, 1)), (6, (0, 1)), (3, (0,)), (1, (0,))]:
    el = S[key]; print("graph/colouring", key, "terms", el.nterms(), "den c0^%d W^%d" % (el.a, el.b), "Lrange", el.Ldeg())
