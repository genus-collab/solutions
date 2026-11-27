# -*- coding: utf-8 -*-
"""chk7.py -- referee NOTE: the two special points at b = +zeta^2 with the roots paired the other way round between the two primes."""
import sys, os, json
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from genus3_qzeta_rec import reconstruct, crt, fmt
R1 = FpRing(P1); R2 = FpRing(P2); z = crt(R1.zeta, R2.zeta)
c1 = json.load(open(os.path.join(PROJ, "g3_cubic_61.json"))); c2 = json.load(open(os.path.join(PROJ, "g3_cubic_89.json")))
names = c1['vis_names'] + ['c5_0', 'c5_1', 'c5_2', 'c6_1', 'c6_2']
def norm(e, p): i = pow(e[1], p - 2, p); return [x * i % p for x in e]
for a, b in ((0, 1), (1, 0)):
    v1 = norm(c1['circ_b_extra'][a], P1); v2 = norm(c2['circ_b_extra'][b], P2)
    print("== extra[%d] at 2^61-1 paired with extra[%d] at 2^89-1" % (a, b))
    for j in range(3, 10):
        d, n, bits = reconstruct(v1[j], v2[j], z); ok = bits <= 22
        print("   %-5s %s   [%d bits]" % (names[j], fmt(d, n) if ok else "not reconstructed", bits))
# symmetric functions of the two roots tau = c2_1 at b = +zeta^2 (rational even if the two points are conjugate)
t1 = [norm(c1['circ_b_extra'][k], P1)[3] for k in (0, 1)]; t2 = [norm(c2['circ_b_extra'][k], P2)[3] for k in (0, 1)]
for nm, s1, s2 in (("sum of the two c2_1", sum(t1) % P1, sum(t2) % P2), ("product of the two c2_1", t1[0] * t1[1] % P1, t2[0] * t2[1] % P2)):
    d, n, bits = reconstruct(s1, s2, z); print("   %s: %s   [%d bits]" % (nm, fmt(d, n) if bits <= 22 else "not reconstructed", bits))
