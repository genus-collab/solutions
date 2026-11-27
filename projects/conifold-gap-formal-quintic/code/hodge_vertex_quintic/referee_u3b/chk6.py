# -*- coding: utf-8 -*-
"""chk6.py -- referee NOTE: lattice reconstruction in Q(zeta) (project routine genus3_qzeta_rec.reconstruct, read-only import) of the
isolated solution Delta* (key circ_star) and of the special points circ_b_extra from their residues at the two primes."""
import sys, os, json
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq_fp import FpRing, P1, P2
from genus3_qzeta_rec import reconstruct, crt, fmt
R1 = FpRing(P1); R2 = FpRing(P2); z = crt(R1.zeta, R2.zeta)
c1 = json.load(open(os.path.join(PROJ, "g3_cubic_61.json"))); c2 = json.load(open(os.path.join(PROJ, "g3_cubic_89.json")))
names = c1['vis_names'] + ['c5_0', 'c5_1', 'c5_2', 'c6_1', 'c6_2']
items = [("Delta* (circ_star)", c1['circ_star'], c2['circ_star'])]
for n in range(4):
    e1 = c1['circ_b_extra'][n]; e2 = c2['circ_b_extra'][n]
    i1 = pow(e1[1], P1 - 2, P1); i2 = pow(e2[1], P2 - 2, P2)
    items.append(("circ_b_extra[%d] normalised c1_1 = 1" % n, [x * i1 % P1 for x in e1], [x * i2 % P2 for x in e2]))
for key, v1, v2 in items:
    print("== %s" % key); nrec = 0; tot = 0
    for j in range(15):
        if v1[j] == 0 and v2[j] == 0: print("   %-5s 0" % names[j]); continue
        d, n, bits = reconstruct(v1[j], v2[j], z); ok = bits <= 22; tot += 1; nrec += ok
        print("   %-5s %s   [%d bits: %s]" % (names[j], fmt(d, n) if ok else "not reconstructed", bits, "strong" if ok else "generic length, rejected"))
    print("   reconstructed %d of %d non-zero coordinates" % (nrec, tot))
