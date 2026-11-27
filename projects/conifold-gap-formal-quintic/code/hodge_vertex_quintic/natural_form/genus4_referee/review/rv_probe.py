# -*- coding: utf-8 -*-
"""rv_probe.py -- read-only probe: edge counts of the genus-4 graphs, rows present in the stored nat table, valuations of E_out."""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, PROJ)
from fractions import Fraction as Fr
from fq_fp import P2
from fq4_split import Split
sp = Split(P2, 4, (Fr(0), Fr(0)), 6, 'nat')
Eof = [G['E'] for G in sp.inner_graphs]
print("graphs by E:", {E: Eof.count(E) for E in range(10)})
d = json.load(open(os.path.join(PROJ, "g4_split_g4_89_nat.json")))
rows = sorted(int(k) for k in d['table'])
print("rows present by E:", {E: sum(1 for i in rows if Eof[i] == E) for E in range(10)})
print("valO table (a,b): val")
for a in range(sp.KZ):
    print("  a=%d: " % a + " ".join("%d" % sp.valO[(a, b)] for b in range(sp.KZ - a)))
print("frozen_ok", sp.frozen_ok)
# valuation per colour pair for (0,0)
print("val of EO(0,0) entries:", [[sp.val(x) for x in row] for row in sp.EO[(0, 0)]])
print("val of EO(1,0) entries:", [[sp.val(x) for x in row] for row in sp.EO[(1, 0)]])
