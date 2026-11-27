# -*- coding: utf-8 -*-
"""rv_cmp.py -- the table of the full baseline run of the COPY (out/g4mut_baseline_nat_89_full.json) against the table of record
g4_split_g4_89_nat.json (project folder, read only): same set of rows, same residues."""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
a = json.load(open(os.path.join(HERE, "out", "g4mut_baseline_nat_89_full.json")))['table']
b = json.load(open(os.path.join(PROJ, "g4_split_g4_89_nat.json")))['table']
p = (1 << 89) - 1
print("rows: copy %d, record %d; same keys: %s; all residues equal: %s" % (len(a), len(b), set(a) == set(b), all([x % p for x in a[k]] == [x % p for x in b[k]] for k in a if k in b)))
