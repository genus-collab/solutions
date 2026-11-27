# -*- coding: utf-8 -*-
"""c9_eout.py -- is E_out really non-degenerate?  The collect check "every entry of E_out = E_R - E_M is O(u)"
is  min(valO) >= 1  and would PASS if E_out were identically zero (val(None) is the full precision), in which
case every outer graph with edges would be regular for free and the vertex graph would carry everything, i.e.
exactly the claimed output.  Here the valuation profile of E_out is printed explicitly for the natural M.
Run:  py -I -B c9_eout.py
"""
import os
import sys
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
from g5_engine import make_split
from fq_fp import P2

for variant in ("nat", "M0"):
    sp = make_split(P2, 5, ("0", "0"), 8, variant)
    tot = 0
    nz = 0
    allnone = 0
    vals = {}
    for key, Eo in sp.EO.items():
        n = 0
        for row in Eo:
            for x in row:
                tot += 1
                if x is not None and any(x):
                    nz += 1
                    n += 1
        if n == 0:
            allnone += 1
        vals[sp.valO[key]] = vals.get(sp.valO[key], 0) + 1
    print("variant %-4s: %d (a,b) keys; %d of %d colour entries of E_out are non-zero; %d keys are entirely zero"
          % (variant, len(sp.EO), nz, tot, allnone))
    print("              valuation of E_out per key:", dict(sorted(vals.items())),
          "(8 = the entry is zero to the stored precision)")
    print("              min over keys =", min(sp.valO.values()), "-- the check only asks for >= 1")
