# -*- coding: utf-8 -*-
"""genus3_cross.py -- cross-prime comparison of the genus-3 state files g3_state_89.json (reporting prime 2^89-1) and g3_state_61.json:
every residue of the 2^89-1 state that reconstructs to a rational (bound 1e10, fallback 1e13) is reduced mod 2^61-1 and compared with the
residue at the same position in the 2^61-1 state; structural data (ranks, dimensions, consistency, regular sets) are compared directly.
Run:  py -I genus3_cross.py  -> logs/genus3_cross.log
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_rec import cross_states, Rec

if __name__ == "__main__":
    dA = json.load(open(os.path.join(HERE, "g3_state_89.json"))); dB = json.load(open(os.path.join(HERE, "g3_state_61.json")))
    RA = FpRing(P2); out = []
    def log(s): print(s); out.append(s)
    n, mism, norec, missing, big = cross_states(dA, dB, RA, P1)
    log("cross-prime: %d values reconstructed at 2^89-1 with the default bound 1e10 compared, %d mismatches, %d not reconstructible at that bound (compared by nothing), %d missing in the 2^61-1 state" % (n, len(mism), len(norec), len(missing)))
    nb_ok = sum(1 for x in big if x[2]); log("   1e13-fallback reconstructions: %d, of which %d confirmed at 2^61-1 and %d NOT confirmed (spurious reconstructions; residues are the values of record there)" % (len(big), nb_ok, len(big) - nb_ok))
    for m in mism[:20]: log("   MISMATCH " + str(m))
    for m in missing[:20]: log("   missing " + str(m))
    struct_keys = [k for k in dA if k.startswith("sys4_")] + ["M0_regular", "M0_nonregular", "J_rank", "J_maxdeg"]
    bad = []
    for k in struct_keys:
        a = dA.get(k); b = dB.get(k)
        if isinstance(a, dict):
            for f in ("rank", "consistent", "dim", "names"):
                if a.get(f) != (b or {}).get(f): bad.append((k, f, a.get(f), (b or {}).get(f)))
        elif a != b: bad.append((k, a, b))
    log("structural data (ranks, consistency, dimensions, regular sets, J rank/degree) identical at both primes: %s %s" % (not bad, bad[:5]))
    # the Q-valued headline numbers, reconstructed at 2^89-1 and confirmed above
    rec = Rec(RA)
    log("headline (2^89-1 reconstructions): total F_3 u-coefficients %s; total in t_c %s; g1 graphs 2, 5 = %s, %s; vertex [u^0] = %s" % (
        [rec.fmt(x) for x in dA['total_id']], [rec.fmt(x) for x in dA['total_tc'][:4]], rec.fmt(dA['A_M0']['2'][1]), rec.fmt(dA['A_M0']['5'][1]), rec.fmt(dA['A_M0']['1'][0])))
    log("ALL CROSS-PRIME CHECKS PASSED" if not mism and not missing and not bad else "CROSS-PRIME FAILURES")
    open(os.path.join(HERE, "logs", "genus3_cross.log"), "w", encoding="utf-8").write("\n".join(out) + "\n")
