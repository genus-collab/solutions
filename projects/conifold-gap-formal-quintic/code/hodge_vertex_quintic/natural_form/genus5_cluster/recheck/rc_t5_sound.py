# -*- coding: utf-8 -*-
"""rc_t5_sound.py -- part (4) of the re-check: what a missing inner graph could hide.

The natural run's TOTAL is the sum over all inner graphs t and all outer-edge subsets S of the pieces.
If per-inner-graph additivity holds, the total equals the sum of the 4555 point-frame terms, which the
shipped parts/g5_<p>_id/local.json give graph by graph.  So an inner graph entirely absent from the
cluster part files would change the total unless its point-frame term vanishes at both primes.  This
script counts how many point-frame terms vanish, i.e. how large the blind spot of the total check is.
"""
import os
import sys
import json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rc_rec

PKG = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "cluster_genus5_quintic"))
K = 8


def idterms(pr):
    with open(os.path.join(PKG, "parts", "g5_%s_id" % pr, "local.json")) as f:
        pieces = json.load(f)['pieces']
    out = {}
    for t, r in pieces.items():
        out[int(t)] = r.get(t, [0] * K)
    return out


def main():
    with open(os.path.join(PKG, "results", "g5_89_nat.json")) as f:
        nat = json.load(f)
    out = nat['outer']
    a = idterms("89")
    b = idterms("61")
    z89 = sorted(t for t, v in a.items() if not any(v))
    z61 = sorted(t for t, v in b.items() if not any(v))
    zboth = sorted(set(z89) & set(z61))
    print("point-frame terms of the 4555 inner graphs (from the shipped parts)")
    print("   zero at p89: %d; zero at p61: %d; zero at both: %d" % (len(z89), len(z61), len(zboth)))
    print("   non-zero at both primes: %d" % (4555 - len(set(z89) | set(z61))))
    import collections
    cE = collections.Counter(out[t]['E'] for t in zboth)
    print("   the %d inner graphs with a vanishing point-frame term, by number of edges: %s"
          % (len(zboth), sorted(cE.items())))
    print()
    print("   -> an inner graph whose point-frame term is non-zero at either prime cannot be missing")
    print("      from the part files without breaking 'natural total == point-frame total', which this")
    print("      re-check verified independently at both primes.  The blind spot of that argument is the")
    print("      %d graphs with a term that vanishes at both primes." % len(zboth))
    print()
    v89 = [sum(a[t][k] for t in a) % rc_rec.P89 for k in range(K)]
    v61 = [sum(b[t][k] for t in b) % rc_rec.P61 for k in range(K)]
    print("   sum of the 4555 point-frame terms, recomputed here:")
    print("      p89 == natural total:", v89 == nat['info']['total'])
    with open(os.path.join(PKG, "results", "g5_61_nat.json")) as f:
        n61 = json.load(f)
    print("      p61 == natural total:", v61 == n61['info']['total'])


if __name__ == "__main__":
    main()
