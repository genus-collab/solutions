# -*- coding: utf-8 -*-
"""rc_t3_census.py -- part (3) of the re-check: regularity exactly as the stored tables state it,
the census of outer graphs (by number of edges and by genus vector) and the genus-4 analogue.
"""
import os
import sys
import json
import collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

RES = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..",
                                   "cluster_genus5_quintic", "results"))
HVQ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
NAT = ["g5_89_nat", "g5_61_nat", "g5_89_nat_ab1o3_m2o7"]


def load(path):
    with open(path) as f:
        return json.load(f)


def census(out, tab, cap, K, label):
    nE = collections.Counter(g['E'] for g in out)
    stored = set(int(k) for k in tab)
    sE = collections.Counter(out[i]['E'] for i in stored)
    nzE = collections.Counter(out[i]['E'] for i in stored if any(tab[str(i)]))
    print("--- %s" % label)
    print("    graphs: %d; K = %d coefficients per row; cap on outer edges = %d" % (len(out), K, cap))
    print("     E   graphs   stored rows   non-zero rows   no row (E <= cap)")
    for E in sorted(nE):
        nor = (nE[E] - sE[E]) if E <= cap else "-"
        print("    %2d   %6d   %11d   %13d   %s" % (E, nE[E], sE[E], nzE[E], nor))
    tot_le = sum(v for k, v in nE.items() if k <= cap)
    print("    graphs with E <= %d: %d (of which the vertex graph 1, edge graphs %d)" % (cap, tot_le, tot_le - 1))
    print("    graphs with E >  %d: %d  -- O(u^%d) for free: every outer edge carries one power of u"
          % (cap, len(out) - tot_le, K))
    print("    stored rows: %d; non-zero rows: %d  -> %s"
          % (len(stored), sum(nzE.values()), sorted(i for i in stored if any(tab[str(i)]))))
    print("    edge graphs (E >= 1) with a stored row: %d; identically zero among them: %d"
          % (sum(v for k, v in sE.items() if k >= 1),
             sum(1 for i in stored if out[i]['E'] >= 1 and not any(tab[str(i)]))))
    print("    edge graphs with E <= %d and NO stored row (no inner graph gave a non-zero piece): %d"
          % (cap, tot_le - 1 - sum(v for k, v in sE.items() if k >= 1)))
    return nE, sE


def gv_census(out, tab, cap):
    stored = set(int(k) for k in tab)
    byn = collections.Counter()
    bynS = collections.Counter()
    zero_gv = collections.Counter()
    zero_gvS = collections.Counter()
    for i, g in enumerate(out):
        if g['E'] > cap:
            continue
        n = len(g['gv'])
        byn[n] += 1
        z = sum(1 for x in g['gv'] if x == 0)
        zero_gv[z] += 1
        if i in stored:
            bynS[n] += 1
            zero_gvS[z] += 1
    print("    distribution over the number of vertices V (graphs with E <= %d):" % cap)
    print("      V         :", " ".join("%5d" % v for v in sorted(byn)))
    print("      graphs    :", " ".join("%5d" % byn[v] for v in sorted(byn)))
    print("      stored row:", " ".join("%5d" % bynS[v] for v in sorted(byn)))
    print("    distribution over the number of genus-0 vertices:")
    print("      #(g_v=0)  :", " ".join("%5d" % v for v in sorted(zero_gv)))
    print("      graphs    :", " ".join("%5d" % zero_gv[v] for v in sorted(zero_gv)))
    print("      stored row:", " ".join("%5d" % zero_gvS[v] for v in sorted(zero_gv)))


def main():
    K = 8
    cap = K - 1
    base = None
    for tag in NAT:
        d = load(os.path.join(RES, tag + ".json"))
        out = d['outer']
        tab = d['table']
        print()
        nE, sE = census(out, tab, cap, K, "genus 5, %s   failed=%s   nonregular=%s"
                        % (tag, d['failed'], d['info']['nonregular']))
        bad = [i for i in (int(k) for k in tab) if out[i]['E'] >= 1 and any(tab[str(i)])]
        print("    CHECK every stored row with E >= 1 is identically zero:", not bad, bad[:10])
        if base is None:
            base = (sorted(int(k) for k in tab), nE, sE)
        else:
            print("    same set of stored outer graph indices as %s: %s" % (NAT[0], sorted(int(k) for k in tab) == base[0]))
    print()
    gv_census(load(os.path.join(RES, NAT[0] + ".json"))['outer'],
              load(os.path.join(RES, NAT[0] + ".json"))['table'], cap)
    print()
    g4 = os.path.join(HVQ, "g4_split_g4_89_nat.json")
    if os.path.exists(g4):
        d = load(g4)
        census(d['outer'], d['table'], 5, 6, "genus 4, %s (for comparison)   failed=%s" % (d['tag'], d['failed']))
        bad = [i for i in (int(k) for k in d['table']) if d['outer'][i]['E'] >= 1 and any(d['table'][str(i)])]
        print("    CHECK every stored row with E >= 1 is identically zero:", not bad, bad[:10])
    else:
        print("genus-4 reference %s not found" % g4)


if __name__ == "__main__":
    main()
