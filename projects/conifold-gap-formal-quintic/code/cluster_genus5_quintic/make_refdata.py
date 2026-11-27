# -*- coding: utf-8 -*-
"""make_refdata.py -- run ONCE on the workstation (needs the project folder hodge_vertex_quintic next to this package); it is shipped only
as a record of how data/ref_*.json were made.

Writes, keyed by the FAST canonical key of the outer graph (json.dumps of [gv, edges], see key_str):
  data/graphs_2_0_3.json, data/graphs_3_0_6.json        fast graph lists (genus 4, 5 lists are the ones of the feasibility study)
  data/ref_g23.json        genus 3: stored per-graph tables of the OLD vertex-tensor code (g3_state_<p>.json: A_id, A_M0 through u^5;
                           g3_nat_<p>_nat.json: A3 through u^3), genus 2: values of OU_EXTENSION_genus2/3.md (nat, M0), both primes
  data/ref_g4.json         genus 4: tables of record g4_split_g4_<p>_<variant>.json (variants nat, id, M0 (89 only), nat_t5 (89 only),
                           nat in frame (1/3,-2/7)), through u^5
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.join(os.path.dirname(HERE), "hodge_vertex_quintic")
sys.path.insert(0, HERE)
from g5_graphs import canon_fast, graphs, levels, DATA
from fq_fp import P1, P2


def key_str(gv, ed):
    k = canon_fast(tuple(gv), [tuple(e) for e in ed], ())
    return json.dumps([list(k[0]), [list(e) for e in k[1]]])


def main():
    for g, E in ((2, 3), (3, 6)):
        graphs(g, E)
    sys.path.insert(1, PROJ)
    from fq3_graphs import outer_graphs
    out = {}
    for prime, p in (("89", P2), ("61", P1)):
        g3 = json.load(open(os.path.join(PROJ, "g3_state_%s.json" % prime))); nat = json.load(open(os.path.join(PROJ, "g3_nat_%s_nat.json" % prime)))
        for variant, ref, kcmp in (("id", g3['A_id'], 6), ("M0", g3['A_M0'], 6), ("nat", nat['A3'], 4)):
            out["g3_%s_%s" % (prime, variant)] = dict(K=kcmp, table={key_str(G['gv'], G['ed']): [x % p for x in ref[str(G['idx'])][:kcmp]] for G in outer_graphs(3)})
        for variant, exp in (("nat", {1: [Fr(1, 1200), Fr(-7, 6000), Fr(287, 864000)], 2: [0, 0, Fr(479, 756000)], 3: [0, 0, Fr(1, 2250)],
                                       4: [0, 0, Fr(-127, 2016000)], 5: [0, 0, 0], 6: [0, 0, 0], 7: [0, 0, 0]}),
                             ("M0", {2: [0, Fr(1, 750)], 4: [0, Fr(-13, 8000)], 3: [0, 0], 5: [0, 0], 6: [0, 0], 7: [0, 0]})):
            tab = {}
            for G in outer_graphs(2):
                if G['idx'] in exp:
                    tab[key_str(G['gv'], G['ed'])] = [Fr(x).numerator * pow(Fr(x).denominator, -1, p) % p for x in exp[G['idx']]]
            out["g2_%s_%s" % (prime, variant)] = dict(K=len(next(iter(exp.values()))), table=tab, partial=True)
    json.dump(out, open(os.path.join(DATA, "ref_g23.json"), "w"))
    g4 = {}
    for name in ("g4_89_nat", "g4_61_nat", "g4_89_id", "g4_61_id", "g4_89_M0", "g4_89_nat_t5", "g4_89_nat_ab1o3_m2o7", "g4_61_nat_ab1o3_m2o7"):
        d = json.load(open(os.path.join(PROJ, "g4_split_%s.json" % name)))
        tab = {key_str(o['gv'], o['ed']): d['table'].get(str(i), [0] * 6) for i, o in enumerate(d['outer'])}
        assert len(tab) == 379
        g4[name] = dict(K=6, table=tab, total=d['total'])
    json.dump(g4, open(os.path.join(DATA, "ref_g4.json"), "w"))
    print("written:", sorted(out), sorted(g4))


if __name__ == "__main__":
    main()
