# -*- coding: utf-8 -*-
"""cen_ab_add.py -- the per-inner-graph additivity of the SECOND FRAME run, which the cluster could not check.

The collect step of g5_89_nat_ab1o3_m2o7 skipped additivity because no point-frame ('id') run existed in the frame
(alpha, beta) = (1/3, -2/7): parts/g5_89_id_ab1o3_m2o7/ held only a stub log.  That run was made here on 2026-11-13
(genus5_run.py --variant id --ab 1/3,-2/7 --run --nproc 20; 4555 inner graphs, wall 650 s, CPU 2.0 h, Python 3.14.6), so the
4555 identities can now be checked.  Same engine as the cluster (g5_engine), but the other code path (graph_id, not graph_dp).

  py -I -B cen_ab_add.py <out.json>
"""
import os, sys, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L

AB = "g5_89_nat_ab1o3_m2o7"
IDAB = "g5_89_id_ab1o3_m2o7"


def main(outp):
    p = L.PRIME["89"]
    with open(os.path.join(L.PKG, "parts", IDAB, "local.json")) as f:
        d = json.load(f)
    assert d["tag"] == IDAB
    idt = {int(t): v[str(t)] for t, v in d["pieces"].items()}
    A = L.load_all_pieces(AB)
    bad = []
    nz_id = 0
    for t in range(L.NIN):
        s = [0] * L.K
        for q, c in A["pieces"][t].items():
            for k in range(L.K):
                s[k] = (s[k] + c[k]) % p
        ref = [x % p for x in idt[t]]
        if any(ref):
            nz_id += 1
        if s != ref:
            bad.append(t)
    rep = dict(frame="(alpha,beta) = (1/3,-2/7)", prime="2^89-1", nat=AB, id=IDAB,
               id_run=dict(made="2026-11-13 on the workstation, genus5_run.py --variant id --ab 1/3,-2/7 --run --nproc 20",
                           inner_graphs=len(idt), cpu_hours=round(sum(d["times"].values()) / 3600.0, 2),
                           wall_seconds=650, python="3.14.6",
                           nonzero_point_frame_terms=nz_id),
               inner_graphs_checked=L.NIN, failures=len(bad), failing=bad[:20], all_pass=(not bad))
    json.dump(rep, open(outp, "w"), indent=1)
    print(json.dumps({k: v for k, v in rep.items() if k != "failing"}, indent=1))
    print("written", outp)


if __name__ == "__main__":
    main(sys.argv[1])
