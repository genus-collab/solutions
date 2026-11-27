# -*- coding: utf-8 -*-
"""f2_add.py -- PER-INNER-GRAPH ADDITIVITY in the second conifold frame (alpha, beta) = (1/3, -2/7), p = 2^89-1.

For every one of the 4555 leg-free stable genus-5 graphs t:
    sum over the outer graphs q of the stored piece (t -> q) of the natural run   ==   the point-frame term of t.
Left side: the cluster's parts/g5_89_nat_ab1o3_m2o7/task*of100.json (copied back, read-only).
Right side: the local point-frame run parts/g5_89_id_ab1o3_m2o7 (variant 'id', same prime, same frame).
The summation is this file's own code; g5_checks.checks is not called.

usage:  py -I -B f2_add.py <nat parts dir> <id parts dir> [--prime 89] [--json OUT]
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_lib import K, NIN, AB, prime, read_parts_dir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("natdir")
    ap.add_argument("iddir")
    ap.add_argument("--prime", default="89")
    ap.add_argument("--json", default="")
    ap.add_argument("--natTag", default="")
    ap.add_argument("--idTag", default="")
    a = ap.parse_args()
    p = prime(a.prime)
    nat_tag = a.natTag or "g5_%s_nat_ab1o3_m2o7" % a.prime
    id_tag = a.idTag or "g5_%s_id_ab1o3_m2o7" % a.prime

    nnames, npieces, ntimes = read_parts_dir(a.natdir, nat_tag)
    inames, ipieces, itimes = read_parts_dir(a.iddir, id_tag)
    print("natural run  %s: %d part files, %d inner graphs, CPU %.1f h" % (nat_tag, len(nnames), len(npieces), sum(ntimes.values()) / 3600))
    print("point frame  %s: %d part files, %d inner graphs, CPU %.1f h" % (id_tag, len(inames), len(ipieces), sum(itimes.values()) / 3600))

    out = {"nat_tag": nat_tag, "id_tag": id_tag, "ab": AB, "prime": a.prime}
    miss_nat = sorted(set(range(NIN)) - set(int(t) for t in npieces))
    miss_id = sorted(set(range(NIN)) - set(int(t) for t in ipieces))
    out["missing_nat"] = miss_nat
    out["missing_id"] = miss_id
    print("coverage: natural %d/%d (missing %s), point frame %d/%d (missing %s)"
          % (len(npieces), NIN, miss_nat[:10], len(ipieces), NIN, miss_id[:10]))

    # the 'id' run must store exactly the diagonal piece q = t and nothing else
    offdiag = sorted(int(t) for t, res in ipieces.items() if any(q != t for q in res))
    out["id_offdiagonal"] = offdiag
    print("point-frame part files storing a piece q != t: %d %s" % (len(offdiag), offdiag[:10]))

    bad = []
    checked = 0
    nzero_id = 0
    npieces_tot = 0
    for t in range(NIN):
        res = npieces.get(str(t), {})
        s = [0] * K
        for q, c in res.items():
            assert len(c) == K, "piece (%d -> %s) has %d residues" % (t, q, len(c))
            for k in range(K):
                s[k] = (s[k] + c[k]) % p
        npieces_tot += len(res)
        ref = ipieces.get(str(t), {}).get(str(t), [0] * K)
        ref = [x % p for x in ref]
        if not any(ref):
            nzero_id += 1
        if s != ref:
            bad.append({"t": t, "sum_of_pieces": s, "point_frame": ref})
        checked += 1
    out["identities_checked"] = checked
    out["failures"] = bad
    out["stored_pieces_nat"] = npieces_tot
    out["inner_graphs_with_zero_point_frame_term"] = nzero_id
    print("stored (t -> q) pieces in the natural run: %d" % npieces_tot)
    print("inner graphs whose point-frame term is identically zero: %d" % nzero_id)
    print("ADDITIVITY: %d / %d identities verified; %d failures%s"
          % (checked - len(bad), checked, len(bad), "" if not bad else ": " + json.dumps(bad[:5])))

    # global corollary, checked separately against the cluster's own merged table
    tn = [0] * K
    ti = [0] * K
    for t, res in npieces.items():
        for q, c in res.items():
            for k in range(K):
                tn[k] = (tn[k] + c[k]) % p
    for t, res in ipieces.items():
        for q, c in res.items():
            for k in range(K):
                ti[k] = (ti[k] + c[k]) % p
    out["total_nat"] = tn
    out["total_id"] = ti
    out["totals_agree"] = (tn == ti)
    print("total of the natural pieces == total of the point-frame terms: %s" % (tn == ti))
    rp = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))), "cluster_genus5_quintic", "results", "%s.json" % nat_tag)
    if os.path.exists(rp):
        with open(rp) as f:
            info = json.load(f)["info"]
        out["total_matches_cluster_result"] = (info.get("total") == tn)
        print("total matches info['total'] of the cluster's results/%s.json: %s" % (nat_tag, info.get("total") == tn))
    print("RESULT: " + ("PASS" if not bad and not miss_nat and not miss_id and not offdiag else "FAIL"))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1)
        print("written %s" % a.json)
    return 0 if not bad and not miss_nat and not miss_id and not offdiag else 1


if __name__ == "__main__":
    sys.exit(main())
