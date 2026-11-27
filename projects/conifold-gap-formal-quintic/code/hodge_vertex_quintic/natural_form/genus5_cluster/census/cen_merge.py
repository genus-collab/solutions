# -*- coding: utf-8 -*-
"""cen_merge.py -- step 1 of the census: re-merge the cluster part files with my own code and compare every row and every
residue with cluster_genus5_quintic/results/<tag>.json.  Also checks that the five result files agree, index by index, on the
graph list (E, gv, aut), and that the vertex row and the total of each full run are what the note records.

  py -I -B cen_merge.py <out.json>
"""
import os, sys, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L


def compare_tables(mine, theirs):
    """mine: {int q -> [K]}; theirs: {str q -> [K]} from the result file.  Returns a report dict."""
    t = {int(q): list(v) for q, v in theirs.items()}
    only_mine = sorted(set(mine) - set(t))
    only_theirs = sorted(set(t) - set(mine))
    diff = []
    for q in sorted(set(mine) & set(t)):
        if mine[q] != t[q]:
            diff.append(dict(row=q, mine=[str(x) for x in mine[q]], stored=[str(x) for x in t[q]],
                             coeffs=[k for k in range(L.K) if mine[q][k] != t[q][k]]))
    return dict(rows_mine=len(mine), rows_stored=len(t), only_mine=only_mine, only_stored=only_theirs,
                differing_rows=diff, equal=(not only_mine and not only_theirs and not diff))


def main(outpath):
    meta = L.graph_meta()
    rep = dict(tags={}, graph_meta_agreement={}, nin=len(meta))
    base = [(g["E"], list(g["gv"]), g["aut"]) for g in meta]
    for tag in L.TAGS + L.IDTAGS:
        r = L.load_result(tag)
        o = [(g["E"], list(g["gv"]), g["aut"]) for g in r["outer"]]
        rep["graph_meta_agreement"][tag] = dict(n=len(o), equal=(o == base))
    for tag in L.TAGS + L.IDTAGS:
        p = L.prime_of(tag)
        A = L.load_all_pieces(tag)
        mine = L.merge_pieces(A["pieces"], p)
        res = L.load_result(tag)
        cmp = compare_tables(mine, res["table"])
        nz_pieces = sum(len(v) for v in A["pieces"].values())
        cpu = sum(A["times"].values())
        slowest = sorted(A["times"].items(), key=lambda kv: -kv[1])[:5]
        tot = [0] * L.K
        for q, row in mine.items():
            for k in range(L.K):
                tot[k] = (tot[k] + row[k]) % p
        stored_tot = res["info"].get("total")
        rep["tags"][tag] = dict(
            prime=("2^89-1" if p == L.PRIME["89"] else "2^61-1"), part_files=len(A["per_file"]),
            rows_nonzero=sorted(q for q, r in mine.items() if any(r)) if sum(1 for r in mine.values() if any(r)) <= 40
            else sum(1 for r in mine.values() if any(r)),
            n_rows_nonzero=sum(1 for r in mine.values() if any(r)),
            n_rows_allzero=sum(1 for r in mine.values() if not any(r)),
            inner_covered=len(A["pieces"]), missing=sorted(set(range(L.NIN)) - set(A["pieces"])),
            duplicates=A["dup"][:20], n_duplicates=len(A["dup"]), conflicts=A["conflict"],
            nonzero_pieces=nz_pieces, cpu_seconds=cpu, cpu_hours=cpu / 3600.0,
            slowest_inner=[[t, round(s, 1)] for t, s in slowest],
            table_comparison=cmp, failed_in_result=res["failed"],
            total_recomputed=[str(x) for x in tot],
            total_matches_stored_info=(stored_tot is None or [int(x) for x in stored_tot] == tot),
            vertex_row=[str(x) for x in mine.get(0, [0] * L.K)],
            vertex_row_matches=(0 in mine and [int(x) for x in res["table"]["0"]] == mine[0]),
        )
        print("%-24s files=%3d inner=%4d missing=%d dup=%d rows mine/stored=%d/%d nz=%d equal=%s pieces=%d CPU=%.1f h"
              % (tag, len(A["per_file"]), len(A["pieces"]), len(rep["tags"][tag]["missing"]), len(A["dup"]),
                 cmp["rows_mine"], cmp["rows_stored"], rep["tags"][tag]["n_rows_nonzero"], cmp["equal"],
                 nz_pieces, cpu / 3600.0))
        del A, mine
    ok = all(v["table_comparison"]["equal"] and not v["missing"] and not v["conflicts"] for v in rep["tags"].values())
    rep["all_tables_reproduced"] = ok
    with open(outpath, "w") as f:
        json.dump(rep, f, indent=1)
    print("all tables reproduced:", ok)
    print("written", outpath)


if __name__ == "__main__":
    main(sys.argv[1])
