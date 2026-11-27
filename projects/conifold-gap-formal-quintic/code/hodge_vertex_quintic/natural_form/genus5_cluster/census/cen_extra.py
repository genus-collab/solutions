# -*- coding: utf-8 -*-
"""cen_extra.py -- small aggregates over the per-row census of cen_census.py, for the tables of CENSUS.md.

For the 1655 coefficients at stake inside the 627 cancelling rows: how many non-zero pieces contribute to each.  A cell with
exactly one contributing piece would be impossible (one non-zero summand cannot cancel), so the count of such cells is a check.

  py -I -B cen_extra.py <census_core.json> [out.json]
"""
import os, sys, json, collections, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L


def main(path, outp=None):
    core = json.load(open(path))
    rep = {}
    for tag in core["tags"]:
        pr = core["tags"][tag]["per_row"]
        cells = collections.Counter()      # number of contributing non-zero pieces -> number of (row, coefficient) cells
        at_stake = 0
        rowsize = []
        wide = []
        for q, r in pr.items():
            if r["E"] == 0:
                continue
            rowsize.append(r["pieces"])
            for k in range(r["E"], L.K):
                at_stake += 1
                n = r["per_coeff"][k]
                cells[n if n <= 2 else (10 if n <= 100 else 100)] += 1
                if n >= 2:
                    wide.append((n, int(q), k))
        wide.sort(reverse=True)
        rep[tag] = dict(
            cancelling_rows=len(rowsize), coefficients_at_stake_in_those_rows=at_stake,
            cells_with_0_pieces=cells[0], cells_with_1_piece=cells[1], cells_with_2_pieces=cells[2],
            cells_with_3_to_100=cells[10], cells_with_more_than_100=cells[100],
            pieces_per_row=dict(total=sum(rowsize), min=min(rowsize), median=statistics.median(rowsize),
                                mean=round(statistics.mean(rowsize), 1), max=max(rowsize)),
            widest_cells=[[n, q, k] for n, q, k in wide[:10]])
        print(tag, json.dumps(rep[tag], indent=1))
    if outp:
        json.dump(rep, open(outp, "w"), indent=1)
        print("written", outp)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
