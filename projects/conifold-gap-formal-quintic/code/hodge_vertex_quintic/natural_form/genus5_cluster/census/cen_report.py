# -*- coding: utf-8 -*-
"""cen_report.py -- assemble census_g5.json (machine-readable summary of the whole census) and print the tables that go into
CENSUS.md.  Inputs: the JSON files written by cen_merge.py, cen_prov.py and cen_census.py (two runs of the latter).

  py -I -B cen_report.py <merge.json> <prov.json> <census_core.json> <census_ab.json> <out census_g5.json> [<ab_add.json>]
"""
import os, sys, json, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L

G4 = dict(genus=4, graphs=379, K=6, cap=5,
          E_hist={0: 1, 1: 3, 2: 7, 3: 21, 4: 43, 5: 75, 6: 89, 7: 81, 8: 42, 9: 17},
          stored_rows_by_E={1: 3, 2: 7, 3: 18, 4: 24, 5: 16}, nonzero_pieces=6767, cpu_minutes=54,
          source="NATURAL_FORM.md s.5 (6767 pieces, 229/81/68) and logs/validate.log of the genus-5 package, where the "
                 "genus-4 run of record is re-checked outer graph by outer graph")


def per_E_table(census, Eof, K):
    rows = census["per_row"]
    byE = collections.Counter()
    pieces_byE = collections.Counter()
    for q, r in rows.items():
        if r["E"] >= 1:
            byE[r["E"]] += 1
            pieces_byE[r["E"]] += r["pieces"]
    hist = collections.Counter(Eof)
    out = []
    for e in range(1, K):
        tot = hist[e]
        got = byE[e]
        out.append(dict(E=e, graphs=tot, rows_with_pieces=got, no_piece=tot - got,
                        coeffs_at_stake=(K - e) * tot, coeffs_in_cancelling_rows=(K - e) * got,
                        pieces=pieces_byE[e]))
    return out


def main(mp, pp, cp, abp, outp, abadd=None):
    merge = json.load(open(mp))
    prov = json.load(open(pp))
    core = json.load(open(cp))
    ab = json.load(open(abp))
    Eof = [g["E"] for g in L.graph_meta()]
    A = core["tags"]["g5_89_nat"]
    B = core["tags"]["g5_61_nat"]
    C = ab["tags"]["g5_89_nat_ab1o3_m2o7"]
    rep = dict(
        what="census and audit of the genus-5 cluster evidence (formal quintic, natural factorisation, edge splitting)",
        date="2026-11-13",
        source=dict(tarball="formalquintic/g5_evidence.tar.gz",
                    extracted=L.EV, results=os.path.join(L.PKG, "results"),
                    code=os.path.dirname(os.path.abspath(__file__))),
        step1_remerge=dict(
            statement="every row and every residue of all seven result files reproduced from the part files by independent code",
            all_tables_reproduced=merge["all_tables_reproduced"],
            graph_meta_agreement=merge["graph_meta_agreement"],
            per_tag={t: dict(part_files=v["part_files"], inner_covered=v["inner_covered"], missing=len(v["missing"]),
                             duplicate_inner=v["n_duplicates"], conflicts=v["conflicts"],
                             stored_rows=v["table_comparison"]["rows_stored"], rows_nonzero=v["n_rows_nonzero"],
                             rows_allzero=v["n_rows_allzero"], nonzero_pieces=v["nonzero_pieces"],
                             cpu_hours=round(v["cpu_hours"], 2), tables_equal=v["table_comparison"]["equal"],
                             differing_rows=v["table_comparison"]["differing_rows"],
                             total_matches_stored_info=v["total_matches_stored_info"], failed_in_result=v["failed_in_result"])
                     for t, v in merge["tags"].items()}),
        step2_provenance=dict(
            statement="task decomposition recomputed from data/cost_table_g5.json with independent LPT code; logs read for "
                      "hosts, Python, wall times, restarts and SLURM ids",
            validate=prov["validate"], slurm=prov["slurm"],
            per_run={t: dict(mode=v["mode"], variant=v["variant"], cost_column=v.get("cost_column"),
                             ntasks=v.get("ntasks"), assignment=v["assignment_check"], coverage=v["coverage"],
                             cpu_hours=round(v["cpu_hours"], 2), logs=v.get("logs"))
                     for t, v in prov["runs"].items()}),
        step3_census=dict(
            statement="cancellation census from the cluster's own per-graph pieces",
            pieces_total={"g5_89_nat": A["nonzero_pieces"], "g5_61_nat": B["nonzero_pieces"],
                          "g5_89_nat_ab1o3_m2o7": C["nonzero_pieces"]},
            pieces_by_S={"g5_89_nat": A["pieces_by_S"], "g5_61_nat": B["pieces_by_S"],
                         "g5_89_nat_ab1o3_m2o7": C["pieces_by_S"]},
            pieces_by_E_inner={"g5_89_nat": A["pieces_by_E_inner"], "g5_61_nat": B["pieces_by_E_inner"]},
            split=A["split"], coefficients=A["coefficients"],
            split_same_in_other_frame=(C["split"] == A["split"] and C["coefficients"] == A["coefficients"]),
            per_E=per_E_table(A, Eof, L.K),
            largest_cancelling_sets=A["largest_cancelling_sets"],
            cancelling_set_size_hist=A["cancelling_set_size_hist"],
            per_row={"g5_89_nat": A["per_row"]},
            two_primes=core["prime_comparison"],
            frame_comparison=dict(note="same comparison between the two frames at p89 (not a prime comparison)",
                                  **ab["prime_comparison"])),
        step4_valuation=dict(
            statement="a stored piece landing on an outer graph with |S| edges has its first |S| coefficients zero",
            violations={"g5_89_nat": A["valuation"]["violations"], "g5_61_nat": B["valuation"]["violations"],
                        "g5_89_nat_ab1o3_m2o7": C["valuation"]["violations"]},
            floor_attained={"g5_89_nat": A["valuation"]["floor_attained"],
                            "g5_61_nat": B["valuation"]["floor_attained"],
                            "g5_89_nat_ab1o3_m2o7": C["valuation"]["floor_attained"]}),
        step5_additivity=dict(
            statement="for every inner graph the sum over the outer graphs of its pieces equals its point-frame term",
            g5_89_nat=A["additivity"], g5_61_nat=B["additivity"],
            g5_89_nat_ab1o3_m2o7=C["additivity"]),
        step6_genus4_vs_genus5=g4g5(A, Eof))
    if abadd and os.path.exists(abadd):
        rep["step5_additivity"]["g5_89_nat_ab1o3_m2o7_closed"] = json.load(open(abadd))
    with open(outp, "w") as f:
        json.dump(rep, f, indent=1)
    print_tables(rep)
    print("written", outp)


def g4g5(A, Eof):
    hist5 = collections.Counter(Eof)
    g4 = dict(genus=4, graphs=G4["graphs"], polar_coefficients=G4["K"], cap_on_S=G4["K"] - 1,
              free_by_selection=sum(v for e, v in G4["E_hist"].items() if e >= G4["K"]),
              edge_graphs_E1_to_cap=sum(v for e, v in G4["E_hist"].items() if 1 <= e <= G4["K"] - 1),
              stored_rows=1 + sum(G4["stored_rows_by_E"].values()),
              cancelling_rows=sum(G4["stored_rows_by_E"].values()),
              coeffs_at_stake=sum((G4["K"] - e) * v for e, v in G4["E_hist"].items() if 1 <= e <= G4["K"] - 1),
              coeffs_in_cancelling_rows=sum((G4["K"] - e) * v for e, v in G4["stored_rows_by_E"].items()),
              nonzero_pieces=G4["nonzero_pieces"], cpu_minutes=G4["cpu_minutes"], source=G4["source"])
    g4["no_nonzero_piece"] = g4["edge_graphs_E1_to_cap"] - g4["cancelling_rows"]
    g5 = dict(genus=5, graphs=L.NIN, polar_coefficients=L.K, cap_on_S=L.K - 1,
              free_by_selection=A["split"]["free_by_selection"],
              edge_graphs_E1_to_cap=A["split"]["small_E1_to_7"], stored_rows=A["stored_rows"],
              cancelling_rows=A["split"]["cancelling"], no_nonzero_piece=A["split"]["no_nonzero_piece"],
              coeffs_at_stake=A["coefficients"]["at_stake_E1_to_7"],
              coeffs_in_cancelling_rows=A["coefficients"]["in_cancelling_rows"],
              nonzero_pieces=A["nonzero_pieces"], source="this census")
    return dict(genus4=g4, genus5=g5, g5_E_hist={str(k): v for k, v in sorted(hist5.items())},
                g4_E_hist={str(k): v for k, v in sorted(G4["E_hist"].items())})


def print_tables(rep):
    s3 = rep["step3_census"]
    print("\nper |S|:", s3["pieces_by_S"]["g5_89_nat"])
    print("per E table:")
    for r in s3["per_E"]:
        print("  E=%d graphs=%4d rows=%3d nopiece=%4d stake=%4d cancel=%4d pieces=%6d" %
              (r["E"], r["graphs"], r["rows_with_pieces"], r["no_piece"], r["coeffs_at_stake"],
               r["coeffs_in_cancelling_rows"], r["pieces"]))
    g = rep["step6_genus4_vs_genus5"]
    print("genus 4:", {k: v for k, v in g["genus4"].items() if k != "source"})
    print("genus 5:", {k: v for k, v in g["genus5"].items() if k != "source"})


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], a[3], a[4], a[5] if len(a) > 5 else None)
