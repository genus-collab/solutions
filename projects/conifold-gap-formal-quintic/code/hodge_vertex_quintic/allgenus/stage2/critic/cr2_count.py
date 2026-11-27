"""Stage-2 critic: recount the non-zero regular coefficients of the genus-3 flat / hatted tables of the flat worker
(S15_PROOF s.0 item 2 says 112 for the edge graphs; its referee says 108 + 4 for the vertex graph).
Input: JSON tables {"p": P, "A": {graph index: [coefficients of u^-4 .. u^3]}} as written by flat/s15_flat.py, s15_hat.py
(kept in the session scratch folder allgenus_s2_flat/).  Graph '1' is the vertex graph (Mathematica order).
Usage: py -I -B cr2_count.py TABLE.json [TABLE.json ...]      (stdlib only; writes nothing)
"""
import json
import sys


def main():
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        P, A = d["p"], d["A"]
        edge_reg = sum(1 for g, row in A.items() if g != "1" for c in row[4:8] if c % P)
        vert_reg = sum(1 for c in A["1"][4:8] if c % P)
        polar = sorted((g for g, row in A.items() if g != "1" and any(c % P for c in row[0:4])), key=int)
        print("%s: %d graphs; non-zero regular coefficients (u^0..u^3): edge graphs %d, vertex graph %d; "
              "edge graphs with a polar term: %s" % (fn.replace("\\", "/").split("/")[-1], len(A), edge_reg, vert_reg, polar))


if __name__ == "__main__":
    main()
