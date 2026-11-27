# -*- coding: utf-8 -*-
"""f2_sel.py -- is the SELECTION RULE frame dependent?  Compares which pieces / rows / coefficients are non-zero in the
second conifold frame (alpha, beta) = (1/3, -2/7) with the frame (0, 0), both natural runs at p = 2^89-1.

Three levels, all computed here from the part files (nothing of g5_checks is used except nothing at all):
  (a) PIECES      : the set of (inner t -> outer q) pairs carrying a stored (i.e. somewhere non-zero) piece;
  (b) COEFFICIENTS: for each common piece, which of the K = 8 coefficients of u^0..u^7 are non-zero (a bitmask);
  (c) ROWS        : the set of outer graphs q whose merged row (sum over t) is non-zero, and its u-support.
"inner graph" and "outer graph" are indices into the same list of 4555 graphs, so the two runs are directly comparable.

usage:  py -I -B f2_sel.py <parts dir frame (0,0)> <parts dir frame2> [--prime 89] [--json OUT]
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_lib import K, NIN, prime, read_parts_dir


def mask(c, p):
    return sum(1 << k for k in range(K) if c[k] % p)


def collect(pieces, p):
    """-> (set of (t, q), {(t, q): u-support mask}, merged rows {q: [K residues]})"""
    keys = set()
    masks = {}
    rows = {}
    for t, res in pieces.items():
        for q, c in res.items():
            key = (int(t), int(q))
            keys.add(key)
            masks[key] = mask(c, p)
            r = rows.setdefault(int(q), [0] * K)
            for k in range(K):
                r[k] = (r[k] + c[k]) % p
    return keys, masks, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirA", help="parts dir of the frame (0,0) natural run")
    ap.add_argument("dirB", help="parts dir of the second-frame natural run")
    ap.add_argument("--prime", default="89")
    ap.add_argument("--json", default="")
    ap.add_argument("--tagA", default="")
    ap.add_argument("--tagB", default="")
    ap.add_argument("--primeA", default="")
    a = ap.parse_args()
    p = prime(a.prime)
    tagA = a.tagA or "g5_%s_nat" % a.prime
    tagB = a.tagB or "g5_%s_nat_ab1o3_m2o7" % a.prime
    _, pA, _ = read_parts_dir(a.dirA, tagA)
    _, pB, _ = read_parts_dir(a.dirB, tagB)
    print("%s: %d inner graphs;  %s: %d inner graphs" % (tagA, len(pA), tagB, len(pB)))
    assert len(pA) == NIN and len(pB) == NIN, "both runs must be complete"

    pa = prime(a.primeA) if a.primeA else p
    kA, mA, rA = collect(pA, pa)
    kB, mB, rB = collect(pB, p)
    out = {"tagA": tagA, "tagB": tagB, "prime": a.prime}

    onlyA = sorted(kA - kB)
    onlyB = sorted(kB - kA)
    out["pieces_A"] = len(kA)
    out["pieces_B"] = len(kB)
    out["pieces_only_A"] = [list(x) for x in onlyA[:50]]
    out["pieces_only_B"] = [list(x) for x in onlyB[:50]]
    out["n_pieces_only_A"] = len(onlyA)
    out["n_pieces_only_B"] = len(onlyB)
    print("(a) PIECES: %d in frame (0,0), %d in the second frame; only in (0,0): %d; only in the second frame: %d"
          % (len(kA), len(kB), len(onlyA), len(onlyB)))

    common = kA & kB
    diffm = [(t, q, mA[(t, q)], mB[(t, q)]) for (t, q) in sorted(common) if mA[(t, q)] != mB[(t, q)]]
    out["n_common_pieces"] = len(common)
    out["n_pieces_with_different_u_support"] = len(diffm)
    out["pieces_with_different_u_support"] = [[t, q, x, y] for t, q, x, y in diffm[:50]]
    print("(b) COEFFICIENTS: %d common pieces, %d with a different u-support bitmask" % (len(common), len(diffm)))
    # the whole distribution of masks, as a structural fingerprint
    from collections import Counter
    cA = Counter(mA.values())
    cB = Counter(mB.values())
    out["mask_hist_A"] = {str(k): v for k, v in sorted(cA.items())}
    out["mask_hist_B"] = {str(k): v for k, v in sorted(cB.items())}
    print("    u-support histogram identical: %s" % (cA == cB))

    nzA = sorted(q for q, r in rA.items() if any(r))
    nzB = sorted(q for q, r in rB.items() if any(r))
    out["nonzero_rows_A"] = len(nzA)
    out["nonzero_rows_B"] = len(nzB)
    out["rows_only_A"] = [q for q in nzA if q not in set(nzB)][:50]
    out["rows_only_B"] = [q for q in nzB if q not in set(nzA)][:50]
    rmask = {"A": {}, "B": {}}
    for lab, rows, nz, pp in (("A", rA, nzA, pa), ("B", rB, nzB, p)):
        for q in nz:
            rmask[lab][q] = mask(rows[q], pp)
    diffr = [(q, rmask["A"][q], rmask["B"][q]) for q in nzA if q in rmask["B"] and rmask["A"][q] != rmask["B"][q]]
    out["rows_with_different_u_support"] = [[q, x, y] for q, x, y in diffr[:50]]
    out["n_rows_with_different_u_support"] = len(diffr)
    print("(c) ROWS: %d non-zero merged rows in frame (0,0), %d in the second frame; sets equal: %s; "
          "rows with a different u-support: %d" % (len(nzA), len(nzB), set(nzA) == set(nzB), len(diffr)))
    # rows that are zero although some piece feeding them is not: the cancellations
    canA = sorted(q for q, r in rA.items() if not any(r))
    canB = sorted(q for q, r in rB.items() if not any(r))
    out["cancelling_rows_A"] = len(canA)
    out["cancelling_rows_B"] = len(canB)
    out["cancelling_sets_equal"] = (set(canA) == set(canB))
    print("    outer graphs that receive a piece but whose row cancels to 0: %d (0,0) / %d second frame; sets equal: %s"
          % (len(canA), len(canB), set(canA) == set(canB)))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1)
        print("written %s" % a.json)
    identical = (not onlyA and not onlyB and not diffm and set(nzA) == set(nzB) and not diffr)
    print("VERDICT: the non-zero pattern is %s across the two frames"
          % ("IDENTICAL" if identical else "DIFFERENT"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
