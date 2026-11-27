# -*- coding: utf-8 -*-
"""g5l_cancel.py -- item (4) of the local audit: the test is not vacuous.

For every outer graph whose row is COMPLETE in the local --maxout m run (E(q) <= m) the row is a sum of pieces coming from many
different inner graphs.  This script exhibits those pieces: how many are non-zero, what they are, and that they cancel.  It also
counts all non-zero pieces of the run, so that the "all rows zero" statement cannot be mistaken for "nothing was computed".

  py -I -B g5l_cancel.py [--trunc g5_89_nat_out1] [--maxout 1] [--prime 89] [--show 12]
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g5l_lib as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trunc", default="g5_89_nat_out1")
    ap.add_argument("--maxout", type=int, default=1)
    ap.add_argument("--prime", default="89")
    ap.add_argument("--show", type=int, default=12)
    a = ap.parse_args()
    log, fh = L.log_opener(os.path.join(OUT, "cancel_%s.log" % a.trunc))
    p = L.prime(a.prime)
    m = a.maxout
    meta = L.outer_meta()
    Eof = [G["E"] for G in meta]
    _, tp, _ = L.load_parts(a.trunc)

    npieces = sum(len(v) for v in tp.values())
    log("local run %s: %d inner graphs, %d non-zero pieces in all" % (a.trunc, len(tp), npieces))
    nz_by_E = {}
    for t, r in tp.items():
        for q in r:
            e = Eof[int(q)]
            nz_by_E[e] = nz_by_E.get(e, 0) + 1
    log("   non-zero pieces by the edge number of the outer graph: %s" % sorted(nz_by_E.items()))

    # contributions per outer graph
    contrib = {}
    for t, r in tp.items():
        for q, c in r.items():
            contrib.setdefault(int(q), []).append((int(t), [x % p for x in c]))
    rows = {}
    summary = []
    for q in sorted(contrib):
        if Eof[q] > m:
            continue                      # incomplete row, cancellation would be meaningless
        lst = contrib[q]
        s = [0] * L.K
        for _, c in lst:
            for k in range(L.K):
                s[k] = (s[k] + c[k]) % p
        rows[q] = s
        summary.append((q, Eof[q], len(lst), not any(s)))

    log("")
    log("COMPLETE rows of the local run (E(q) <= %d), with the number of non-zero contributions:" % m)
    for q, E, n, zero in summary:
        log("   outer graph %3d  E=%d  gv=%s  aut=%d : %5d non-zero pieces, row %s" % (
            q, E, meta[q]["gv"], meta[q]["aut"], n, "= 0" if zero else "NON-ZERO"))
    edged = [s for s in summary if s[1] >= 1]
    worst = max((n for _, _, n, _ in edged), default=0)
    allzero = all(z for _, E, _, z in edged for z in [z])
    log("")
    log("%s: every complete row with edges is exactly zero in all %d coefficients" % (
        "PASS" if allzero and edged else "FAIL", L.K))
    log("largest number of non-zero contributions cancelling for a single outer graph: %d" % worst)

    # concrete sample: the first --show non-zero contributions of each outer graph with edges, and running partial sums.
    # An outer edge carries a factor E_out = O(u), so for E(q) = 1 the coefficient u^0 of every piece vanishes identically and
    # the first coefficient that says anything is u^1; residues mod p are printed (a single prime rarely reconstructs a piece).
    for q, E, n, zero in edged:
        lst = sorted(contrib[q], key=lambda x: x[0])
        kk = next((k for k in range(L.K) if any(c[k] for _, c in lst)), 1)
        log("")
        log("outer graph %d (E=%d, gv=%s): %d non-zero pieces; the first %d of them, coefficient u^%d (residues mod 2^%s-1):" % (
            q, E, meta[q]["gv"], n, min(a.show, n), kk, a.prime))
        run = 0
        for t, c in lst[:a.show]:
            run = (run + c[kk]) % p
            log("   inner graph %4d (E=%2d, gv=%-22s): %28d   partial sum %28d" % (
                t, meta[t]["E"], str(meta[t]["gv"]), c[kk], run))
        tot = [0] * L.K
        for t, c in lst:
            for k in range(L.K):
                tot[k] = (tot[k] + c[k]) % p
        log("   ... sum over ALL %d pieces: %s" % (n, "ZERO in every one of the %d coefficients" % L.K if not any(tot) else tot))
        nzk = [k for k in range(L.K) if sum(1 for _, c in lst if c[k]) > 1]
        log("   coefficients u^k with more than one non-zero contribution: %s (counts %s)" % (
            nzk, [sum(1 for _, c in lst if c[k]) for k in nzk]))
        drop = lst[0]
        d = [(tot[k] - drop[1][k]) % p for k in range(L.K)]
        log("   the sum of the other %d pieces alone (piece of inner graph %d removed) is %s" % (
            n - 1, drop[0], "NON-ZERO: no piece is redundant" if any(d) else "still zero"))

    json.dump(dict(trunc=a.trunc, maxout=m, prime=a.prime, pieces_total=npieces, pieces_by_E=nz_by_E,
                   complete_rows=[dict(q=q, E=E, contributions=n, row_zero=z) for q, E, n, z in summary],
                   largest_cancelling=worst, all_edged_rows_zero=bool(allzero and edged)),
              open(os.path.join(OUT, "cancel_%s.json" % a.trunc), "w"), indent=1)
    log("")
    log("G5L_CANCEL done")
    fh.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
