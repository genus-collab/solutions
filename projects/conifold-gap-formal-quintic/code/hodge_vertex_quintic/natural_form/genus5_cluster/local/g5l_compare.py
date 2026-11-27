# -*- coding: utf-8 -*-
"""g5l_compare.py -- item (2) of the local audit: compare the rows of the LOCAL --maxout m run, entry by entry, with the rows
of the cluster run results/g5_89_nat.json.

With --maxout m the engine keeps exactly the outer-edge subsets S with |S| <= m, and the outer graph of a subset S has E = |S|
edges; so the local rows of the outer graphs with E <= m are COMPLETE and must equal the cluster's rows exactly, while every
other outer graph must receive nothing at all.  Both facts are checked.  The local table is re-merged here from the per-inner-graph
part files, so the package's own collect step is not trusted; the collect result (if present) is compared with it as well.

  py -I -B g5l_compare.py [--local g5_89_nat_out1] [--ref g5_89_nat] [--maxout 1] [--prime 89]
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g5l_lib as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", default="g5_89_nat_out1")
    ap.add_argument("--ref", default="g5_89_nat")
    ap.add_argument("--maxout", type=int, default=1)
    ap.add_argument("--prime", default="89")
    ap.add_argument("--variant", default="nat")
    ap.add_argument("--control", action="store_true",
                    help="the run is a FALSIFICATION control (bare M_0 or a shifted tau): the complete rows with edges must then "
                         "be NON-zero and the vertex graph must NOT be c t_c^-8 + 0 + ... + 0")
    a = ap.parse_args()
    log, fh = L.log_opener(os.path.join(OUT, "compare_%s.log" % a.local))
    p = L.prime(a.prime)
    m = a.maxout
    FAIL = []

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok:
            FAIL.append(name)

    ref = L.load_result(a.ref)
    meta = ref["outer"]
    Eof = [G["E"] for G in meta]
    names, pieces, times = L.load_parts(a.local)
    log("local part files %s: %d inner graphs, CPU %.2f h (wall was logged by the run)" % (names, len(pieces), sum(times.values()) / 3600))
    log("reference %s: tag %s, failed %s, %d stored rows, %d outer graphs" % (
        L.result_path(a.ref), ref["tag"], ref["failed"], len(ref["table"]), len(meta)))

    covered = set(int(t) for t in pieces)
    check("the local run covers all %d inner graphs" % L.NIN, covered == set(range(L.NIN)),
          "%d covered, %d missing" % (len(covered), L.NIN - len(covered)))

    tab = L.merge(pieces, p)
    log("local merged table: %d rows with a stored piece" % len(tab))

    # structural: with maxout m only the outer graphs with E <= m may receive anything
    wrongE = sorted(q for q in tab if Eof[q] > m)
    check("every outer graph with a local piece has E <= %d" % m, not wrongE, wrongE[:10])
    testable = [q for q in range(len(meta)) if Eof[q] <= m]
    log("outer graphs with E <= %d (exact in the local run): %s" % (m, testable))

    # entry-by-entry comparison on the complete rows
    nentry = 0
    diffs = []
    for q in testable:
        lo = L.row(tab, q)
        rf = L.row(ref["table"], q)
        for k in range(L.K):
            nentry += 1
            if lo[k] != rf[k]:
                diffs.append((q, k, lo[k], rf[k]))
    check("all %d entries of the %d complete rows agree with the cluster run" % (nentry, len(testable)),
          not diffs, diffs[:8])
    for q in testable:
        lo = L.row(tab, q)
        log("   outer graph %3d (E=%d, gv=%s, aut=%d): %s" % (
            q, Eof[q], meta[q]["gv"], meta[q]["aut"],
            "ZERO (all %d coefficients)" % L.K if not any(lo) else [L.fmt_rat(x, p) for x in lo]))
    nz = [q for q in testable if any(L.row(tab, q))]
    if a.control:
        check("CONTROL: at least one complete row with edges is non-zero (the test discriminates)",
              [q for q in nz if Eof[q] >= 1], nz)
    else:
        check("exactly one complete row is non-zero, the vertex graph q = 0", nz == [0], nz)

    # the vertex row carries the whole polar part: in t_c it must be c_5 t_c^-8 + 0 ... + 0
    try:
        from fq_hat import tc_series
        from fq3_msum import tc_expansion3
        from fq3_rec import Rec
        from g5_engine import make_split
        sp = make_split(p, 5, ("0", "0"), L.K, a.variant)
        R = sp.R
        tc = tc_series(sp.C, R, L.K + 4)
        T = tc_expansion3(L.row(tab, 0), tc, R, 5)[:8]
        log("   local vertex graph in t_c: %s" % " + ".join(
            "(%s) t_c^%d" % (L.fmt_rat(T[k], p), k - 8) for k in range(8)))
        if a.control:
            check("CONTROL: the local vertex graph is NOT (c) t_c^-8 + 0 + ... + 0 t_c^-1", any(T[1:]))
        else:
            check("local vertex graph = (c_5) t_c^-8 + 0 + ... + 0 t_c^-1", not any(T[1:]))
        c5 = R.of(__import__("g5_checks").c_conj(5))
        check("local c_5 equals |B_10|/(10*8*5^4) = 1/660000", T[0] == c5, L.fmt_rat(T[0], p))
        check("local vertex t_c vector equals the cluster's info['vertex']", T == ref["info"]["vertex"])
    except Exception as e:
        log("   (t_c check skipped: %r)" % (e,))
        FAIL.append("t_c check raised")

    # discrimination: the same complete rows must DIFFER from the two shipped controls (bare M_0, tau with its u^7 coefficient
    # shifted by 1), which were run on the cluster with --maxout 2 and are therefore exact on these rows too
    if not a.control:
        for ctl in ("g5_89_M0_out2", "g5_89_nat_t7_out2"):
            if not os.path.exists(L.result_path(ctl)):
                continue
            c = L.load_result(ctl)
            differ = [q for q in testable if L.row(c["table"], q) != L.row(tab, q)]
            nzc = [q for q in testable if Eof[q] >= 1 and any(L.row(c["table"], q))]
            check("the local rows DIFFER from the control %s (it is a real test, not an identity)" % ctl,
                  bool(differ), differ)
            log("   control %s: complete rows differing from the local run: %s; its non-zero rows with edges: %s" % (
                ctl, differ, nzc))

    # the local collect step, if it ran, must agree with this independent merge
    if os.path.exists(L.result_path(a.local)):
        mine = L.load_result(a.local)
        bad = [q for q in testable if L.row(mine["table"], q) != L.row(tab, q)]
        check("the local collect result agrees with this independent re-merge on the complete rows", not bad, bad[:8])
        log("   local collect reported failed = %s; nonregular = %s" % (mine["failed"], mine["info"].get("nonregular")))

    json.dump(dict(local=a.local, ref=a.ref, maxout=m, prime=a.prime, testable=testable,
                   entries_compared=nentry, mismatches=len(diffs), inner_covered=len(covered),
                   rows_with_pieces=len(tab), failed=FAIL),
              open(os.path.join(OUT, "compare_%s.json" % a.local), "w"), indent=1)
    log("G5L_COMPARE %s" % ("ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL))
    fh.close()
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
