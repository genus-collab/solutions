# -*- coding: utf-8 -*-
"""g5l_additivity.py -- item (3) of the local audit: per-inner-graph additivity, checked on this machine.

The claim of the method (fq4_split.py): for EVERY inner graph t, the sum over all outer-edge subsets S of its pieces equals the
point-frame term of t (variant 'id', M = 1, where the outer graph is t itself).  The cluster's collect step checks this against the
shipped parts/g5_<p>_id/local.json, but its per-graph pieces were not copied back, so the check is redone here.

WHAT --maxout m CAN AND CANNOT GIVE.  With --maxout m the pieces of t are only those of the subsets with |S| <= m, i.e. a PARTIAL
sum, so the sum of the stored pieces equals the point-frame term only for the inner graphs with E(t) <= m (then every subset is
kept).  For the other inner graphs the honest statement is weaker: the truncated sum is the restriction of the full piece family to
the outer graphs with E <= m.  This script therefore checks
  (A) FULL additivity, sum over all pieces = point-frame term, for every inner graph computed with maxout = None locally
      (tag --full, here the low-E inner graphs that fit in the local budget) -- the real check;
  (B) the same for the inner graphs with E(t) <= m of the --maxout m run, where the truncation is vacuous;
  (C) RESTRICTION consistency: for every inner graph present in both runs, the maxout-m pieces equal the pieces of the full run
      restricted to the outer graphs with E <= m (so the cheap run really is a sub-sum of the full one).

  py -I -B g5l_additivity.py [--trunc g5_89_nat_out1] [--full g5_89_nat] [--id g5_89_id] [--maxout 1] [--prime 89]
"""
import os, sys, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g5l_lib as L

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trunc", default="g5_89_nat_out1")
    ap.add_argument("--full", default="g5_89_nat")
    ap.add_argument("--id", default="g5_89_id")
    ap.add_argument("--maxout", type=int, default=1)
    ap.add_argument("--prime", default="89")
    a = ap.parse_args()
    log, fh = L.log_opener(os.path.join(OUT, "additivity_%s.log" % a.trunc))
    p = L.prime(a.prime)
    m = a.maxout
    FAIL = []
    res = {}

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok:
            FAIL.append(name)

    Eof = [G["E"] for G in L.outer_meta()]
    _, idp, _ = L.load_parts(a.id)
    log("point-frame reference parts/%s: %d inner graphs" % (a.id, len(idp)))
    check("the shipped point-frame run covers all %d inner graphs" % L.NIN, len(idp) == L.NIN, len(idp))
    nz_id = sum(1 for t in idp if any(idp[t].get(t, [0] * L.K)))
    log("   of them %d have a non-zero point-frame term" % nz_id)
    # the reference parts must be the ones behind the shipped point-frame result
    if os.path.exists(L.result_path(a.id)):
        idr = L.load_result(a.id)
        tab = L.merge(idp, p)
        badid = [q for q in set(list(tab) + [int(x) for x in idr["table"]]) if L.row(tab, q) != L.row(idr["table"], q)]
        tot = [0] * L.K
        for q in tab:
            for k in range(L.K):
                tot[k] = (tot[k] + tab[q][k]) % p
        check("the point-frame parts re-merge to the shipped results/%s.json table" % a.id, not badid, badid[:8])
        check("their sum equals its info['total']", tot == [x % p for x in idr["info"]["total"]],
              [L.fmt_rat(x, p) for x in tot[:3]])

    def ref_term(t):
        return [x % p for x in idp.get(str(t), {}).get(str(t), [0] * L.K)]

    # (A) full additivity on the locally recomputed full-maxout graphs
    nf, bad = 0, []
    byE = {}
    if a.full:
        names, fullp, ftimes = L.load_parts(a.full)
        log("local full-maxout part files %s: %d inner graphs, CPU %.2f h" % (names, len(fullp), sum(ftimes.values()) / 3600))
        for t, pieces_t in sorted(fullp.items(), key=lambda kv: int(kv[0])):
            s = L.piece_sum(pieces_t, p)
            nf += 1
            E = Eof[int(t)]
            byE[E] = byE.get(E, 0) + 1
            if s != ref_term(int(t)):
                bad.append(int(t))
        check("(A) FULL additivity: for each of the %d locally recomputed inner graphs the pieces add up to its point-frame term" % nf,
              nf > 0 and not bad, bad[:20])
        log("   inner graphs checked by edge number E(t): %s" % sorted(byE.items()))
        nzt = sum(1 for t in fullp if any(ref_term(int(t))))
        log("   of them %d have a non-zero point-frame term (so the identity is not 0 = 0)" % nzt)
        res.update(full_checked=nf, full_bad=bad, full_byE=byE, full_nonzero_target=nzt)

    # (B) the inner graphs of the truncated run whose truncation is vacuous
    _, tp, ttimes = L.load_parts(a.trunc)
    log("local maxout-%d part files: %d inner graphs, CPU %.2f h" % (m, len(tp), sum(ttimes.values()) / 3600))
    small = [t for t in range(L.NIN) if Eof[t] <= m and str(t) in tp]
    badB = [t for t in small if L.piece_sum(tp[str(t)], p) != ref_term(t)]
    check("(B) additivity for the %d inner graphs with E(t) <= %d of the maxout-%d run (truncation vacuous)" % (len(small), m, m),
          not badB, badB)

    # (C) restriction consistency between the two local runs
    both = sorted(set(tp) & set(fullp), key=int) if a.full else []
    badC = []
    for t in both:
        rtrunc = {q: [x % p for x in c] for q, c in tp[t].items()}
        rfull = {q: [x % p for x in c] for q, c in fullp[t].items() if Eof[int(q)] <= m and any(c)}
        if rtrunc != rfull:
            badC.append(int(t))
    check("(C) RESTRICTION: on the %d inner graphs in both runs, the maxout-%d pieces are exactly the full pieces with E(q) <= %d" % (
        len(both), m, m), bool(both) and not badC, badC[:20])
    res.update(trunc_small=len(small), trunc_small_bad=badB, both=len(both), restriction_bad=badC)

    # how much of the full sum the truncation sees, on the graphs we have both for (quantifies how weak (B)/(C) are alone)
    if both:
        nsub_tr = sum(len(tp[t]) for t in both)
        nsub_fu = sum(len(fullp[t]) for t in both)
        log("   non-zero pieces on those %d inner graphs: %d in the maxout-%d run, %d in the full run" % (
            len(both), nsub_tr, m, nsub_fu))
        res.update(pieces_trunc_on_both=nsub_tr, pieces_full_on_both=nsub_fu)

    json.dump(dict(prime=a.prime, maxout=m, failed=FAIL, **res),
              open(os.path.join(OUT, "additivity_%s.json" % a.trunc), "w"), indent=1)
    log("G5L_ADDITIVITY %s" % ("ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL))
    fh.close()
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
