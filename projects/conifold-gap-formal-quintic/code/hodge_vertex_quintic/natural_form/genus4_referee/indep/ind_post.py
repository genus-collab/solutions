# -*- coding: utf-8 -*-
"""ind_post.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): post-processing of finished ind_run runs (no heavy computation).

  py -I ind_post.py --prime 89 --variant nat
 (1) PIECEWISE comparison with the project's per-inner-graph part file <project>/parts/g4_<prime>_<variant>_run.json (read only):
     for every inner graph t and every outer index the six coefficients of the piece (inner graph t -> outer graph) are compared.
     These pieces are non-zero numbers, so this is a much finer test than the comparison of the final tables (whose entries for
     the outer graphs with edges are all zero for the natural M).
 (2) Cancellation statistics: for every outer graph, the number of inner graphs with a non-zero piece through u^(2g-3).
 (3) The vertex graph and the total rewritten in t_c (coefficients of t_c^-(2g-2), .., t_c^-1 from [u^0..u^(2g-3)]).
 (4) out/<tag>_verified.txt: the list of outer graphs compared (index, E, genera, edges, |Aut|, number of non-zero pieces).
"""
import os, sys, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "referee")))
from ref_fp import P89, P61, of, inv, fmt, smul, sinv
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))


def tc_polar(coefs, tcu, p, npol):
    """u^-npol Sum a_k u^k = Sum_j A_j t_c^(j-npol): A_0..A_(npol-1), with t_c/u = tcu (series)"""
    K = npol; r = sinv(tcu[:K], p)                                   # u/t_c
    rp = [[1] + [0] * (K - 1)]
    for _ in range(npol): rp.append(smul(rp[-1], r, p, K))
    tp = [[1] + [0] * (K - 1)]                                       # (t_c/u)^j
    for _ in range(npol): tp.append(smul(tp[-1], tcu[:K], p, K))
    # G(u) = Sum_j A_j t_c^j (u/t_c)^npol ... solve triangularly: [u^k]: Sum_j A_j [u^(k-j)] (t_c/u)^j (u/t_c)^npol
    A = []
    for k in range(npol):
        s = coefs[k]
        for j in range(k): s -= A[j] * smul(tp[j], rp[npol], p, K)[k - j]
        A.append(s % p)                                              # the j = k term has coefficient 1
    return A


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--variant", default="nat")
    ap.add_argument("--genus", type=int, default=4); a = ap.parse_args(); p = {"89": P89, "61": P61}[a.prime]; g = a.genus; npol = 2 * g - 2
    tag = "g%d_%s_%s" % (g, a.prime, a.variant); F = lambda x: fmt(x, p)
    LOG = open(os.path.join(HERE, "out", tag + "_post.log"), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 100); log("ind_post %s" % tag)
    mine = json.load(open(os.path.join(HERE, "out", tag + "_parts.json")))['res']
    fin = json.load(open(os.path.join(HERE, "out", tag + ".json"))); outer = fin['outer']; nin = len(outer); K = len(next(iter(fin['table'].values())))
    zero = [0] * K
    # ---------------------------------------------------------------- (1)
    ppath = os.path.join(PROJ, "parts", "g%d_%s_%s_run.json" % (g, a.prime, a.variant))
    if os.path.exists(ppath):
        theirs = json.load(open(ppath))['res']; npieces = 0; nonzero = 0; bad = []
        for t in range(nin):
            A = mine.get(str(t), {}); B = theirs.get(str(t), {})
            for pos in set(A) | set(B):
                x = A.get(pos, zero); y = [c % p for c in B.get(pos, zero)]
                if any(x) or any(y): npieces += 1; nonzero += any(x)
                if x != y: bad.append((t, int(pos)))
        log("%s (1) piecewise: %d (inner graph, outer graph) pieces compared with parts/%s, %d non-zero; mismatches: %s"
            % ("PASS" if not bad else "FAIL", npieces, os.path.basename(ppath), nonzero, bad[:20]))
    else: log("INFO (1) no project part file %s" % os.path.basename(ppath))
    # ---------------------------------------------------------------- (2)
    cnt = [0] * nin; cnt_any = [0] * nin
    for t, r in mine.items():
        for pos, c in r.items():
            if any(c[:npol]): cnt[int(pos)] += 1
            if any(c): cnt_any[int(pos)] += 1
    tab = {int(k): v for k, v in fin['table'].items()}
    byE = {}
    for i in range(nin):
        E = outer[i]['E']; d = byE.setdefault(E, [0, 0, 0, 0]); d[0] += 1; d[1] += cnt[i]; d[2] += (1 if any(tab.get(i, zero)[:npol]) else 0); d[3] += (1 if cnt[i] else 0)
    log("(2) by number of outer edges E: [outer graphs, non-zero pieces through u^%d, outer graphs with a non-zero sum, outer graphs with at least one non-zero piece]" % (npol - 1))
    for E in sorted(byE): log("      E = %d: %s" % (E, byE[E]))
    for i in range(min(nin, 11)):
        log("      outer graph %2d (E=%d, gv=%s, ed=%s): %d non-zero pieces, sum = %s" % (i, outer[i]['E'], outer[i]['gv'], outer[i]['ed'], cnt[i], [F(x) for x in tab.get(i, zero)]))
    # ---------------------------------------------------------------- (3)
    from ref_hat import HatRows
    H = HatRows(p, 1, K + 2); tcu = [H.tc[k + 1] for k in range(K + 1)]
    for name, c in (("vertex graph", tab.get(0, zero)), ("total", fin['total'])):
        if c is None: continue
        A = tc_polar(c, tcu, p, min(npol, K))
        log("(3) %s in t_c: %s" % (name, " + ".join("(%s) t_c^%d" % (F(A[j]), j - npol) for j in range(len(A)))))
    # ---------------------------------------------------------------- (4)
    with open(os.path.join(HERE, "out", tag + "_verified.txt"), "w", encoding="utf-8") as f:
        f.write("# outer graphs of %s compared with the project table (index, E, genera, edges, |Aut|, non-zero pieces, sum through u^%d is zero)\n" % (tag, npol - 1))
        for i in fin['complete']:
            f.write("%3d  E=%d  gv=%s  ed=%s  aut=%d  pieces=%d  zero=%s\n" % (i, outer[i]['E'], outer[i]['gv'], outer[i]['ed'], outer[i]['aut'], cnt[i], not any(tab.get(i, zero)[:npol])))
    log("(4) %d complete outer graphs listed in out/%s_verified.txt" % (len(fin['complete']), tag))


if __name__ == "__main__":
    main()
