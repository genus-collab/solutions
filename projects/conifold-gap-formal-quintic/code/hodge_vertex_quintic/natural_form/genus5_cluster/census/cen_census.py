# -*- coding: utf-8 -*-
"""cen_census.py -- steps 3, 4 and 5 of the census, from the cluster's own per-graph pieces.

One streaming pass over parts/<tag>/task*.json per run (the part files are read one at a time, so the memory stays small):
  (3) the cancellation census: pieces by |S| = E(outer) and by E(inner), per-row piece and per-coefficient counts, the
      627/1036/2891 split, the 3142/1655 coefficient counts, and the non-zero pattern of every piece (an 8-bit mask) so the
      two primes can be compared piece by piece and coefficient by coefficient;
  (4) the valuation law: a piece landing on an outer graph with |S| edges must have its first |S| coefficients zero, and the
      floor |S| is attained when coefficient |S| is non-zero;
  (5) additivity: for every inner graph the sum over the outer graphs of its pieces equals its point-frame term (the id run).

  py -I -B cen_census.py <out.json> [tag ...]
"""
import os, sys, json, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L


def scan(tag, Eof):
    p = L.prime_of(tag)
    rows = {}          # outer q -> [n_pieces, per-coefficient non-zero counts (8), summed row]
    byS = collections.Counter()
    byEin = collections.Counter()
    byS_Ein = collections.Counter()
    attained = collections.Counter()
    viol = []
    masks = {}         # inner t -> {outer q -> 8-bit non-zero mask}
    addsum = {}        # inner t -> summed row over all its outer graphs
    npieces = 0
    inner_seen = set()
    for name in L.part_files(tag):
        d = L.load_part(tag, name)
        for ts, res in d["pieces"].items():
            t = int(ts)
            inner_seen.add(t)
            s = [0] * L.K
            mk = {}
            for qs, c in res.items():
                q = int(qs)
                if not any(c):
                    continue                      # never stored, but be safe
                npieces += 1
                e = Eof[q]
                byS[e] += 1
                byEin[Eof[t]] += 1
                byS_Ein[(e, Eof[t])] += 1
                if any(c[:e]):
                    viol.append([t, q, e, [k for k in range(e) if c[k]]])
                if e < L.K and c[e]:
                    attained[e] += 1
                r = rows.get(q)
                if r is None:
                    r = rows[q] = [0, [0] * L.K, [0] * L.K]
                r[0] += 1
                for k in range(L.K):
                    if c[k]:
                        r[1][k] += 1
                    r[2][k] = (r[2][k] + c[k]) % p
                    s[k] = (s[k] + c[k]) % p
                mk[q] = sum(1 << k for k in range(L.K) if c[k])
            masks[t] = mk
            addsum[t] = s
    return dict(p=p, rows=rows, byS=byS, byEin=byEin, byS_Ein=byS_Ein, attained=attained, viol=viol,
                masks=masks, addsum=addsum, npieces=npieces, inner=inner_seen)


def id_terms(tag):
    """point-frame term of every inner graph, from the id run of the same prime: pieces[t][t]"""
    idtag = "g5_%s_id" % tag.split("_")[1]
    A = L.load_all_pieces(idtag)
    return idtag, {t: v.get(t, [0] * L.K) for t, v in A["pieces"].items()}


def summarise(tag, S, Eof):
    p = S["p"]
    rows = S["rows"]
    nz_rows = sorted(q for q, r in rows.items() if any(r[2]))
    zero_rows = sorted(q for q, r in rows.items() if not any(r[2]))
    edge_rows = [q for q in rows if Eof[q] >= 1]
    # the three classes among the 4554 outer graphs with at least one edge
    free = [i for i in range(1, L.NIN) if Eof[i] >= L.K]                       # regular a priori (|S| <= K-1 = 7)
    small = [i for i in range(1, L.NIN) if 1 <= Eof[i] <= L.K - 1]
    nopiece = [i for i in small if i not in rows]
    cancel = [i for i in small if i in rows]
    at_stake = sum(L.K - Eof[i] for i in small)
    cancel_coeffs = sum(L.K - Eof[i] for i in cancel)
    # the per-coefficient cancellation structure of the rows that do receive pieces
    multi = collections.Counter()
    biggest = []
    for q in sorted(rows):
        if Eof[q] == 0:
            continue
        r = rows[q]
        m = max(r[1][Eof[q]:]) if Eof[q] < L.K else 0
        multi[m] += 1
        biggest.append([q, Eof[q], r[0], max(r[1]), r[1]])
    biggest.sort(key=lambda x: (-x[3], -x[2], x[0]))
    out = dict(
        tag=tag, prime=("2^89-1" if p == L.PRIME["89"] else "2^61-1"), nonzero_pieces=S["npieces"],
        pieces_by_S={str(k): v for k, v in sorted(S["byS"].items())},
        pieces_by_E_inner={str(k): v for k, v in sorted(S["byEin"].items())},
        stored_rows=len(rows), rows_nonzero=nz_rows, rows_zero_count=len(zero_rows),
        rows_with_edges=len(edge_rows),
        split=dict(edge_graphs=L.NIN - 1, free_by_selection=len(free), small_E1_to_7=len(small),
                   no_nonzero_piece=len(nopiece), cancelling=len(cancel)),
        coefficients=dict(at_stake_E1_to_7=at_stake, in_cancelling_rows=cancel_coeffs,
                          in_rows_without_pieces=at_stake - cancel_coeffs),
        valuation=dict(violations=len(S["viol"]), examples=S["viol"][:10],
                       floor_attained={str(k): [S["attained"].get(k, 0), S["byS"].get(k, 0)] for k in sorted(S["byS"])}),
        largest_cancelling_sets=biggest[:25],
        cancelling_set_size_hist={str(k): v for k, v in sorted(multi.items())},
        per_row={str(q): dict(E=Eof[q], pieces=rows[q][0], per_coeff=rows[q][1],
                              nonzero_total=[k for k in range(L.K) if rows[q][2][k]]) for q in sorted(rows)},
    )
    idtag, idt = id_terms(tag)
    bad = [t for t in range(L.NIN) if S["addsum"].get(t, [0] * L.K) != [x % p for x in idt.get(t, [0] * L.K)]]
    if "_ab" in tag:
        out["additivity"] = dict(
            reference=None, inner_graphs_checked=len(S["addsum"]),
            skipped="no point-frame (id) run exists in the frame (alpha,beta) = (1/3,-2/7): the sum over the outer graphs "
                    "of a graph's pieces is its term in that frame, and no independent reference value is stored",
            differing_from_frame_00_point_term=len(bad))
    else:
        out["additivity"] = dict(reference=idtag, inner_graphs_checked=len(S["addsum"]), failures=len(bad),
                                 failing=bad[:20], all_pass=(not bad and len(S["addsum"]) == L.NIN))
    return out


def main(outpath, tags):
    meta = L.graph_meta()
    Eof = [g["E"] for g in meta]
    rep = dict(tags={}, Eof_hist={str(k): v for k, v in sorted(collections.Counter(Eof).items())})
    masks = {}
    for tag in tags:
        S = scan(tag, Eof)
        rep["tags"][tag] = summarise(tag, S, Eof)
        masks[tag] = S["masks"]
        r = rep["tags"][tag]
        print("%-22s pieces=%d rows=%d nz_rows=%s split=%s viol=%d add_fail=%d"
              % (tag, r["nonzero_pieces"], r["stored_rows"], r["rows_nonzero"], r["split"],
                 r["valuation"]["violations"], r["additivity"].get("failures", -1)))
        del S
    if len(tags) >= 2:
        a, b = tags[0], tags[1]
        ka = set((t, q) for t, m in masks[a].items() for q in m)
        kb = set((t, q) for t, m in masks[b].items() for q in m)
        only_a = sorted(ka - kb)
        only_b = sorted(kb - ka)
        cdiff = []
        for t, q in sorted(ka & kb):
            x = masks[a][t][q] ^ masks[b][t][q]
            if x:
                cdiff.append([t, q, Eof[q], [k for k in range(L.K) if x >> k & 1],
                              masks[a][t][q], masks[b][t][q]])
        rep["prime_comparison"] = dict(
            tags=[a, b], pieces=[len(ka), len(kb)], same_piece_set=(not only_a and not only_b),
            only_in_first=only_a[:50], only_in_second=only_b[:50],
            n_only_first=len(only_a), n_only_second=len(only_b),
            coefficientwise_differences=cdiff[:50], n_coefficientwise_differences=len(cdiff))
        print("prime comparison: same piece set = %s (%d / %d), coefficient-wise differences = %d"
              % (rep["prime_comparison"]["same_piece_set"], len(only_a), len(only_b), len(cdiff)))
    with open(outpath, "w") as f:
        json.dump(rep, f, indent=1)
    print("written", outpath)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or ["g5_89_nat", "g5_61_nat"])
