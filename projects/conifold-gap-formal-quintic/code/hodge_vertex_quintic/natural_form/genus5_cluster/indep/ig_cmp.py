# -*- coding: utf-8 -*-
"""ig_cmp.py -- compare the independent runs (out/*.jsonl) with the cluster's own stored pieces, piece by piece.

  py -I -B ig_cmp.py --prime 89 [--tag nat_89] [--idtag id_89] [--save cmp_89.json]

For every inner graph t in the independent run:
  (1) the whole family {outer graph q -> [8 residues]} is compared entry by entry with the cluster's
      parts/g5_<p>_nat/task*of*.json; a piece stored by only one of the two codes counts as a mismatch
      (the cluster stores a piece only when it is non-zero, so an omission is a claim too);
  (2) if an independent point-frame run (variant id) is present, its term is compared with the shipped
      parts/g5_<p>_id/local.json pieces[t][t];
  (3) additivity inside the independent numbers: Sum_q (independent nat piece)[q] == (independent id term).
Shards of the same graph are added before the comparison.
"""
import os, sys, json, glob, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
EV = os.environ.get("G5EV", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "g5ev"))
PRIME = {"89": (1 << 89) - 1, "61": (1 << 61) - 1}


def load_jsonl(tag):
    path = os.path.join(OUT, "%s.jsonl" % tag)
    rows = {}
    if not os.path.exists(path):
        return rows
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line:
            d = json.loads(line)
            rows.setdefault(d["t"], []).append(d)
    return rows


def merged(rows, p):
    """{t: ({q: [res]}, total seconds, n shards, aut)}"""
    out = {}
    for t, lst in rows.items():
        acc = {}
        for d in lst:
            for q, v in d["res"].items():
                a = acc.setdefault(int(q), [0] * len(v))
                for i, x in enumerate(v):
                    a[i] = (a[i] + x) % p
        out[t] = ({q: v for q, v in acc.items() if any(v)}, sum(d["sec"] for d in lst), len(lst), lst[0]["aut"])
    return out


def cluster_pieces(tag, want):
    out = {}
    for fn in sorted(glob.glob(os.path.join(EV, "parts", tag, "*.json"))):
        d = json.load(open(fn))
        assert d["tag"] == tag
        for t, fam in d["pieces"].items():
            if int(t) in want:
                out[int(t)] = {int(q): v for q, v in fam.items() if any(v)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", default="89")
    ap.add_argument("--tag", default="")
    ap.add_argument("--idtag", default="")
    ap.add_argument("--save", default="")
    a = ap.parse_args()
    p = PRIME[a.prime]
    tag = a.tag or ("nat_%s" % a.prime)
    idtag = a.idtag or ("id_%s" % a.prime)
    nat = merged(load_jsonl(tag), p)
    idd = merged(load_jsonl(idtag), p)
    ix = json.load(open(os.path.join(OUT, "index.json")))
    E = {int(t): v["E"] for t, v in ix["graphs"].items()}
    cancel = set(ix["cancel"])
    log = open(os.path.join(OUT, "cmp_%s.log" % a.prime), "a", encoding="utf-8")

    def say(*x):
        s = " ".join(str(y) for y in x)
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    say("=" * 110)
    say("ig_cmp prime %s: independent nat graphs %d, independent id graphs %d" % (a.prime, len(nat), len(idd)))
    ref = cluster_pieces("g5_%s_nat" % a.prime, set(nat))
    say("cluster families loaded for %d of them" % len(ref))
    bad = []
    npieces = 0
    nentries = 0
    onlyind = 0
    onlypkg = 0
    per_E = {}
    cancel_hit = set()
    for t in sorted(nat):
        mine, sec, nsh, aut = nat[t]
        theirs = ref[t]
        qs = set(mine) | set(theirs)
        npieces += len(qs)
        nentries += 8 * len(qs)
        for q in qs:
            if mine.get(q) != theirs.get(q):
                bad.append((t, q, "only in independent" if q not in theirs else ("only in cluster" if q not in mine else "values differ")))
                if q not in theirs:
                    onlyind += 1
                elif q not in mine:
                    onlypkg += 1
            if q in cancel:
                cancel_hit.add(q)
        d = per_E.setdefault(E[t], [0, 0, 0.0])
        d[0] += 1
        d[1] += len(qs)
        d[2] += sec
    say("pieces compared: %d (%d individual residues); pieces present in only one of the two codes: independent %d, cluster %d"
        % (npieces, nentries, onlyind, onlypkg))
    say("MISMATCHES: %d%s" % (len(bad), "" if not bad else " -> " + str(bad[:20])))
    say("distinct outer graphs reached: %d, of which in the cancellation set (1 <= E <= 7, non-zero pieces): %d of %d"
        % (len(set(q for t in nat for q in (set(nat[t][0]) | set(ref[t])))), len(cancel_hit), len(cancel)))
    say("by number of edges of the inner graph:")
    say("   %3s %7s %9s %12s" % ("E", "graphs", "pieces", "CPU s"))
    for e in sorted(per_E):
        say("   %3d %7d %9d %12.0f" % (e, per_E[e][0], per_E[e][1], per_E[e][2]))
    say("   total %5d %9d %12.0f" % (sum(v[0] for v in per_E.values()), npieces, sum(v[2] for v in per_E.values())))
    # ---- point frame and additivity
    add_bad = []
    idp_bad = []
    nid = 0
    if idd:
        shipped = json.load(open(os.path.join(EV, "parts", "g5_%s_id" % a.prime, "local.json")))["pieces"]
        for t in sorted(idd):
            mine_id = idd[t][0]
            rows = {q: v for q, v in mine_id.items()}
            exp = {int(q): [x % p for x in v] for q, v in shipped[str(t)].items() if any(x % p for x in v)}
            if rows != exp:
                idp_bad.append(t)
            nid += 1
            if t in nat:
                s = [0] * 8
                for q, v in nat[t][0].items():
                    for i, x in enumerate(v):
                        s[i] = (s[i] + x) % p
                own = rows.get(t, [0] * 8)
                if s != own:
                    add_bad.append(t)
        say("point frame: %d independent id terms compared with parts/g5_%s_id/local.json; mismatches %d%s"
            % (nid, a.prime, len(idp_bad), "" if not idp_bad else " -> %s" % idp_bad[:20]))
        say("additivity inside the independent code (Sum_q nat pieces == own id term): %d graphs checked, %d failures%s"
            % (sum(1 for t in idd if t in nat), len(add_bad), "" if not add_bad else " -> %s" % add_bad[:20]))
    brutes = [d for lst in load_jsonl(tag).values() for d in lst if d.get("brute_equal") is not None]
    if brutes:
        say("colour-space brute force inside the independent code: %d graphs, all equal: %s"
            % (len(brutes), all(d["brute_equal"] for d in brutes)))
    ok = not bad and not idp_bad and not add_bad
    say("VERDICT prime %s: %s" % (a.prime, "ALL EQUAL" if ok else "FAILURES PRESENT"))
    if a.save:
        json.dump(dict(prime=a.prime, graphs=sorted(nat), npieces=npieces, nentries=nentries, per_E={str(k): v for k, v in per_E.items()},
                       cancel_hit=sorted(cancel_hit), bad=bad, idp_bad=idp_bad, add_bad=add_bad, nid=nid),
                  open(os.path.join(OUT, a.save), "w"))
        say("wrote %s" % os.path.join(OUT, a.save))


if __name__ == "__main__":
    main()
