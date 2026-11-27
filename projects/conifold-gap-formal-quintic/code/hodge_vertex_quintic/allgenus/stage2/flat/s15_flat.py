# -*- coding: utf-8 -*-
"""s15_flat.py -- driver: the flat P-action graph sum (s15_core.Flat) for every leg-free stable graph of genus g (fq3_graphs.outer_graphs,
Mathematica order = the order of the project's hatted tables), modulo p, written to a JSON table {idx: [u^-(2g-2) .. u^UMAX]} (residues).

Graph value = (1/|Aut|) Sum_{psi-decorations of the half-edges} Sum_{flat colours} prod_v vertex_v prod_e E_e   (s15_core docstring).
Run:  py -I -B s15_flat.py <89|61> <g> <UMAX> [variant|-] [graph indices comma-separated|all] [alpha,beta (conifold frame)]
Output: <scratch>/flat_g<g>_<p>_<variant>[_ab..].json ; log logs/s15_flat_g<g>_<p>_<variant>[_ab..].log (overwritten per run).
"""
import os, sys, json, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from s15_core import Flat, LS, outer_graphs
from fq_fp import P1, P2
SCRATCH = os.environ.get("S15_SCRATCH", os.path.join(HERE, "scratch"))


def tuples_le(n, s):
    if n == 0:
        yield (); return
    for k in range(s + 1):
        for rest in tuples_le(n - 1, s - k): yield (k,) + rest


def loop_contract(F, T, lv, ks, pairs):
    """contract self-loops: pairs = [(i0, i1)] positions in lv of the half-edges 2e, 2e+1"""
    while pairs:
        (i0, i1) = pairs[0]; new = {}
        for key, val in T.items():
            Ee = F.E.get(((ks[i0], key[i0]), (ks[i1], key[i1])))
            if Ee is None: continue
            nk = tuple(x for j, x in enumerate(key) if j not in (i0, i1))
            new[nk] = F.add(new.get(nk), F.mul(val, Ee))
        T = new
        rem = [j for j in range(len(lv)) if j not in (i0, i1)]
        # re-index the remaining pairs
        mp = {j: r for r, j in enumerate(rem)}
        pairs = [(mp[a], mp[b]) for (a, b) in pairs[1:]]
        lv = [lv[j] for j in rem]; ks = [ks[j] for j in rem]
    return T, lv, ks


def graph_value(F, G):
    V = G['V']; gv = G['gv']; ed = G['ed']; hsAt = G['hsAt']
    hv = {}
    for e, (a, b) in enumerate(ed): hv[2 * e] = a; hv[2 * e + 1] = b
    order = [0]; seen = {0}
    while len(order) < V:
        for e, (a, b) in enumerate(ed):
            for x, y in ((a, b), (b, a)):
                if x in seen and y not in seen: order.append(y); seen.add(y)
    opts = [list(tuples_le(len(hsAt[v]), 3 * gv[v] - 3 + len(hsAt[v]))) for v in range(V)]
    total = None
    for dec in itertools.product(*opts):
        kh = {}
        for v in range(V):
            for r, h in enumerate(hsAt[v]): kh[h] = dec[v][r]
        cur = None; opn = []
        for v in order:
            lv = list(hsAt[v]); ks = list(dec[v])
            T = F.vertex(gv[v], tuple(ks))
            if not T: cur = {}; break
            pairs = [(lv.index(2 * e), lv.index(2 * e + 1)) for e in range(len(ed)) if ed[e][0] == v and ed[e][1] == v]
            T, lv, ks = loop_contract(F, T, lv, ks, pairs)
            if cur is None: cur, opn = T, lv; continue
            conn = []          # (position in opn, position in lv, edge key order)
            for i, h in enumerate(opn):
                o = h ^ 1
                if o in lv: conn.append((i, lv.index(o), h % 2 == 0))
            ci = set(i for i, _, _ in conn); cj = set(j for _, j, _ in conn)
            new = {}
            for ka, va in cur.items():
                for kb, vb in T.items():
                    val = F.mul(va, vb); ok = True
                    for (i, j, first) in conn:
                        hi = opn[i]; hj = lv[j]
                        key = ((kh[hi], ka[i]), (kh[hj], kb[j])) if first else ((kh[hj], kb[j]), (kh[hi], ka[i]))
                        Ee = F.E.get(key)
                        if Ee is None: ok = False; break
                        val = F.mul(val, Ee)
                    if not ok: continue
                    nk = tuple(x for i, x in enumerate(ka) if i not in ci) + tuple(x for j, x in enumerate(kb) if j not in cj)
                    new[nk] = F.add(new.get(nk), val)
            cur = new; opn = [h for i, h in enumerate(opn) if i not in ci] + [h for j, h in enumerate(lv) if j not in cj]
            if not cur: break
        if not cur: continue
        assert list(cur) == [()] and not opn, (G['idx'], list(cur)[:3], opn)
        total = F.add(total, cur[()])
    if total is None: return None
    return F.scal(total, pow(G['aut'], -1, F.p))


def main():
    pr = sys.argv[1]; g = int(sys.argv[2]); umax = int(sys.argv[3])
    var = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] != "-" else None
    only = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 and sys.argv[5] != "all" else None
    ab = sys.argv[6].split(",") if len(sys.argv) > 6 else ["0", "0"]
    ftag = "" if ab == ["0", "0"] else "_ab" + "_".join(x.replace("/", "o").replace("-", "m") for x in ab)
    p = P2 if pr == "89" else P1; N = umax + 2 * g - 2 + 1
    os.makedirs(SCRATCH, exist_ok=True)
    LOG = open(os.path.join(HERE, "logs", "s15_flat_g%d_%s_%s%s.log" % (g, pr, var or "can", ftag)), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 100); log("s15_flat: genus %d, p = 2^%s-1, window u^%d..u^%d, variant %s, frame %s [%s]" % (g, pr, 2 - 2 * g, umax, var, ab, time.strftime("%Y-%m-%dT%H:%M:%S")))
    t0 = time.time(); F = Flat(p, N, K=8, a=ab[0], b=ab[1], variant=var, log=log)
    log("  setup %.0f s" % (time.time() - t0))
    out = {}; R = F.R
    for G in outer_graphs(g):
        if only and G['idx'] not in only: continue
        t = time.time(); val = graph_value(F, G)
        co = F.coeffs(val, 2 - 2 * g, umax); out[G['idx']] = co
        fm = [str(R.ratrec(x)) if R.ratrec(x) is not None else "[%d]" % x for x in co]
        log("  graph %2d (E=%d, gv=%s, |Aut|=%d): %s  (%.1f s)" % (G['idx'], len(G['ed']), tuple(G['gv']), G['aut'], fm, time.time() - t))
    tot = [sum(out[i][k] for i in out) % p for k in range(len(next(iter(out.values()))))]
    log("  TOTAL: %s" % [str(R.ratrec(x)) for x in tot])
    log("  teff = T_P (u-series): %s ; t_c: %s" % ([str(x) for x in F.teff[:5]], [str(x) for x in F.tc[:5]]))
    fn = os.path.join(SCRATCH, "flat_g%d_%s_%s%s.json" % (g, pr, var or "can", ftag))
    json.dump(dict(p=p, g=g, umin=2 - 2 * g, umax=umax, variant=var, A={str(k): v for k, v in out.items()}), open(fn, "w"))
    log("  written %s ; total time %.0f s" % (fn, time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
