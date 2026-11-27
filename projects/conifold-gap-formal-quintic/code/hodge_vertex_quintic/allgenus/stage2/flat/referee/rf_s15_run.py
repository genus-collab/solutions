# -*- coding: utf-8 -*-
"""rf_s15_run.py -- referee driver: per-graph flat P-action graph sum at genus g (own engine rf_s15_sum + rf_s15_outer).
Run: py -I -B rf_s15_run.py <89|61> <g> <N> <variant> [a,b]
  N = series length: coefficients u^(2-2g) .. u^(N+1-2g) are exact.  variant in can, stored, nodil, noT, nosign, zminus
  (frozen R(-z)), dress (P -> exp(u z N) P, N = E_04 + E_43: symplectic, graded, flat-regular, NOT the gauge; the E3-type
  dressing of theory/referee/rf_g2flat.py, used here only as a failing control).
Graph labels: genus 3 in the order of fq3_graphs.outer_graphs3 (read-only import, used ONLY for the index and |Aut| cross-check);
genus 2 in the order of theory/referee/rf_g2flat.log (fixed by hand below).
Output: <scratch>/rf_g<g>_<p>_<variant>[_ab].json, log referee/logs/rf_s15_run_g<g>_<p>_<variant>[_ab].log.
"""
import os, sys, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
from rf_s15_base import P61, P89, Frozen, ratrec, canon2, smul
from rf_s15_sum import Engine
from rf_s15_outer import Outer, load_gauge
SCR = os.path.join(HERE, "scratch")
G2 = [((2,), (), ()), ((1,), ((0, 0),), ()), ((1, 1), ((0, 1),), ()), ((0,), ((0, 0), (0, 0)), ()),
      ((0, 1), ((0, 0), (0, 1)), ()), ((0, 0), ((0, 1), (0, 1), (0, 1)), ()), ((0, 0), ((0, 0), (0, 1), (1, 1)), ())]
G2NAME = ["vertex g=2", "loop on g=1", "bridge g1-g1", "two loops on g=0", "g0(loop)-g1", "g0=g0 triple edge", "g0(loop)-g0(loop)"]


def labels(g):
    if g == 2: return [(i, canon2(h, list(e), ()), G2NAME[i], None) for i, (h, e, _) in enumerate(G2)]
    sys.path.insert(0, PROJ)
    from fq3_graphs import outer_graphs3
    return [(G['idx'], canon2(tuple(G['gv']), [tuple(sorted(x)) for x in G['ed']], ()), "E=%d gv=%s" % (len(G['ed']), tuple(G['gv'])), G['aut'])
            for G in outer_graphs3()]


def main():
    pr, g, N, var = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    ab = sys.argv[5].split(",") if len(sys.argv) > 5 else ["0", "0"]
    tag = "" if ab == ["0", "0"] else "_ab" + "_".join(x.replace("/", "o").replace("-", "m") for x in ab)
    p = P89 if pr == "89" else P61
    os.makedirs(SCR, exist_ok=True); os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    LOG = open(os.path.join(HERE, "logs", "rf_s15_run_g%d_%s_%s%s.log" % (g, pr, var, tag)), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    t0 = time.time()
    log("rf_s15_run: genus %d, p = 2^%s-1, N = %d, variant %s, frame %s [%s]" % (g, pr, N, var, ab, time.strftime("%Y-%m-%dT%H:%M:%S")))
    K = 9
    Fz = Frozen(p, K, -1 if var == 'zminus' else 1)
    P, tc = load_gauge(K, N + 2, Fr(ab[0]), Fr(ab[1]))
    if var == 'dress':      # P -> exp(u z N) P, N = E_04 + E_43 (lowered; graded, eta N symmetric, N^3 = 0): symplectic NON-gauge dressing
        def nmul(Mx):       # (N M)[0] = M[4], (N M)[4] = M[3], other rows 0
            out = [[[Fr(0)] * (N + 2) for _ in range(5)] for _ in range(5)]
            out[0] = [list(r) for r in Mx[4]]; out[4] = [list(r) for r in Mx[3]]
            return out
        def ush(Mx, c):     # c u Mx
            return [[[Fr(0)] + [c * x for x in r[:-1]] for r in row] for row in Mx]
        P2 = {}
        for k in range(K + 1):
            M = [[list(r) for r in row] for row in P[k]]
            if k >= 1:
                A = ush(nmul(P[k - 1]), Fr(1)); M = [[[a + b for a, b in zip(M[i][j], A[i][j])] for j in range(5)] for i in range(5)]
            if k >= 2:
                A = ush(ush(nmul(nmul(P[k - 2])), Fr(1, 2)), Fr(1))
                M = [[[a + b for a, b in zip(M[i][j], A[i][j])] for j in range(5)] for i in range(5)]
            P2[k] = M
        P = P2
    O = Outer(Fz, P, N, 'can' if var in ('zminus', 'dress') else var)
    if hasattr(O, 'cz'): log("  c(z) = P_11(0,z) of the stored gauge: %s" % {k: str(v) for k, v in O.cz.items()})
    log("  d = %s ..; v = u (%s ..)" % ([str(ratrec(x, p)) for x in O.d[:4]], [str(ratrec(x, p)) for x in O.V[:4]]))
    Eo, XO = O.to_idempotent(Fz, 7)
    log("  outer edges %d entries, outer leaves %d entries; setup %.0f s" % (len(Eo), len(XO), time.time() - t0))
    E = Engine(Fz, N, XO, Eo)
    res, info = E.outer_sums(g)
    off, gf = O.global_factor(g)
    out = {}
    for idx, key, name, aut in labels(g):
        s = smul(res.get(key, [0] * N), gf, p, N); out[idx] = s
        log("  graph %2s (%s; %d composite classes): %s" % (idx, name, info.get(key, 0), [str(ratrec(x, p)) for x in s[:8]]))
    assert set(res) <= set(k for _, k, _, _ in labels(g)), "unlabelled outer graph"
    full = {}
    for idx, key, name, aut in labels(g):
        full[idx] = smul(E.full.get(key, [0] * N), gf, p, N)
    log("  point-frame regrouping (by uncontracted graph), first graphs: %s" % {i: [str(ratrec(x, p)) for x in full[i][:3]] for i in list(full)[:3]})
    tot = [sum(out[i][k] for i in out) % p for k in range(N)]
    log("  TOTAL: %s" % [str(ratrec(x, p)) for x in tot])
    tcm = [x.numerator % p * pow(x.denominator, p - 2, p) % p for x in tc[:N + 1]]
    T = tcm[1:N + 1]                                     # t_c = u T
    Ti = [1] + [0] * (N - 1)
    from rf_s15_base import sinv
    Tm = sinv(T, p, N); Tg = [1] + [0] * (N - 1)
    for _ in range(2 * g - 2): Tg = smul(Tg, Tm, p, N)
    log("  t_c^(2-2g) / u^(2-2g): %s" % [str(ratrec(x, p)) for x in Tg[:6]])
    json.dump(dict(p=p, g=g, N=N, offset=off, variant=var, frame=ab, A={str(k): v for k, v in out.items()}, tcpow=Tg,
                   Afull={str(k): v for k, v in full.items()}),
              open(os.path.join(SCR, "rf_g%d_%s_%s%s.json" % (g, pr, var, tag)), "w"))
    log("  done, %.0f s" % (time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
