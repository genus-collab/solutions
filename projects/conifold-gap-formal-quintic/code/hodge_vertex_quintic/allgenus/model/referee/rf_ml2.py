# -*- coding: utf-8 -*-
"""rf_ml2.py -- (ML) test for large n (beyond the worker's range): generator tuples with mH copies of tau_0(H) and at most QMAX other
generators.  The amputated amplitude is first contracted with tau_0(H) on mH legs (legs are symmetric), then the remaining legs are
contracted with all flat insertions (rf_tables method).  Same generators, zero test at both primes, same wrong-sign control as rf_ml.py.
Run: py -I -B rf_ml2.py h n QMAX [--one]   (--one: zero test at 2^61-1 only, when the 2^89-1 table was not computed)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat as rf
import rf_gsum as gsm
import rf_ml as ml


def table_after_H(m, h, n, mH):
    p = m.p; A = gsm.amp(m, h, n)
    L = [[[sum(m.Rinv[j][c][b] * m.flat[a][b] for b in range(5)) % p for c in range(5)] for j in range(m.K + 1)] for a in range(5)]
    cur = {}
    for (cols, exps), val in A.items():
        x = val
        for l in range(mH):
            x = x * L[1][exps[l]][cols[l]] % p
            if not x: break
        if x:
            key = ((), cols[mH:], exps[mH:]); cur[key] = (cur.get(key, 0) + x) % p
    for l in range(n - mH):
        nxt = {}
        for (done, cols, exps), val in cur.items():
            c = cols[0]; d = exps[0]
            for a in range(5):
                for j in range(d + 1):
                    x = L[a][j][c]
                    if x:
                        k2 = (done + ((a, d - j),), cols[1:], exps[1:]); nxt[k2] = (nxt.get(k2, 0) + val * x) % p
        cur = nxt
    return {k[0]: v for k, v in cur.items() if v}


def main():
    h, n, QMAX = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); t0 = time.time()
    m1, m2 = rf.models(); MS = (m1,) if "--one" in sys.argv else (m1, m2)
    for variant in ("kappa", "wrongsign"):
        G = ml.gens(3 * h + n + 8, variant); tabs = {}
        nt = nb = 0; ex = []
        for combo in ml.singular_tuples(G, h, n):
            gl = [G[i] for i in combo]; lead = [g[2] for g in gl]
            if not rf.charge_ok(lead): continue
            others = [g for g in gl if g[0] != "H"]; mH = n - len(others)
            if not others or len(others) > QMAX: continue
            if mH not in tabs: tabs[mH] = [(m, table_after_H(m, h, n, mH)) for m in MS]
            # value: the H legs are already contracted; expand the others (dimension pruning inside ml.value uses D of the full type)
            vals = []
            for m, tab in tabs[mH]:
                D = 3 * h - 3 + n; comps = [[(m.fr(cf), ins) for cf, ins in g[1]] for g in others]
                tot = 0; stack = [(0, 1, (), 0)]
                while stack:
                    i, c, ins, deg = stack.pop()
                    if i == len(comps):
                        x = tab.get(ins)
                        if x: tot = (tot + c * x) % m.p
                        continue
                    for cf, x in comps[i]:
                        if deg + x[1] <= D: stack.append((i + 1, c * cf % m.p, ins + (x,), deg + x[1]))
                vals.append(tot)
            nt += 1
            if any(vals):
                nb += 1
                if len(ex) < 2: ex.append(" ".join("H" if g[0] == "H" else "k(H%d,%d)" % (g[0][1], g[0][2]) for g in gl))
        print("(%d,%d) %-9s: %5d singular generator tuples with 1..%d kappas, %5d NON-ZERO %s" % (h, n, variant, nt, QMAX, nb,
              ("e.g. " + "; ".join(ex)) if ex else ""), flush=True)
    print("rf_ml2 (%d,%d) done, primes: %s (%.0f s)" % (h, n, "2^61-1 only" if len(MS) == 1 else "2^61-1 and 2^89-1", time.time() - t0))


if __name__ == "__main__":
    main()
