# -*- coding: utf-8 -*-
"""gi_diag_old.py -- reproduces the (G*) failure of the FIRST version of gi_relations.py with its own completeness bound (e + 2n + 2 <= K)
and its sign (-1) ** m, and shows that it disappears with an integer sign (GAUGE_IDENTITY.md s.8).  Frame (0, 0), K = 20, n <= 8.
Run: py -I -B gi_diag_old.py      log: gi_diag_old.log
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import gauge_frame, model_vanishing, quintic_vanishing, p11_at_zero, inv_poly, normalise
PAIR = {0: 3, 1: 2, 2: 1, 3: 0, 4: 4}
K, NU, NMAX = 20, 66, 8
OUT = open(os.path.join(HERE, "gi_diag_old.log"), "w", encoding="ascii")
def log(s): print(s, flush=True); OUT.write(s + "\n")
res = gauge_frame(K, NU, 0, 0); p = res['p']; v = res['v']; bq = res['bq']
pc = normalise(p, inv_poly(p11_at_zero(p), K))
P_ = []
for n in range(NMAX + 1):
    out = {}
    for k in range(K + 1):
        for b in range(5):
            c = pc[k][b].coef(n)
            if c != 0: out.setdefault(((b + k) % 5, b), {})[k] = c
    P_.append(out)
YC = model_vanishing(v, NMAX + 1, 2 * NMAX + 4); YQ = quintic_vanishing(bq, NMAX + 1, 2 * NMAX + 4)
def ucoef(Y, n):
    out = [dict() for _ in range(5)]
    for m, vec in Y.items():
        for a in range(5):
            c = vec[a].coef(n) if n < vec[a].prec else None
            if c: out[a][m] = c
    return out
eta = [ucoef(YC, n) for n in range(NMAX + 1)]; theta = [ucoef(YQ, n) for n in range(NMAX + 1)]
ok_all = True
for label, sign in (("first version: (-1) ** m (a float for m < 0)", lambda m: (-1) ** m), ("integer sign", lambda m: 1 if m % 2 == 0 else -1)):
    star = lambda vec: [{m: 5 * c * sign(m) for m, c in vec[PAIR[a]].items()} for a in range(5)]
    tot = 0; bad = []
    for n in range(NMAX + 1):
        acc = [dict() for _ in range(5)]
        for j in range(n + 1):
            ts = star(theta[n - j])
            for (a, b), pz in P_[j].items():
                for k, c in pz.items():
                    for m, y in ts[a].items(): acc[b][k + m] = acc[b].get(k + m, Fr(0)) + c * y
        es = star(eta[n]); mneg = 2 * n + 2
        for b in range(5):
            for e in set(acc[b]) | set(es[b]):
                if e + mneg > K: continue
                tot += 1; l = acc[b].get(e, 0); r = es[b].get(e, 0)
                if l != r: bad.append((n, b, e, type(l).__name__, abs(float(l) - float(r))))
    res_ = [x[4] for x in bad]
    log("%s: %d identities tested, %d failing; types %s; |lhs - rhs| in [%s, %s]; first four %s" % (label, tot, len(bad), sorted(set(x[3] for x in bad)),
        "%.1e" % min(res_) if res_ else "-", "%.1e" % max(res_) if res_ else "-", [x[:3] for x in bad[:4]]))
    if label.startswith("integer"): ok_all = ok_all and not bad
    else: ok_all = ok_all and bool(bad) and all(x[3] == 'float' for x in bad) and max(res_) < 1e-12
log("GI_DIAG_OLD: %s" % ("DIAGNOSIS CONFIRMED (float rounding only; exact identities hold)" if ok_all else "NOT CONFIRMED"))
OUT.close()
