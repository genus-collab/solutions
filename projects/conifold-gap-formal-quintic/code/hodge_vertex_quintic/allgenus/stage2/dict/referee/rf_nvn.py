"""Referee: direct check of the target of Lemma 3 (DICT_PROOF.md s.5) by linear algebra in the polynomial ring Q[p, q].
Usage: py -I -B rf_nvn.py [WMAX]
  V_w := Q[v,x]_w + (pq Q[p,q])_w inside Q[p,q]_w   (v = pq, x = p^4/4 + q; monomial basis p^a q^c, a + 4c = w).
  Checks: dim Q[p,q]_w / V_w = [w >= 1]  (so (N/vN)_w = that (+) log p Q[x]_w has dimension [w >= 1] + [4|w]);
          p^w is not in V_w (w >= 1);  q^n + p^(4n)/4^n lies in V_(4n)  (q = -p^4/4 mod V);
          Control (must fail): the wrong sign, q^n - p^(4n)/4^n in V_(4n).
"""
import os
import sys
from fractions import Fraction as Fr
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf_core import Ar, rank  # noqa: E402

ar = Ar(None)
WMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200


def vx_elements(w, sx=1):
    """v^m x^n (5m + 4n = w) as dicts {(a, c): coeff} in Q[p,q]; x = p^4/4 + sx q."""
    out = []
    for n in range(w // 4 + 1):
        if (w - 4 * n) % 5:
            continue
        m = (w - 4 * n) // 5
        out.append({(4 * (n - t) + m, t + m): Fr(comb(n, t) * sx ** t, 4 ** (n - t)) for t in range(n + 1)})
    return out


def V_rows(w):
    rows = vx_elements(w)
    rows += [{(a + 1, c + 1): Fr(1)} for c in range((w - 5) // 4 + 1) for a in [w - 5 - 4 * c] if a >= 0] if w >= 5 else []
    return rows


if __name__ == "__main__":
    bad = []
    for w in range(WMAX + 1):
        V = V_rows(w)
        rV = rank(ar, V)
        dimQpq = w // 4 + 1
        ok = dimQpq - rV == (1 if w >= 1 else 0)
        if w >= 1:
            ok = ok and rank(ar, V + [{(w, 0): Fr(1)}]) == rV + 1
        if w % 4 == 0 and w >= 4:
            n = w // 4
            ok = ok and rank(ar, V + [{(0, n): Fr(1), (w, 0): Fr(1, 4 ** n)}]) == rV
        if not ok:
            bad.append(w)
    print("rf_nvn.py: dim Q[p,q]_w/V_w = [w>=1], p^w not in V_w, q^n = -p^(4n)/4^n mod V: w <= %d, failing weights %s"
          % (WMAX, bad))
    ctrl = [w for w in range(4, WMAX + 1, 4)
            if rank(ar, V_rows(w) + [{(0, w // 4): Fr(1), (w, 0): Fr(-1, 4 ** (w // 4))}]) == rank(ar, V_rows(w))]
    print("control (must fail): the wrong sign q^n = +p^(4n)/4^n mod V holds at %d of %d weights 4n <= %d (expected 0)"
          % (len(ctrl), WMAX // 4, WMAX))
    if ctrl:
        bad.append("control")
    print("RESULT:", "ALL CHECKS PASSED" if not bad else "FAILURES at %s" % bad[:10])
