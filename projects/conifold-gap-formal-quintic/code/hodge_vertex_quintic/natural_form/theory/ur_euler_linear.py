# -*- coding: utf-8 -*-
"""ur_euler_linear.py -- UNIFORM_REDUCTION.md, check (D'): FIRST-ORDER response of the frozen edge graphs to the Euler-dilaton rank-one
element n(zhat) = zhat (H - zhat)(x)(H + zhat) (left dressing M = exp(eps n) M_0), by exact interpolation in eps (degree <= 3 at genus 2,
<= 6 at genus 3).  Prediction (string + dilaton + homogeneity of the primary sector F_g = c_g t^(2-2g) of the frozen vertex theory):
the eps^1 coefficient vanishes on a one-edge graph whenever a vertex of genus >= 2 carries a half-edge of the edge (genus 3: graphs 2 and
5), and does NOT vanish when only genus-1 vertices do (genus 2: graphs 2 and 4; the genus-1 anomaly <H,H>_1 - <tau_1(1)>_1 != 0).
Control: the opposite relative sign zhat (H + zhat)(x)(H - zhat).
Run: py -I ur_euler_linear.py [89|61]      (single process, about a minute)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_par import setup, make_M
from fq3_msum import frozen_sums6
from fq_msum import outer_graphs
from fq_hat import vander_inv, interp_fp
from fq3_rec import Rec


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "89"
    OUT = open(os.path.join(HERE, "ur_euler_linear_%s.log" % tag), "w", encoding="utf-8")

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    p = P2 if tag == "89" else P1
    st = setup((p, "0", "0", 34, 6)); R = st['R']; rec = Rec(R); z = R.zeta

    def from_fourier(comp):
        return [[sum(c * pow(z, (-m * i - mp * j) % 5, p) for (m, mp), c in comp.items()) % p for j in range(5)] for i in range(5)]

    Z = [[0] * 5 for _ in range(5)]
    HH = from_fourier({(3, 3): 1}); H1 = from_fourier({(3, 4): 1, (4, 3): -1}); II = from_fourier({(4, 4): -1})
    H1m = from_fourier({(3, 4): -1, (4, 3): 1})
    cases = {"n = zhat (H - zhat)(x)(H + zhat)": [Z, HH, H1, II, Z, Z, Z],
             "CONTROL zhat (H + zhat)(x)(H - zhat)": [Z, HH, H1m, II, Z, Z, Z],
             "CONTROL bare zhat H(x)H": [Z, HH, Z, Z, Z, Z, Z]}
    g2 = outer_graphs(); g3 = [st['graphs'][i] for i in (2, 5)]
    nodes = list(range(7)); Vinv = vander_inv(nodes, p)
    log("p = 2^%s-1; eps-nodes %s" % (tag, nodes)); ok = True
    for name, N in cases.items():
        t0 = time.time(); v2 = {G['idx']: [] for G in g2}; v3 = {G['idx']: [] for G in g3}
        for eps in nodes:
            M, kind = make_M(st, ('eps_M0', N, eps))
            r2 = frozen_sums6(st['frozen'], M, g2); r3 = frozen_sums6(st['frozen'], M, g3)
            for i in v2: v2[i].append(r2[i])
            for i in v3: v3[i].append(r3[i])
        c2 = {i: interp_fp(v2[i], Vinv, p) for i in v2}; c3 = {i: interp_fp(v3[i], Vinv, p) for i in v3}
        log("-" * 100); log("%s   (%.0f s)" % (name, time.time() - t0))
        log("   genus 2, eps^1 coefficients of graphs 1..7: %s" % {i: rec.fmt(c2[i][1]) for i in sorted(c2)})
        log("   genus 2, eps-degree check (coefficients of eps^4..eps^6 vanish): %s" % all(c2[i][k] == 0 for i in c2 for k in (4, 5, 6)))
        log("   genus 3, eps^1 coefficients of the one-edge graphs 2 (loop on g = 2) and 5 (bridge g1 - g2): %s" % {i: rec.fmt(c3[i][1]) for i in sorted(c3)})
        log("   genus 3, eps^2 coefficients: %s" % {i: rec.fmt(c3[i][2]) for i in sorted(c3)})
        if name.startswith("n ="):
            good = c3[2][1] == 0 and c3[5][1] == 0 and (c2[2][1] != 0 or c2[4][1] != 0)
            log(("PASS " if good else "FAIL ") + "first-order response of n: zero on the genus-3 one-edge graphs, non-zero on genus-2 one-edge graphs (genus-1 anomaly)")
            ok = ok and good
        else:
            vis = c3[2][1] != 0 or c3[5][1] != 0
            log(("PASS " if vis else "FAIL ") + "control has a non-zero first-order response on a genus-3 one-edge graph")
            ok = ok and vis
    log("UR_EULER_LINEAR (%s) %s" % (tag, "ALL CHECKS PASSED" if ok else "FAILURES")); OUT.close()


if __name__ == "__main__":
    main()
