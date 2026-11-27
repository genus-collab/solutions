# -*- coding: utf-8 -*-
"""ur_euler_check.py -- UNIFORM_REDUCTION.md, check (D): the Euler-dilaton rank-one elements are invisible to the frozen vertex theory.

In hatted variables the 'H-part' of the regular gauge is, to leading order, alpha u^2 n(zhat) with
        n(zhat) = zhat (H - zhat 1) (x) (H + zhat 1) = zhat H(x)H + zhat^2 (H(x)1 - 1(x)H) - zhat^3 1(x)1
(Fourier components Xhat(3,3) = 1 at zhat; Xhat(3,4) = 1, Xhat(4,3) = -1 at zhat^2; Xhat(4,4) = -1 at zhat^3; this is the pattern of
[u^2] G with alpha = -127/15750).  Claim tested here, in the FROZEN theory (u = 0, rows R^, vertex CohFT M.omega^):
   for M = exp(eps f(zhat) n_0(zhat)) M_0  (LEFT dressing; n_0 = (H - zhat)(x)(H + zhat), f odd: f = zhat or zhat^3) every stable graph
   with edges vanishes for every eps and the vertex graph keeps its value (1/1200 at genus 2, 1/25200 at genus 3),
i.e. R' = Psi_0 exp(-eps f n_0) is invisible: f n_0 is a product of two insertions of the Euler-dilaton vector e(z) = z 1 -+ H, which
annihilates the homogeneous model theory (dilaton equation + homogeneity in the one variable t).  Controls: the opposite relative sign
n_0' = (H + zhat)(x)(H - zhat), and the bare propagator zhat H(x)H, must NOT be invisible.
OUTCOME (2026-11-10, 2^89-1): the claim is REFUTED for the covariant element f = zhat: n is NOT an exact symmetry (all genus-2 edge graphs
respond; the genus-1 vertices carry the conformal anomaly <H,H>_1 - <tau_1(1)>_1 = 1/12 - 5/24).  The case f = zhat^3 passes only
trivially (it is not Z/5-covariant, so the colour sums vanish by the Z/5 selection rule).  What survives is the FIRST-ORDER statement on
one-edge graphs with a vertex of genus >= 2, see ur_euler_linear.py.  The two FAIL lines of the log are this negative result.
Graphs: genus 2: all 7; genus 3: the vertex graph is skipped (heavy), the one-edge graphs 2 and 5 and (optionally) two-edge graphs.
Run: py -I ur_euler_check.py [89|61] [g3 graph indices ...]      (single process)
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_par import setup, make_M
from fq3_m import in_sp
from fq3_msum import frozen_sums6
from fq_msum import outer_graphs
from fq3_rec import Rec
FAILED = []


def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else "89"
    g3 = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else [2, 5]
    OUT = open(os.path.join(HERE, "ur_euler_check_%s.log" % tag), "w", encoding="utf-8")

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()

    def check(name, ok, detail=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
        if not ok: FAILED.append(name)

    p = P2 if tag == "89" else P1
    st = setup((p, "0", "0", 34, 6)); R = st['R']; rec = Rec(R); z = R.zeta

    def from_fourier(comp):
        """X_ij = Sum Xhat(m, m') zeta^(-m i) zeta^(-m' j)"""
        return [[sum(c * pow(z, (-m * i - mp * j) % 5, p) for (m, mp), c in comp.items()) % p for j in range(5)] for i in range(5)]

    Z = [[0] * 5 for _ in range(5)]
    HH = from_fourier({(3, 3): 1}); H1 = from_fourier({(3, 4): 1, (4, 3): -1}); II = from_fourier({(4, 4): -1})
    H1m = from_fourier({(3, 4): -1, (4, 3): 1})
    cases = {
        "n = zhat (H - zhat)(x)(H + zhat)": [Z, HH, H1, II, Z, Z, Z],
        "zhat^3 (H - zhat)(x)(H + zhat)": [Z, Z, Z, HH, H1, II, Z],
        "CONTROL opposite sign zhat (H + zhat)(x)(H - zhat)": [Z, HH, H1m, II, Z, Z, Z],
        "CONTROL bare zhat H(x)H": [Z, HH, Z, Z, Z, Z, Z],
    }
    g2 = outer_graphs(); g3graphs = [st['graphs'][i] for i in g3]
    ref2 = frozen_sums6(st['frozen'], st['M0'], g2)
    log("p = 2^%s-1; frozen genus-2 sums with M_0: %s" % (tag, {i: rec.fmt(v) for i, v in ref2.items()}))
    t0 = time.time(); ref3 = frozen_sums6(st['frozen'], st['M0'], g3graphs)
    log("frozen genus-3 sums with M_0, graphs %s: %s (%.0f s)" % (g3, {i: rec.fmt(v) for i, v in ref3.items()}, time.time() - t0))
    for name, N in cases.items():
        log("-" * 100); log(name + ("   [in sp: %s]" % in_sp(N, p)))
        ok2 = True; ok3 = True
        for eps in (1, 3, R.of(Fr(-2, 7))):
            M, kind = make_M(st, ('eps_M0', N, eps))
            r2 = frozen_sums6(st['frozen'], M, g2)
            good2 = all(r2[i] == 0 for i in range(2, 8)) and r2[1] == R.of(Fr(1, 1200))
            ok2 = ok2 and good2
            t0 = time.time(); r3 = frozen_sums6(st['frozen'], M, g3graphs)
            good3 = all(v == 0 for v in r3.values())
            ok3 = ok3 and good3
            log("   eps = %s: genus 2 %s; genus 3 %s (%.0f s)" % (rec.fmt(eps), {i: rec.fmt(v) for i, v in r2.items()}, {i: rec.fmt(v) for i, v in r3.items()}, time.time() - t0))
        if name.startswith("CONTROL"):
            check("control is visible (some edge graph non-zero): " + name, (not ok2) or (not ok3))
        else:
            check("INVISIBLE at genus 2 (six edge graphs = 0, vertex graph = 1/1200, three values of eps): " + name, ok2)
            check("INVISIBLE on the genus-3 graphs %s (three values of eps): " % g3 + name, ok3)
    log("UR_EULER_CHECK (%s) " % tag + ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED)); OUT.close()


if __name__ == "__main__":
    main()
