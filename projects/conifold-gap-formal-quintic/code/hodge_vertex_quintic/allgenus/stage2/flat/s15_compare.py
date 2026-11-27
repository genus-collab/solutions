# -*- coding: utf-8 -*-
"""s15_compare.py -- graph-by-graph comparison of the FLAT P-action sum (s15_flat.py) with the HATTED outer sum (stored project tables
and s15_hat.py), genus 3 (and genus 2 against rf_g2flat.log), modulo one prime.  Pure table arithmetic (no graph sums).

Checks (index k of a genus-3 table = coefficient of u^(k-4)):
  H1  s15_hat M0-mode  = stored g3_state_<p>.json A_M0           (all stored orders)   -- validates the hatted driver
  H2  s15_hat nat-mode = stored g3_nat_<p>_nat.json A3, k <= 3  (the valid part: the stored pass used M mod u^4)
  H3  s15_hat nat K=6 and K=8 agree on k <= 5 (precision of the rows)
  F1  flat = stored A3, k <= 3, all 42 graphs                   -- s.1.5 against the record (polar part)
  F2  flat = s15_hat nat, all 42 graphs, all common orders      -- s.1.5 beyond the polar part (the main test)
  F3  sum over graphs of flat = stored point-frame total total_id (M-independence of the total; k <= 5)
  G2  flat genus 2 = rf_g2flat.log canonical frame (0,0) (7 graphs, u^-2..u^1, rationals reduced mod p)
Controls (MUST mismatch): X1 flat vs stored A_M0 (the split of the constant-M_0 outer sum; totals agree, graphs do not);
  X2-X5 flat variants nodil / noT / nosign / stored vs hatted nat.
Run: py -I -B s15_compare.py <89|61> [frame tag, e.g. _ab1o3_m2o7]     Log: logs/s15_compare_<p><tag>.log (overwritten per run).
"""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
from s15_flat import SCRATCH

G2_RECORD = {1: ['1/1200', '-7/6000', '287/864000', '936083/193914000000'], 2: ['0', '0', '479/756000', '-601319/4847850000'],
             3: ['0', '0', '1/2250', '571/2126250'], 4: ['0', '0', '-127/2016000', '979/115425000'], 5: ['0', '0', '0', '-1397/5670000'],
             6: ['0', '0', '0', '-44/151875'], 7: ['0', '0', '0', '0']}       # rf_g2flat.log lines 5-11 (frame (0,0), canonical gauge)


def load(fn):
    fn = os.path.join(SCRATCH, fn)
    if not os.path.exists(fn): return None
    return {int(k): v for k, v in json.load(open(fn))['A'].items()}


def main():
    pr = sys.argv[1]; p = (1 << 89) - 1 if pr == "89" else (1 << 61) - 1
    ft = sys.argv[2] if len(sys.argv) > 2 else ""          # conifold-frame tag, e.g. _ab1o3_m2o7 (stored tables exist at 2^89-1 only)
    LOG = open(os.path.join(HERE, "logs", "s15_compare_%s%s.log" % (pr, ft)), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    res = {}
    def check(name, ok, detail="", control=False):
        tag = ("PASS " if ok else "FAIL ") if not control else ("CONTROL FAILS AS IT MUST " if not ok else "CONTROL DID NOT FAIL ")
        log(tag + name + ((": " + detail) if detail else "")); res[name] = ok
    def cmp(A, B, ks, idxs=None):
        """(number of compared coefficients, list of (graph, k) mismatches)"""
        idxs = sorted(set(A) & set(B)) if idxs is None else idxs; n = 0; bad = []
        for i in idxs:
            for k in ks:
                a, b = A[i][k], B[i][k]
                if a is None or b is None: continue
                n += 1
                if (a - b) % p: bad.append((i, k))
        return n, bad
    log("s15_compare: p = 2^%s-1, conifold frame tag '%s', genus-3 index k = coefficient of u^(k-4)" % (pr, ft))
    nat = json.load(open(os.path.join(PROJ, "g3_nat_%s_nat%s.json" % (pr, ft)))); A3 = {int(k): v for k, v in nat['A3'].items()}
    st = json.load(open(os.path.join(PROJ, "g3_state_%s%s.json" % (pr, ft)))); AM0 = {int(k): v for k, v in st['A_M0'].items()}
    flat = load("flat_g3_%s_can%s.json" % (pr, ft)); hM0 = load("hat_g3_%s_M0_K6%s.json" % (pr, ft))
    h6 = load("hat_g3_%s_nat_K6%s.json" % (pr, ft)); h8 = load("hat_g3_%s_nat_K8%s.json" % (pr, ft))
    log("tables: flat %s, hat M0 K6 %s, hat nat K6 %s, hat nat K8 %s (graphs / orders)" % tuple(
        ("%d/%d" % (len(T), len(next(iter(T.values())))) if T else "missing") for T in (flat, hM0, h6, h8)))
    if hM0:
        n, bad = cmp(hM0, AM0, range(6)); check("H1 hatted driver (M0 mode) = stored A_M0", not bad, "%d coefficients, mismatches %s" % (n, bad[:8]))
    hat = h8 or h6
    for nm, H in (("K6", h6), ("K8", h8)):
        if H:
            n, bad = cmp(H, A3, range(4)); check("H2 hatted nat %s = stored A3 on k <= 3" % nm, not bad, "%d coefficients, mismatches %s" % (n, bad[:8]))
    if h6 and h8:
        n, bad = cmp(h6, h8, range(6)); check("H3 hatted nat K6 = K8 on k <= 5", not bad, "%d coefficients, mismatches %s" % (n, bad[:8]))
    if flat:
        n, bad = cmp(flat, A3, range(4)); check("F1 flat = stored hatted record A3 (k <= 3, all 42 graphs)", not bad, "%d coefficients, mismatches %s" % (n, bad[:8]))
        if hat:
            K = min(len(next(iter(flat.values()))), len(next(iter(hat.values()))))
            n, bad = cmp(flat, hat, range(K))
            nz = sum(1 for i in flat for k in range(4, K) if flat[i][k] % p)
            check("F2 flat = hatted nat-form outer sum, all 42 graphs, k = 0..%d (u^-4..u^%d)" % (K - 1, K - 5), not bad,
                  "%d coefficients (%d non-zero regular coefficients k >= 4), mismatches %s" % (n, nz, bad[:8]))
        tot = [sum(flat[i][k] for i in flat) % p for k in range(6)]
        check("F3 flat total = stored point-frame total total_id (k <= 5)", tot == [x % p for x in st['total_id']][:6])
        n, bad = cmp(flat, AM0, range(6))
        check("X1 flat vs the constant-M_0 split A_M0 (per graph)", not bad, "%d coefficients, %d mismatches, e.g. %s" % (n, len(bad), bad[:6]), control=True)
        totM0 = [sum(AM0[i][k] for i in AM0) % p for k in range(6)]
        check("X1' totals of the two splits agree (k <= 5)", totM0 == tot)
    for var, nm in (("nodil", "X2 dilaton factor dropped"), ("noT", "X3 leaves removed"), ("nosign", "X4 edge without (-1)^t"),
                     ("stored", "X5 stored gauge normalisation P_can c(z) instead of P_can")):
        V = load("flat_g3_%s_%s%s.json" % (pr, var, ft))
        if V and hat:
            K = min(len(next(iter(V.values()))), len(next(iter(hat.values()))))
            n, bad = cmp(V, hat, range(K)); gs = sorted(set(i for i, k in bad))
            check(nm + " vs hatted nat", not bad, "%d coefficients, %d mismatches in graphs %s" % (n, len(bad), gs), control=True)
    g2 = load("flat_g2_%s_can.json" % pr if not ft else "none")
    if g2:
        R = {i: [Fr(x) for x in v] for i, v in G2_RECORD.items()}
        bad = [(i, k) for i in R for k in range(4) if (g2[i][k] - R[i][k].numerator * pow(R[i][k].denominator, -1, p)) % p]
        check("G2 flat genus 2 (own code) = rf_g2flat.log record, 7 graphs x 4 orders", not bad, str(bad))
    log("SUMMARY: %s" % ("ALL CHECKS PASSED, ALL CONTROLS FAILED AS THEY MUST" if all(
        v if not k.startswith("X") or k.startswith("X1'") else not v for k, v in res.items()) else "SEE ABOVE"))
    LOG.close()


if __name__ == "__main__":
    main()
