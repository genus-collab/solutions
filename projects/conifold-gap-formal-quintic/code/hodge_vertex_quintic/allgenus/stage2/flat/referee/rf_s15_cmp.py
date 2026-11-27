# -*- coding: utf-8 -*-
"""rf_s15_cmp.py -- referee comparisons (reads JSON only).  Index k = coefficient of u^(k-4) (genus 3), as in the stored tables.
Sources: own engine tables (scratch allgenus_s2_flat_ref/rf_g*.json); project record g3_nat_<p>_nat.json (A3) and g3_state_<p>.json
(total_id, A_M0); worker tables (scratch allgenus_s2_flat/flat_g3_*.json, hat_g3_*_nat_K8*.json); theory/referee/rf_g2flat.log.
Run: py -I -B rf_s15_cmp.py <89|61> [ab]
"""
import os, sys, json, re
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_s15_base import P61, P89, ratrec
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
MINE = os.path.join(HERE, "scratch"); WORK = os.path.join(HERE, "..", "scratch")


def ld(path):
    return json.load(open(path)) if os.path.exists(path) else None


def mine(g, pr, var, tag=""):
    d = ld(os.path.join(MINE, "rf_g%d_%s_%s%s.json" % (g, pr, var, tag)))
    return None if d is None else {int(k): v for k, v in d['A'].items()}


def cmp(A, B, ks, idx=None):
    bad = []
    for i in (idx or sorted(A)):
        for k in ks:
            a = A[i][k]; b = B[i][k] if B[i] is not None and k < len(B[i]) else None
            if b is None: continue
            if a % P != b % P: bad.append((i, k))
    return bad


def main():
    global P
    pr = sys.argv[1]; P = P89 if pr == "89" else P61
    tag = "_ab1o3_m2o7" if len(sys.argv) > 2 else ""
    print("rf_s15_cmp: p = 2^%s-1, frame tag '%s'; k = coefficient of u^(k-4)" % (pr, tag))
    M = mine(3, pr, "can", tag); assert M is not None
    nat = ld(os.path.join(PROJ, "g3_nat_%s_nat%s.json" % (pr, tag))); st = ld(os.path.join(PROJ, "g3_state_%s%s.json" % (pr, tag)))
    A3 = {int(k): v for k, v in nat['A3'].items()}
    b = cmp(M, A3, range(4)); print("%s R1 own flat = stored A3 (production record), 42 graphs, k <= 3 (168): mismatches %s" % ("PASS" if not b else "FAIL", b[:8]))
    b45 = cmp(M, A3, (4, 5)); print("INFO R1' own flat vs stored A3 at k = 4, 5 (the stored pass used M mod u^4): %d of 84 differ, graphs %s"
                                    % (len(b45), sorted(set(i for i, _ in b45))))
    H = ld(os.path.join(WORK, "hat_g3_%s_nat_K8%s.json" % (pr, tag)))
    if H:
        H = {int(k): v for k, v in H['A'].items()}
        b = cmp(M, H, range(8)); nz = sum(1 for i in M for k in range(4, 8) if M[i][k] % P and i != 1)
        print("%s R2 own flat = worker hatted (M_nat exact, K8), k = 0..7 (336; %d non-zero regular edge-graph coefficients): mismatches %s"
              % ("PASS" if not b else "FAIL", nz, b[:8]))
    W = ld(os.path.join(WORK, "flat_g3_%s_can%s.json" % (pr, tag)))
    if W:
        W = {int(k): v for k, v in W['A'].items()}
        b = cmp(M, W, range(8)); print("%s R3 own flat = worker flat, k = 0..7: mismatches %s" % ("PASS" if not b else "FAIL", b[:8]))
    tot = [sum(M[i][k] for i in M) % P for k in range(6)]
    b = [k for k in range(6) if tot[k] != st['total_id'][k] % P]
    print("%s R4 sum over graphs = stored point-frame total total_id, k <= 5: mismatches %s" % ("PASS" if not b else "FAIL", b))
    ed = [(i, k) for i in M if i != 1 for k in range(4) if M[i][k] % P]
    print("%s R5a all 41 edge graphs regular (k <= 3 vanish): non-zero %s" % ("PASS" if not ed else "FAIL", ed[:8]))
    d = ld(os.path.join(MINE, "rf_g3_%s_can%s.json" % (pr, tag))); c3 = Fr(1, 25200)
    vt = [c3.numerator * pow(c3.denominator, P - 2, P) * x % P for x in d['tcpow'][:4]]
    b = [k for k in range(4) if vt[k] != M[1][k] % P]
    print("%s R5b vertex graph polar part = (1/25200) t_c^-4 polar part: mismatches %s; vertex k<=3 = %s" % ("PASS" if not b else "FAIL", b, [str(ratrec(x, P)) for x in M[1][:4]]))
    if 'A_id' in st and 'Afull' in d:      # point frame (M = 1): the full R-action graph by graph, production record
        Aid = {int(k): v for k, v in st['A_id'].items()}; Fu = {int(k): v for k, v in d['Afull'].items()}
        b = cmp(Fu, Aid, range(6)); nz = sum(1 for i in Fu for k in range(6) if Fu[i][k] % P)
        print("%s R7 own composite regrouped by uncontracted graph = stored point-frame table A_id, 42 graphs, k <= 5 (252, %d non-zero): mismatches %s"
              % ("PASS" if not b else "FAIL", nz, b[:8]))
        for var in ("stored", "dress"):
            dv = ld(os.path.join(MINE, "rf_g3_%s_%s%s.json" % (pr, var, tag)))
            if dv and 'Afull' in dv:
                b = cmp({int(k): v for k, v in dv['Afull'].items()}, Aid, range(6))
                print("CONTROL %s R7[%s] vs A_id: %d of 252 differ" % ("FAILS AS IT MUST" if b else "DID NOT FAIL", var, len(b)))
    if 'A_M0' in st:
        A0 = {int(k): v for k, v in st['A_M0'].items()}
        b = cmp(M, A0, range(6)); t0 = [sum(A0[i][k] for i in A0) % P for k in range(6)]
        print("CONTROL %s X1 own flat vs constant-M_0 split A_M0: %d of 252 differ; totals agree: %s"
              % ("FAILS AS IT MUST" if b else "DID NOT FAIL", len(b), t0 == tot))
    for var, wv in (("stored", "stored"), ("nodil", "nodil"), ("noT", "noT"), ("nosign", "nosign"), ("zminus", None), ("dress", None)):
        X = mine(3, pr, var, tag)
        if X is None: continue
        ref = H if H else A3; ks = range(8) if H else range(4)
        b = cmp(X, ref, ks)
        line = "CONTROL %s %s vs hatted: %d of %d differ, graphs %s" % ("FAILS AS IT MUST" if b else "DID NOT FAIL", var, len(b), 42 * len(ks), sorted(set(i for i, _ in b)))
        if wv:
            Wv = ld(os.path.join(WORK, "flat_g3_%s_%s%s.json" % (pr, wv, tag)))
            if Wv:
                Wv = {int(k): v for k, v in Wv['A'].items()}; bb = cmp(X, Wv, range(8))
                line += "; = worker's '%s' table: %s" % (wv, "PASS" if not bb else "FAIL %s" % bb[:6])
        print(line)
    # genus 2 against rf_g2flat.log (first block = frame (0,0), canonical gauge)
    G2 = mine(2, pr, "can", tag)
    if G2:
        txt = open(os.path.join(PROJ, "allgenus", "theory", "referee", "rf_g2flat.log")).read()
        if tag: txt = txt.split("rf_g2flat: frame (1/3, -2/7)")[1]
        txt = txt.split("--- control X1")[0]
        rows = re.findall(r"^\s+(\d) .*?G = \[(.*?)\]", txt, re.M)[:7]; bad = []
        for i, s in rows:
            vals = [Fr(x.strip().strip("'")) for x in s.split(",")]
            for k, x in enumerate(vals):
                if (x.numerator * pow(x.denominator, P - 2, P) - G2[int(i)][k]) % P: bad.append((int(i), k))
        print("%s R6 own genus-2 flat = rf_g2flat.log (7 graphs x 4 orders): mismatches %s" % ("PASS" if not bad and len(rows) == 7 else "FAIL", bad))


if __name__ == "__main__":
    main()
