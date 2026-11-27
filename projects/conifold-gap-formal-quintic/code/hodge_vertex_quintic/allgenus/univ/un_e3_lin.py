# -*- coding: utf-8 -*-
"""un_e3_lin.py -- E3, first order: does annihilation of the vanishing covector (A) + flat regularity + symplecticity force the
first-order regularity of every edge graph?  Linear algebra mod p on the complete response map of un_e3_resp (un_resp_<p>.json).

Variables (un_lin.variables): the Taylor coefficients p_{k,b,j} of a flat-regular, Z/5-graded, infinitesimally symplectic left
dressing p(u, z) = Sum p_{k,b}(u) H(k, b) of R' that are visible at genus g (j + mu(b, k) <= W = 2g - 3), k = 0..3g-3.
Matrices: Resp (polar coefficients [u^0..u^W] of every outer graph with edges), Vert (the vertex graph), Ann (C_m[u^i] = 0, the
conditions (p Y^C_Q)_+ = 0 in the window, un_lin.annihilation_rows).  Reported:
  dim V, rank Ann, dim ker Ann, rank Resp, rank(Resp on ker Ann)   [0 <=> (A) + flat regularity suffice at first order]
  rank(Ann on ker Resp)                                           [0 <=> regularity forces (A) at first order]
  the vertex-graph law on ker Ann: delta(vertex) = c_g (2-2g) (u/t_c)^(2g-1) (delta T/u), delta T = -[z^-1 e_0] p Y^C_Q.
Run: py -I -B un_e3_lin.py [89]
"""
import os, sys, json, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, os.path.join(PROJ, "natural_form", "theory")); sys.path.insert(0, HERE)
from un_core import PRIMES, mu
from un_lin import variables, response_rows, vertex_rows, yc_series, annihilation_rows, teff_row, rank, nullspace, rref
from fq_fp import FpRing
from fq3_rec import Rec

LOGF = os.path.join(HERE, "un_e3_lin.log"); FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True)
    with open(LOGF, "a", encoding="utf-8") as f: f.write(s + "\n")


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def matmul(A, B, p):
    """A (rows x n) times B given as a list of column vectors (each of length n): rows x len(B)"""
    return [[sum(r[c] * v[c] for c in range(len(r)) if r[c] and v[c]) % p for v in B] for r in A]


def hat_variables(g):
    """all hatted Taylor coefficients x_{k,b,j} (orbit representatives, k = 0..3g-3, j = 0..2g-3)"""
    W = 2 * g - 3; out = []
    for k in range(0, 3 * g - 2):
        seen = set()
        for b in range(5):
            if b in seen: continue
            bs = (3 - b - k) % 5; seen.add(b); seen.add(bs)
            if bs == b and k % 2 == 0: continue
            for j in range(W + 1): out.append((k, b, j))
    return out


def hat_response_rows(resp_g, g, HV, p):
    W = 2 * g - 3; rows = []; labels = []
    idxs = sorted(set(int(i) for d in resp_g.values() for i in d))
    for i in idxs:
        if i == 0: continue
        for n in range(W + 1):
            row = []
            for (k, b, j) in HV:
                c = resp_g.get("%d,%d" % (k, b), {}).get(str(i)); d = n - j
                row.append(c[d] % p if (c is not None and 0 <= d < len(c)) else 0)
            if any(row): rows.append(row); labels.append((i, n))
    return rows, labels


def hat_annihilation_rows(g, HV, ys, R):
    """delta_m[u^i] = 0 (hatted), m = 1..3g-3, i = 0..2g-3; yhat_n = y_n/u^(n - n//5)"""
    from un_lin import orbit_of
    W = 2 * g - 3; KZ = 3 * g - 3; pos = {x: n for n, x in enumerate(HV)}; rows = []; labels = []
    yh = []
    for n, y in enumerate(ys):
        nu = n - n // 5; yh.append([R.of(c) for c in y[nu:nu + W + 1]])
    for m in range(1, KZ + 1):
        for i in range(W + 1):
            row = [0] * len(HV)
            for n in range(len(yh)):
                k = m + n; col = (1 - n) % 5
                if k > KZ: break
                ob = orbit_of(k, col)
                if ob is None: continue
                b, sg = ob
                for j in range(i + 1):
                    c = yh[n][i - j]
                    if c and (k, b, j) in pos: row[pos[(k, b, j)]] = (row[pos[(k, b, j)]] + sg * c) % R.p
            if any(row): rows.append(row); labels.append((m, i))
    return rows, labels


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"; p = PRIMES[prime]; R = FpRing(p); rec = Rec(R)
    import ur_gauge
    d = json.load(open(os.path.join(HERE, "un_resp_%s.json" % prime)))
    log("=" * 110); log("un_e3_lin: p = 2^%s-1   [%s]" % (prime, datetime.datetime.now().isoformat(timespec='seconds')))
    bq, bt, tau, v, C = ur_gauge.coefficients(30)
    ys = yc_series(bq, 14, 24)
    log("   y_n (n = 0..6) of Y^C_Q: %s" % [[str(c) for c in y[:4]] for y in ys[:7]])
    vals = [next((i for i, c in enumerate(y) if c != 0), None) for y in ys[:12]]
    check("   y_n has u-valuation n - floor(n/5) (n = 0..11)", vals == [n - n // 5 for n in range(12)], vals)
    out = {}
    for g in (2, 3):
        W = 2 * g - 3; VARS = variables(g); nv = len(VARS); resp = d['resp'][str(g)]; outer = d['outer'][str(g)]
        Rr, Rl = response_rows(resp, outer, g, VARS, p); Vr = vertex_rows(resp, g, VARS, p)
        Ar, Al = annihilation_rows(g, VARS, ys, R)
        N = nullspace(Ar, p, nv); RN = matmul(Rr, N, p); rRN = rank(RN, p, len(N)) if N else 0
        KR = nullspace(Rr, p, nv); AK = matmul(Ar, KR, p); rAK = rank(AK, p, len(KR)) if KR else 0
        log("-" * 110); log("genus %d (window u^0..u^%d): %d variables %s" % (g, W, nv, VARS))
        log("   rank Ann = %d (%d conditions), dim ker Ann = %d;  rank Resp = %d (%d non-zero rows), dim ker Resp = %d" % (rank(Ar, p, nv), len(Ar), len(N), rank(Rr, p, nv), len(Rr), len(KR)))
        log("   rank(Resp restricted to ker Ann) = %d;   rank(Ann restricted to ker Resp) = %d" % (rRN, rAK))
        check("   genus %d: (A) + flat regularity + symplecticity => every edge graph regular at FIRST order" % g, rRN == 0)
        log("   %s" % ("ker Resp = ker Ann (regularity <=> annihilation at first order)" if (rRN == 0 and rAK == 0) else
                     "ker Ann %s ker Resp; ker Resp %s ker Ann" % ("is in" if rRN == 0 else "is NOT in", "is in" if rAK == 0 else "is NOT in")))
        if rRN:
            piv, Mr = rref([list(r) for r in RN], p, len(N))
            bad = sorted(set(Rl[i][0] for i, r in enumerate(RN) if any(r)))
            log("   edge graphs responding on ker Ann: %s (outer index: (E, gv))" % {i: (outer[i]['E'], outer[i]['gv']) for i in bad})
        if rAK:
            badA = sorted(set(Al[i] for i, r in enumerate(AK) if any(r)))
            log("   annihilation conditions (m, u^i) violated somewhere on ker Resp: %s" % badA)
        # vertex law on ker Ann:  [u^n] delta(vertex G) = [u^n] c_g (2-2g) (u/t_c)^(2g-1) (delta T/u)
        cg = {2: Fr(1, 1200), 3: Fr(1, 25200)}[g]; tc = [C.tc.coef(k) for k in range(1, W + 3)]          # t_c/u
        from fq_series import Ser
        from fq_fp import QRing
        Qr = QRing(); tcu = Ser.from_list(Qr, tc, 0, W + 2); fac = tcu.inv().pow(2 * g - 1).scale(cg * (2 - 2 * g))
        TT = [teff_row(g, VARS, ys, R, n) for n in range(0, W + 2)]           # rows of [u^n] delta T
        bad_v = []; VL = []
        for n in range(W + 1):
            # [u^n] fac * (delta T/u) = Sum_j fac_j [u^(n - j + 1)] delta T
            row = [0] * nv
            for j in range(n + 1):
                f = R.of(fac.coef(j))
                if f:
                    for c in range(nv): row[c] = (row[c] + f * TT[n - j + 1][c]) % p
            diff = [(Vr[n][c] - row[c]) % p for c in range(nv)]; VL.append(diff)
            if any(sum(diff[c] * w[c] for c in range(nv)) % p for w in N): bad_v.append(n)
        check("   genus %d: on ker Ann the vertex graph moves as c_g T^(2-2g) with T -> t_c + delta T (orders u^0..u^%d)" % (g, W), not bad_v, bad_v)
        KRV = nullspace(Rr + VL, p, nv); rAKV = rank(matmul(Ar, KRV, p), p, len(KRV)) if KRV else 0
        log("   edge graphs regular AND vertex law (vertex = c_g T_eff^(2-2g) to first order): dim %d; rank(Ann on it) = %d -> %s" % (len(KRV), rAKV,
            "equivalent to annihilation at first order" if (rAKV == 0 and len(KRV) == len(N)) else "NOT equivalent to annihilation"))
        out[g] = dict(nvars=nv, rankAnn=rank(Ar, p, nv), kerAnn=len(N), rankResp=rank(Rr, p, nv), RN=rRN, AK=rAK, regular_and_vertex_law=len(KRV), AKV=rAKV)
        # ---- HATTED variables (flat regularity dropped): x_{k,b,j}, j <= W, any j >= 0
        HV = hat_variables(g); nh = len(HV); pos = {x: i for i, x in enumerate(HV)}
        Rh, _ = hat_response_rows(resp, g, HV, p); Ah, Ahl = hat_annihilation_rows(g, HV, ys, R)
        Nh = nullspace(Ah, p, nh); RNh = matmul(Rh, Nh, p); rh = rank(RNh, p, len(Nh)) if Nh else 0
        flat = [i for i, (k, b, j) in enumerate(HV) if j >= mu(b, k)]
        KRA = nullspace(Rh + Ah, p, nh)                                   # regular AND annihilating, hatted
        nonflat = [x for i, x in enumerate(HV) if i not in flat]
        proj = [[v[i] for i in range(nh) if i not in flat] for v in KRA]
        rirr = rank(proj, p, len(nonflat)) if KRA else 0
        log("   HATTED variables (flat regularity dropped): %d (%d of them flat-irregular); rank Ann = %d, dim ker Ann = %d; rank(Resp on ker Ann) = %d"
            % (nh, len(nonflat), rank(Ah, p, nh), len(Nh), rh))
        log("   dim(ker Resp cap ker Ann) = %d = %d flat-regular + %d with a flat-irregular component" % (len(KRA), len(KRA) - rirr, rirr))
        check("   genus %d: annihilation + symplecticity WITHOUT flat regularity is NOT sufficient at first order (rank > 0)" % g, rh > 0)
        if rirr:
            piv, Mr = rref(proj, p, len(nonflat))
            log("   flat-irregular directions invisible at first order (annihilating and regular): pivots %s" % [nonflat[c] for c in piv])
        out[g].update(hat_vars=nh, hat_rank_resp_on_ann=rh, hat_reg_ann=len(KRA), hat_irregular_invisible=rirr)
    log("UN_E3_LIN %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
    json.dump(out, open(os.path.join(HERE, "un_e3_lin_%s.json" % prime), "w"))


if __name__ == "__main__":
    main()
