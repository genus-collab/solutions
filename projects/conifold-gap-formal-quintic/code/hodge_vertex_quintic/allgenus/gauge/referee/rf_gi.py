# -*- coding: utf-8 -*-
"""rf_gi.py -- REFEREE: Theorem G, (G*), Proposition S, (R1)-(R3), Lemma Q, the diagonal, and the float diagnosis, exact over Q,
with OWN periods/systems/vanishing solutions (rf_core) and the PROJECT gauge ur_gauge (rf_ur), first certified by direct substitution
into the gauge equation built from the OWN coefficients.  Both the canonical gauge Pc = P c(z)^-1 and the STORED gauge P (c = P_11(0,z))
are tested: Theorem G claims P Y^C = c(z) Y^C_Q for EVERY regular gauge, (G*) claims Y^C_Q(-z)^T eta P(z) = c(z) Y^C(-z)^T eta.
Run: py -I -B rf_gi.py K NU al be      log: rf_gi_<frame>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import systems, model_vc, quintic_vc, mul, deriv, powf, inv, CEXP
from rf_ur import stored_gauge, p11_zero, normalise


def conv(a, b, n): return sum(a[i] * b[n - i] for i in range(n + 1) if i < len(a) and n - i < len(b))


def main():
    K = int(sys.argv[1]); NU = int(sys.argv[2]); al = Fr(sys.argv[3]); be = Fr(sys.argv[4])
    tag = "ab_%s_%s" % (str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "rf_gi_%s.log" % tag), "w", encoding="ascii"); FAIL = []
    def log(s): print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(name, ok, det=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(name)
    t0 = time.time()
    p, prec = stored_gauge(K, NU, al, be); c = p11_zero(p); pc, ci = normalise(p, c)
    N = prec - 1; S = systems(N + 2, al, be); b, bt = S['b'], S['bt']
    log("rf_gi: frame (%s, %s), project gauge z^0..z^%d, common precision O(u^%d); own coefficients; %.0f s" % (al, be, K, prec, time.time() - t0))
    log("   stored c(z) = P_11(0, z): %s" % ", ".join("z^%d: %s" % (k, c[k]) for k in sorted(c)))
    ub = lambda B, a: B[1] if a == 1 else [Fr(0)] + B[a]
    for name, g in (("stored", p), ("canonical", pc)):
        bad = 0; tested = 0
        for k in range(K + 1):
            for cc in range(5):
                lhs = [Fr(0)] + deriv(g[k - 1][cc]) if k else [Fr(0)] * N
                rhs1 = mul(ub(b, (cc + k - 1) % 5), g[k][cc], N); rhs2 = mul(ub(bt, (cc - 1) % 5), g[k][(cc - 1) % 5], N)
                for n in range(N - 1):
                    tested += 1
                    if lhs[n] != rhs1[n] - rhs2[n]: bad += 1
        check("direct substitution: the %s gauge solves u p_{k-1,c}' = (u b_{c+k-1}) p_{k,c} - (u bt_{c-1}) p_{k,c-1} with OWN b, bt (%d coefficients)" % (name, tested), bad == 0, bad)
    r1 = all(pc[k][1][0] == 0 for k in range(K + 1) if k % 5) and all(pc[k][(2 - k) % 5][0] == 0 for k in range(K + 1) if k % 5)
    r1 = r1 and all(pc[k][2][0] == pc[k][1][0] for k in range(0, K + 1, 5)) and all(p[k][2][0] == p[k][1][0] for k in range(0, K + 1, 5))
    check("(R1) column 1 and row 2 of P(0, z) diagonal, P_22(0, z) = P_11(0, z) (stored and canonical)", r1)
    r23 = all(pc[k][1][0] == 0 and pc[k][1][1] == 0 and (k + 2 > K or pc[k][1][2] == -pc[k + 2][0][0]) for k in range(1, K + 1, 5))
    r23 = r23 and all(k + 1 > K or (pc[k][0][1] == pc[k + 1][0][0] and pc[k][1][1] == -pc[k + 1][0][0]) for k in range(2, K + 1, 5))
    check("(R2) P_21 = O(u^2); (R3) P_20 = u pi/z, P_31 = -u pi/z, P_21 = -u^2 pi/z^2 + .. (pi = P_30(0, z)), canonical gauge", r23)
    log("   pi(z) = %s" % ", ".join("z^%d: %s" % (k, pc[k][0][0]) for k in range(3, K + 1, 5)))
    # vanishing solutions
    nmax = (5 * N) // 4 + 2
    YC = model_vc(S['v'][:N], N, nmax); YQ = quintic_vc(b, N, nmax)
    val_ok = True; nval = 0
    for n in range(0, 46):
        for a in range(5):
            w5 = 4 * n + CEXP[a] + 1
            yq = YQ[n][a] if n in YQ else None; ym = YC[n][a] if n in YC else None
            nz_q = yq is not None and any(yq); nz_m = ym is not None and any(ym)
            if (n - 1 + a) % 5:
                if nz_q or nz_m: val_ok = False
                continue
            w = w5 // 5
            if w >= N: continue
            nval += 1
            if not (nz_q and nz_m): val_ok = False; continue
            vq = min(i for i, x in enumerate(yq) if x); vm = min(i for i, x in enumerate(ym) if x)
            if vq != w or vm != w or yq[w] != ym[w]: val_ok = False
    check("Lemma Q (V): ord_u [z^-n] Y^C_Q,a = (4n + c_a + 1)/5 exactly, zero off n = 1 - a mod 5, same leading term as the model (%d slots, n <= 45)" % nval, val_ok)
    check("Lemma Q (d): [z^-1] Y^C_Q = -t_c e_0", all(YQ[1][0][n] == -S['tc'][n] for n in range(N)) and not any(any(YQ[1][a]) for a in range(1, 5)))
    # Theorem G for both gauges: P Y^C = c(z) Y^C_Q
    for name, g, cz in (("canonical", pc, {0: Fr(1)}), ("stored", p, c)):
        tested = 0; bad = []; pos_nz = 0
        for a in range(5):
            for m in range(-nmax, K + 1):
                if (m - (a - 1)) % 5: continue
                for n in range(N):
                    if m + (5 * n) // 4 > K: continue
                    lhs = sum(conv(g[k][(a - k) % 5], YC[k - m][(a - k) % 5], n) for k in range(max(m, 0), K + 1) if k - m in YC)
                    rhs = sum(cz[j] * YQ[j - m][a][n] for j in cz if j - m in YQ and j <= K)
                    tested += 1
                    if m > 0 and lhs != 0: pos_nz += 1
                    if lhs != rhs: bad.append((a, m, n))
        check("Theorem G, %s gauge: P Y^C = c(z) Y^C_Q on %d complete graded slots (positive-z slots with non-zero P Y^C: %d)" % (name, tested, pos_nz), not bad, bad[:4])
    # (G*) and the float diagnosis
    def gstar(g, cz, sgn, bound):
        tested = 0; bad = []; fl = 0
        for cc in range(5):
            for e in range(-nmax, K + 1):
                for n in range(N):
                    if e + bound(n) > K: continue
                    lhs = 0
                    for k in range(max(e, 0), K + 1):
                        a = (cc + k) % 5; m = k - e
                        if m in YQ: lhs += 5 * sgn(-m) * conv(YQ[m][(3 - a) % 5], g[k][cc], n)      # sign of the z-POWER -m, as in gi_*
                    rhs = sum(cz[j] * 5 * sgn(e - j) * YC[j - e][(3 - cc) % 5][n] for j in cz if j - e in YC and j - e >= 0)
                    if lhs == 0 and rhs == 0: continue
                    tested += 1
                    if isinstance(lhs, float) or isinstance(rhs, float): fl += 1
                    if lhs != rhs: bad.append((n, cc, e, type(lhs).__name__, abs(float(lhs) - float(rhs))))
        return tested, bad, fl
    isg = lambda m: 1 if m % 2 == 0 else -1
    sharp = lambda n: (5 * n) // 4; old = lambda n: 2 * n + 2
    for name, g, cz, expect in (("canonical, c = 1", pc, {0: Fr(1)}, True), ("stored, twisted by its c(z)", p, c, True),
                                ("CONTROL stored, untwisted", p, {0: Fr(1)}, False), ("CONTROL canonical, wrong twist c(z)", pc, c, False)):
        t, bad, fl = gstar(g, cz, isg, sharp)
        check("(G*) %s: %d non-trivial complete identities, %d failures (expected %s)" % (name, t, len(bad), "none" if expect else "some"),
              (not bad) == expect, bad[:3] if expect else "")
    for bname, bound in (("old bound 2n+2", old), ("sharp bound floor(5n/4)", sharp)):
        t, bad, fl = gstar(pc, {0: Fr(1)}, lambda m: (-1) ** m, bound)
        t2, bad2, fl2 = gstar(pc, {0: Fr(1)}, isg, bound)
        log("   float diagnosis (%s): with (-1)**m: %d tested, %d float-valued, %d fail, max |lhs-rhs| %s; integer sign: %d tested, %d fail"
            % (bname, t, fl, len(bad), ("%.1e" % max(x[4] for x in bad)) if bad else "-", t2, len(bad2)))
        check("   float failures only in float-valued slots, all of them exact identities ([%s])" % bname,
              all(x[3] == 'float' for x in bad) and not bad2 and len(bad) > 0)
    # Proposition S
    def sdef(g):
        out = {}
        for k in range(K + 1):
            for k1 in range(k + 1):
                k2 = k - k1
                for b1 in range(5):
                    a1 = (b1 + k1) % 5; a2 = (3 - a1) % 5; b2 = (a2 - k2) % 5
                    t = mul(g[k1][b1], g[k2][b2], N)
                    key = (k, b1, b2); acc = out.setdefault(key, [Fr(0)] * N)
                    for n in range(N): acc[n] += 5 * isg(k1) * t[n]
        return out
    for name, g, cz in (("canonical", pc, {0: Fr(1)}), ("stored", p, c)):
        tt = {}
        for k1, x in cz.items():
            for k2, y in cz.items(): tt[k1 + k2] = tt.get(k1 + k2, 0) + x * y * isg(k2)
        D = sdef(g); bad = 0
        for (k, b1, b2), s in D.items():
            want = 5 * tt.get(k, 0) if (b1 + b2) % 5 == 3 else 0
            if s[0] != want or any(s[1:]): bad += 1
        check("Prop S, %s gauge: P(-z)^T eta P(z) = c(z) c(-z) eta, all entries z^0..z^%d mod u^%d" % (name, K, N), bad == 0, bad)
    lt = S['tau']; tq = powf(lt[:N], Fr(-1, 4), N); target = mul(tq, inv(S['varpi0'][:N], N), N)
    check("diagonal: p_{0,0} = tau^(-1/4)/varpi_0 (canonical), and prod_a p_{0,a} = 1", pc[0][0][:N] == target[:N] and
          mul(mul(mul(pc[0][0], pc[0][1], N), mul(pc[0][2], pc[0][3], N), N), pc[0][4], N) == [Fr(1)] + [Fr(0)] * (N - 1))
    log("RF_GI %s: %s (%.0f s)" % (tag, "ALL PASS" if not FAIL else "FAIL %s" % FAIL, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
