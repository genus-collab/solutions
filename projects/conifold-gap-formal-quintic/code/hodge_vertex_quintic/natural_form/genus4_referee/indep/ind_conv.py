# -*- coding: utf-8 -*-
"""ind_conv.py -- INDEPENDENT RE-IMPLEMENTATION (genus-4 referee): CONVENTION CHECK ONLY.

Compares the tables built by ind_tables (from the referee rows, own formulas) with those inside the project's fq4_split.Split
(from fq_rows / fq_conifold / fq_sser / fq3_msum / fq3_mvertex): E_R(a, b)[i][j], T_R(e)[i], E_M(a, b)[i][j], the u-valuations of
E_out, xhat, and the Witten-Kontsevich numbers used by the genus-g vertex tables.  This is the only place where the project's
genus-4 code is imported; nothing computed here enters ind_run.  The comparison is entrywise (the two codes happen to choose the
same primitive fifth root of unity; if they did not, the colours would be permuted and the check would be done up to that
permutation).   py -I ind_conv.py [--genus 4] [--prime 89]
"""
import os, sys, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--genus", type=int, default=4); ap.add_argument("--prime", default="89")
    ap.add_argument("--variant", default="nat"); a = ap.parse_args(); g = a.genus; K = max(2 * g - 2, 4)
    from ind_tables import Tables
    T = Tables(a.prime, g, K, a.variant); p = T.p
    sys.path.insert(0, PROJ)
    from fq4_split import Split
    sp = Split(p, g, (Fr(0), Fr(0)), K, a.variant)
    # colour permutation c -> m c relating the two primitive roots (m = 1 when they coincide)
    m = next(k for k in range(1, 5) if pow(T.zeta, 4 * k, p) == sp.zph[1]); perm = [(m * c) % 5 for c in range(5)]
    res = {}
    res["primitive roots coincide"] = (m == 1)
    res["E_R: %d tables x 25 entries x %d coefficients" % (len(T.ER), K)] = set(T.ER) == set(sp.ER) and all(
        T.ER[k][perm[i]][perm[j]] == sp.ER[k][i][j] for k in T.ER for i in range(5) for j in range(5))
    res["T_R: leaves e = 2..%d" % (T.NZ + 1)] = set(T.TR) == set(sp.TR) and all(T.TR[e][perm[c]] == sp.TR[e][c] for e in T.TR for c in range(5))
    res["E_M"] = all(T.EM[k][perm[i]][perm[j]] == sp.EM[k][i][j] for k in T.EM for i in range(5) for j in range(5))
    res["E_out"] = all(T.EO[k][perm[i]][perm[j]] == sp.EO[k][i][j] for k in T.EO for i in range(5) for j in range(5))
    res["u-valuations of E_out"] = (T.valO == sp.valO)
    res["xhat powers"] = all(T.xp[n] == sp.xpow[n] for n in range(min(len(T.xp), len(sp.xpow))))
    # vertex tables: all (g_v, hs) occurring in the genus-g graphs
    from ind_graphs import graphs
    keys = set((G['gv'][v], len(G['hs'][v])) for G in graphs(g) for v in range(G['V'])); nent = 0; ok = True
    for (gv, hs) in sorted(keys):
        A = T.vertex_table(gv, hs); B = sp.vertex_table(gv, hs)
        if set(A) != set(B): ok = False; continue
        for hp in A:
            nent += 1
            if any(A[hp][2][perm[c]] != B[hp][c] for c in range(5)): ok = False
    res["vertex vectors: %d types (g_v, valence), %d exponent patterns" % (len(keys), nent)] = ok
    with open(os.path.join(HERE, "out", "conv_g%d_%s_%s.log" % (g, a.prime, a.variant)), "a", encoding="utf-8") as f:
        for name, v in res.items():
            s = "%s %s" % ("EQUAL   " if v else "DIFFERENT", name); print(s); f.write(s + "\n")
        s = "IND_CONV (g%d, 2^%s-1, %s): %s" % (g, a.prime, a.variant, "ALL TABLES EQUAL" if all(v for k, v in res.items() if k != "primitive roots coincide") else "DIFFERENCES"); print(s); f.write(s + "\n")


if __name__ == "__main__":
    main()
