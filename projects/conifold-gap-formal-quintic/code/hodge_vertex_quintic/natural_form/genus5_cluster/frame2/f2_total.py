# -*- coding: utf-8 -*-
"""f2_total.py -- the TWO-PRIME total of the second conifold frame (alpha, beta) = (1/3, -2/7), from the point-frame runs.

Sums, with this file's own code, the point-frame term of all 4555 inner graphs of the two 'id' runs (p = 2^89-1 and
p = 2^61-1) in the second frame, CRT-combines the 8 coefficients of u^0..u^7 and reconstructs them as rationals.  Then
rewrites u^-8 Sum_k a_k u^k in the conifold coordinate t_c OVER THE RATIONALS (QRing, no prime) and checks the polar part
    (1/660000) t_c^-8 + 0 t_c^-7 + ... + 0 t_c^-1 ,     1/660000 = |B_10| / (10 . 8 . 5^4).
Optionally the same for the frame (0, 0) 'id' runs, to display that the u-series totals of the two frames differ while the
t_c polar parts agree.

usage:  py -I -B f2_total.py --dir89 <parts> --dir61 <parts> [--dir89b <parts (0,0)> --dir61b <parts (0,0)>] [--json OUT]
"""
import os, sys, json, argparse
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_lib import K, NIN, P89, P61, read_parts_dir, crt2, ratrec

C5 = Fr(1, 660000)


def total_of(d, tag, p):
    """sum over all inner graphs of the stored point-frame term; asserts completeness and the diagonal-only shape"""
    _, pieces, times = read_parts_dir(d, tag)
    assert len(pieces) == NIN, "%s: %d inner graphs, expected %d" % (tag, len(pieces), NIN)
    tot = [0] * K
    for t, res in pieces.items():
        for q, c in res.items():
            assert q == t, "%s: inner graph %s stores a piece for q = %s" % (tag, t, q)
            for k in range(K):
                tot[k] = (tot[k] + c[k]) % p
    return tot, sum(times.values())


def crt_series(a89, a61):
    M = P89 * P61
    out = []
    for k in range(K):
        x = crt2(a89[k], P89, a61[k], P61)
        r = ratrec(x, M)
        out.append(r)
    return out, M


def polar(acoefs, ab):
    """the 8 coefficients of t_c^-8 .. t_c^-1 of u^-8 Sum a_k u^k, over the rationals"""
    from fq_fp import QRing
    from fq_conifold import Conifold
    from fq_hat import tc_series
    from fq3_msum import tc_expansion3
    Q = QRing()
    C = Conifold(Q, Nu=K + 4, alpha=Fr(ab[0]), beta=Fr(ab[1]))
    tc = tc_series(C, Q, K + 4)
    T = tc_expansion3([Q.of(x) for x in acoefs], tc, Q, 5)[:2 * 5 - 2]
    return T, tc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir89", required=True)
    ap.add_argument("--dir61", required=True)
    ap.add_argument("--dir89b", default="")
    ap.add_argument("--dir61b", default="")
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    out = {}
    ok = True

    for lab, ab, d89, d61, t89, t61 in (
            ("second frame (1/3, -2/7)", (Fr(1, 3), Fr(-2, 7)), a.dir89, a.dir61, "g5_89_id_ab1o3_m2o7", "g5_61_id_ab1o3_m2o7"),
            ("frame (0, 0)", (Fr(0), Fr(0)), a.dir89b, a.dir61b, "g5_89_id", "g5_61_id")):
        if not d89 or not d61:
            continue
        print("=" * 90)
        print(lab)
        s89, c89 = total_of(d89, t89, P89)
        s61, c61 = total_of(d61, t61, P61)
        print("   point-frame totals: p89 %s" % s89)
        print("                       p61 %s" % s61)
        print("   CPU of the two runs: %.2f h + %.2f h" % (c89 / 3600, c61 / 3600))
        rat, M = crt_series(s89, s61)
        print("   CRT modulus 2^89-1 times 2^61-1 ~ 2^150; reconstruction bound ~2^71")
        print("   total u-series  = u^-8 ( %s )" % " + ".join("(%s) u^%d" % (rat[k], k) for k in range(K)))
        bad = [k for k in range(K) if rat[k] is None]
        if bad:
            ok = False
            print("   FAIL: coefficients %s did not reconstruct" % bad)
        else:
            # confirm each rational against both primes separately
            for k in range(K):
                assert (rat[k].numerator * pow(rat[k].denominator, -1, P89) - s89[k]) % P89 == 0
                assert (rat[k].numerator * pow(rat[k].denominator, -1, P61) - s61[k]) % P61 == 0
            print("   each reconstructed rational re-reduces to the stored residue at both primes: OK")
        T, tc = polar(rat, ab)
        print("   t_c(u) = %s + ..." % " + ".join("(%s) u^%d" % (tc.coef(k), k) for k in range(tc.v, min(tc.v + 4, tc.prec))))
        print("   polar part: " + " + ".join("(%s) t_c^%d" % (T[k], k - 8) for k in range(8)))
        good = (T[0] == C5) and not any(T[1:])
        print("   %s: polar part = (1/660000) t_c^-8 + 0 + ... + 0" % ("PASS" if good else "FAIL"))
        ok = ok and good
        out[lab] = dict(res89=s89, res61=s61, u_series=[str(x) for x in rat], polar=[str(x) for x in T],
                        leading_is_c5=(T[0] == C5), higher_all_zero=not any(T[1:]))

    if "second frame (1/3, -2/7)" in out and "frame (0, 0)" in out:
        A = out["second frame (1/3, -2/7)"]["u_series"]
        B = out["frame (0, 0)"]["u_series"]
        diff = [k for k in range(K) if A[k] != B[k]]
        out["u_series_differ_at"] = diff
        out["polar_parts_equal"] = out["second frame (1/3, -2/7)"]["polar"] == out["frame (0, 0)"]["polar"]
        print("=" * 90)
        print("the two frames' u-series totals differ at the orders %s (equal at u^0)" % diff)
        print("their t_c polar parts are equal: %s" % out["polar_parts_equal"])
        ok = ok and bool(diff) and out["polar_parts_equal"]
    print("RESULT: " + ("PASS" if ok else "FAIL"))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1)
        print("written %s" % a.json)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
