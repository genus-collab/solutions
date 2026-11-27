# -*- coding: utf-8 -*-
"""ctl_cross.py -- the two genus-5 controls measure THE SAME four numbers.  Decisive cross-check, from the copied-back JSONs only.

   Floor-order response (ctl_floor.py, confirmed at genus 3): for the vertex graph and the one-edge graphs,
       [u^k] (row(tau + c u^k) - row(tau)) = c . A_Gamma      with A_Gamma independent of k,
   because the perturbed edge factor contributes m c Eh_0(a,b) u^k (m = a+b+1) and every other factor of the row
   must then be taken at u^0.  The bare symbol M_0 is tau' = 1 = tau + (1 - tau), whose lowest term is
   c_1 = -tau_1 = 4/9, so
       [u^1] (row(M0) - row(nat)) = (4/9) A_Gamma   and   [u^7] (row(nat_t7) - row(nat)) = A_Gamma
   for Gamma = the vertex graph and the three one-edge graphs.  The two controls were run independently (different
   variants, different tasks) and the natural run is a third run, so the four identities below are a real test.

   py -I -B ctl_cross.py
"""
import os, sys, json
from fractions import Fraction as Fr
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
RES = os.path.join(PKG, "results")
sys.path.insert(0, PKG)


def main():
    from fq_fp import FpRing, P1, P2
    for prime, p in (("89", P2),):
        R = FpRing(p)
        tabs = {}
        for tag in ("g5_%s_nat" % prime, "g5_%s_M0_out2" % prime, "g5_%s_nat_t7_out2" % prime):
            d = json.load(open(os.path.join(RES, tag + ".json")))
            tabs[tag.split(prime + "_")[1]] = {int(k): v for k, v in d['table'].items()}
            outer = d['outer']
        nat = tabs['nat']; m0 = tabs['M0_out2']; t7 = tabs['nat_t7_out2']
        f = R.of(Fr(4, 9))
        print("=" * 100)
        print("2^%s-1:  [u^1](row(M0) - row(nat))  ==  (4/9) . [u^7](row(nat_t7) - row(nat))   for E <= 1" % prime)
        print("=" * 100)
        ok = True
        for i in (0, 1, 2, 3):
            A = (t7.get(i, [0] * 8)[7] - nat.get(i, [0] * 8)[7]) % p
            lhs = (m0.get(i, [0] * 8)[1] - nat.get(i, [0] * 8)[1]) % p
            rhs = f * A % p
            good = lhs == rhs and A != 0
            ok = ok and good
            print("   row %d (E=%d, gv=%s): A = %d" % (i, outer[i]['E'], list(outer[i]['gv']), A))
            print("      [u^1] M0 - nat = %d" % lhs)
            print("      (4/9) A        = %d   -> %s" % (rhs, "EQUAL, A non-zero" if good else "MISMATCH" if lhs != rhs else "A = 0"))
        print()
        print("   verdict: %s" % ("the two controls agree on all four numbers and all four are non-zero"
                                  if ok else "DISAGREEMENT"))
        print()
        print("   consequence: for ANY tau' = tau + sum_{k >= k0} c_k u^k with c_k0 != 0 and k0 <= 7, the three one-edge")
        print("   rows have [u^k0] = c_k0 A != 0, so the regularity test fails.  The natural tau is the unique tau mod u^8")
        print("   for which the one-edge graphs are regular, and the test is sensitive at every order u^1..u^7.")
    # the same identity on the full genus-4 tables of record (window u^0..u^5, so k = 5 against k = 1)
    print()
    print("=" * 100)
    print("genus 4, 2^89-1, full tables of record:  [u^1](M0 - nat) == (4/9) [u^5](nat_t5 - nat)  for E <= 1")
    print("=" * 100)
    ref = json.load(open(os.path.join(PKG, "data", "ref_g4.json")))
    nat4 = ref['g4_89_nat']['table']; m04 = ref['g4_89_M0']['table']; t54 = ref['g4_89_nat_t5']['table']
    R = FpRing(P2); f = R.of(Fr(4, 9)); ok = True
    for key in sorted(nat4, key=lambda k: (len(json.loads(k)[1]), k)):
        gv, ed = json.loads(key)
        if len(ed) > 1: continue
        A = (t54[key][5] - nat4[key][5]) % P2
        lhs = (m04[key][1] - nat4[key][1]) % P2
        good = (lhs == f * A % P2) and A != 0
        ok = ok and good
        print("   gv=%-10s ed=%-12s A = %d ; [u^1] M0 - nat = %d ; (4/9) A = %d -> %s"
              % (gv, ed, A, lhs, f * A % P2, "EQUAL, A non-zero" if good else "MISMATCH"))
    print("   verdict: %s" % ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
