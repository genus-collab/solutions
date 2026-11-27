# -*- coding: utf-8 -*-
"""rf_s15_frozen.py -- calibration and self-tests of the referee engine (own code only; no project import).
  (1) WK spot values; graph counts (genus 2: 7, genus 3: 42, genus 4: 379);
  (2) the frozen R-matrix (Frobenius recursion along the v-line, rf_s15_base.Frozen) is symplectic for eta_c = diag(1/Delta), and its
      unit row at colour 0 equals F_0 of UNIFORM_REDUCTION s.2 (stationary phase of the LG integral) through t^7;
  (3) frozen free energies F_h(v = 1) = c_h = |B_2h|/(2h(2h-2)5^(h-1)) for h = 2, 3 (and 4 with argument 'g4'); R(-z) fails.
Run: py -I -B rf_s15_frozen.py [g4]
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_s15_base import P61, P89, Frozen, wk, stable_graphs, fr
from rf_s15_sum import Engine

F0 = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000),
      Fr(52989974037, 25600000000), Fr(-7749624763881, 512000000000)]          # UNIFORM_REDUCTION s.2
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}


def main():
    ok = True
    def rep(flag, s):
        nonlocal ok; ok = ok and flag; print(("PASS " if flag else "FAIL ") + s, flush=True)
    rep(wk(1, (1,)) == Fr(1, 24) and wk(2, (4,)) == Fr(1, 1152) and wk(3, (7,)) == Fr(1, 82944) and wk(2, (2, 3)) == Fr(29, 5760)
        and wk(3, (2,) * 6) == Fr(1225, 144) and wk(4, (10,)) == Fr(1, 7962624), "WK spot values (incl. <tau_2^6>_3 = 1225/144, <tau_10>_4 = 1/7962624)")
    gs = {g: len(stable_graphs(g)) for g in (2, 3, 4)}
    rep(gs == {2: 7, 3: 42, 4: 379}, "stable leg-free graph counts %s" % gs)
    hmax = 4 if len(sys.argv) > 1 and sys.argv[1] == 'g4' else 3
    for p in (P89, P61):
        Fz = Frozen(p, 10)
        rep(not Fz.symplectic_defect(), "p = 2^%d-1: frozen R symplectic through z^10" % (p.bit_length()))
        rep(Fz.unit_row(0, 7) == [fr(x, p) for x in F0], "p = 2^%d-1: unit row of colour 0 = F_0 (UR s.2) through t^7" % p.bit_length())
        for zs in (1, -1):
            Fz = Frozen(p, 3 * hmax, zs); E = Engine(Fz, 1)
            for h in range(2, hmax + 1):
                t = time.time(); r, _ = E.outer_sums(h); tot = sum(v[0] for v in r.values()) % p
                good = tot == fr(CG[h], p)
                if zs == 1: rep(good, "p = 2^%d-1: frozen F_%d(1) = %s (%.0f s)" % (p.bit_length(), h, CG[h], time.time() - t))
                elif h == 2: rep(not good, "control: R(-z) gives F_2(1) != 1/1200 (as it must)")
                else: print("INFO R(-z): F_%d(1) %s c_%d (the z-sign is fixed by F_0 and F_2, not by F_3)" % (h, "=" if good else "!=", h))
    print("ALL PASS" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
