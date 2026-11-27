# -*- coding: utf-8 -*-
"""rf_mla_scope.py -- scope of the ML_Lambda test (th_mla.py / rf_mla.py) without computing any correlator: which singular lists are
skipped by the psi-deficit cut (<= 6 for MVertex6), and whether any TESTED singular list contains a frame that depends on the
i >= 1 terms of the vanishing series Y^C (I_m = Sum_i (-1/(4z))^i/i! ..., i.e. the z^-5 and deeper corrections).  A frame "depends on
i >= 1" if it changes when the factor (-1)^i of those terms is flipped (control C4 of rf_mla.py).
Run: py -I -B rf_mla_scope.py   Log: rf_mla_scope.log (overwritten).
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_mla_core import frames_upto, singular_lists, w
from rf_flat_core import yc_monomials


def main():
    LOG = open(os.path.join(HERE, "rf_mla_scope.log"), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    Ybad = [(b, j, n, cf * (-1) ** (n - j)) for (b, j, n, cf) in yc_monomials(60)]
    for h, NT in ((0, 5), (1, 5), (1, 6), (2, 4), (2, 5), (3, 3), (3, 4), (4, 2), (4, 3)):
        wmax = 2 * h - 2 + NT + 2
        FR = frames_upto(wmax); FB = frames_upto(wmax, Y=Ybad)
        sens = sorted(c for c in FR if FR[c] != FB[c])
        st = dict(sing=0, skip=0, sens=0, sens_skip=0, sens_tested=0)
        mind = None
        for combo, m, e in singular_lists(h, NT, FR):
            if e >= 0: continue
            st['sing'] += 1
            n = len(combo) + m; dim = 3 * h - 3 + n
            # psi-degree of the lowest flat term of the expansion (each frame contributes its least z-power)
            low = sum(min(k for (k, c) in FR[cb]) for cb in combo)
            dmax = dim - low
            skipped = dmax > 6
            st['skip'] += skipped
            if any(c in sens for c in combo):
                st['sens'] += 1; st['sens_skip'] += skipped; st['sens_tested'] += (not skipped)
                if not skipped and mind is None: mind = (combo, m, e)
        log("genus %d, <= %d insertions: %d singular lists, %d skipped by the deficit cut; C4-sensitive frames (weight <= %d): %s; "
            "singular lists containing one: %d (%d skipped, %d testable), first testable: %s"
            % (h, NT, st['sing'], st['skip'], wmax, sens[:6], st['sens'], st['sens_skip'], st['sens_tested'], mind))
    for k, a in ((5, 0), (6, 0), (5, 4), (6, 4), (7, 3)):
        log("frame(%d,%d) = %s" % (k, a, frames_upto(w(k, a))[(k, a)]))
    LOG.close()


if __name__ == "__main__":
    main()
