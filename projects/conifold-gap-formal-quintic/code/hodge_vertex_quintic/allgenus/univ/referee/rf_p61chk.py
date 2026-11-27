# -*- coding: utf-8 -*-
"""rf_p61chk.py -- REFEREE (univ): the worker's normalisation-response rationals ((b - a)/m5, (d - c1)/m5 tables of UNIVERSALITY.md
s.2.3, reconstructed by the worker from 2^89-1 residues) reduced mod 2^61-1 and compared with the referee's own (shift - canon)
per unit m at 2^61-1 (rf_e4.py with RF_PRIME=61; tables in the scratch folder).  Also the vanishing of the sums over graphs.
Run: py -I -B rf_p61chk.py  -> rf_p61chk.log
"""
import os, sys, json
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_e4 import SCR

P61 = 2 ** 61 - 1; P89 = 2 ** 89 - 1
WORKER = {   # outer graph: ([u^4], [u^5]) per unit m5, UNIVERSALITY.md s.2.3
    "q": {0: (F(-1, 8400000), F(-187, 680400000)), 1: (F(-1, 1008000), F(13, 6804000)), 2: (F(1, 900000), F(-53, 32400000)),
          3: (0, F(-41, 68040000)), 4: (0, F(-11, 34020000)), 7: (0, F(49, 16200000)), 8: (0, F(-17, 8100000))},
    "c1": {0: (F(-1, 8400000), F(3673, 567000000)), 1: (F(-1, 1008000), F(1, 280000)), 2: (F(1, 900000), F(-407, 40500000)),
           3: (0, F(169, 17010000)), 4: (0, F(1, 850500)), 7: (0, F(-2, 253125)), 8: (0, F(-13, 4050000))},
}


def red(x, p):
    x = F(x); return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def main():
    LOG = open(os.path.join(HERE, "rf_p61chk.log"), "a", encoding="ascii")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 100)
    for case, tabw in WORKER.items():
        s4 = sum(v[0] for v in tabw.values()); s5 = sum(v[1] for v in tabw.values())
        log("%s: sum over graphs of the worker's per-unit responses: [u^4] %s, [u^5] %s (exact rationals)" % (case, s4, s5))
        for p, suf in ((P61, "_p61"), (P89, "")):
            fs = os.path.join(SCR, "tab_%s_shift_g4_E2%s.json" % (case, suf)); fc = os.path.join(SCR, "tab_%s_canon_g4_E2%s.json" % (case, suf))
            if not (os.path.exists(fs) and os.path.exists(fc)): log("   p = %d: tables missing" % p.bit_length()); continue
            S = {int(k): v for k, v in json.load(open(fs)).items()}; C = {int(k): v for k, v in json.load(open(fc)).items()}
            bad = []; nz = []
            for i in set(S) | set(C):
                d = [(S.get(i, [0] * 6)[k] - C.get(i, [0] * 6)[k]) % p for k in range(6)]
                w = tabw.get(i, (0, 0)); want = [0, 0, 0, 0, red(w[0], p), red(w[1], p)]
                if d != want: bad.append(i)
                if any(d): nz.append(i)
            missing = [i for i in tabw if i not in nz and any(tabw[i])]
            log("   p = 2^%d-1: own (shift - canon) per unit m equals the worker's rationals on %d graphs: %s; graphs that differ %s, missing %s"
                % (p.bit_length(), len(nz), "PASS" if not bad and not missing else "FAIL", bad, missing))


if __name__ == "__main__":
    main()
