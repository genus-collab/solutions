"""Synthesis cross-checks on the STORED genus-4 tables (read only; a few seconds).

Run:  py -I syn_check.py        (writes syn_check.log next to this file)

Nothing is recomputed from the rows.  The script only compares stored result
and part files of genus4_split.py with each other:
 (1) the non-zero rows of each result table;
 (2) stored / non-zero piece counts of the part files;
 (3) control consistency (2^89-1): for the three one-edge outer graphs
     [u^5] (tau_5 shifted by 1) = (9/4) [u^1] (bare M_0); both are the first-order
     response to a change of tau (delta tau_1 = +4/9 for M_0, since
     tau = 1 - 4u/9 + ...), so they must be the same linear functional;
 (4) the totals of all variants agree, and the vertex row of the tau_5 control
     differs from the row of record by minus the three one-edge terms;
 (5) second frame: rationals printed at 2^89-1 reduce to the residues at 2^61-1.
"""
import sys
sys.dont_write_bytecode = True
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, "..", ".."))
P = {"89": 2 ** 89 - 1, "61": 2 ** 61 - 1}
OUT = []


def say(s):
    print(s)
    OUT.append(s)


def load(name):
    with open(os.path.join(BASE, name)) as f:
        return json.load(f)


def rows(tab, p):
    return {int(k): [int(x) % p for x in r] for k, r in tab.items()}


def red(fr, p):
    return fr.numerator % p * pow(fr.denominator, -1, p) % p


def main():
    ok = True
    res = {}
    for tag in ("89_nat", "61_nat", "89_id", "61_id", "89_nat_ab1o3_m2o7",
                "61_nat_ab1o3_m2o7", "89_M0", "89_nat_t5"):
        p = P[tag[:2]]
        d = load("g4_split_g4_%s.json" % tag)
        t = rows(d["table"], p)
        res[tag] = (t, [int(x) % p for x in d["total"]])
        nz = sorted(k for k, r in t.items() if any(r))
        low = {}
        for k in nz:
            o = next(i for i, x in enumerate(t[k]) if x)
            low[o] = low.get(o, 0) + 1
        say("(1) %-20s rows stored %3d, non-zero rows %3d %s" % (
            tag, len(t), len(nz), nz if len(nz) <= 6 else "lowest order histogram %s" % low))
        pt = load(os.path.join("parts", "g4_%s_run.json" % tag))["res"]
        st = sum(len(v) for v in pt.values())
        nzp = sum(1 for v in pt.values() for r in v.values() if any(int(x) % p for x in r))
        say("(2) %-20s inner graphs %d, stored pieces %d, non-zero pieces %d" % (tag, len(pt), st, nzp))
    p = P["89"]
    m0, t5, nat = res["89_M0"][0], res["89_nat_t5"][0], res["89_nat"][0]
    for k in (1, 2, 3):
        c = t5[k][5] == m0[k][1] * 9 * pow(4, -1, p) % p
        ok &= c
        say("(3) one-edge graph %d: [u^5](nat_t5) = (9/4) [u^1](M0): %s" % (k, "PASS" if c else "FAIL"))
    c = res["89_M0"][1] == res["89_nat_t5"][1] == res["89_nat"][1] == res["89_id"][1]
    ok &= c
    say("(4) totals of nat, id, M0, nat_t5 equal at 2^89-1: %s" % ("PASS" if c else "FAIL"))
    c = res["61_nat"][1] == res["61_id"][1]
    ok &= c
    say("(4) totals of nat, id equal at 2^61-1: %s" % ("PASS" if c else "FAIL"))
    c = all((t5[0][i] - nat[0][i]) % p == 0 for i in range(5)) and \
        (t5[0][5] - nat[0][5] + t5[1][5] + t5[2][5] + t5[3][5]) % p == 0
    ok &= c
    say("(4) vertex row of nat_t5 = row of record - one-edge [u^5] terms: %s" % ("PASS" if c else "FAIL"))
    ab89, ab61 = res["89_nat_ab1o3_m2o7"], res["61_nat_ab1o3_m2o7"]
    want = [F(1, 180000), F(-11, 900000), F(-79, 378000000), F(20081, 2126250000)]
    c = all(red(w, P["89"]) == ab89[1][i] and red(w, P["61"]) == ab61[1][i] for i, w in enumerate(want))
    ok &= c
    say("(5) second frame, total [u^0..u^3] = 1/180000, -11/900000, -79/378000000, "
        "20081/2126250000 at both primes: %s" % ("PASS" if c else "FAIL"))
    say("SYN_CHECK %s" % ("ALL PASS" if ok else "FAILURES"))
    with open(os.path.join(HERE, "syn_check.log"), "w") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
