# -*- coding: utf-8 -*-
"""rv_tc.py -- independent form of the GAP / VERTEX checks of genus4_split.checks, and a read-only demonstration of the part-file
prefix collision of the --collect step.

(1) u^-6 G(u) = c t_c^-6 + O(t_c^0)  <=>  G(u) = c (u/t_c(u))^6 mod u^6.  With own series arithmetic (no fq_series compose/reverse):
    the stored totals and vertex-graph rows (g4_split_g4_<p>_<variant>.json, both primes; nat and id; also every other g4_split_g4_*.json
    present) are compared with (1/180000) (u/t_c)^6, t_c from fq_conifold.Conifold (the frame is read from the tag).
(2) t_c(u) in the frame (0, 0) printed as rationals (for comparison with the literature).
(3) the file filter of genus4_split.py --collect ( f.startswith(tag + "_") and f.endswith(".json") ) applied to the present content
    of parts/ for each tag: lists the tags whose collect step would merge part files of a DIFFERENT tag.
Usage: py -I rv_tc.py   -> out/tc.log
"""
import os, sys, json, re
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT
from fractions import Fraction as Fr
from fq_fp import FpRing, P1, P2
from fq_conifold import Conifold
LOG = open(os.path.join(OUT, "tc.log"), "w", encoding="utf-8"); FAILED = []
def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)
K = 6
def gap_series(R, al, be):
    """coefficients u^0..u^5 of (1/180000) (u/t_c)^6"""
    p = R.p; C = Conifold(R, Nu=K + 4, alpha=al, beta=be); tc = C.tc
    assert tc.v == 1 and tc.prec >= K + 1, (tc.v, tc.prec)
    f = [R.of(tc.coef(k + 1)) for k in range(K)]                      # t_c/u
    inv = [pow(f[0], p - 2, p)]
    for k in range(1, K): inv.append((-inv[0] * sum(f[j] * inv[k - j] for j in range(1, k + 1))) % p)
    out = [1] + [0] * (K - 1)
    for _ in range(6): out = [sum(out[i] * inv[k - i] for i in range(k + 1)) % p for k in range(K)]
    c = R.of(Fr(1, 180000)); return [x * c % p for x in out], f
for fn in sorted(os.listdir(PROJ)):
    m = re.match(r"g4_split_g4_(89|61)_(.+)\.json$", fn)
    if not m: continue
    prime, rest = m.group(1), m.group(2); R = FpRing(P2 if prime == "89" else P1); p = R.p
    ab = re.search(r"_ab(.+?)_(.+?)(_K\d+)?$", rest)
    def fr(s): return Fr(s.replace("m", "-").replace("o", "/"))
    al, be = (fr(ab.group(1)), fr(ab.group(2))) if ab else (Fr(0), Fr(0))
    variant = rest.split("_ab")[0]
    d = json.load(open(os.path.join(PROJ, fn))); gap, f = gap_series(R, al, be)
    tot = [x % p for x in d['total']]; v0 = [x % p for x in d['table'].get('0', [0] * K)]
    okT = tot == gap; okV = v0 == gap
    nzrows = [i for i, c in d['table'].items() if i != '0' and any(x % p for x in c)]
    log("%s: variant %s, frame (%s, %s): total = gap series: %s; vertex row = gap series: %s; non-zero rows with edges: %d; failed list of the run: %s" % (
        fn, variant, al, be, okT, okV, len(nzrows), d['failed']))
    check("%s: total = (1/180000)(u/t_c)^6 mod u^6" % fn, okT)
    if variant == 'nat': check("%s: vertex graph row = (1/180000)(u/t_c)^6 mod u^6 and no other non-zero row" % fn, okV and not nzrows)
    elif variant != 'id': log("   (control variant %s: vertex row equal to the gap series: %s -- expected False)" % (variant, okV))
R = FpRing(P2); C = Conifold(R, Nu=K + 4, alpha=Fr(0), beta=Fr(0))
log("t_c(u), frame (0,0): " + " + ".join("(%s) u^%d" % (C.tc.coef(k), k) for k in range(1, K + 2)))
# (3) prefix collision
parts = sorted(f for f in os.listdir(os.path.join(PROJ, "parts")) if f.endswith(".json"))
tags = sorted(set(re.sub(r"_(run|task\d+of\d+)\.json$", "", f) for f in parts))
log("part files present: %s" % parts)
coll = {}
for tag in tags:
    got = [f for f in parts if f.startswith(tag + "_") and f.endswith(".json")]
    foreign = [f for f in got if re.sub(r"_(run|task\d+of\d+)\.json$", "", f) != tag]
    if foreign: coll[tag] = foreign
log("tags whose --collect would merge part files of another tag with the filter of the version of 08:44 (BUG B1: prefix filter): %s" % coll)
coll2 = {}
for tag in tags:
    pat = re.compile("^" + re.escape(tag) + r"_(run|task\d+of\d+)\.json$")
    foreign = [f for f in parts if pat.match(f) and re.sub(r"_(run|task\d+of\d+)\.json$", "", f) != tag]
    if foreign: coll2[tag] = foreign
check("with the exact pattern of the version of 09:40 no tag picks up a foreign part file", not coll2, coll2)
log("RV_TC %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
