# -*- coding: utf-8 -*-
"""rv_table.py -- prints the mutation table of REVIEW.md from out/mut_g3.json and out/g4mut_*.json (read-only summary)."""
import os, sys, json, re
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
d = json.load(open(os.path.join(OUT, "mut_g3.json")))
g4 = {}
for f in sorted(os.listdir(OUT)):
    m = re.match(r"g4mut_(.+)_(nat|id)_(89|61)_(full|pop\d)\.json$", f)
    if m:
        r = json.load(open(os.path.join(OUT, f))); short = []
        for x in r['failed']:
            short.append("REG" if x.startswith("REGULARITY") else "VERTEX" if x.startswith("VERTEX") else "GAP" if x.startswith("GAP") else "rows" if x.startswith("rows") else x[:25])
        g4.setdefault(m.group(1), []).append("%s/%s: %s" % (m.group(2), m.group(4), ",".join(short) if short else "all pass"))
print("| mutant | change | genus-2/3 regression (2^89-1 and 2^61-1) | genus 4 (2^89-1) |"); print("|---|---|---|---|")
for name in sorted(d):
    r = d[name]; cols = []
    if r['crash']: cols.append("crash: " + r['crash'][:60])
    for key, lab in (("genus 3, 2^89-1, variant id", "g3 id"), ("genus 3, 2^89-1, variant M0", "g3 M0"), ("genus 3, 2^89-1, variant nat", "g3 nat"),
                     ("genus 2, 2^89-1, variant nat", "g2 nat"), ("genus 2, 2^89-1, variant M0", "g2 M0")):
        if any(l.startswith("FAIL " + key) for l in r['lines']): cols.append(lab)
    same = all(any(l.startswith("FAIL") and l[5:].startswith(k.replace("89", "61")) for l in r['lines']) == any(l.startswith("FAIL") and l[5:].startswith(k) for l in r['lines'])
               for k in ("genus 3, 2^89-1, variant id", "genus 3, 2^89-1, variant M0", "genus 3, 2^89-1, variant nat", "genus 2, 2^89-1, variant nat", "genus 2, 2^89-1, variant M0"))
    verdict = ("FAIL: " + ", ".join(cols)) if cols else "**all 10 pass**"
    if not same: verdict += " (primes differ!)"
    print("| %s | %s | %s | %s |" % (name, r['desc'], verdict, "; ".join(g4.get(name, ["-"]))))
