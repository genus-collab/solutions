# -*- coding: utf-8 -*-
"""true u-dependent vertex tensors (UVertex, NN nodes) at the random full-class point of tA_frozen.py (same seed): one configuration per process
usage: py -I -B tB_run.py <tag> <cfg a|b|c> <NN> <graphs comma separated> <outfile>"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import json, random, time
from evlib import *
tag, cfg, NN, gl, outf = sys.argv[1], sys.argv[2], int(sys.argv[3]), [int(x) for x in sys.argv[4].split(",")], sys.argv[5]
p = P2 if tag == "89" else P1; rng = random.Random(4242)
n1 = [rng.randrange(p) for _ in range(NSP)]; n2 = [rng.randrange(p) for _ in range(NSP)]; n2b = [rng.randrange(p) for _ in range(NSP)]; n3 = [rng.randrange(p) for _ in range(NSP)]
Ns = {'a': [n1], 'b': [n1, n2], 'c': [n1, n2, n3]}[cfg]
t = time.time(); out, md = evaluate(tag, Ns, gl, mode='u', NN=NN)
json.dump(dict(cfg=cfg, NN=NN, maxdeg=md, vals={str(i): v for i, v in out.items()}, time=time.time() - t), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), outf), "w"))
print("cfg %s NN %d graphs %s: %.0f s, maxdeg %s" % (cfg, NN, gl, time.time() - t, md))
