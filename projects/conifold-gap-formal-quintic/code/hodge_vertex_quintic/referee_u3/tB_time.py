import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
import random, time
from evlib import *
p = P1; rng = random.Random(1)
n1 = [rng.randrange(p) for _ in range(NSP)]
t = time.time(); a, md = evaluate("61", [n1], [2, 5, 3, 6, 7, 8, 17, 9, 10, 11, 12, 18, 19, 20, 32], mode='u', NN=2, verbose=True); print("NN=2: %.0f s, maxdeg %s" % (time.time() - t, md))
