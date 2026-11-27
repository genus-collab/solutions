# -*- coding: utf-8 -*-
"""rf_run_amp.py -- compute and cache referee amplitudes amp(g, n) at both primes.  Run: py -I -B rf_run_amp.py g,n [g,n ...]"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat as rf
import rf_gsum as gsm

if __name__ == "__main__":
    m1, m2 = rf.models()
    for arg in sys.argv[1:]:
        g, n = map(int, arg.split(","))
        for m in (m1, m2):
            t0 = time.time(); A = gsm.amp(m, g, n, verbose=True)
            print("DONE amp(%d,%d) p=2^%d-1: %d entries, %.0f s" % (g, n, m.p.bit_length(), len(A), time.time() - t0), flush=True)
