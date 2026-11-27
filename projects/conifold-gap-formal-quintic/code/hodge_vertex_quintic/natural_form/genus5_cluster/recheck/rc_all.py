# -*- coding: utf-8 -*-
"""rc_all.py -- run the whole re-check in order:  py -I -B rc_all.py

  rc_t0_pf.py          the independent Picard-Fuchs / t_c solution (smoke test)
  rc_t1_recon.py       (1) CRT + rational reconstruction, margins; writes rc_total_rational.json
  rc_t2_tc.py          (2) the gap in t_c, the Jacobian, the perturbation test, the closed form
  rc_t3_census.py      (3) regularity and the census of outer graphs, genus-4 comparison
  rc_t4_coherence.py   (4)/(5) cross-run coherence, the shipped id parts, the controls
  rc_t5_sound.py       (4) the blind spot of the total check

Nothing is imported from the project; the only inputs are cluster_genus5_quintic/results/*.json,
cluster_genus5_quintic/parts/g5_{89,61}_id/local.json and hodge_vertex_quintic/g4_split_g4_89_nat.json.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = ["rc_t0_pf.py", "rc_t1_recon.py", "rc_t2_tc.py", "rc_t3_census.py",
         "rc_t4_coherence.py", "rc_t5_sound.py"]


def main():
    bad = []
    for s in STEPS:
        print("=" * 100)
        print("### %s" % s)
        print("=" * 100)
        r = subprocess.run([sys.executable, "-I", "-B", os.path.join(HERE, s)])
        if r.returncode:
            bad.append(s)
    print("=" * 100)
    print("RC_ALL: %s" % ("ALL STEPS RAN" if not bad else "FAILURES: %s" % bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
