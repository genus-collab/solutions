import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.dont_write_bytecode = True
from evlib import *
from fq3_steps import frame_of
for tag, p in (("61", P1), ("89", P2)):
    st = setup(frame_of(p, (0, 0), 34, 4)); H = st['H']; R = st['R']
    nz = {k: sum(1 for a in range(5) for i in range(5) if H.Q[a][i][0].coef(k)) for k in range(4)}
    # Psi(u) diagonal between Fourier modes and flat index: sum_i Q[a][i][0] zeta^(-m i) non-zero only for one m per a
    pat = {}
    for a in range(5):
        pat[a] = sorted(set(m for m in range(5) for k in range(4) if sum(H.Q[a][i][0].coef(k) * R.z(-m * i) for i in range(5)) % p))
    print("p=2^%s-1: number of non-zero entries of [u^k] Psi(u) (zhat^0 rows), k = 0..3: %s; Fourier support of row a (flat index a -> modes m): %s" % (tag, nz, pat))
