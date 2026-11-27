# Independent sanity check (referee unit U04) of the uniqueness mechanism of Theorem U (b) /
# lem:gau-lead at n = 0, by brute-force linear algebra on truncated series:
#   unknowns P_0, P_1 in gl_5(Q[u]/u^N);
#   equations  A P_0 - P_0 At = 0   and   P_0' = A P_1 - P_1 At   at orders u^{-1}, ..., u^{N-2}.
# Prediction of the lemma: the projection of the solution space to (P_0 mod u^M) is a line
# (spanned by d*Theta_U).  Without the second equation the projection has dimension 5*M-ish
# (the integral commutant), so the test is not vacuous.
from fractions import Fraction as Fr
import random, sys

def rank(rows, ncols):
    """rank of a list of sparse rows (dict col->Fr) by fraction-exact elimination"""
    piv = {}
    r = 0
    for row in rows:
        row = dict(row)
        while row:
            c = min(row)
            if c in piv:
                p = piv[c]; f = row[c]
                for k, v in p.items():
                    nv = row.get(k, 0) - f * v
                    if nv == 0: row.pop(k, None)
                    else: row[k] = nv
            else:
                f = row[c]
                piv[c] = {k: v / f for k, v in row.items()}
                r += 1
                break
    return r

def series_inv(f, n):
    g = [Fr(0)] * n; g[0] = 1 / f[0]
    for i in range(1, n):
        g[i] = -sum(f[k] * g[i - k] for k in range(1, min(i, len(f) - 1) + 1)) / f[0]
    return g
def series_mul(f, g, n):
    h = [Fr(0)] * n
    for i, a in enumerate(f[:n]):
        if a:
            for j, b in enumerate(g[:n - i]): h[i + j] += a * b
    return h

def run(N, M, seed):
    random.seed(seed)
    def runit():
        return [Fr(random.choice([1, -1, 2, -3]))] + [Fr(random.randint(-3, 3), random.randint(1, 3)) for _ in range(3)]
    # unit parts; index 1 carries the factor u^{-1}
    bt = [runit() for _ in range(5)]
    b = [runit() for _ in range(4)]
    n = N + 3
    num = [Fr(1)]
    for a in range(5): num = series_mul(num, bt[a] + [Fr(0)] * n, n)
    den = [Fr(1)]
    for a in range(4): den = series_mul(den, b[a] + [Fr(0)] * n, n)
    b.append(series_mul(num, series_inv(den, n), n))      # (EP): b_4 = Bt/(b_0 b_1 b_2 b_3)
    def co(x, a, k):   # coefficient of u^k of b_a (x = b or bt), pole at a = 1
        kk = k + 1 if a == 1 else k
        return x[a][kk] if 0 <= kk < len(x[a]) else Fr(0)
    def var(k, a, c, nn): return ((k * 5 + a) * 5 + c) * N + nn
    nv = 2 * 25 * N
    rows0, rows1 = [], []
    for a in range(5):
        for c in range(5):
            for o in range(-1, N - 1):
                for k, rows in ((0, rows0), (1, rows1)):
                    row = {}
                    # (A P_k)_{a,c} = b_a (P_k)_{a+1,c};  (P_k At)_{a,c} = (P_k)_{a,c-1} bt_{c-1}
                    for nn in range(N):
                        x = co(b, a, o - nn)
                        if x: row[var(k, (a + 1) % 5, c, nn)] = row.get(var(k, (a + 1) % 5, c, nn), 0) + x
                        y = co(bt, (c - 1) % 5, o - nn)
                        if y: row[var(k, a, (c - 1) % 5, nn)] = row.get(var(k, a, (c - 1) % 5, nn), 0) - y
                    if k == 1 and o >= 0:
                        row[var(0, a, c, o + 1)] = row.get(var(0, a, c, o + 1), 0) - (o + 1)   # - (P_0')_o
                    row = {kx: v for kx, v in row.items() if v != 0}
                    if row: rows.append(row)
    cons = [{var(0, a, c, nn): Fr(1)} for a in range(5) for c in range(5) for nn in range(M)]
    def projdim(rows):
        return (nv - rank(rows, nv)) - (nv - rank(rows + cons, nv))
    # (i) commutant only: restrict to P_0 (set P_1 free -> irrelevant): use rows0 only
    d_comm = projdim(rows0)
    d_full = projdim(rows0 + rows1)
    return d_comm, d_full

for N in (7, 9, 11):
    M = 2
    dc, df = run(N, M, 777)
    print("N=%d, M=%d: dim proj (commutant equation only) = %d ; with solvability of the z^1 equation = %d" % (N, M, dc, df))
