"""cr5_cyc_run.py -- Stage 5 critic: checks of the critic's proof of (CYC) (STAGE5_SYNTHESIS_2 s.9), exact over Q.

Uses only cr5_cyc.py (critic's own module; stdlib).  Frobenius basis of PF4 at the conifold u = 0 (exponents 0, 1, 1, 2):
analytic f0, f1, f2 (leading 1, u, u^2) and ylog = g log u + h, g = f1 + gamma f2, gamma FOUND from the log obstruction.
Checks: (a) PF4 is formally self-adjoint; (b) the concomitant B is constant, alternating and non-degenerate on the solutions, and
B(f_i, g) = 0 for every analytic f_i; (c) T_k := Sum_ij Binv_ij y_i D^k y_j = 0 (k = 0, 1, 2) and = 1/u (k = 3); (d) B(varpi0, g) = 0
in a frame, the control pairs not isotropic; (e) (CYC): M y = 0
for all four solutions and M u^n = PF4 u^n (n = 0..4) in three conifold frames.  Controls (must fail): non-self-adjoint perturbation;
top part of the concomitant only; t = f1/varpi0 (gamma = 0, non-isotropic); g with gamma = 1/2.
Usage: py -I -B cr5_cyc_run.py N   (log logs/cr5_cyc_run_N.log)
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from cr5_cyc import LS, LF, lf, lf_th, ls_th, is_zero, A, apply_op, adjoint, concomitant, INF

N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
LOGREL = 'logs/cr5_cyc_run_%d.log' % N
OUT, NFAIL = [], [0]


def log(s):
    OUT.append(s); print(s)


def check(name, ok, control=False):
    if control:
        log(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        if ok: NFAIL[0] += 1
    else:
        log(('PASS ' if ok else 'FAIL ') + name)
        if not ok: NFAIL[0] += 1


t0 = time.time()
I = lambda n: n * (n - 1) ** 2 * (n - 2)
# columns of PF4 on monomials: PF4 u^n = Sum_m c[n][m] u^m, m in n-3..n+1
c = {n: apply_op(A, lf(LS({n: Fr(1)}))).F[0].c for n in range(0, N + 6)}
check('PF4 maps power series to power series (no negative powers from u^n, n >= 0)',
      all(m >= 0 for n in c for m in c[n]))
check('indicial: [u^(n-3)] PF4 u^n = n(n-1)^2(n-2) for n <= N+5', all(c[n].get(n - 3, 0) == I(n) for n in range(3, N + 6)))


def solve(rhs, a0):  # PF4 (Sum a_n u^n) = -rhs  (rhs a dict of Laurent coefficients), a_0..a_2 given
    a = dict(a0)
    for m in range(-3, 0):
        assert sum(c[n].get(m, 0) * a.get(n, 0) for n in range(0, 3)) + rhs.get(m, 0) == 0, ('obstruction', m)
    for n in range(3, N):
        m = n - 3
        s = sum(c[k].get(m, 0) * a.get(k, 0) for k in range(max(0, n - 4), n)) + rhs.get(m, 0)
        a[n] = -s / I(n)
    return LS({k: v for k, v in a.items()}, N)


f = [solve({}, {0: Fr(int(i == 0)), 1: Fr(int(i == 1)), 2: Fr(int(i == 2))}) for i in range(3)]
for i in range(3):
    check('f%d solves PF4 to u^%d' % (i, N - 6), is_zero(apply_op(A, lf(f[i])), N - 6))
log('f0 = 1 + %s u^3 + ..., f1 = u + %s u^3 + ..., f2 = u^2 + %s u^3 + ...' % (f[0].c.get(3), f[1].c.get(3), f[2].c.get(3)))
R = [apply_op(A, LF({1: f[i]})) for i in (1, 2)]  # PF4(f_i log u) = log-free part (PF4 f_i = 0)
for i in (0, 1):
    check('PF4(f%d log u) has no log part to u^%d' % (i + 1, N - 6), is_zero(LF({1: R[i].F[1]}), N - 6))
r1, r2 = R[0].F[0], R[1].F[0]
neg = [m for m in range(-4, 0) if r1.c.get(m, 0) or r2.c.get(m, 0)]
m0 = neg[-1]
gam = -r1.c.get(m0, 0) / r2.c.get(m0, 0)
log('log obstruction at u^%d gives gamma = %s (UR, POLY: 7/10)' % (m0, gam))
check('gamma = 7/10, and all obstructions at u^-3, u^-2, u^-1 vanish',
      gam == Fr(7, 10) and all(r1.c.get(m, 0) + gam * r2.c.get(m, 0) == 0 for m in (-3, -2, -1)))
g = f[1] + f[2].sc(gam)
rhs = (r1 + r2.sc(gam)).c
h = solve(rhs, {0: Fr(0), 1: Fr(0), 2: Fr(0)})
ylog = LF({1: g, 0: h})
check('ylog = g log u + h solves PF4 to u^%d' % (N - 7), is_zero(apply_op(A, ylog), N - 7))
Y = [lf(f[0]), lf(f[1]), lf(f[2]), ylog]

# (a) self-adjointness
Aad = adjoint(A)
check('(a) PF4* = PF4 (formal adjoint for D, measure dq/q), exact',
      all(not (Aad.get(k, LS()) - A[k]).c for k in A) and all(not v.c for k, v in Aad.items() if k not in A))
Ap = dict(A); Ap[1] = A[1] + LS({1: Fr(1)})
Apad = adjoint(Ap)
check('(a) PF4 + u D is self-adjoint', all(not (Apad.get(k, LS()) - Ap[k]).c for k in Ap), control=True)

# (b) concomitant
P = N - 8
Bm = [[None] * 4 for _ in range(4)]
const_ok, logfree_ok = True, True
for i in range(4):
    for j in range(4):
        b = concomitant(A, Y[i], Y[j])
        for k, v in b.F.items():
            if k and not is_zero(LF({0: v}), P): logfree_ok = False
        v0 = b.F[0]
        if any(e != 0 and e < P for e in v0.c): const_ok = False
        Bm[i][j] = v0.c.get(0, Fr(0))
check('(b) B(y_i, y_j) is log-free to u^%d' % P, logfree_ok)
check('(b) B(y_i, y_j) is a constant (Lagrange identity) to u^%d' % P, const_ok)
check('(b) B is alternating', all(Bm[i][j] == -Bm[j][i] for i in range(4) for j in range(4)))
Bg = [Bm[i][1] + gam * Bm[i][2] for i in range(4)]
check('(b) B(f_i, g) = 0 for i = 0, 1, 2 (analytic solutions are B-orthogonal to the vanishing period)', Bg[:3] == [0, 0, 0])
check('(b) B(ylog, g) != 0 and B(f0, f2) != 0', Bg[3] != 0 and Bm[0][2] != 0)
log('B matrix (basis f0, f1, f2, ylog): %s' % [[str(x) for x in row] for row in Bm])


def minv(M):  # Gauss-Jordan over Q; returns (det, inverse)
    n = len(M); a = [list(map(Fr, r)) + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(M)]; det = Fr(1)
    for col in range(n):
        piv = next(r for r in range(col, n) if a[r][col] != 0)
        if piv != col: a[col], a[piv] = a[piv], a[col]; det = -det
        det *= a[col][col]; pv = a[col][col]; a[col] = [x / pv for x in a[col]]
        for r in range(n):
            if r != col and a[r][col] != 0:
                fac = a[r][col]; a[r] = [x - fac * y for x, y in zip(a[r], a[col])]
    return det, [row[n:] for row in a]


det, Bi = minv(Bm)
check('(b) B is non-degenerate (det = %s)' % det, det != 0)
Bt = concomitant({4: A[4]}, Y[0], Y[2]).F[0]
check('(b) top-order part of the concomitant alone is constant on f0, f2', not any(e != 0 and e < P for e in Bt.c), control=True)

# (c) Griffiths transversality T_k = (P^-1)_(0,k)
thY = [Y]
for k in range(3):
    thY.append([lf_th(y) for y in thY[-1]])
Q = N - 9
for k in range(4):
    T = None
    for i in range(4):
        for j in range(4):
            if Bi[i][j]:
                t = Y[i].mul(thY[k][j]).sc(Bi[i][j]); T = t if T is None else T + t
    target = LF({0: LS({-1: Fr(1)})}) if k == 3 else LF({0: LS()})
    diff = T + target.sc(-1)
    check('(c) T_%d = Sum Binv_ij y_i D^%d y_j = %s to u^%d' % (k, k, '1/u' if k == 3 else '0', Q), is_zero(diff, Q))


# (e) (CYC) in three conifold frames, varpi0 = f0 + a f1 + b f2, t = g/varpi0
def M_apply(y, w0, tt):
    iw = w0.inv(); Dt = ls_th(tt); iD = Dt.inv(); mid = LS({1: Fr(1)}) * w0 * w0 * Dt * Dt
    s = y.mls(iw); s = lf_th(s); s = s.mls(iD); s = lf_th(s); s = s.mls(mid); s = lf_th(s); s = s.mls(iD); s = lf_th(s)
    return s.mls(iw)


Z = N - 10
for (a, b) in [(Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5))]:
    w0 = f[0] + f[1].sc(a) + f[2].sc(b)
    tt = g * w0.inv()
    ok_sol = all(is_zero(M_apply(y, w0, tt), Z) for y in Y)
    check('(e) frame (%s, %s): M y = 0 for f0, f1, f2, ylog to u^%d' % (a, b, Z), ok_sol)
    ok_mon = True
    for n in range(5):
        un = lf(LS({n: Fr(1)}))
        d = M_apply(un, w0, tt) + apply_op(A, un).sc(-1)
        if not is_zero(d, Z): ok_mon = False
    check('(e) frame (%s, %s): M u^n = PF4 u^n for n = 0..4 to u^%d' % (a, b, Z), ok_mon)
w0 = f[0] + f[1].sc(Fr(1, 3)) + f[2].sc(Fr(-2, 7))
check('(e) t = f1/varpi0 (gamma = 0, non-isotropic pair): M y = 0 for all four solutions',
      all(is_zero(M_apply(y, w0, f[1] * w0.inv()), Z) for y in Y), control=True)
g2 = f[1] + f[2].sc(Fr(1, 2))
check('(e) t = (f1 + f2/2)/varpi0 (gamma = 1/2): M y = 0 for all four solutions',
      all(is_zero(M_apply(y, w0, g2 * w0.inv()), Z) for y in Y), control=True)
Bw = lambda j: Bm[0][j] + Fr(1, 3) * Bm[1][j] + Fr(-2, 7) * Bm[2][j]  # B(varpi0, f_j), varpi0 = f0 + f1/3 - 2f2/7
log('control pairs: B(varpi0, f1) = %s, B(varpi0, f1 + f2/2) = %s (non-zero: not isotropic); B(varpi0, g) = %s'
    % (Bw(1), Bw(1) + Fr(1, 2) * Bw(2), Bw(1) + gam * Bw(2)))
check('(d) B(varpi0, g) = 0 in the frame (1/3, -2/7), and the two control pairs are not isotropic',
      Bw(1) + gam * Bw(2) == 0 and Bw(1) != 0 and Bw(1) + Fr(1, 2) * Bw(2) != 0)
log('%s   (%d s)' % ('CR5_CYC_RUN %d: ALL CHECKS PASSED' % N if NFAIL[0] == 0 else 'CR5_CYC_RUN %d: %d FAILURES' % (N, NFAIL[0]),
                    time.time() - t0))
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
with open(os.path.join(HERE, LOGREL), 'w') as fh:
    fh.write('# %s (py -I -B cr5_cyc_run.py %d)\n' % (LOGREL, N) + '\n'.join(OUT) + '\n')
