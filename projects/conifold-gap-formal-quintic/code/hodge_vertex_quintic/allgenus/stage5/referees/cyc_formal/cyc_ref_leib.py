# -*- coding: utf-8 -*-
"""cyc_ref_leib.py N -- referee cyc_formal: the pure-algebra identities of CYC s.6-s.7 on RANDOM data (exact over Q).

These identities hold for arbitrary series, so they are tested where they are non-trivial (tau_1, tau_2 != 0):
  L1  s_0 = tau_0 = 0 for an antisymmetric matrix (Leibniz on y_j = y_0 f_j):
      T_1 = y_0^2 s_1, T_2 = 2 y_0 Dy_0 s_1 + y_0^2 s_2, T_3 = 3 y_0 D^2y_0 s_1 + 3 y_0 Dy_0 s_2 + y_0^2 s_3;
  L2  with e = Dt, D = e th: s_1 = e tau_1, s_2 = De tau_1 + e^2 tau_2, s_3 = D^2e tau_1 + 3 e De tau_2 + e^3 tau_3;
  L3  for Binv = -J (adapted basis, f_0 = 1, f_1 = t): tau_k = -th^k f_3 + f_3 th^k f_0 - t th^k f_2 + f_2 th^k f_1, and for
      k >= 2 (th^k 1 = th^k t = 0): tau_2 = -th^2 f_3 - t th^2 f_2, tau_3 = -th^3 f_3 - t th^3 f_2;
  L4  operator identities D (1/e) D = e th^2 and M = (e/y_0) th^2 (u y_0^2 e^3) th^2 (1/y_0) (any unit y_0, any t, e = Dt a unit);
      leading coefficient of M is u.
  Controls: T_3 with 2 y_0 D^2y_0 s_1 (wrong binomial); s_3 without the factor 3; M with (u y_0^2 e^2) in the th-form.
Usage: py -I -B cyc_ref_leib.py N    (stored as logs/cyc_ref_leib_N.log)
"""
import os, sys, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U
from cyc_ref_ops import Op, REl

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NF = [0]
DOP = Op.Dk(1)


def check(name, ok, control=False):
    if control:
        print(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        NF[0] += bool(ok)
    else:
        print(('PASS ' if ok else 'FAIL ') + name)
        NF[0] += (not ok)


def rser(v0, N, lead=None):
    c = [Fr(random.randint(-9, 9), random.randint(1, 9)) for _ in range(N - v0)]
    if lead is not None:
        c[0] = Fr(lead)
    return Ser(v0, c, N)


def z(x, M):
    return x.zero_below(M)


t0 = time.time()
random.seed(31415)
y0 = rser(0, N, lead=3)
t = rser(1, N + 1, lead=-2)                    # t = -2u + ..., Dt a unit
fs = [Ser(0, [Fr(1)]), t, rser(0, N), rser(0, N)]
A = [[Fr(0)] * 4 for _ in range(4)]
for i in range(4):
    for j in range(i + 1, 4):
        A[i][j] = Fr(random.randint(-9, 9), random.randint(1, 9)); A[j][i] = -A[i][j]
e = t.D(); ie = e.inv()
D = lambda x: x.D()
th = lambda x: x.D() * ie


def contr(op_list_i, op_list_j, B):            # Sum B_ij X_i Y_j
    acc = Ser(0, [])
    for i in range(4):
        for j in range(4):
            if B[i][j]:
                acc = acc + (op_list_i[i] * op_list_j[j]).sc(B[i][j])
    return acc


def iterate(F, k, fn):
    out = list(F)
    for _ in range(k):
        out = [fn(x) for x in out]
    return out


ys = [y0 * x for x in fs]
s = [contr(fs, iterate(fs, r, D), A) for r in range(4)]
T = [contr(ys, iterate(ys, k, D), A) for k in range(4)]
tau = [contr(fs, iterate(fs, k, th), A) for k in range(4)]
M0 = N - 6
check('L1 s_0 = tau_0 = T_0 = 0 (antisymmetry)', z(s[0], M0) and z(tau[0], M0) and z(T[0], M0))
check('L1 T_1 = y_0^2 s_1; T_2 = 2 y_0 Dy_0 s_1 + y_0^2 s_2; T_3 = 3 y_0 D^2y_0 s_1 + 3 y_0 Dy_0 s_2 + y_0^2 s_3',
      z(T[1] - y0 * y0 * s[1], M0) and z(T[2] - (y0 * D(y0) * s[1]).sc(2) - y0 * y0 * s[2], M0)
      and z(T[3] - (y0 * D(D(y0)) * s[1]).sc(3) - (y0 * D(y0) * s[2]).sc(3) - y0 * y0 * s[3], M0))
check('L1 T_3 with 2 y_0 D^2y_0 s_1', z(T[3] - (y0 * D(D(y0)) * s[1]).sc(2) - (y0 * D(y0) * s[2]).sc(3) - y0 * y0 * s[3], M0),
      control=True)
check('L2 s_1 = e tau_1; s_2 = De tau_1 + e^2 tau_2; s_3 = D^2e tau_1 + 3 e De tau_2 + e^3 tau_3',
      z(s[1] - e * tau[1], M0) and z(s[2] - D(e) * tau[1] - e * e * tau[2], M0)
      and z(s[3] - D(D(e)) * tau[1] - (e * D(e) * tau[2]).sc(3) - e * e * e * tau[3], M0))
check('L2 s_3 without the factor 3', z(s[3] - D(D(e)) * tau[1] - e * D(e) * tau[2] - e * e * e * tau[3], M0), control=True)
check('L2 (random data: tau_1 and tau_2 are non-zero, so L1-L2 are tested non-trivially)', not z(tau[1], M0) and not z(tau[2], M0))
J = [[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]
Bi = [[-x for x in r] for r in J]
tauJ = [contr(fs, iterate(fs, k, th), Bi) for k in range(4)]
ok3 = True
for k in range(4):
    th_k = iterate(fs, k, th)
    rhs = -th_k[3] + fs[3] * th_k[0] - t * th_k[2] + fs[2] * th_k[1]
    ok3 = ok3 and z(tauJ[k] - rhs, M0)
check('L3 tau_k = -th^k f_3 + f_3 th^k f_0 - t th^k f_2 + f_2 th^k f_1 (k = 0..3, Binv = -J)', ok3)
t2, t3 = iterate(fs, 2, th), iterate(fs, 3, th)
check('L3 th^k 1 = th^k t = 0 (k = 2, 3); tau_2 = -th^2 f_3 - t th^2 f_2; tau_3 = -th^3 f_3 - t th^3 f_2',
      all(z(x, M0) for x in (t2[0], t2[1], t3[0], t3[1])) and z(tauJ[2] + t2[3] + t * t2[2], M0)
      and z(tauJ[3] + t3[3] + t * t3[2], M0))
THop = Op({1: ie})
check('L4 D o (1/e) o D = e th^2 as operators', (DOP * Op.fn(ie) * DOP).diff(Op.fn(e) * THop * THop, N - 6))
iy = y0.inv()
Mfac = Op.fn(iy) * (DOP * (Op.fn(ie) * (DOP * (Op.fn(U * y0 * y0 * e * e) * (DOP * (Op.fn(ie) * (DOP * Op.fn(iy))))))))
Mth = Op.fn(e * iy) * THop * THop * Op.fn(U * y0 * y0 * e * e * e) * THop * THop * Op.fn(iy)
p = min(x.N for x in Mfac.A.values())
check('L4 (1/y_0) D (1/e) D (u y_0^2 e^2) D (1/e) D (1/y_0) = (e/y_0) th^2 (u y_0^2 e^3) th^2 (1/y_0)', Mfac.diff(Mth, p - 2))
check('L4 leading coefficient of M is u (random y_0, t)', z(Mfac.A[4] - U, p) and max(Mfac.A) == 4)
Mbad = Op.fn(e * iy) * THop * THop * Op.fn(U * y0 * y0 * e * e) * THop * THop * Op.fn(iy)
check('L4 the th-form with (u y_0^2 e^2)', Mfac.diff(Mbad, p - 2), control=True)
print('CYC_REF_LEIB %d: %s   (%.1f s)' % (N, 'ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
