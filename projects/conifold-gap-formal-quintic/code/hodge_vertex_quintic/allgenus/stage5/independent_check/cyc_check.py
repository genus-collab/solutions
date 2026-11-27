# cyc_check.py N -- independent check (main session, 2026-11-14) of Theorem CYC (stage5/critic/CYC_PROOF.md):
#   PF4 = (1/y0) D (1/Dt) D (u y0^2 Dt^2) D (1/Dt) D (1/y0)   in every conifold frame,
# with D = q d/dq = -(1-u) d/du, u = 1 - 5^5 q, PF4 = D^4 - 5^5 q prod_(k=1..4)(D + k/5), t = g/y0, g the vanishing period.
# Own code: truncated power series in u over Q (no project, worker, referee or critic module is imported or read). Method:
# Frobenius solutions from the recursion, the vanishing period g found from the log obstruction (Laurent bookkeeping of
# PF4(g log u)), then both operators applied to u^0..u^5 and to the three power-series solutions.
# Controls (must fail): non-isotropic pairs t = f1/y0, t = (f1 + f2/2)/y0, and g -> g + f2/7.
import sys, time, itertools
from fractions import Fraction as Fr
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
P = N + 12                      # working precision: coefficients u^0..u^(P-1); each D loses one order
T0 = time.time()


def add(a, b): return [x + y for x, y in zip(a, b)]
def scal(c, a): return [c * x for x in a]


def mul(a, b):
    n = min(len(a), len(b)); r = [Fr(0)] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j in range(n - i):
                r[i + j] += x * b[j]
    return r


def inv(a):
    assert a[0] != 0
    n = len(a); r = [Fr(0)] * n; r[0] = 1 / a[0]
    for k in range(1, n):
        r[k] = -sum(a[j] * r[k - j] for j in range(1, k + 1)) / a[0]
    return r


def ddu(a): return [(k + 1) * a[k + 1] for k in range(len(a) - 1)]      # loses one order


def D(a):                                                               # D = -(1-u) d/du
    d = ddu(a); r = [-x for x in d]
    for k in range(1, len(d)):
        r[k] += d[k - 1]
    return r


def mono(n, p): return [Fr(1) if k == n else Fr(0) for k in range(p)]
def U(a): return [Fr(0)] + a[:-1]                                       # times u (same length; valid part grows by one)


def PF4(f):
    d1 = D(f); d2 = D(d1); d3 = D(d2); d4 = D(d3)
    n = len(d4); f, d1, d2, d3 = f[:n], d1[:n], d2[:n], d3[:n]
    inner = add(add(add(scal(Fr(2), d3), scal(Fr(7, 5), d2)), scal(Fr(2, 5), d1)), scal(Fr(24, 625), f))
    omu = [Fr(1), Fr(-1)] + [Fr(0)] * (n - 2)
    return add(U(d4), scal(Fr(-1), mul(omu, inner)))


# ---- Frobenius basis: [u^(n-3)] PF4 u^n = n (n-1)^2 (n-2); u^0, u^1, u^2 free ----
def solve_ps(a0, a1, a2, p):
    f = [Fr(a0), Fr(a1), Fr(a2)] + [Fr(0)] * (p - 3)
    for n in range(3, p):
        cur = PF4(f[:n + 5] if n + 5 <= p else f + [Fr(0)] * (n + 5 - p))
        f[n] = -cur[n - 3] / (n * (n - 1) ** 2 * (n - 2))
    return f


t1 = time.time()
f0 = solve_ps(1, 0, 0, P); f1 = solve_ps(0, 1, 0, P); f2 = solve_ps(0, 0, 1, P)
print('Frobenius basis to u^%d (%.1f s): f0 = 1 + (%s) u^3 + .., f1 = u + (%s) u^3 + .., f2 = u^2 + (%s) u^3 + ..'
      % (P - 1, time.time() - t1, f0[3], f1[3], f2[3]))
fails = []
for nm, f in (('f0', f0), ('f1', f1), ('f2', f2)):
    r = PF4(f); bad = [k for k, x in enumerate(r) if x]
    print('%s PF4 %s = 0 through u^%d' % ('PASS' if not bad else 'FAIL', nm, len(r) - 1))
    if bad: fails.append('PF4 ' + nm)

# ---- the log obstruction. Elements A + B l with l = log u, D l = 1 - 1/u; A, B Laurent, stored with shift S ----
S = 4


def Dlaur(X):                       # X = sum_k c_k u^(k-S); D u^e = -e u^(e-1) + e u^e
    r = [Fr(0)] * len(X)
    for k, c in enumerate(X):
        e = k - S
        if c and e:
            if k == 0: raise SystemExit('pole beyond the Laurent window')
            r[k - 1] -= e * c; r[k] += e * c
    return r


def Dpair(pr):
    A, B = pr; nA = Dlaur(A); nB = Dlaur(B)
    for k, c in enumerate(B):       # + B (1 - 1/u)
        if c:
            if k == 0: raise SystemExit('pole beyond the Laurent window')
            nA[k] += c; nA[k - 1] -= c
    return nA, nB


def PF4_pair(pr):
    d = [pr]
    for _ in range(4): d.append(Dpair(d[-1]))
    out = []
    for i in (0, 1):
        f, a1, a2, a3, a4 = d[0][i], d[1][i], d[2][i], d[3][i], d[4][i]
        n = len(f); o = [Fr(0)] * n
        for k in range(n):
            inner = 2 * a3[k] + Fr(7, 5) * a2[k] + Fr(2, 5) * a1[k] + Fr(24, 625) * f[k]
            o[k] -= inner
            if k >= 1:
                o[k] += a4[k - 1]           # u D^4
                o[k] += 2 * a3[k - 1] + Fr(7, 5) * a2[k - 1] + Fr(2, 5) * a1[k - 1] + Fr(24, 625) * f[k - 1]
        out.append(o)
    return out


def obstruction(g):
    W = 30
    A, B = PF4_pair(([Fr(0)] * (S + W), [Fr(0)] * S + g[:W]))
    return [A[S + e] for e in (-3, -2, -1)]


a, b, c = obstruction(f0), obstruction(f1), obstruction(f2)
print('log obstruction [u^-3, u^-2, u^-1] PF4(f_i log u): f0 %s, f1 %s, f2 %s' % ([str(x) for x in a], [str(x) for x in b], [str(x) for x in c]))
sol = None
for i, j in itertools.combinations(range(3), 2):
    det = a[i] * c[j] - a[j] * c[i]
    if det:
        sol = ((-b[i] * c[j] + b[j] * c[i]) / det, (-a[i] * b[j] + a[j] * b[i]) / det); break
c0, gam = sol
ok = all(c0 * a[k] + b[k] + gam * c[k] == 0 for k in range(3))
print('%s vanishing period g = (%s) f0 + f1 + (%s) f2 (all three obstructions vanish)' % ('PASS' if ok else 'FAIL', c0, gam))
if not ok or c0 != 0 or gam != Fr(7, 10): fails.append('gamma')
g = add(add(scal(c0, f0), f1), scal(gam, f2))


# ---- the factorised operator in a frame ----
def M_op(fr, y0, tt):
    Dt = D(tt); n = len(Dt)
    iy0 = inv(y0[:n]); iDt = inv(Dt); y0n = y0[:n]
    mid = U(mul(mul(y0n, y0n), mul(Dt, Dt)))
    x = mul(fr[:n], iy0); x = D(x); x = mul(x, iDt[:len(x)]); x = D(x); x = mul(x, mid[:len(x)]); x = D(x)
    x = mul(x, iDt[:len(x)]); x = D(x); x = mul(x, iy0[:len(x)])
    return x


def check_frame(label, y0, gg, expect):
    tt = mul(gg, inv(y0)); ok = True; first = None
    for n in range(6):
        m = M_op(mono(n, P), y0, tt); p4 = PF4(mono(n, P)); k = min(len(m), len(p4), N + 1)
        d = [e for e in range(k) if m[e] != p4[e]]
        if d: ok = False; first = 'u^%d on u^%d' % (d[0], n); break
    if ok:
        for nm, f in (('f0', f0), ('f1', f1), ('f2', f2)):
            m = M_op(f, y0, tt); d = [e for e in range(min(len(m), N + 1)) if m[e] != 0]
            if d: ok = False; first = 'M %s at u^%d' % (nm, d[0]); break
    good = (ok == expect)
    if not good: fails.append(label)
    print('%s %-58s M = PF4 on u^0..u^5, M f_i = 0, through u^%d: %s%s' % (
        'PASS' if good else 'FAIL', label, N, 'holds' if ok else 'fails (first %s)' % first, '' if expect else '  [control]'))


t2 = time.time()
for (a_, b_) in ((0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5)), (Fr(-5, 2), Fr(9, 4))):
    y0 = add(add(f0, scal(Fr(a_), f1)), scal(Fr(b_), f2))
    check_frame('frame y0 = f0 + (%s) f1 + (%s) f2, t = g/y0' % (a_, b_), y0, g, True)
    check_frame('  same frame, g -> 3 g', y0, scal(Fr(3), g), True)
y0 = add(add(f0, scal(Fr(1, 3), f1)), scal(Fr(-2, 7), f2))
check_frame('control: t = f1/y0 (non-isotropic pair)', y0, f1, False)
check_frame('control: t = (f1 + f2/2)/y0', y0, add(f1, scal(Fr(1, 2), f2)), False)
check_frame('control: g -> g + f2/7', y0, add(g, scal(Fr(1, 7), f2)), False)
print('operator checks %.1f s' % (time.time() - t2))
print('CYC_CHECK %d: %s   (%.0f s)' % (N, 'ALL CHECKS PASSED' if not fails else 'FAILED: %s' % fails, time.time() - T0))
