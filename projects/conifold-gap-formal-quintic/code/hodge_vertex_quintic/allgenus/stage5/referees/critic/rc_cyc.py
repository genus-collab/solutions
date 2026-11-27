# -*- coding: utf-8 -*-
"""rc_cyc.py N -- critic of the Stage-5 referees: decisive re-check of (CYC) and of the points the two CYC referees raise.

Exact over Q, own code (rc_core). Method new relative to cr5_cyc / pb_cyc / cyc_ref_* / cy_*: operator identities are tested on
u^x for six NON-INTEGER rational x (twisted derivation Dx), which determines an order-<=4 operator over Q((u)).
  A  PF4 from its large-radius form; C1 (self-adjoint); PF5 identity D^5 - q prod_(1..5)(5D+k) = D o PF4 (operator algebra).
  B  Frobenius basis; full 3x3 log-obstruction matrix (cyc_formal E1); gamma; u^-1 row = B(f_0, .) on span(f_1, f_2) (cyc_adversary).
  C  B from the definition of the concomitant: log-free constants; hand values b01, b02, b03, b12, b13, b23 (cyc_formal E2);
     Pfaffian 2, det 4; det Wr = -2/u^2 (Abel).
  E  frame relabelling (E3 = K1); M(y0, y1) = PF4 in programme frames, critic-convention frames, y0(0) = 0 frames (K2), rescaled;
     the Corollary in operator form; unit shapes val(c) (K2); controls that must fail.
Usage: py -I -B rc_cyc.py 40   (log: logs/rc_cyc_N.log next to this script; relative paths only)."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from rc_core import LS, RL, ONE, U, BIG, apply_op, twisted

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
PREC = N + 14
XS = [Fr(1, 3), Fr(-2, 7), Fr(5, 4), Fr(11, 13), Fr(-7, 3), Fr(3, 5)]
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
LOGREL = 'logs/rc_cyc_%d.log' % N
OUT, NFAIL, T0 = [], [0], time.time()


def say(s):
    print(s)
    OUT.append(s)


def check(name, ok, control=False):
    good = (not ok) if control else ok
    say(('PASS ' if good else 'FAIL ') + name + ('  [control: must fail; it fails]' if control and good else
                                                 '  [control: must fail; it HOLDS]' if control else ''))
    NFAIL[0] += (not good)


def lsp(d):                      # exact polynomial in u
    return LS(d, BIG)


def compose(P, Q):               # operators as lists of left coefficients (LS): sum P[i] D^i
    out = {}
    for i, a in enumerate(P):
        for j, b in enumerate(Q):
            bm, binom = b, 1
            for m in range(i + 1):
                k = i - m + j
                out[k] = out.get(k, LS({}, BIG)) + a * bm * binom
                bm, binom = bm.D(), binom * (i - m) // (m + 1)
    return [out.get(k, LS({}, BIG)) for k in range(max(out) + 1)]


def op_eq(P, Q):
    n = max(len(P), len(Q))
    z = LS({}, BIG)
    return all(((P[k] if k < len(P) else z) - (Q[k] if k < len(Q) else z)).c == {} for k in range(n))


q = lsp({0: Fr(1, 3125), 1: Fr(-1, 3125)})
A = [lsp({0: Fr(-24, 625), 1: Fr(24, 625)}), lsp({0: Fr(-2, 5), 1: Fr(2, 5)}), lsp({0: Fr(-7, 5), 1: Fr(7, 5)}),
     lsp({0: -2, 1: 2}), lsp({1: 1})]
say('rc_cyc N = %d (critic of the Stage-5 referees; exact over Q; x-values %s)' % (N, ', '.join(str(x) for x in XS)))
# ---- A: operator algebra
prod4 = [lsp({0: 1})]
for k in range(1, 5):
    prod4 = compose(prod4, [lsp({0: k}), lsp({0: 5})])
lr = compose([q * 5], prod4)
pf4_lr = [-c for c in lr] + [LS({}, BIG)] * (5 - len(lr))
pf4_lr[4] = pf4_lr[4] + lsp({0: 1})
check('A1 D^4 - 5q prod_(k=1..4)(5D+k) has the coefficients A_k = (u, -2(1-u), -7/5(1-u), -2/5(1-u), -24/625(1-u))',
      op_eq(pf4_lr, A))
prod5 = compose(prod4, [lsp({0: 5}), lsp({0: 5})])
lhs5 = [-c for c in compose([q], prod5)] + [LS({}, BIG)]
lhs5[5] = lhs5[5] + lsp({0: 1})
check('A2 D^5 - q prod_(k=1..5)(5D+k) = D o PF4 (so PF5 = z^5 D PF4 - 1)', op_eq(lhs5, compose([LS({}, BIG), ONE], A)))
ctl5 = [-c for c in compose(prod5, [q])]
ctl5[5] = ctl5[5] + ONE
check('A2c control: with q on the right of the product', op_eq(ctl5, compose([LS({}, BIG), ONE], A)), control=True)


def adjoint(P):
    out = [LS({}, BIG)] * len(P)
    for k, a in enumerate(P):
        Dk = [LS({}, BIG)] * k + [lsp({0: (-1) ** k})]
        t = compose(Dk, [a])
        for i in range(len(t)):
            out[i] = out[i] + t[i]
    return out


check('A3 C1: PF4* = PF4 (formal adjoint, D* = -D)', op_eq(adjoint(A), A))
check('A3c control: PF4 + u D is self-adjoint', op_eq(adjoint([A[0], A[1] + U] + A[2:]), [A[0], A[1] + U] + A[2:]), True)
# ---- B: Frobenius basis, obstruction matrix
Pj = {}
for n in range(PREC + 4):
    img = apply_op(A, RL.of(LS.mono(n))).a[0]
    for j in range(-3, 2):
        Pj[(n, j)] = img.coef(n + j)
check('B1 PF4 maps Q[[u]] to itself: [u^m] PF4 u^n = 0 for m < 0 (n = 0, 1, 2)',
      all(Pj[(n, j)] == 0 for n in range(3) for j in range(-3, 2) if n + j < 0))
check('B1 pivot [u^(n-3)] PF4 u^n = n (n-1)^2 (n-2) for n < %d' % (PREC + 4),
      all(Pj[(n, -3)] == n * (n - 1) ** 2 * (n - 2) for n in range(PREC + 4)))


def solve(start, rhs=None):      # PF4 f = rhs (rhs an LS), f = sum a_n u^n with a_0..a_2 = start
    a = list(start) + [Fr(0)] * PREC
    top = PREC if rhs is None else min(PREC, rhs.prec + 3)
    for m in range(top - 3):
        s = sum(a[n] * Pj[(n, m - n)] for n in range(max(0, m - 1), m + 3))
        r = rhs.coef(m) if rhs is not None else 0
        a[m + 3] = (r - s) / Pj[(m + 3, -3)]
    return LS({k: a[k] for k in range(top)}, top)


f = [solve([1, 0, 0]), solve([0, 1, 0]), solve([0, 0, 1])]
check('B2 f_0, f_1, f_2 solve PF4 below u^%d; first terms 1 + (2/625)u^3, u - (19/60)u^3, u^2 + (37/30)u^3' % (PREC - 4),
      all(apply_op(A, RL.of(x)).a[0].zero_below(PREC - 4) for x in f) and
      [f[0].coef(3), f[1].coef(3), f[2].coef(3)] == [Fr(2, 625), Fr(-19, 60), Fr(37, 30)])


def r_of(x):                      # r(x) = PF4(x l) - l PF4(x), the l^0 part for a solution x
    return apply_op(A, RL([LS({}, BIG), x])).a[0]


Ob = [[r_of(f[i]).coef(m) for i in range(3)] for m in (-3, -2, -1)]
say('   obstruction matrix (rows u^-3, u^-2, u^-1; columns f_0, f_1, f_2): %s' % [[str(x) for x in row] for row in Ob])
check('B3 obstruction matrix = [[-2,0,0],[17/5,0,0],[-7/5,-7/5,2]] (cyc_formal E1): kernel Q(0, 1, 7/10)',
      Ob == [[-2, 0, 0], [Fr(17, 5), 0, 0], [Fr(-7, 5), Fr(-7, 5), 2]])
gamma = -Ob[2][1] / Ob[2][2]
g = f[1] + f[2] * gamma
check('B4 gamma = 7/10, g = u + (7/10)u^2 + (41/75)u^3 + ...', gamma == Fr(7, 10) and
      [g.coef(k) for k in range(4)] == [0, 1, Fr(7, 10), Fr(41, 75)])
h = solve([0, 0, 0], -r_of(g))
y3 = RL([h, g])
res3 = apply_op(A, y3)
check('B5 y_3 = g l + h solves PF4 (l^0 and l^1 parts) below u^%d' % (PREC - 6),
      all(x.zero_below(PREC - 6) for x in res3.a))
check('B5c control: g + f_0 (in place of g) is unobstructed', all(r_of(g + f[0]).coef(m) == 0 for m in (-3, -2, -1)), True)
check('B5c control: gamma = 1/2 is unobstructed', all(r_of(f[1] + f[2] * Fr(1, 2)).coef(m) == 0 for m in (-3, -2, -1)), True)


# ---- C: the concomitant from its definition
def conc(x, y):                   # P(x, y) = sum_k sum_j (-1)^j D^j(A_k y) D^(k-1-j) x
    tot = RL([LS({}, BIG)])
    for k in range(1, 5):
        for j in range(k):
            t = y.scale(A[k])
            for _ in range(j):
                t = t.D()
            s = x
            for _ in range(k - 1 - j):
                s = s.D()
            tot = tot + (t * s).scale(Fr((-1) ** j))
    return tot


basis = [RL.of(f[0]), RL.of(f[1]), RL.of(f[2]), y3]
Bm, okc = [[None] * 4 for _ in range(4)], True
for i in range(4):
    for j in range(4):
        b = conc(basis[i], basis[j])
        cut = PREC - 10
        okc &= b.logfree(cut) and all(b.a[0].coef(k) == 0 for k in range(b.a[0].val(), cut) if k != 0)
        Bm[i][j] = b.a[0].coef(0)
say('   Bmat (f_0, f_1, f_2, y_3) = %s' % [[str(x) for x in row] for row in Bm])
check('C1 B(y_i, y_j) is log-free and constant below u^%d, alternating' % (PREC - 10),
      okc and all(Bm[i][j] == -Bm[j][i] for i in range(4) for j in range(4)))
check('C2 hand values (cyc_formal E2): b01 = -7/5, b02 = 2, b03 = -9/10, b12 = 0, b13 = -1, b23 = 0',
      [Bm[0][1], Bm[0][2], Bm[0][3], Bm[1][2], Bm[1][3], Bm[2][3]] == [Fr(-7, 5), 2, Fr(-9, 10), 0, -1, 0])
pf = Bm[0][1] * Bm[2][3] - Bm[0][2] * Bm[1][3] + Bm[0][3] * Bm[1][2]
check('C3 Pf(Bmat) = -B(f_0,f_2) B(f_1,y_3) = 2, det Bmat = 4', pf == 2 and pf * pf == 4)
check('C4 u^-1 obstruction row = B(f_0, .) on (f_1, f_2): (-7/5, 2) (cyc_adversary s.3)', [Ob[2][1], Ob[2][2]] == [Bm[0][1], Bm[0][2]])
check('C5 C4 (isotropy): B(f_i, g) = 0 for i = 0, 1, 2 and B(y_3, g) = 1',
      all(Bm[i][1] + gamma * Bm[i][2] == 0 for i in range(3)) and Bm[3][1] + gamma * Bm[3][2] == 1)
W = [[basis[i] for i in range(4)]]
for a in range(1, 4):
    W.append([w.D() for w in W[-1]])


def det4(Mx):
    from itertools import permutations
    tot = RL([LS({}, BIG)])
    for p in permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if p[i] > p[j]:
                    sg = -sg
        t = RL.of(ONE).scale(Fr(sg))
        for a in range(4):
            t = t * Mx[a][p[a]]
        tot = tot + t
    return tot


dw = det4(W)
cut = PREC - 12
check('C6 det Wr is log-free and equals -2/u^2 below u^%d (Abel: D(u^2 det Wr) = 0; cyc_formal E2)' % cut,
      dw.logfree(cut) and dw.a[0].coef(-2) == -2 and all(dw.a[0].coef(k) == 0 for k in range(dw.a[0].val(), cut) if k != -2))
exec(open(os.path.join(HERE, 'rc_cyc_frames.py'), encoding='utf-8').read())
