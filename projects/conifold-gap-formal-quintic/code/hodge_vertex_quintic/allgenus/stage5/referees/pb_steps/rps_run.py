"""rps_run.py -- referee pb_steps: driver.  Usage: py -I -B rps_run.py E NH NPI M
  E   : original u-order of V (exact over Q)          NH : hatted orders 0..NH of the direct hatted residue H
  NPI : u-order of the reduction coefficients pi_m     M  : t-order of the (S2) / Prop. S checks
Tests (all exact over Q):
  T1  relabelling V = s Res_X[...]:  [u^N w^j] H  =  [u^(N+(1-4j)/5) z^j] V  for every (N, j) in range, both directions;
  T2  hatted FL3 + reduction:  H = Sum_m pi_m Res_Y[Y^m e^(-f/w) dY]  at every hatted order <= NPI;
  T3  key lemma:  Sum_m pi_m s^(-m) I^C_m(v(u), z) = V  on admissible keys E' + (4b-1)/5 <= NPI; truncation layer;
  T4  (S2) and Prop. S: pivots, eps_b, SP rows satisfy the model system.
Every test has failing controls.  The log is written to logs/ next to this script (relative path printed)."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import rps_core as rc
import rps_mech as rm
import rps_s2 as rs

E, NH, NPI, M = (int(a) for a in sys.argv[1:5])
NH = max(NH, NPI + 1)          # T2 compares at hatted orders <= NPI + 1
HERE = os.path.dirname(os.path.abspath(__file__))
LOGREL = 'logs/rps_run_v3_%d_%d_%d_%d.log' % (E, NH, NPI, M)   # v3: T2 degree diagnostic as weight checks (v2 had a wrong prediction)
lines, fails = [], []
def log(s=''):
    lines.append(s); print(s, flush=True)
def check(name, ok, info=''):
    log('  %-74s %s %s' % (name, 'PASS' if ok else 'FAIL', info))
    if not ok:
        fails.append(name)
def control(name, failed, info=''):
    log('  control %-66s %s %s' % (name, 'PASS (fails as it must)' if failed else 'FAIL (did not fail)', info))
    if not failed:
        fails.append('control ' + name)

def compare_HV(H, V, Nmax, Elim):
    """returns (n_compared, n_nonzero, first_mismatch, offgrade_ok, touched)."""
    ncmp = nnz = 0; first = None; offok = True; touched = {}
    for N in range(Nmax + 1):
        for j in range(rc.jmin_for(N, Elim), 2 * N + 6):
            h = H[N].get(j, 0)
            if (1 - 4 * j) % 5:
                offok = offok and h == 0
                continue
            Ep = N + (1 - 4 * j) // 5
            if Ep > Elim:
                continue
            v = V.get((Ep, j), 0)
            ncmp += 1; nnz += (v != 0)
            if v:
                touched.setdefault(N, set()).add(Ep)
            if h != v and first is None:
                first = (N, j)
    return ncmp, nnz, first, offok, touched

t0 = time.time()
log('rps_run.py %d %d %d %d   (exact over Q; referee pb_steps)' % (E, NH, NPI, M))
# ------------------------------------------------------------------ T1
log('T1  relabelling: original-variable V versus the hatted residue H computed directly')
t = time.time()
V = rc.orig_V(E)
okV = all(b % 5 == 4 and b <= -1 and n >= 1 for (n, b) in V)
check('V: support in u^>=1, z^beta with beta = 4 mod 5 and beta <= -1', okV, '(%d terms through u^%d)' % (len(V), E))
check('V at u^1 is -1/z', V.get((1, -1)) == -1 and all(b == -1 for (n, b) in V if n == 1))
H, P = rc.hat_H(NH, E)
l0 = [H[0].get(-5 * l - 1) == Fr((-1) ** (5 * l + 1), rc.fact(4 * l + 1) * 4 ** l * rc.fact(l)) for l in range(4)]
check('H at hatted u^0: Sum_l (-1)^(5l+1) w^-(5l+1)/((4l+1)! 4^l l!)  (l < 4)', all(l0))
ncmp, nnz, first, offok, touched = compare_HV(H, V, NH, E)
check('H = s^-1 V under u^N w^j <-> u^(N+(1-4j)/5) z^j, hatted orders 0..%d' % NH, first is None,
      '(%d coefficients, %d non-zero; first mismatch %s)' % (ncmp, nnz, first))
check('H has no w-power outside 4 mod 5 in the window', offok)
pos = sum(1 for N in range(NH + 1) for j in H[N] if j > 0)
check('H has no positive w-power (V has no positive z-power)', pos == 0)
cov = sum(1 for (n, b) in V if n + (4 * b - 1) // 5 <= NH)
log('    V-terms covered by hatted orders <= %d: %d of %d (each matched exactly once: the map is injective)' % (NH, cov, len(V)))
for N in sorted(touched)[:4]:
    log('    hatted order u^%d collects original u-orders %s' % (N, sorted(touched[N])[:8] + (['...'] if len(touched[N]) > 8 else [])))
log('    (%.1f s)' % (time.time() - t))
NHc = min(NH, 5)
for name, kw in (('Lam/u -> 1 in the hatted phase (no regrouping factor)', dict(noregroup=True)),
                 ('hatted amplitude without (1-uX^5)^(-1/2)', dict(nosqrt=True)),
                 ('hatted amplitude with b_k -> -b_k', dict(flip=True)),
                 ('hatted amplitude with an extra u per Stirling term', dict(wrongw=True))):
    Hc, _ = rc.hat_H(NHc, E, **kw)
    r = compare_HV(Hc, V, NHc, E)
    control(name, r[2] is not None, '(first mismatch %s)' % (r[2],))
Vc = rc.orig_V(E, gfirst=True)
r = compare_HV(H, Vc, NHc, E)
control('original G truncated to x^4/4', r[2] is not None, '(first mismatch %s)' % (r[2],))
# ------------------------------------------------------------------ T2
log('T2  CFU, reduction, and the hatted identity H = Sum_m pi_m Res_Y[Y^m e^(-f/w) dY]')
t = time.time()
NC = NPI + 1
y, Vc_ = rm.cfu(NC)
tinv = [Fr(4, 4 + 5 * n) for n in range(NC + 1)]
lhs = [Fr(1)] + [Fr(0)] * NC
for _ in range(4):
    lhs = rc.mul_u(lhs, Vc_, NC)
rhs = [Fr(1)] + [Fr(0)] * NC
for _ in range(5):
    rhs = rc.mul_u(rhs, tinv, NC)
check('CFU: V^4 = tau^-5 (V = tau^(-5/4), 1/tau = 2F1(1,4/5;9/5;u)) through u^%d' % NC, lhs == rhs)
check('CFU: y_1 = (1 + T)/18, y_k polynomial', [c for c in y[1] if True][:2] == [Fr(1, 18), Fr(1, 18)])
x = rm.cfu_inverse(y, NC)
Tx5 = rm.smul([[Fr(0), Fr(1)]], rm.spow(x, 5, NC), NC)
idt = rm.smul(x, rm.scompose(y, Tx5, NC), NC)
check('inverse: x(T) y(T x(T)^5) = 1 through u^%d' % NC,
      all(all(c == 0 for c in idt[i][(1 if i == 0 else 0):]) for i in range(NC + 1)) and idt[0][0] == 1)
B = rc.hat_B(NC)
beta = rm.beta_Y(B, x, NC)
pi = rm.reduce_forms(beta, Vc_, NC)
grade = all((j - m) % 5 == 0 for m in pi for (i, j) in pi[m])
check('grading: pi_m contains only w^j with j = m mod 5', grade)
degok = True; sharp = 0; apri = 0
for m in pi:
    for (i, j) in pi[m]:
        if 4 * j > 5 * i - m:
            degok = False
        if 4 * j > 5 * i - m - 4:
            sharp += 1
for (i, n, j) in beta:
    apri = max(apri, j - i)
check('degree bound deg_w [u^i] pi_m <= (5i - m)/4 (regularity of P^SP), i <= %d' % NC, degok,
      '(%d coefficients in the top band)' % sharp)
log('    beta(Y) has %d terms; max (w-degree - u-order) in beta: %d' % (len(beta), apri))
apr = [max((j + (n + 1) // 4 for (i0, n, j) in beta if i0 <= i), default=0) for i in range(NC + 1)]
act = [max((j for m in pi for (i0, j) in pi[m] if i0 == i), default=None) for i in range(NC + 1)]
bnd = [(5 * i + 1) // 4 for i in range(NC + 1)]
log('    w-degree of pi at u^i, i = 0..%d: actual %s' % (NC, act))
log('      bound floor((5i+1)/4) from regularity %s;  a-priori bound of the reduction %s' % (bnd, apr))
check('weight: every beta monomial u^i Y^n w^j has n + 4j <= 5i (an original-variable object)',
      all(n + 4 * j <= 5 * i for (i, n, j) in beta))
check('hence the reduction alone gives deg_w [u^i] pi <= floor((5i+1)/4) (a-priori = bound)', apr == bnd)
check('weight of pi: every u^i w^j in pi_m has m + 4j <= 5i', all(m + 4 * j <= 5 * i for m in pi for (i, j) in pi[m]))
emax = 4 * (2 * NC + 6 - rc.jmin_for(0, E)) // 5 + 8
Vp = rm.vpowers(Vc_, emax, NC)
first = None; ncmp = 0
for N in range(NC + 1):
    jm = rc.jmin_for(N, E)
    R = rm.resid_hat(pi, Vp, N, jm, 2 * N + 5)
    for j in range(jm, 2 * N + 6):
        ncmp += 1
        if R.get(j, 0) != H[N].get(j, 0) and first is None:
            first = (N, j)
check('H = Sum_m pi_m R_m at hatted orders 0..%d (all w-powers in the window)' % NC, first is None,
      '(%d coefficients; first mismatch %s)' % (ncmp, first))
log('    (%.1f s)' % (time.time() - t))
NCc = min(NC, 5)
xc = [[Fr(1)]] + [[] for _ in range(NCc)]
for name, pic in (('pull-back without X\'(Y)', rm.reduce_forms(rm.beta_Y(rc.hat_B(NCc), x[:NCc + 1], NCc, noJac=True), Vc_, NCc)),
                  ('reduction with V -> 1', rm.reduce_forms(rm.beta_Y(rc.hat_B(NCc), x[:NCc + 1], NCc), [Fr(1)] + [Fr(0)] * NCc, NCc)),
                  ('no CFU (X = Y)', rm.reduce_forms(rm.beta_Y(rc.hat_B(NCc), xc, NCc), Vc_, NCc))):
    first = None
    for N in range(NCc + 1):
        jm = rc.jmin_for(N, E)
        R = rm.resid_hat(pic, Vp, N, jm, 2 * N + 5)
        for j in range(jm, 2 * N + 6):
            if R.get(j, 0) != H[N].get(j, 0) and first is None:
                first = (N, j)
    control(name, first is not None, '(first mismatch %s)' % (first,))
exec(open(os.path.join(HERE, 'rps_run_b.py')).read())
