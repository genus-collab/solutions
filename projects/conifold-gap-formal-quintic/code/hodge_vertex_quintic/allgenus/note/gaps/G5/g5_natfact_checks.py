# g5_natfact_checks.py -- cluster G5, item 3 (cor:mai-natfact): evidence only, exact over Q(zeta), zeta^5 = 1.
# (1) Psi_0 = Rhat^mod(0) (flat index b, colour i; entries eps_b zeta^(i(b+1))) satisfies Psi_0^T eta^-1 Psi_0 = 1 (orthonormality,
#     eta^-1_ab = (1/5) delta_(a+b = 3 mod 5)); the relation "Psi eta^-1 Psi^T = eta^-1" printed in OU3 s.3 FAILS for this Psi_0 (erratum).
# (2) Row 0 of Psi_0 is the unit 1 = sum_j zeta^(-4j) e~_j (hatted, at u = 0).
# (3) The exponents mu(b,k) = (4k + nu_b - nu_(b+k))/5 of the hatted gauge: integers >= max(0, (4k-4)/5); mu = 0 only for k = 0 and
#     (b,k) = (1,1)  (the entry (2,1), killed at u = 0 by (R1)).
# (4) The model rows from the T-system (prop:set-modrows): Fs_b(0) = eps_b, Fs_0 = 1 - (3/20) t + (9/800) t^2 + ...; M_0 = Psi_0^-1 Rhat^mod
#     is symplectic, M_0(-t)^T M_0(t) = 1, through t^K (input of cor:mai-natfact (c); def:set-M0).  Control: a perturbed Fs_0 breaks it.
# Usage: py -I -B g5_natfact_checks.py [K]
import sys
from fractions import Fraction as Fr

def zred(v):                          # v: list of 5 coefficients of 1, z, .., z^4 -> reduce to z^4 coefficient 0
    c4 = v[4]
    return [v[0] - c4, v[1] - c4, v[2] - c4, v[3] - c4, Fr(0)]
def zmul(a, b):
    r = [Fr(0)] * 5
    for i in range(5):
        if a[i]:
            for j in range(5):
                if b[j]:
                    r[(i + j) % 5] += a[i] * b[j]
    return zred(r)
def zadd(a, b, c=1): return [x + c * y for x, y in zip(a, b)]
def zc(q): return [Fr(q), Fr(0), Fr(0), Fr(0), Fr(0)]
def zpow(k):
    v = [Fr(0)] * 5; v[k % 5] = Fr(1); return zred(v)
ZERO = zc(0); ONE = zc(1)
NU = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)

def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum_z(zmul(A[i][k], B[k][j]) for k in range(m)) for j in range(p)] for i in range(n)]
def sum_z(it):
    r = list(ZERO)
    for x in it:
        r = zadd(r, x)
    return r
def transpose(A): return [list(r) for r in zip(*A)]
ETAINV = [[zc(Fr(1, 5)) if (a + b) % 5 == 3 else list(ZERO) for b in range(5)] for a in range(5)]
ID = [[list(ONE) if a == b else list(ZERO) for b in range(5)] for a in range(5)]

def model_rows(K, perturb=None):      # colour 0, T_k = -1 + (t/5)(k - 4 t d/dt); returns Fs_b as coefficient lists (over Q)
    def T(k, f):
        r = [-x for x in f] + [Fr(0)]
        for n, x in enumerate(f):
            r[n + 1] += Fr(k - 4 * n, 5) * x
        return r
    def chain(f0):
        f1 = T(-1, f0); f2 = [-x for x in T(-2, f1)]; f3 = T(2, f2); f4 = [-x for x in T(1, f3)]; back = [-x for x in T(0, f4)]
        return [f0, f1, f2, f3, f4], back
    f0 = [Fr(1)]
    for n in range(1, K + 1):
        trial = f0 + [Fr(0)]
        _, back = chain(trial)
        res = back[n + 1] - (trial[n + 1] if n + 1 < len(trial) else 0)   # [t^(n+1)] of (O - 1) Fs_0 with f_n = 0
        trial[n] = Fr(1)
        _, back1 = chain(trial)
        res1 = back1[n + 1] - (trial[n + 1] if n + 1 < len(trial) else 0)
        piv = res1 - res
        f0 = trial[:n] + [-res / piv]
    if perturb:
        f0 = list(f0); f0[perturb[0]] += perturb[1]
    rows, back = chain(f0)
    closing = all(back[n] == (f0[n] if n < len(f0) else 0) for n in range(K + 1))
    return [r[:K + 1] for r in rows], closing, piv

def rhat(rows, K):                    # Rhat^mod_b^(i)(t) = zeta^(i(b+1)) Fs_b(zeta^-i t): list over t-orders of 5x5 matrices (b, i)
    out = []
    for n in range(K + 1):
        out.append([[zmul(zpow(i * (b + 1) - i * n), zc(rows[b][n] if n < len(rows[b]) else 0)) for i in range(5)] for b in range(5)])
    return out

def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print('g5_natfact_checks: K = %d, exact over Q(zeta)' % K)
    Psi0 = [[zmul(zc(EPS[b]), zpow(i * (b + 1))) for i in range(5)] for b in range(5)]
    G = mat_mul(mat_mul(transpose(Psi0), ETAINV), Psi0)
    print('(1) Psi_0^T eta^-1 Psi_0 = 1:', G == ID)
    G2 = mat_mul(mat_mul(Psi0, ETAINV), transpose(Psi0))
    print('(1) control/erratum: Psi_0 eta^-1 Psi_0^T = eta^-1 (OU3 s.3 wording):', G2 == ETAINV, '(False expected)')
    print('(2) row 0 of Psi_0 equals (zeta^(-4j))_j:', all(Psi0[0][j] == zpow(-4 * j) for j in range(5)))
    ok3 = True; zero_mu = []
    for b in range(5):
        for k in range(0, 41):
            num = 4 * k + NU[b] - NU[(b + k) % 5]
            if num % 5:
                ok3 = False; continue
            mu = num // 5
            if mu < 0 or 5 * mu < 4 * k - 4:
                ok3 = False
            if mu == 0:
                zero_mu.append((b, k))
    print('(3) mu(b,k) integral, >= max(0,(4k-4)/5) for k <= 40:', ok3, '; mu = 0 exactly at (b,k) =', zero_mu)
    rows, closing, piv = model_rows(K)
    print('(4) model rows: closing relation through t^%d:' % K, closing, '; Fs_b(0) =', [str(r[0]) for r in rows],
          '; Fs_0 =', [str(x) for x in rows[0][:4]], '...')
    R = rhat(rows, K)
    print('(4) Rhat^mod(0) = Psi_0:', R[0] == Psi0)
    PsiInv = mat_mul(transpose(Psi0), ETAINV)
    M = [mat_mul(PsiInv, R[n]) for n in range(K + 1)]
    def symp(M):
        ok = True
        for n in range(K + 1):
            acc = [[list(ZERO) for _ in range(5)] for _ in range(5)]
            for a in range(n + 1):
                Ma = transpose(M[a]); sgn = (-1) ** a
                P = mat_mul(Ma, M[n - a])
                acc = [[zadd(acc[i][j], P[i][j], sgn) for j in range(5)] for i in range(5)]
            ok &= (acc == (ID if n == 0 else [[list(ZERO) for _ in range(5)] for _ in range(5)]))
        return ok
    print('(4) M_0(-t)^T M_0(t) = 1 through t^%d:' % K, symp(M))
    rowsX, _, _ = model_rows(K, perturb=(3, Fr(1, 7)))
    RX = rhat(rowsX, K); MX = [mat_mul(PsiInv, RX[n]) for n in range(K + 1)]
    print('(4) control: Fs_0 perturbed at t^3 (closing relation then fails) -> symplectic:', symp(MX), '(False expected)')

if __name__ == '__main__':
    main()
