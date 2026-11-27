# ref5_item4.py -- referee of gaps/G5, item 4 (prop:con-deg step 5).  Evidence only; exact over Q (fractions).
# Independent of g5_degree_bound.py: the operator P is built as a NORMAL-ORDERED differential operator
#   sum c z^j x^a theta^b   (theta = x d/dx),  composition  theta^b x^c = x^c (theta + c)^b,
# and the identities (L0)-(L3) of GAPS_G5_3 s.4.2 are checked as OPERATOR identities (not on monomials).  New: (L4) every term
# z^j x^a theta^b of P has a <= 4j (so P_j raises the x-degree by at most 4j); this is the degree bound that the existence of
# A_n in x^(4n+1) Q[[x^-5]] needs and that GAPS_G5_3 s.4.3 leaves implicit.
# Recursion: E A_n = rho_n (flat), checked against an independent hatted HPF built from A_e (my own chain rule, no POLY coefficients).
# Controls (must fail): P + z^2 x^3 (keeps (L0),(L2),(L4), the grading; breaks (L1)); P + z^2 (1 - x^5) x^-2 (keeps (L1); breaks (L2)).
# Usage: py -I -B ref5_item4.py [N]
import sys
from fractions import Fraction as Fr
from math import comb

def oadd(A, B, c=1):
    R = dict(A)
    for k, v in B.items():
        R[k] = R.get(k, 0) + c * v
    return {k: v for k, v in R.items() if v != 0}

def omul(A, B):
    R = {}
    for (j, a, b), u in A.items():
        for (k, c, d), w in B.items():
            for i in range(b + 1):                      # theta^b x^c = x^c sum_i C(b,i) c^(b-i) theta^i
                t = comb(b, i) * c ** (b - i)
                if t:
                    key = (j + k, a + c, i + d)
                    R[key] = R.get(key, 0) + u * w * t
    return {k: v for k, v in R.items() if v != 0}

def build_P(KS=(1, 2, 3, 4, 5)):
    M = {(0, 1, 0): Fr(1), (1, 0, 1): Fr(-1, 5), (1, 5, 1): Fr(1, 5)}
    M5 = {(0, 0, 0): Fr(1)}
    for _ in range(5):
        M5 = omul(M5, M)
    Pi = {(0, 0, 0): Fr(1, 3125)}
    for k in KS:
        Pi = omul(Pi, {(0, 1, 0): Fr(5), (1, 0, 1): Fr(-1), (1, 5, 1): Fr(1), (1, 0, 0): Fr(k)})
    P = oadd(M5, {(0, 0, 0): Fr(1)}, -1)
    return oadd(P, omul({(0, -5, 0): Fr(1), (0, 0, 0): Fr(-1)}, Pi))

def part(P, j):
    return {(a, b): c for (jj, a, b), c in P.items() if jj == j}

def apply(Pj, f):                                       # Pj: {(a,b): c}; f: {m: c}
    R = {}
    for (a, b), c in Pj.items():
        for m, v in f.items():
            t = c * v * m ** b
            if t:
                R[a + m] = R.get(a + m, 0) + t
    return {k: v for k, v in R.items() if v != 0}

def div_1mx5(f):                                        # exact division of a Laurent polynomial by (1 - x^5); None if not exact
    g = dict(f); q = {}
    while g:
        m = min(g); c = g.pop(m); q[m] = c
        g[m + 5] = g.get(m + 5, 0) + c
        if g[m + 5] == 0:
            del g[m + 5]
        if g and min(g) > max(f) + 5:
            return None
        if m > max(f):
            return None
    return q

def div_1mu(f):                                         # exact division of a polynomial in u by (1 - u); None if not exact
    if not f:
        return {}
    top = max(f); q = {}; acc = 0
    for k in range(0, top + 1):
        acc += f.get(k, 0); q[k] = acc                  # f/(1-u) = sum_k (sum_{i<=k} f_i) u^k
    if acc != 0:
        return None
    return {k: v for k, v in q.items() if v != 0}

def operator_checks(P, label):
    out = {}
    out['L0'] = not part(P, 0)
    L1 = True
    for j in range(1, 6):
        for b in {bb for (_, bb) in part(P, j)}:
            col = {a: c for (a, bb), c in part(P, j).items() if bb == b}
            L1 &= all(sum(c for a, c in col.items() if a % 5 == r) == 0 for r in range(5))
    out['L1'] = L1
    want = {(-1, 0): Fr(1), (-1, 1): Fr(-1), (4, 0): Fr(-1), (4, 1): Fr(1)}   # (1-x^5) x^-1 (1 - theta)
    out['L3'] = part(P, 1) == want
    out['L4'] = all(a <= 4 * j for (j, a, b) in P)
    out['grade'] = all((a + j) % 5 == 0 for (j, a, b) in P)
    out['L2'] = all(min(apply(part(P, j), {m: Fr(1)}), default=99) >= 1 for j in range(1, 6) for m in range(1, 81))
    out['L2_fails_at_m0'] = any(min(apply(part(P, j), {0: Fr(1)}), default=99) < 1 for j in range(1, 6))
    out['maxj'] = max(j for (j, a, b) in P)
    out['top'] = {j: max(a for (jj, a, b) in P if jj == j) for j in range(1, 6)}
    print('[%s] (L0) P_0 = 0: %s; (L1) left factor 1-x^5 in every z^j, theta^b coefficient: %s; (L3) P_1 = (1-x^5) x^-1 (1-theta): %s'
          % (label, out['L0'], out['L1'], out['L3']))
    print('[%s] (L4) a <= 4j for every term: %s (top a per j: %s); grade a + j = 0 mod 5: %s; (L2) P_j(x^m) in xQ[x], 1<=m<=80: %s;'
          ' P_j(1) leaves xQ[x] (expected): %s' % (label, out['L4'], out['top'], out['grade'], out['L2'], out['L2_fails_at_m0']))
    return out

def flat_recursion(P, N, free):
    Pj = {j: part(P, j) for j in range(1, 6)}
    A = [{1: Fr(1)}]; ok = True; rhos = []
    for n in range(1, N + 1):
        rhs = {}
        for j in range(2, 6):
            if n + 1 - j >= 0:
                rhs = oadd(rhs, apply(Pj[j], A[n + 1 - j]), -1)
        rho = div_1mx5(rhs) if rhs else {}
        if rho is None:
            print('  n = %d: RHS_n not divisible by 1 - x^5' % n); return A, rhos, False
        rhos.append(rho)
        good = (not rho or (min(rho) >= 1 and max(rho) <= 4 * n)) and all((m - 4 * n) % 5 == 0 for m in rho)
        An = {m + 1: c / (1 - (m + 1)) for m, c in rho.items() if m + 1 != 1}
        if rho.get(0, 0) != 0:
            good = False
        if n % 5 == 0 and free.get(n):
            An[1] = free[n]
        tot = {}
        A.append(An)
        for j in range(1, 6):
            if n + 1 - j >= 0:
                tot = oadd(tot, apply(Pj[j], A[n + 1 - j]))
        good &= (not tot) and min(An, default=1) >= 1 and max(An, default=0) <= 4 * n + 1
        ok &= good
    return A, rhos, ok

def Ae(F, e):                                           # F: {(k,n): c} for u^k zhat^n;  A_e = 1 + zhat(1-u)((4/5)th_zhat - e/5 - th_u)
    R = dict(F)
    for (k, n), c in F.items():
        w = c * (Fr(4 * n, 5) - Fr(e, 5) - k)
        if w:
            R[(k, n + 1)] = R.get((k, n + 1), 0) + w
            R[(k + 1, n + 1)] = R.get((k + 1, n + 1), 0) - w
    return {k: v for k, v in R.items() if v != 0}

def hpf_full(F):                                        # P(s^-1 Phi) = s^-6 hpf_full(Phi)  (my chain rule)
    T = F
    for e in (-1, -2, -3, -4, -5):
        T = Ae(T, e)
    Q = F
    for e, k in ((-1, 1), (-2, 2), (-3, 3), (-4, 4), (-5, 5)):
        Q5 = {key: 5 * v for key, v in Ae(Q, e).items()}
        Q = oadd(Q5, {(kk + 1, n + 1): k * c for (kk, n), c in Q.items()})
    Q = oadd({kn: c / 3125 for kn, c in Q.items()}, {(kk + 1, n): c / 3125 for (kk, n), c in Q.items()}, -1)   # (1-u)/3125
    R = oadd(T, {(k + 1, n): c for (k, n), c in F.items()}, -1)
    return oadd(R, Q, -1)

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    print('ref5_item4: N = %d, exact over Q, colour 0' % N)
    P = build_P(); operator_checks(P, 'P')
    free = {5 * m: Fr((-1) ** m * (3 * m + 1), 11) for m in range(1, N // 5 + 1)}
    for lab, fr in (('free = 0', {}), ('free = (-1)^m (3m+1)/11', free)):
        A, rhos, ok = flat_recursion(P, N, fr)
        print('[%s] n <= %d: RHS_n = (1-x^5) rho_n, rho_n in x Q[x] of degree <= 4n, [x^0] rho_n = 0, A_n in x Q[x] of degree <= 4n+1,'
              ' [z^(n+1)] P psi = 0: %s' % (lab, N, ok))
        Phi = {(k, n): A[n][4 * n + 1 - 5 * k] for n in range(len(A)) for k in range(0, n + 1) if A[n].get(4 * n + 1 - 5 * k)}
        H = hpf_full(Phi)
        print('[%s] independent hatted HPF (from A_e): [zhat^m] hpf_full(Phi) = 0 for m <= %d: %s' % (lab, N + 1,
              all(c == 0 for (k, n), c in H.items() if n <= N + 1)))
        dic = True
        for n in range(1, min(N, 30) + 1):
            Phib = {kn: c for kn, c in Phi.items() if kn[1] < n}
            Hn = {k: c for (k, m), c in hpf_full(Phib).items() if m == n + 1}
            rn = div_1mu({k - 1: c for k, c in Hn.items()}) if Hn.get(0, 0) == 0 else None   # / u, then / (1 - u), exact
            dic &= rn is not None and all(rn.get(k, 0) == rhos[n - 1].get(4 * n - 5 * k, 0) for k in range(0, n + 1)) \
                and all((4 * n - m) % 5 == 0 and 0 <= (4 * n - m) // 5 <= n for m in rhos[n - 1])
        print('[%s] dictionary rho_n = x^(4n) r_n(x^-5), r_n := [zhat^(n+1)] hpf_full(Phi|ph_n=0)/(u(1-u)), n <= %d: %s'
              % (lab, min(N, 30), dic))
        if lab == 'free = 0':
            print('  ph_5 =', [str(A[5].get(21 - 5 * k, 0)) for k in range(5)])
    print('[P] hatted solver, true operator, n <= 10: obstructions =', resonance_scan({}, 10), '(must be empty)')
    for lab, pert in (('P + z^2 x^3', {(2, 3, 0): Fr(1)}), ('P + z^2 (1-x^5) x^-2', {(2, -2, 0): Fr(1), (2, 3, 0): Fr(-1)})):
        Pc = oadd(P, pert); operator_checks(Pc, lab)
        _, _, okc = flat_recursion(Pc, min(N, 12), {})
        print('[%s] flat recursion with trivial resonance succeeds (must be False): %s' % (lab, okc))
        hatted = {(2, 2): Fr(1)} if lab == 'P + z^2 x^3' else {(2, 2): Fr(-1), (3, 2): Fr(1)}   # u^2 zh^2, resp. -u^2 (1-u) zh^2
        res = resonance_scan(hatted, 10)
        print('[%s] hatted solver: first obstruction (n, [u^(4n/5)] coefficient) =' % lab, res, '(must be non-empty)')

def resonance_scan(pert, NM):                           # solve hpf_full(Phi) + pert * Phi = 0 order by order in Q[[u]] (truncated)
    U = 4 * NM
    ph = [{0: Fr(1)}]
    for n in range(0, NM):
        Phi = {(k, m): c for m in range(len(ph)) for k, c in ph[m].items()}
        H = hpf_full(Phi)
        for (a, b), c in pert.items():
            for (k, m), v in Phi.items():
                H[(k + a, m + b)] = H.get((k + a, m + b), 0) + c * v
        Y = {k: c for (k, m), c in H.items() if m == n + 2 and c != 0 and k <= U}      # [zhat^(n+2)], unknown ph_(n+1)
        # [zhat^(n+2)] = u(1-u)(4(n+1) - 5 theta) ph_(n+1) + Y  ->  divide Y by u(1-u) as a power series
        if Y.get(0, 0):
            return [(n + 1, 'no u-divisibility')]
        y = {}; acc = 0
        for k in range(1, U + 1):
            acc += Y.get(k, 0); y[k - 1] = acc                 # Y/(u(1-u)) = sum_k (sum_{i<=k+1} Y_i) u^k
        new = {}
        for k in range(0, U - 1):
            piv = 4 * (n + 1) - 5 * k
            if piv == 0:
                if y.get(k, 0) != 0:
                    return [(n + 1, str(y[k]))]
                continue
            if y.get(k, 0):
                new[k] = -y[k] / piv
        ph.append(new)
    return []

if __name__ == '__main__':
    main()
