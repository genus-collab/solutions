# g5_degree_bound.py -- cluster G5, item 4 (prop:con-deg, step 5): evidence only, exact over Q.
# The flat WKB recursion of PF5 in x = L (colour 0; colour i is colour 0 at zeta^i x, since M_i(x) = M_0(zeta^i x)).
#   M = x - (z/5)(1 - x^5) theta_x,  Pi = (1/3125) prod_{k in KS} (5M + kz),  P = M^5 - 1 + (x^-5 - 1) Pi,  P = sum_j z^j P_j.
# Checks: (A) P_0 = 0; (B) P_j(x^m) divisible by (1 - x^5) for j >= 1; (C) P_j(x^m) in x Q[x] for m >= 1; (D) Pi(x^m) in x^6 Q[x][z];
#   (T) T_1 := P_1 = (1 - x^5) E, E x^m = (1 - m) x^(m-1); (R) the recursion T_1 psi_n = RHS_n: RHS_n = (1 - x^5) rho_n with
#   rho_n in x Q[x] (no x^0 term: the resonance is trivial) for EVERY choice of the free coefficients; psi_n in x Q[x], degree <= 4n+1;
#   (H) psi_n = x^(4n+1) ph_n(x^-5) with ph_n from the independent HPF recursion in u (POLY s.1) and the same free coefficients.
# Controls (must fail): (X1) a lower amplitude with a term of x-valuation <= 0 (psi_0 = x + x^-4, same Z/5 grade);
#   (X2) twist factors KS = (1,2,3,4,10) (leading exponent 2): property (C) fails at m = 2.
# Usage: py -I -B g5_degree_bound.py [N]      (prints relative information only)
import sys
from fractions import Fraction as Fr

def ladd(a, b, c=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + c * v
    return {k: v for k, v in r.items() if v != 0}

def zadd(F, G, c=1):
    n = max(len(F), len(G)); R = []
    for j in range(n):
        R.append(ladd(F[j] if j < len(F) else {}, G[j] if j < len(G) else {}, c))
    return R

def apply_M(F):                       # F: list over z-powers of Laurent polys in x
    R = [dict() for _ in range(len(F) + 1)]
    for j, f in enumerate(F):
        for m, c in f.items():
            R[j][m + 1] = R[j].get(m + 1, 0) + c
            if m != 0:
                R[j + 1][m] = R[j + 1].get(m, 0) - Fr(m, 5) * c
                R[j + 1][m + 5] = R[j + 1].get(m + 5, 0) + Fr(m, 5) * c
    return [{k: v for k, v in r.items() if v != 0} for r in R]

def apply_factor(k, F):               # 5M + k z
    MF = apply_M(F)
    zF = [dict()] + [dict(f) for f in F]
    R = zadd([{m: 5 * c for m, c in f.items()} for f in MF], [{m: k * c for m, c in f.items()} for f in zF])
    return R

def apply_Pi(F, KS):
    G = F
    for k in KS:
        G = apply_factor(k, G)
    return [{m: c / 3125 for m, c in g.items()} for g in G]

def apply_P(F, KS):
    G = F
    for _ in range(5):
        G = apply_M(G)
    G = zadd(G, F, -1)
    Pi = apply_Pi(F, KS)
    shifted = [ladd({m - 5: c for m, c in p.items()}, p, -1) for p in Pi]   # (x^-5 - 1) Pi
    return zadd(G, shifted)

def P_parts(m, KS):
    R = apply_P([{m: Fr(1)}], KS)
    return R + [dict() for _ in range(6 - len(R))]

def div_1mx5(f):                      # exact division of a Laurent polynomial by (1 - x^5); returns (q, ok)
    if not f:
        return {}, True
    g = dict(f); q = {}
    lo = min(g)
    while g:
        m = min(g)
        c = g.pop(m)
        q[m] = c
        g[m + 5] = g.get(m + 5, 0) + c
        if g[m + 5] == 0:
            del g[m + 5]
        if g and min(g) > max(f) + 5:
            return q, False
        if m > max(f) + 10:
            return q, False
    return q, True

def flat_recursion(N, KS, psi0, free, e_of_m):
    psi = [psi0]
    parts = {}
    def Pj(j, m):
        if m not in parts:
            parts[m] = P_parts(m, KS)
        return parts[m][j]
    out = []
    for n in range(1, N + 1):
        rhs = {}
        for j in range(2, 6):
            k = n + 1 - j
            if k < 0:
                continue
            for m, c in psi[k].items():
                rhs = ladd(rhs, Pj(j, m), -c)
        rho, ok = div_1mx5(rhs)
        info = {'n': n, 'div_ok': ok, 'rho_min': min(rho) if rho else None}
        new = {}
        res_ok = True
        for m1, c in rho.items():                     # E psi = rho: e(m) [x^m] psi = [x^(m-1)] rho
            m = m1 + 1
            e = e_of_m(m)
            if e == 0:
                if c != 0:
                    res_ok = False
                continue
            new[m] = c / e
        mk = [m for m in range(-60, 4 * n + 30) if e_of_m(m) == 0]
        for m in mk:
            if (m - (4 * n + 1)) % 5 == 0 and n in free:
                new[m] = free[n]
        info['res_ok'] = res_ok
        psi.append({k: v for k, v in new.items() if v != 0})
        # direct check of [z^(n+1)] P psi = 0
        tot = {}
        for j in range(1, 6):
            k = n + 1 - j
            if k < 0:
                continue
            for m, c in psi[k].items():
                tot = ladd(tot, Pj(j, m), c)
        info['eq_ok'] = (len(tot) == 0)
        p = psi[-1]
        info['psi_min'] = min(p) if p else None
        info['psi_max'] = max(p) if p else None
        info['grade_ok'] = all((m - (4 * n + 1)) % 5 == 0 for m in p)
        out.append(info)
    return psi, out

def hpf_polys(N, free):               # POLY s.1 recursion in u; returns ph_n as coefficient lists
    ph = [[Fr(1)]]
    def theta(p): return [k * c for k, c in enumerate(p)]
    def padd(a, b, c=1):
        n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + c * (b[i] if i < len(b) else 0) for i in range(n)]
    def pmul_1mu(p): return padd(p + [0], [0] + p, -1)
    def pmul_u(p, e): return [Fr(0)] * e + p
    def phis(phl, nmax):
        Phi = [[list(phl[m]) if m < len(phl) else [] for m in range(nmax + 1)]]
        for j in range(1, 6):
            prev = Phi[-1]; cur = []
            for m in range(nmax + 1):
                t = list(prev[m])
                if m >= 1:
                    q = padd([Fr(4 * (m - 1) + j, 5) * c for c in prev[m - 1]], theta(prev[m - 1]), -1)
                    t = padd(t, pmul_1mu(q))
                cur.append(t)
            Phi.append(cur)
        return Phi
    coef = [Fr(3), Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625)]
    for n in range(1, N + 1):
        trial = ph + [[]]
        Phi = phis(trial, n + 1)
        m = n + 1
        X = padd(Phi[5][m], Phi[0][m], -1)
        br = []
        for i, c in enumerate(coef):        # 3 w Phi_4 + (17/5) u w^2 Phi_3 + ... + (24/625) u^4 w^5 Phi_0
            jj = 4 - i; mm = m - (i + 1)
            if mm >= 0:
                br = padd(br, pmul_u([c * x for x in Phi[jj][mm]], i))
        X = padd(X, pmul_1mu(br), -1)
        # X = (1 - u) r_n ; divide
        r = []; rem = list(X)
        for k in range(len(rem) - 1):
            r.append(rem[k]); rem[k + 1] += rem[k]
        if rem and rem[-1] != 0:
            raise RuntimeError('HPF: [w^(n+1)] not divisible by (1-u) at n=%d' % n)
        new = []
        for k, c in enumerate(r):
            if 5 * k == 4 * n:
                if c != 0:
                    raise RuntimeError('HPF resonance obstructed at n=%d' % n)
                new.append(free.get(n, Fr(0)))
            else:
                new.append(c / (5 * k - 4 * n))
        if 4 * n % 5 == 0 and 4 * n // 5 >= len(new):
            new += [Fr(0)] * (4 * n // 5 - len(new)) + [free.get(n, Fr(0))]
        while new and new[-1] == 0:
            new.pop()
        ph.append(new)
    return ph

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    KS = (1, 2, 3, 4, 5)
    print('g5_degree_bound: N = %d, exact over Q, colour 0' % N)
    okA = okB = okC = okD = okT = True
    for m in range(-25, 41):
        Pp = P_parts(m, KS)
        if Pp[0]:
            okA = False
        for j in range(1, 6):
            q, ok = div_1mx5(Pp[j])
            okB &= ok
            if m >= 1 and Pp[j] and min(Pp[j]) < 1:
                okC = False
        t = Pp[1]; want = {} if m == 1 else {m - 1: Fr(1 - m), m + 4: Fr(m - 1)}
        okT &= (t == want)
        if m >= 1:
            Pi = apply_Pi([{m: Fr(1)}], KS)
            okD &= all((not p) or min(p) >= 6 for p in Pi)
    print('(A) P_0 = 0 on x^m, -25 <= m <= 40:', okA)
    print('(B) P_j(x^m) divisible by (1 - x^5), j = 1..5:', okB)
    print('(C) P_j(x^m) in x Q[x] for 1 <= m <= 40:', okC)
    print('(D) Pi(x^m) in x^6 Q[x][z] for 1 <= m <= 40:', okD)
    print('(T) P_1 x^m = (1 - m)(x^(m-1) - x^(m+4)):', okT)
    for label, free in (('free coefficients 0', {}),
                        ('free coefficients f_m = (-1)^m (2m+3)/7', {5 * m: Fr((-1) ** m * (2 * m + 3), 7) for m in range(1, N // 5 + 1)})):
        psi, info = flat_recursion(N, KS, {1: Fr(1)}, free, lambda m: Fr(1 - m))
        allok = all(i['div_ok'] and i['res_ok'] and i['eq_ok'] and i['grade_ok'] and (i['rho_min'] is None or i['rho_min'] >= 1)
                    and (i['psi_min'] is None or i['psi_min'] >= 1) and (i['psi_max'] is None or i['psi_max'] <= 4 * i['n'] + 1)
                    for i in info)
        ph = hpf_polys(N, free)
        match = True
        for n in range(N + 1):
            want = {4 * n + 1 - 5 * k: c for k, c in enumerate(ph[n]) if c != 0}
            match &= (want == psi[n])
        degs = [len(ph[n]) - 1 for n in range(N + 1)]
        bound = all(degs[n] <= (4 * n) // 5 for n in range(N + 1))
        print('(R) %s: rho_n in x Q[x], resonance trivial, psi_n in x Q[x], deg <= 4n+1, n <= %d:' % (label, N), allok)
        print('(H) %s: psi_n = x^(4n+1) ph_n(x^-5) with HPF ph_n, deg ph_n <= floor(4n/5):' % label, match, bound)
    print('first HPF polynomial ph_5 (free coefficient 0):', hpf_polys(5, {})[5])
    # control X1: corrupt psi_0 by x^-4 (same grade 1 mod 5)
    psi, info = flat_recursion(min(N, 12), KS, {1: Fr(1), -4: Fr(1)}, {}, lambda m: Fr(1 - m))
    bad = [i['n'] for i in info if (not i['res_ok']) or (i['rho_min'] is not None and i['rho_min'] < 1)]
    obst = [i['n'] for i in info if not i['res_ok']]
    print('(X1) control psi_0 = x + x^-4: steps with rho_n not in x Q[x] or obstructed resonance:', bad[:6], '(must be non-empty);',
          'obstructed resonances (x^0 term of rho_n != 0) at n =', obst[:4], '; min x-exponent of rho_1..rho_4:',
          [i['rho_min'] for i in info[:4]])
    # control X2: KS = (1,2,3,4,10); T_1 has kernel x^2
    KS2 = (1, 2, 3, 4, 10)
    okC2 = all((not P_parts(m, KS2)[j]) or min(P_parts(m, KS2)[j]) >= 2 for m in range(2, 12) for j in range(1, 6))
    t1 = P_parts(3, KS2)[1]
    print('(X2) control KS = (1,2,3,4,10): P_1 x^3 =', t1, '; valuation property "P_j(x^m) in x^2 Q[x] for m >= 2":', okC2, '(must be False)')

if __name__ == '__main__':
    main()
