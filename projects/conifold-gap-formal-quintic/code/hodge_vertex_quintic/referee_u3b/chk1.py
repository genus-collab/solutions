# -*- coding: utf-8 -*-
"""chk1.py -- referee NOTE: light checks of statements of s.6d against the stored forms (read-only use of the project data)."""
import sys, json, os, random
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)
from fq3_circsys import CircSys
from fq3_cubic import IDX_3EDGE, CubicForm
from fq3_quad import IDX_2EDGE


def rank(M, p):
    M = [list(r) for r in M]; rk = 0; nc = len(M[0]) if M else 0
    for c in range(nc):
        piv = None
        for r in range(rk, len(M)):
            if M[r][c] % p: piv = r; break
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], p - 2, p); M[rk] = [x * inv % p for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] % p:
                f = M[r][c]; M[r] = [(x - f * y) % p for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def rat(x, p, B=10 ** 6):
    """small rational reconstruction"""
    for d in range(1, 2000):
        n = x * d % p
        if n > p // 2: n -= p
        if abs(n) < B: return "%d/%d" % (n, d)
    return "[%d]" % x


def main(tag):
    S = CircSys(tag); p = S.p; cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % tag)))
    print("=" * 30, "prime", tag)
    A = cs['Pi_A']; B = cs['Pi_B']; star = cs['circ_star']; sol = cs['circ_b_solutions'][0]; extra = cs['circ_b_extra']
    z2 = (-B[2]) % p; z4 = (-2 * A[2]) % p; z = z2 * z4 % p
    print("zeta^5 = 1:", pow(z, 5, p) == 1, " zeta != 1:", z != 1, " z^2 ok:", z * z % p == z2)
    # --- (1) c_32 = kappa m_1^3
    m = S.m1; c32 = S.C[32].coef
    kap1 = c32[(1, 1, 1)] * pow(pow(m[1], 3, p), p - 2, p) % p
    cube = {}
    for i in range(10):
        for j in range(i, 10):
            for k in range(j, 10):
                mult = 6 if (i < j < k) else (1 if i == j == k else 3)
                v = kap1 * mult * m[i] * m[j] * m[k] % p
                if v: cube[(i, j, k)] = v
    print("(1) c_32 == kappa m_1^3 as polynomials:", cube == c32, " kappa != 0:", kap1 != 0, " monomials:", len(c32))
    # --- c_18, c_20, c_32 on m_1 = 0 ; c_11, c_19 divisible by c1_0
    basis = []
    for j in range(10):
        if j == 4: continue
        v = [0] * 10; v[j] = 1; v[4] = (-m[j]) * pow(m[4], p - 2, p) % p; basis.append(v)
    van = [i for i in IDX_3EDGE if not S.C[i].restrict(basis).coef]
    print("    cubics vanishing identically on m_1 = 0:", van)
    print("    cubics all of whose monomials contain c1_0:", [i for i in IDX_3EDGE if all(0 in k for k in S.C[i].coef)])
    # m_1^2 | c_20 ?  (restrict the gradient)
    print("    order-4 dependence (index 8 or 9 in support):", [i for i in IDX_3EDGE if any(8 in k or 9 in k for k in S.C[i].coef)])
    # --- (2) stratum c1 = 0
    nz = {}
    for i in IDX_3EDGE:
        sub = {k: v for k, v in S.C[i].coef.items() if not (set(k) & {0, 1, 2})}
        if sub: nz[i] = sorted(set(x for k in sub for x in k))
    print("(2) c1 = 0: non-zero cubics and their variables:", nz)
    b2 = [v for v in basis if not any(v[j] for j in (0, 1, 2))]
    nz2 = {}
    for i in IDX_3EDGE:
        r = CubicForm(p, 10, {k: v for k, v in S.C[i].coef.items() if not (set(k) & {0, 1, 2})}).restrict(
            [[0, 0, 0, 1, (-m[3]) * pow(m[4], p - 2, p) % p, 0, 0, 0, 0, 0]])
        if r.coef: nz2[i] = r.coef
    print("    c1 = 0 and m_1 = 0 (c2 on its line): non-zero cubics:", sorted(nz2))
    # --- (3) Pi and Delta*
    res = lambda v: [k for k, x in S.residuals(v).items() if x]
    m1v = lambda v: sum(m[j] * v[j] for j in range(15)) % p
    x, y = random.randrange(p), random.randrange(p)
    P = [(x * a + y * b) % p for a, b in zip(A, B)]
    print("(3) residuals at A, B, random xA+yB:", res(A), res(B), res(P), " m_1:", m1v(A), m1v(B))
    print("    Delta*: non-zero residuals", res(star), " m_1 =", m1v(star), " c3_0 =", rat(star[5], p), " c5_0, c6_2 =", star[10], star[14])
    a = star[1]; PiPt = [(u + a * w) % p for u, w in zip(A, B)]
    print("    c1 of Delta* == c1 of A + a B:", star[:3] == PiPt[:3], " c3_0 of that Pi point:", rat(PiPt[5], p),
          " c3_0(A), c3_0(B):", rat(A[5], p), rat(B[5], p))
    print("    Delta* - Pi point, non-zero coordinates:", [S.names[j] for j in range(15) if (star[j] - PiPt[j]) % p])
    print("    rank of (A, B, Delta*) on coordinates 0..9:", rank([A[:10], B[:10], star[:10]], p))
    print("    line of Pi at c1(Delta*): c1_0 + 2 z^3 c1_1 + 2 z c1_2 =", (star[0] + 2 * z * z2 * star[1] + 2 * z * star[2]) % p)
    # tangency of the c1-line of Pi with the conic a_4 = 0
    def a4(c1):
        t = random.randrange(p)
        return (S.cub(4, S.low(c1, t, 1) + [0, 0]) - S.cub(4, S.low(c1, t, 0) + [0, 0])) % p
    i2z = pow(2 * z, p - 2, p)
    line = lambda s: [1, s % p, (-(1 + 2 * z * z2 * s)) * i2z % p]
    f = [a4(line(s)) for s in (0, 1, 2, 3)]
    d3 = (f[3] - 3 * f[2] + 3 * f[1] - f[0]) % p
    q2 = (f[2] - 2 * f[1] + f[0]) * pow(2, p - 2, p) % p; q1 = (f[1] - f[0] - q2) % p; q0 = f[0]
    disc = (q1 * q1 - 4 * q2 * q0) % p
    print("    a_4 on the c1-line of Pi: degree <= 2:", d3 == 0, " leading != 0:", q2 != 0, " discriminant == 0 (tangent):", disc == 0,
          " double root == c1_1(Delta*):", (-q1 * pow(2 * q2, p - 2, p)) % p == a)
    # --- (4) the 17 points on c1_0 = 0 that are stored (1 genuine F_p-rational + 4 special)
    print("(4) stored c1_0 = 0 points: residuals", res(sol), [res(e) for e in extra])
    for n, e in enumerate(extra):
        s = e[1]; lowprop = all((e[j] - s * B[j]) % p == 0 for j in range(10))
        bb = e[2] * pow(e[1], p - 2, p) % p
        print("    extra[%d]: b = %s zeta^2;  coordinates 0..9 proportional to Pi_B: %s" % (n, "+" if bb == z2 else ("-" if bb == (-z2) % p else "?"), lowprop))
    print("    U_b1 degree:", len(cs['U_b1']) - 1, " genuine degree:", cs['U_b1_genuine_degree'], " factor degrees:", cs['U_b1_genuine_factors'],
          " F_p roots of the genuine part:", len(cs['U_b1_genuine_roots']))
    # --- (5) N_2 system rank: rows J3 (2) + grad q_Gamma(Delta) (5) on the 15 circulant coordinates
    def n2rank(v):
        return rank([S.J3[0], S.J3[1]] + [S.quad_grad(i, v) for i in IDX_2EDGE], p)
    hi = [0] * 10 + [random.randrange(p) for _ in range(5)]
    AB2 = [(u + 2 * w) % p for u, w in zip(A, B)]
    print("(5) rank of the 7 x 15 N_2 matrix: random pure-high", n2rank(hi), " A+2B", n2rank(AB2), " random Pi", n2rank(P), " B", n2rank(B),
          " genuine c1_0=0 point", n2rank(sol), " extras", [n2rank(e) for e in extra], " Delta*", n2rank(star),
          " random point", n2rank([random.randrange(p) for _ in range(15)]))
    # scaling: all N_1-conditions are homogeneous in Delta
    s = random.randrange(2, p)
    print("    homogeneity: residuals at s*Delta*, s*sol:", res([s * t % p for t in star]), res([s * t % p for t in sol]))
    # fibre dimension over the low part: linear system for (c4, c5, c6) at fixed low part (8 coordinates)
    def fibre_rank(v):
        rows = []
        for r in range(2): rows.append([S.J3[r][j] for j in range(8, 15)])
        for i in IDX_2EDGE: g = S.quad_grad(i, v); rows.append([g[j] for j in range(8, 15)])
        for i in IDX_3EDGE: g = S.C[i].gradient(v[:10]); rows.append([g[8], g[9], 0, 0, 0, 0, 0])
        return rank(rows, p)
    print("    rank of the conditions on the 7 high coordinates (fibre dimension = 7 - rank): A+2B", fibre_rank(AB2), " B", fibre_rank(B),
          " genuine", fibre_rank(sol), " extras", [fibre_rank(e) for e in extra], " Delta*", fibre_rank(star))


if __name__ == "__main__":
    random.seed(20261110)
    for tag in ("61", "89"): main(tag)
