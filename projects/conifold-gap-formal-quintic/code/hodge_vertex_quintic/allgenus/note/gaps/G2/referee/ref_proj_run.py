# -*- coding: utf-8 -*-
"""ref_proj_run.py -- driver of ref_proj.py (see its docstring).  Usage: py -I -B ref_proj_run.py [npoints]."""
import os, sys, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ref_proj import Du, poly_at, zero, mm, val, der, pr_linalg, check, FAIL


def rnd_poly(rng):
    return [Fr(rng.choice([-1, 1]) * rng.randint(1, 9), rng.randint(1, 7)) for _ in range(3)]


def main():
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    rng = random.Random(1115)
    print("ref_proj_run.py: %d random points; general weights (random quadratic polynomials), dual numbers, exact Q" % npts)
    for t in range(npts):
        u0 = Fr(rng.randint(1, 9), rng.randint(2, 9))
        bt = [poly_at(rnd_poly(rng), u0) for _ in range(5)]
        b = [poly_at(rnd_poly(rng), u0) for _ in range(4)]
        B = bt[0] * bt[1] * bt[2] * bt[3] * bt[4]
        b.append(B / (b[0] * b[1] * b[2] * b[3]))                  # (EP): prod b = prod bt as functions
        if any(x.a == 0 for x in bt + b):
            continue
        At = zero(); Y = zero()
        for a in range(5):
            At[a][(a + 1) % 5] = bt[a]; Y[(a + 1) % 5][a] = bt[a].inv()
        one = [[Du(int(i == j)) for j in range(5)] for i in range(5)]
        Ypw = [one]
        for j in range(1, 5):
            Ypw.append(mm(Ypw[-1], Y))
        Atv = val(At); Yv = [val(M) for M in Ypw]
        d = [Du(1)]
        for a in range(5):
            d.append(d[a] * bt[a] / b[a])
        check("pt %d (u0 = %s): Theta_U closes up under (EP): d_5 = 1 (value and derivative)" % (t, u0), d[5].a == 1 and d[5].b == 0)
        Dl = zero()
        for a in range(5):
            Dl[a][a] = Du(d[a].b / d[a].a)
        trDl = sum(d[a].b / d[a].a for a in range(5))
        dlogB = B.b / B.a
        # (c), (d): rank of ad At and directness, from the linear algebra itself
        X = [[Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(5)] for _ in range(5)]
        cX, rk = pr_linalg(Atv, Yv, X)
        check("pt %d: [Y^0..Y^4 | ad At] has rank 25 (gl_5 = Q[Y] + im ad At, direct)" % t, rk == 25)
        # closed formula of lem:gau-comm (e) versus the linear algebra
        prc = [[Fr(0)] * 5 for _ in range(5)]
        for m in range(5):
            lam = Fr(0)
            for a in range(5):
                cal = Fr(1)
                for i in range(m if m else 5):
                    cal *= bt[(a + i) % 5].a
                lam += X[a][(a + m) % 5] / cal
            Am = one
            for _ in range(m if m else 5):
                Am = mm(Am, At)
            for r in range(5):
                for s in range(5):
                    prc[r][s] += Fr(1, 5) * lam * Am[r][s].a
        prl = [[sum(cX[j] * Yv[j][r][s] for j in range(5)) for s in range(5)] for r in range(5)]
        check("pt %d: closed formula pr(X) = (1/5) lambda_m(X) At^m (each shift) = linear-algebra pr" % t, prc == prl)
        # G3, first fact
        ok1 = True; ok1c = True
        for j in range(5):
            cj, _ = pr_linalg(Atv, Yv, der(Ypw[j]))
            want = [Fr(0)] * 5; want[j] = -Fr(j, 5) * dlogB
            ok1 = ok1 and cj == want
            if j:
                bad = [Fr(0)] * 5; bad[j] = -Fr(j, 4) * dlogB
                ok1c = ok1c and cj != bad
        check("pt %d: pr((Y^j)') = -(j/5)(Bprod'/Bprod) Y^j, j = 0..4" % t, ok1)
        check("pt %d: control: with j/4 the first fact fails (every j >= 1)" % t, ok1c)
        # G3, second fact
        ok2 = True
        for j in range(5):
            cj, _ = pr_linalg(Atv, Yv, val(mm(Dl, Ypw[j])))
            want = [Fr(0)] * 5; want[j] = trDl / 5
            ok2 = ok2 and cj == want
        check("pt %d: pr(dlogTheta_U Y^j) = (tr dlogTheta_U/5) Y^j, j = 0..4" % t, ok2)
        Nd = [[Du(Fr(rng.randint(-5, 5), 3)) for _ in range(5)] for _ in range(5)]
        ok2c = False
        for j in range(1, 5):
            cj, _ = pr_linalg(Atv, Yv, val(mm(Nd, Ypw[j])))
            want = [Fr(0)] * 5; want[j] = sum(Nd[a][a].a for a in range(5)) / 5
            ok2c = ok2c or cj != want
        check("pt %d: control: a non-diagonal matrix in place of dlogTheta_U breaks the trace rule" % t, ok2c)
        # lem:gau-scal (a) on a random element c = Sum c_j Y^j (c_j random polynomials at u0)
        cs = [poly_at(rnd_poly(rng), u0) for _ in range(5)]
        C = zero()
        for j in range(5):
            for r in range(5):
                for s in range(5):
                    C[r][s] = C[r][s] + cs[j] * Ypw[j][r][s]
        lhs = [[der(C)[r][s] + val(mm(Dl, C))[r][s] for s in range(5)] for r in range(5)]
        cl, _ = pr_linalg(Atv, Yv, lhs)
        want = [cs[j].b + (-Fr(j, 5) * dlogB + trDl / 5) * cs[j].a for j in range(5)]
        check("pt %d: pr(c' + dlogTheta c) = Sum_j (c_j' + kappa_j c_j) Y^j" % t, cl == want)
        # control: Frobenius-orthogonal projection onto Q[Y] is NOT pr (the complement matters)
        G = [[sum(Yv[i][r][s] * Yv[j][r][s] for r in range(5) for s in range(5)) for j in range(5)] for i in range(5)]
        V = der(Ypw[2])
        rhs = [sum(Yv[i][r][s] * V[r][s] for r in range(5) for s in range(5)) for i in range(5)]
        from ref_proj import solve
        co, _ = solve(G, rhs)
        want = [Fr(0)] * 5; want[2] = -Fr(2, 5) * dlogB
        check("pt %d: control: Frobenius-orthogonal projection of (Y^2)' differs from pr" % t, co != want)
        # control: (EP) broken -> Theta_U does not close
        bb = list(b); bb[4] = bb[4] * Du(2)
        dd = Du(1)
        for a in range(5):
            dd = dd * bt[a] / bb[a]
        check("pt %d: control: with prod b != prod bt the cyclic definition of Theta_U does not close" % t, dd.a != 1)
    print("RESULT:", "ALL CHECKS PASSED" if not FAIL else "FAILED: %s" % FAIL)


if __name__ == "__main__":
    main()
