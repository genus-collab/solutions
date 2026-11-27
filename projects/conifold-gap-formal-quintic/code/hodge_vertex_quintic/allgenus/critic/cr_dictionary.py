"""Critic check: the dictionary between the model's slice K (vanishing covector tower Pi_i, model worker)
and the TR leg class L (xyswap worker), on the model curve x = W = p^4/4 + v/p, y = p.

TR legs of ancestor insertions (DNOPS (2-35) dualised; TR worker: leg(tau_0(1)) = y = p, leg(tau_1(1)) = -Int y dx):
  primary legs (up to one constant per b): 1 -> p, H -> log p, H^2 -> v/p, H^3 -> p^3, H^4 -> p^2
  (p-primitives of d_alpha W for the flat directions a_4, a_5, a_1, a_2, a_3);  leg(tau_(k+1)) = -Int leg(tau_k) dx.
Test A: every leg is (Laurent polynomial) + log p . f(x), f a polynomial in x; mu_0 := f(0).  Claim tested:
  mu_0(tau_k(H^b)) = c_b eps^k PiC(tau_k(H^b)),   PiC(tau_k(H^b)) = (-1)^(k+1) [z^-k] Y^C_b  (model worker's table),
  Y^C = (I_0, -I_-1, -I_3, I_2, I_1), I^C_m = Sum_i (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!.
  (The shifts mu_i = i! [x^i] f follow mu_i(k+1) = -mu_(i-1)(k) automatically, as Pi_i(tau_(k+1)) = Pi_(i-1)(tau_k).)
Test B: in each weight w (p:1, v:5, x:4) the log-free combinations of v^m tau_k(H^b) (= K_0 if Test A holds) have
  legs in Q[v][p] + Q(v)[x]  (xyswap's class L, modulo functions of x, which pair to zero); their number equals
  #{(b, k, m): b != 1} (dim K_0 in weight w).
Usage: py -I -B cr_dictionary.py [KMAX] [WMAX]
"""
import os
import sys
from fractions import Fraction as Fr
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cr_laurent import add, scal, mul, integ, deriv, split_log, as_poly_in_x, l_reduce  # noqa: E402

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 20
WMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 60
X = {(4, 0, 0): Fr(1, 4), (-1, 0, 1): Fr(1)}
XP = deriv(X)
LEG0 = {0: {(1, 0, 0): Fr(1)}, 1: {(0, 1, 0): Fr(1)}, 2: {(-1, 0, 1): Fr(1)}, 3: {(3, 0, 0): Fr(1)},
        4: {(2, 0, 0): Fr(1)}}
WT0 = {0: 1, 1: 0, 2: 4, 3: 3, 4: 2}
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg)


def leg_tables():
    T = {}
    for b in range(5):
        cur = LEG0[b]
        for k in range(KMAX + 1):
            w = WT0[b] + 4 * k
            for (j, l, m) in cur:
                check(j + 5 * m == w, "weight b=%d k=%d" % (b, k))
            sp = split_log(cur)
            check(max(sp) <= 1, "log power > 1 at b=%d k=%d" % (b, k))
            f = sp.get(1, {})
            fx = as_poly_in_x(f, X, Fr(1, 4)) if f else {}
            check(fx is not None, "log coefficient not a polynomial in x at b=%d k=%d" % (b, k))
            T[(b, k)] = (sp.get(0, {}), fx or {})
            cur = scal(integ(mul(cur, XP)), -1)
    return T


def pic(b, k):
    """PiC(tau_k(H^b)) as (coefficient, v-power) or (0, None)."""
    mm, sgn = [(0, 1), (-1, -1), (3, -1), (2, 1), (1, 1)][b]
    n = k
    if (n - mm - 1) % 5:
        return Fr(0), None
    i = (n - mm - 1) // 5
    e = mm + 4 * i + 1
    if i < 0 or e < 0:
        return Fr(0), None
    c = Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** e, factorial(e)) * sgn
    return Fr((-1) ** (k + 1)) * c, e


def test_a(T):
    print("== Test A: TR monodromy constant term mu_0 versus the vanishing covector PiC")
    ratios = {}
    for b in range(5):
        row = []
        for k in range(KMAX + 1):
            f = T[(b, k)][1]
            mu = [(m, c) for (n, m), c in f.items() if n == 0]
            pc, pv = pic(b, k)
            if pc == 0:
                check(not mu, "mu_0 != 0 where PiC = 0: b=%d k=%d" % (b, k))
                continue
            check(len(mu) == 1 and mu[0][0] == pv, "mu_0 missing or wrong v-power: b=%d k=%d" % (b, k))
            if mu:
                r = mu[0][1] / pc
                ratios[(b, k)] = r
                row.append("k=%d: mu0=%s v^%d PiC=%s ratio=%s" % (k, mu[0][1], pv, pc, r))
        print(" b=%d:" % b, "; ".join(row))
    eps = None
    for b in range(5):
        ks = sorted(k for (bb, k) in ratios if bb == b)
        for k in ks[1:]:
            q = ratios[(b, k)] / ratios[(b, ks[0])]
            e = q if (k - ks[0]) % 2 else None
            if e is not None:
                eps = e if eps is None else eps
                check(q == eps, "ratio not eps^k: b=%d k=%d" % (b, k))
            else:
                check(q == 1, "even shift ratio != 1: b=%d k=%d" % (b, k))
    cb = {b: ratios[(b, min(k for (bb, k) in ratios if bb == b))] for b in range(5)}
    print(" eps =", eps, "; c_b eps^k0 =", cb)
    print(" non-zero comparisons:", len(ratios), "; zero-pattern slots:", 5 * (KMAX + 1) - len(ratios))


def rank_kernel(rows, ncol):
    """Left kernel of a list of rows (dicts col -> Fraction) over Q."""
    n = len(rows)
    M = [dict(r) for r in rows]
    E = [{i: Fr(1)} for i in range(n)]
    piv_rows = []
    used = [False] * n
    for col in range(ncol):
        p = next((i for i in range(n) if not used[i] and M[i].get(col, 0) != 0), None)
        if p is None:
            continue
        used[p] = True
        piv_rows.append(p)
        for i in range(n):
            if i != p and M[i].get(col, 0) != 0:
                f = M[i][col] / M[p][col]
                for c, val in M[p].items():
                    M[i][c] = M[i].get(c, 0) - f * val
                M[i] = {c: val for c, val in M[i].items() if val != 0}
                for c, val in E[p].items():
                    E[i][c] = E[i].get(c, 0) - f * val
    kern = [E[i] for i in range(n) if not used[i]]
    return len(piv_rows), kern


def test_b(T):
    print("== Test B: log-free combinations (= K_0) have legs in L = Q[v][p] + Q(v)[x]")
    tot_ok, tot_vec = 0, 0
    for w in range(0, WMAX + 1):
        S = [(b, k, (w - WT0[b] - 4 * k) // 5) for b in range(5) for k in range(KMAX + 1)
             if w - WT0[b] - 4 * k >= 0 and (w - WT0[b] - 4 * k) % 5 == 0]
        if not S:
            continue
        cols = {}
        rows = []
        for (b, k, m) in S:
            r = {}
            for (n, mm), c in T[(b, k)][1].items():
                key = (n, mm + m)
                if key not in cols:
                    cols[key] = len(cols)
                r[cols[key]] = r.get(cols[key], 0) + c
            rows.append(r)
        rk, kern = rank_kernel(rows, len(cols))
        expect = sum(1 for (b, k, m) in S if b != 1)
        check(len(kern) == expect, "w=%d: dim log-free = %d, expected %d" % (w, len(kern), expect))
        nok = 0
        reps = []
        for vec in kern:
            phi = {}
            for i, c in vec.items():
                b, k, m = S[i]
                for (j, mm), d in T[(b, k)][0].items():
                    phi[(j, mm + m)] = phi.get((j, mm + m), 0) + c * d
            phi = {key: c for key, c in phi.items() if c != 0}
            ok, rem = l_reduce(phi, X, w)
            nok += ok
            check(ok, "w=%d: a log-free combination is NOT in L (remainder %s)" % (w, sorted(rem.items())[:3]))
            reps.append({j: c for (j, l, m), c in rem.items() if j != 0})   # canonical mod functions of x
        cols2 = sorted({j for r in reps for j in r})
        rk2, _ = rank_kernel([{cols2.index(j): c for j, c in r.items()} for r in reps], len(cols2))
        nL = w // 5 + 1 - (1 if w % 5 == 0 else 0)          # dim Q[v][p]_w modulo constants
        check(rk2 == len(kern) == nL, "w=%d: K_0 -> L/(functions of x) not bijective (%d, %d, %d)" % (w, rk2, len(kern), nL))
        tot_ok += nok
        tot_vec += len(kern)
        if w <= 60 or w % 20 == 0:
            print(" w=%2d: |S|=%2d rank(log map)=%2d dim K_0=%2d (expected %2d), in L: %d/%d, rank mod x = %d = dim L_w/Q[x] %d"
                  % (w, len(S), rk, len(kern), expect, nok, len(kern), rk2, nL))
    print(" total log-free basis vectors tested:", tot_vec, "; in L:", tot_ok)


def control(T):
    print("== Control: a combination with a NON-zero log part fails, and z^-5-type legs fail the L test")
    ok, rem = l_reduce({(-5, 2): Fr(1)}, X, 5)
    check(not ok, "control p^-5 v^2 should fail L test")
    print(" p^-5 v^2 in L:", ok, "(expected False)")
    ok2, _ = l_reduce({(-1, 1): Fr(1)}, X, 4)
    print(" v/p in L:", ok2, "(expected True: v/p = x - p^4/4)")
    check(ok2, "v/p should be in L")
    pc = pic(1, 5)
    f = T[(1, 5)][1]
    print(" tau_5(H): PiC =", pc, "; TR log coefficient f(x) =", sorted(f.items()))


if __name__ == "__main__":
    print("cr_dictionary.py KMAX=%d WMAX=%d (exact over Q, v symbolic)" % (KMAX, WMAX))
    T = leg_tables()
    test_a(T)
    test_b(T)
    control(T)
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS))
