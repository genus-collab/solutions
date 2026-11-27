"""Critic check: conventions, signs and normalisations that differ between the worker notes (exact over Q).

(1) vertex constant c_g in the five forms used by NATURAL_FORM / UNIFORM_REDUCTION / model / TR / LG_MODEL s.5.4, and the TR
    free energies on x = W and on x = lambda = p^4 + 4/p (xyswap's -1/3840, 1/258048);
(2) the vanishing covector: model PiC, gauge phi (sign list of GAUGE_IDENTITY_HAT s.2), plan/gauge Y^C, theory frames;
(3) normalisation response at genus 4: univ per-graph [u^4] numbers versus the model's Delta_4(5) (total);
(4) leaf bookkeeping -(g-1) c_g/15 and its frame (in)dependence;
(5) the Euler-dilaton edge legs of XYSWAP s.4 (worker -(2/5) z^9 versus referee -(16/45) z^9).
Usage: py -I -B cr_conv.py
"""
import os
import sys
from fractions import Fraction as Fr
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cr_laurent import add, scal, mul, integ, deriv, split_log  # noqa: E402

FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def bernoulli(N):
    B = [Fr(1)]
    for n in range(1, N + 1):
        B.append(-sum(comb(n + 1, k) * B[k] for k in range(n)) / (n + 1))
    return B


def part1():
    print("== (1) vertex constants")
    B = bernoulli(26)
    rec = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000), 5: Fr(1, 660000)}
    for g in range(2, 13):
        b = B[2 * g]
        fW = b / (2 * g * (2 * g - 2))                   # TR, x = W (TR worker, xyswap: standard EO kernel)
        c_tr = Fr(-5) ** (1 - g) * fW                    # TR worker: frozen = (-5)^(1-g) F_g(x = W)
        c_abs = abs(b) / (2 * g * (2 * g - 2) * 5 ** (g - 1))     # brief / NATURAL_FORM s.6
        c_ur = (-1) ** (g - 1) * b / (2 * g * (2 * g - 2)) / 5 ** (g - 1)   # UNIFORM_REDUCTION s.6.5, model
        c_lg = Fr(20) ** (1 - g) * Fr(-4) ** (3 * g - 3) * Fr(4) ** (2 - 2 * g) * fW   # LG_MODEL s.5.4, a_5 = 4
        ok = c_tr == c_abs == c_ur == c_lg and (g not in rec or c_abs == rec[g])
        check(ok, "g=%d: c_g = %s in all four forms%s" % (g, c_abs, " (= record)" if g in rec else ""))
    lam = {g: Fr(-5, 16) ** (g - 1) * abs(B[2 * g]) / (2 * g * (2 * g - 2) * 5 ** (g - 1)) for g in (2, 3)}
    check(lam[2] == Fr(-1, 3840) and lam[3] == Fr(1, 258048),
          "x = p^4 + 4/p: F_2, F_3 = %s, %s (xyswap direct TR: -1/3840, 1/258048)" % (lam[2], lam[3]))
    check([B[2 * g] / (2 * g * (2 * g - 2)) for g in (2, 3, 4, 5)] == [Fr(-1, 240), Fr(1, 1008), Fr(-1, 1440), Fr(1, 1056)],
          "TR F_2..F_5 on x = W = -1/240, 1/1008, -1/1440, 1/1056")


def pic(b, k):
    mm, sgn = [(0, 1), (-1, -1), (3, -1), (2, 1), (1, 1)][b]
    if (k - mm - 1) % 5:
        return Fr(0)
    i = (k - mm - 1) // 5
    e = mm + 4 * i + 1
    if i < 0 or e < 0:
        return Fr(0)
    yc = Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** e, factorial(e)) * sgn   # [z^-k] Y^C_b at v = 1
    return (-1) ** (k + 1) * yc


def part2():
    print("== (2) vanishing covector conventions (v = 1)")
    model = {(1, 0): 1, (0, 1): -1, (4, 2): Fr(-1, 2), (3, 3): Fr(-1, 6), (2, 4): Fr(1, 24), (1, 5): Fr(1, 96),
             (0, 6): Fr(-1, 480)}
    check(all(pic(b, k) == val for (b, k), val in model.items()),
          "PiC(tau_k(H^b)) = (-1)^(k+1) [z^-k] Y^C_b reproduces the model's printed PiC list (MODEL_SINGULAR s.1)")
    gauge = {(1, 0): -1, (0, 1): 1, (4, 2): Fr(1, 2), (3, 3): Fr(1, 6), (2, 4): Fr(-1, 24), (1, 5): Fr(-1, 96)}
    check(all(-pic(b, k) == val for (b, k), val in gauge.items()),
          "gauge phi(tau_m(H^c)) / (u tau^(-1/2-m)) = -PiC on the printed list (GAUGE_IDENTITY_HAT s.2)")
    ok_sup, ok_mag = True, True
    for m in range(0, 61):
        q = m // 5
        nz = [b for b in range(5) if pic(b, m) != 0]
        ok_sup &= nz == [(1 - m) % 5]
        ok_mag &= all(abs(pic(b, m)) == Fr(1, 4 ** q * factorial(q) * factorial(m - q)) for b in nz)
    check(ok_sup, "support rule: PiC(tau_m(H^b)) != 0 iff b = 1 - m mod 5 (m <= 60)")
    check(ok_mag, "closed form |phi_m| = 1/(4^q q! (m-q)!), q = floor(m/5) (m <= 60)")
    frames = {(1, 0): Fr(1), (2, 4): Fr(1, 2), (3, 3): Fr(1, 6), (4, 2): Fr(-1, 24)}   # frame(k,a) = z^k H^a + f H
    check(all(pic(a, k) + f * pic(1, 0) == 0 for (k, a), f in frames.items()),
          "theory worker's frames frame(1,0), (2,4), (3,3), (4,2) lie in ker PiC")
    c3 = Fr(1, 25200)
    law = pic(1, 5) * (2 - 2 * 3) * c3
    check(law == Fr(-1, 604800) and Fr(53, 1209600) != law,
          "one-point law at <tau_5(H)>_3: PiC F_3' = %s, value on record 53/1209600 (law REFUTED off K)" % law)


def part3():
    print("== (3) z^5 normalisation response at genus 4 (per unit m5, [u^4])")
    vertex, g1, g2, loop = Fr(-1, 8400000), Fr(-1, 1008000), Fr(1, 900000), Fr(0)
    check(g1 + g2 == Fr(1, 8400000), "univ bridges (1,3) + (2,2) = 1/8400000 = model's 'bridges' term")
    check(vertex + g1 + g2 + loop == 0, "sum over graphs = 0 = model's Delta_4(5) (a statement about the TOTAL)")
    check(g1 != 0 and g2 != 0, "per graph the response is non-zero: graphs 1, 2 singular at [u^4] (Conjecture N sees it)")


def part4():
    print("== (4) leaf bookkeeping")
    for name, gt, w1 in (("(0,0)", Fr(7, 10), Fr(0)), ("(1/3,-2/7)", Fr(11, 30), Fr(1, 3))):
        check(gt + w1 == Fr(7, 10) and 2 * gt + 2 * w1 - Fr(4, 3) == Fr(1, 15) and Fr(2, 3) - gt - w1 == Fr(-1, 30),
              "frame %s: gamma_t + w_1 = 7/10, 2 gamma_t + 2 w_1 - 4/3 = 1/15, p_14(0)/2 = -1/30" % name)
    for g, c in ((2, Fr(1, 1200)), (3, Fr(1, 25200))):
        need = -(2 * g - 2) * Fr(7, 10) * c + (3 * g - 3) * Fr(4, 9) * c
        check(need == -(g - 1) * c / 15, "g=%d: gap minus leafless u^1 term = -(g-1) c_g/15 = %s" % (g, need))


def part5():
    print("== (5) Euler-dilaton edge legs, curve x = z^4 + v/z, tau_0(1) <-> -z, tau_0(H) <-> -log z, tau_k = Int tau_(k-1) dx")
    X = {(4, 0, 0): Fr(1), (-1, 0, 1): Fr(1)}
    XP = deriv(X)
    t01, t0H = {(1, 0, 0): Fr(-1)}, {(0, 1, 0): Fr(-1)}
    t11 = integ(mul(t01, XP))
    t1H = integ(mul(t0H, XP))
    t21 = integ(mul(t11, XP))
    e1 = add(t11, mul({(0, 0, 1): Fr(1)}, t0H))                     # psi 1 + v H
    e2 = add(t21, mul({(0, 0, 1): Fr(1)}, t1H))                     # psi (psi 1 + v H)
    check(e1 == {(5, 0, 0): Fr(-4, 5)}, "psi 1 + v H <-> -(4/5) z^5")
    check(e2 == {(9, 0, 0): Fr(-16, 45), (4, 0, 1): Fr(1, 5)},
          "psi(psi 1 + v H) <-> -(16/45) z^9 + (v/5) z^4 (referee right; worker's -(2/5) z^9 is a slip) : %s" % sorted(e2.items()))
    cur, logs = e1, []
    for k in range(1, 7):
        cur = integ(mul(cur, XP))
        logs.append(max(split_log(cur)))
    check(logs[:4] == [0, 0, 0, 0] and logs[4] == 1, "Int^k(-(4/5) z^5) dx log-free for k <= 4, log z at k = 5: %s" % logs)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    part4()
    part5()
    print("RESULT:", "ALL CHECKS PASSED" if not FAILS else "%d FAILURES" % len(FAILS))
