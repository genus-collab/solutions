# -*- coding: utf-8 -*-
"""rf_extra.py -- referee's extra checks, single process (monkeypatches rf_bsys IN MEMORY only; no file is changed):
  branch  S8 = the formal quintic system (B unperturbed) pulled back by the OTHER real branch v = -u - ... of v'^5 = -B v (m = -1;
          real fifth roots).  The system's R-matrix is intrinsic, so the rows must be the project rows and every graph must equal R0.
  nosig   S5 without the constant symplectic scalar sigma (the worker's overwritten first prep run): rows not symplectic.
  count   non-triviality counts of the worker (graphs differing from R0 in the window) from the referee's own json tables.
Usage: py -I -B rf_extra.py branch g prime | nosig g prime | count
"""
import os, sys, json, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_run
import rf_bsys as RB
from rf_run import SCR, PR

_root = RB.root_rat


def root_rat_signed(x, e):
    x = Fr(x)
    if x < 0:
        assert e.denominator % 2 == 1, "even root of a negative number"
        sign = -1 if e.numerator % 2 else 1                    # integer sign: (-1) ** (negative int) would be a FLOAT
        return sign * _root(-x, e)                             # real root: (-|x|)^(a/b) = (-1)^a |x|^(a/b), b odd
    return _root(x, e)


RB.root_rat = root_rat_signed
RB.SYS["S8"] = ((1,), (1,), (0, 0), (1,), (1,), -1)


class A(object):
    pass


def args(sys_, g, prime, out, ctrl="none"):
    a = A(); a.sys = sys_; a.g = g; a.prime = prime; a.norm = "canon"; a.ctrl = ctrl; a.maxout = None; a.kext = 2; a.nproc = 1; a.out = out
    return a


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("mode"); ap.add_argument("g", nargs="?", type=int, default=2)
    ap.add_argument("prime", nargs="?", default="61"); ap.add_argument("--out", default=SCR)
    a = ap.parse_args(); K = 2 * a.g - 2 + 2
    if a.mode == "branch":
        S = RB.build("S8"); print("S8: v = %r ; tau = %r" % (S['v'], S['tau']))
        Qh, nr = rf_run.make_rows(S, rf_run.system_gauge("S8", "canon")[1], a.prime, a.g, K)
        same = all(Qh[x][i][n] == nr['Qh'][x][i][n] for x in range(5) for i in range(5) for n in range(3 * a.g - 2))
        print("S8 rows = project (quintic) rows: %s" % same)
        rf_run.run(args("S8", a.g, a.prime, a.out))
        t8 = json.load(open(os.path.join(a.out, "rf_S8_g%d_%s_canon_none_k%d.json" % (a.g, a.prime, K))))['table']
        t0 = json.load(open(os.path.join(a.out, "rf_R0_g%d_%s_canon_none_k%d.json" % (a.g, a.prime, K))))['table']
        print("RF_EXTRA branch g%d %s: S8 tables = R0 tables on every graph: %s" % (a.g, a.prime, t8 == t0))
    elif a.mode == "nosig":
        orig = RB.normalise
        def nosig(p, mode):
            pn, C, sig, c0, m5 = orig(p, mode)
            return RB.zmul(pn, {0: 1 / sig}), C, Fr(1), c0 / sig, m5
        RB.normalise = nosig
        out = os.path.join(a.out, "nosig"); os.makedirs(out, exist_ok=True)
        rf_run.run(args("S5", a.g, a.prime, out))
    elif a.mode == "staleboth":
        # S1 gauge, but the WHOLE hatted dressing (prefactor tau^((1-c_b)/4), tau^(n-k) on Rh, and E_M) uses the quintic tau
        S, pn, inf = rf_run.system_gauge("S1", "canon"); S0 = RB.build("R0"); S2 = dict(S); S2['w'] = S0['w']; S2['tau'] = S0['tau']
        p = PR[a.prime]; Qh, nr = rf_run.make_rows(S2, pn, a.prime, a.g, K)
        from un_core import Arith, lift_rows, symplectic_check
        ar = Arith(p, K, 1); sb = symplectic_check(ar, lift_rows(ar, Qh), 3 * a.g - 3, nr['eta_inv'])
        taur = [rf_run.fr_mod(S0['tau'].co(n), p) for n in range(K)]
        v, tab = rf_run.graph_sum((a.prime, a.g, K, taur, Qh, None), 1, "sb")
        npol = 2 * a.g - 2; sing = [i for i in range(1, len(v.inner_graphs)) if any(tab.get(i, [0] * K)[:npol])]
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(len(v.inner_graphs))) % p for k in range(K)]
        Tq = S['T']; V = rf_run.in_T(tab.get(0, [0] * K), Tq, a.g, K, p); Vt = rf_run.in_T(tot, Tq, a.g, K, p)
        print("RF_EXTRA staleboth S1 g%d %s: rows symplectic %s; singular edge graphs %s; vertex %s; total %s" % (
            a.g, a.prime, not sb, sing or "none", [rf_run.rat(x, p) for x in V], [rf_run.rat(x, p) for x in Vt]))
    elif a.mode == "norm":
        # (stored - canon)/m5 per outer graph from the referee's genus-4 maxout-2 tables of S2 at a.prime
        p = PR[a.prime]; tg = "rf_S2_g4_%s_%s_none_out2_k6.json"
        ca = json.load(open(os.path.join(a.out, tg % (a.prime, "canon")))); st = json.load(open(os.path.join(a.out, tg % (a.prime, "stored"))))
        m5 = rf_run.system_gauge("S2", "stored")[2]['m5']; im5 = pow(rf_run.fr_mod(m5, p), p - 2, p); outer = ca['outer']
        TA = {int(i): c for i, c in ca['table'].items()}; TB = {int(i): c for i, c in st['table'].items()}; tot = [0] * 6
        print("S2 m5 = %s" % m5)
        for i in range(len(outer)):
            if i > 0 and outer[i][0] > 2: continue
            d = [(TB.get(i, [0] * 6)[k] - TA.get(i, [0] * 6)[k]) * im5 % p for k in range(6)]
            tot = [(x + y) % p for x, y in zip(tot, d)]
            if any(d): print("   graph %d %s: (stored - canon)/m5 = %s" % (i, outer[i], [rf_run.rat(x, p) for x in d]))
        print("RF_EXTRA norm S2 g4 %s: sum over graphs (<= 2 edges) = %s" % (a.prime, [rf_run.rat(x, p) for x in tot]))
    elif a.mode == "rows4":
        # genus-4 R0 regression (rows only) and symplecticity of the genus-4 rows of every system, both primes
        from un_core import Arith, lift_rows, symplectic_check
        for prime in ("89", "61"):
            for s in ("R0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"):
                S_, pn, inf = rf_run.system_gauge(s, "canon"); Qh, nr = rf_run.make_rows(S_, pn, prime, 4, 6)
                ar = Arith(PR[prime], 6, 1); sb = symplectic_check(ar, lift_rows(ar, Qh), 9, nr['eta_inv'])
                same = all(Qh[x][i][n] == nr['Qh'][x][i][n] for x in range(5) for i in range(5) for n in range(10))
                print("RF_EXTRA rows4 %s %s: symplectic mod (zhat^10, u^6) %s; = project rows %s" % (s, prime, not sb, same))
    elif a.mode == "count":
        for g in (2, 3):
            K = 2 * g; t0 = json.load(open(os.path.join(a.out, "rf_R0_g%d_89_canon_none_k%d.json" % (g, K))))['table']
            for s in ("S1", "S2", "S3", "S4", "S5", "S6", "S7"):
                t = json.load(open(os.path.join(a.out, "rf_%s_g%d_89_canon_none_k%d.json" % (s, g, K))))['table']
                keys = set(t) | set(t0); dif = [k for k in keys if t.get(k, [0] * K) != t0.get(k, [0] * K)]
                print("genus %d %s: %d of %d outer graphs differ from R0 in u^0..u^%d" % (g, s, len(dif), len(keys), K - 1))


if __name__ == "__main__":
    main()
