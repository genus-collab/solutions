# -*- coding: utf-8 -*-
"""v3_fourier.py -- the Fourier components of G = Psi(u)^-1 Ghat Psi(0) predicted from the own gauge (v2_gauge.py, periods only),
against the numbers PRINTED by the production run (logs/genus3_natural_89_nat.log and ..._ab1o3_m2o7.log; typed in here by hand).

Derivation (own): Psi(0)[b][i] = eps_b zeta^(i(b+1)), eps = F_b(0) = (1,-1,-1,1,1); Psi(u) = Ghat(u, 0) Psi(0) (Ghat(u,0) diagonal), so
    Xhat(m, m') := (1/25) Sum_ij zeta^(m i) G_ij zeta^(m' j) = (1/5) eps_a eps_b [Ghat(u,0)^-1 Ghat(u,zhat)]_{ab},  a = m - 1, b = -m' - 1 (mod 5),
    [Ghat(u,0)^-1 Ghat]_{ab} at zhat^k (a = b + k) = tau^((c_a - c_b)/4) u^mu(b,k) p_{k,b}(u) / p_{0,a}(u).
In the frame (0,0) the gauge is multiplied by (1 + c_5 z^5), c_5 = 43397/807975000; irrelevant for zhat <= 4.
Run: py -I -B v3_fourier.py
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from v2_gauge import gauge
CEXP = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)
OUT = open(os.path.join(HERE, "v3_fourier.log"), "w", encoding="utf-8")


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def xhat(res, n, k):
    """{(m, m'): value} of [u^n zhat^k] G"""
    out = {}
    for b in range(5):
        a = (b + k) % 5; mu = (4 * k + CEXP[b] - CEXP[a]) // 5
        f = res['tau'].powq(Fr(CEXP[a] - CEXP[b], 4)) * res['p'][k][b] * res['p'][0][a].inv()
        val = f.coef(n - mu) if n >= mu else Fr(0)
        if val: out[((a + 1) % 5, (-b - 1) % 5)] = Fr(EPS[a] * EPS[b], 5) * val
    return out


PROD00 = {
    (1, 1): {(0, 1): '-1/75', (1, 0): '-1/75', (2, 4): '-8/225', (4, 2): '-8/225'},
    (2, 1): {(0, 1): '-1/378', (1, 0): '-1/378', (2, 4): '1/1134', (3, 3): '-127/15750', (4, 2): '1/1134'},
    (2, 2): {(0, 2): '-17/23625', (1, 1): '1/2250', (2, 0): '73/23625', (3, 4): '-127/15750', (4, 3): '127/15750'},
    (2, 3): {(4, 4): '127/15750'},
    (3, 1): {(0, 1): '-43/32319', (1, 0): '-43/32319', (2, 4): '716117/302990625', (3, 3): '-167269/33665625', (4, 2): '716117/302990625'},
    (3, 2): {(0, 2): '71/323190', (1, 1): '1/5670', (2, 0): '31/161595', (3, 4): '-39766/33665625', (4, 3): '29342/11221875'},
    (3, 3): {(0, 3): '-8783/67331250', (1, 2): '-51512/100996875', (2, 1): '-9191/14428125', (3, 0): '1958/4809375', (4, 4): '-13159/11221875'},
    (3, 4): {(0, 4): '-8783/67331250', (4, 0): '-1958/4809375'},
}
PRODAB = {(1, 1): {(0, 1): '-1/75', (1, 0): '-1/75', (2, 4): '7/225', (4, 2): '7/225'}}


def main():
    ok_all = True
    for name, (al, be), table in (("frame (0,0)", (0, 0), PROD00), ("frame (1/3,-2/7)", (Fr(1, 3), Fr(-2, 7)), PRODAB)):
        res = gauge(4, 16, al, be)
        for (n, k), want in sorted(table.items()):
            got = xhat(res, n, k)
            ok = got == {key: Fr(v) for key, v in want.items()}
            ok_all = ok_all and ok
            log(("PASS " if ok else "FAIL ") + "%s [u^%d zhat^%d] G: predicted %s" % (name, n, k, {key: str(v) for key, v in sorted(got.items())}))
        # vanishing beyond the degree: [u^n zhat^k] = 0 for k > deg(n) (k <= 4 here)
        deg = {1: 1, 2: 3, 3: 4}
        z = all(not xhat(res, n, k) for n in deg for k in range(deg[n] + 1, 5))
        z0 = all(not xhat(res, n, 0) for n in (1, 2, 3))
        log(("PASS " if z and z0 else "FAIL ") + "%s: [u^n zhat^k] G = 0 for k > deg_n (k <= 4) and for k = 0, n = 1..3" % name); ok_all = ok_all and z and z0
    log("V3_FOURIER " + ("ALL PASS" if ok_all else "FAILURES")); OUT.close()


if __name__ == "__main__":
    main()
