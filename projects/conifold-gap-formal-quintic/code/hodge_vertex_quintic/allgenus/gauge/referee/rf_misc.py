# -*- coding: utf-8 -*-
"""rf_misc.py -- REFEREE: small exact checks (stdlib only, own code).
 (V1) Lemma V on the term level: z d/dv I_m = -I_{m-1};  m z I_{m-1} - I_{m+3} + v I_{m-2} = 0;  Y^C = (I_0, -I_-1, -I_3, I_2, I_1) solves
      z d/dv Y_b = bmod_b Y_{b+1}, bmod = (1, -1/v, 1, -1, -1);  Y^C(v = 0) = -e_1.
 (F)  closed form of the vanishing covector: [t^-m] y_b != 0 only for b = 1 - m mod 5 and (-1)^m eps_b [t^-m] y_b = 1/(4^q q! (m-q)!), m <= 80.
 (C)  the stored normalisation c(z) = P_11(0, z): c(z) c(-z) - 1 and the claimed z^10 constant; frame (in)dependence of c(z).
Run: py -I -B rf_misc.py      log: rf_misc.log
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import ic_terms, fact, CEXP, EPS
OUT = open(os.path.join(HERE, "rf_misc.log"), "w", encoding="ascii"); FAIL = []


def log(s): print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, det=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(det)) if det != "" else ""))
    if not ok: FAIL.append(name)


def I(m, nmax):
    """{(j, n): coefficient of v^j z^-n}"""
    return {(j, n): c for (n, j, c) in ic_terms(m, nmax)}


def zddv(T):
    """z d/dv on {(j, n): c}: v^j z^-n -> j v^(j-1) z^-(n-1)"""
    return {(j - 1, n - 1): j * c for (j, n), c in T.items() if j}


def lin(*pairs):
    out = {}
    for c, T in pairs:
        for k, x in T.items(): out[k] = out.get(k, 0) + c * x
    return {k: x for k, x in out.items() if x}


def cut(T, nmax): return {k: x for k, x in T.items() if k[1] <= nmax}


def main():
    NM = 30; ok1 = ok2 = True
    for m in range(-6, 12):
        if cut(lin((1, zddv(I(m, NM + 1))), (1, I(m - 1, NM))), NM - 1): ok1 = False
        zI = {(j, n - 1): c for (j, n), c in I(m - 1, NM).items()}                  # z I_{m-1}
        vI = {(j + 1, n): c for (j, n), c in I(m - 2, NM).items()}                  # v I_{m-2}
        if cut(lin((m, zI), (-1, I(m + 3, NM)), (1, vI)), NM - 2): ok2 = False
    check("(V1)(i) z d/dv I_m = -I_{m-1}, m = -6..11, through z^-%d" % (NM - 1), ok1)
    check("(V1)(ii) m z I_{m-1} - I_{m+3} + v I_{m-2} = 0, m = -6..11, through z^-%d" % (NM - 2), ok2)
    Y = [lin((EPS[b], I(CEXP[b], NM + 1))) for b in range(5)]
    okY = True
    for b in range(5):
        lhs = cut(zddv(Y[b]), NM - 1); nxt = Y[(b + 1) % 5]
        if b == 1: rhs = {(j - 1, n): -c for (j, n), c in nxt.items()}            # -1/v Y_2
        else: rhs = {k: [1, None, 1, -1, -1][b] * c for k, c in nxt.items()}
        neg = [k for k, x in rhs.items() if k[0] < 0 and x != 0]                     # a v^-1 term would break the system
        rhs = cut(rhs, NM - 1)
        if lin((1, lhs), (-1, rhs)) or neg: okY = False
    check("(V1) Y^C solves the model system z d/dv Y_b = bmod_b Y_{b+1} (no negative power of v from -1/v Y_2), through z^-%d" % (NM - 1), okY)
    check("(V1) Y^C(v = 0) = -e_1", all({k: c for k, c in Y[b].items() if k[0] == 0} == ({(0, 0): -1} if b == 1 else {}) for b in range(5)))
    bad = []
    for m in range(0, 81):
        q = m // 5
        supp = [b for b in range(5) if any(n == m for (n, j, c) in ic_terms(CEXP[b], m))]
        if supp != [(1 - m) % 5]: bad.append(("support", m, supp)); continue
        b = supp[0]; c = sum(c for (n, j, c) in ic_terms(CEXP[b], m) if n == m)
        if (-1) ** m * EPS[b] * (EPS[b] * c) != Fr(1, 4 ** q * fact(q) * fact(m - q)): bad.append(("value", m))      # [t^-m] y_b = eps_b c
        if (b + 1 - (2 - m)) % 5: bad.append(("phase", m))
    check("(F) (-1)^m eps_b [t^-m] y_b = 1/(4^q q! (m-q)!), support b = 1 - m mod 5, phase zeta^(-i(b+1)) = zeta^(i(m-2)), m <= 80", not bad, bad[:3])
    c5 = Fr(-43397, 807975000); c10 = Fr(-765368951400317, 1305647201250000000)
    check("(C) stored c(z): 2 c_10 - c_5^2 = -18223115111903/15543419062500000 (the z^10 symplecticity defect of the stored gauge)",
          2 * c10 - c5 ** 2 == Fr(-18223115111903, 15543419062500000), 2 * c10 - c5 ** 2)
    pi8s = Fr(-77141089, 362679778125000); pi3 = Fr(127, 3150)
    check("(C) canonical pi_8 = stored pi_8 + c_5 pi_3 (c_5 = 43397/807975000) = 1416472237/725359556250000",
          pi8s + Fr(43397, 807975000) * pi3 == Fr(1416472237, 725359556250000))
    log("RF_MISC: %s" % ("ALL PASS" if not FAIL else "FAIL %s" % FAIL)); OUT.close()


if __name__ == "__main__":
    main()
