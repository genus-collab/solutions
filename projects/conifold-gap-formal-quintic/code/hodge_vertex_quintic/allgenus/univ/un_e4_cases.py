# -*- coding: utf-8 -*-
"""un_e4_cases.py -- E4/U4 inputs: the genus-4 cases (cyclic systems, normalisation of the Theorem-U gauge, hatted rows).

Normalisation of the regular gauge P (Theorem U fixes P up to P -> P c(z^5), c a constant scalar series):
  'stored'  the output of natural_form/theory/ur_gauge.gauge: the free constant of every resonant order k = 0 mod 5 is set to 0,
            i.e. p_{k,0}(0) = 0 (k = 5, 10, ..): P_00(0, z) = 1.  For the quintic P_11(0, z) = 1 - (43397/807975000) z^5 - ...
  'canon'   P -> P c0/P_11(0, z), c0 = P_11(0, 0): then P_11(0, z) = c0 is constant, i.e. (G) P Y^C_model(v(u)) = c0 Y^C_system holds
            with no positive power of z (c0 = 1 for the quintic, c0 = 1/lambda for a system with f(0) = lambda).
In hatted variables z^5 E_{b+5,b} = u^4 zhat^5 E_{b+5,b}, so at genus 4 (window u^0..u^5) only the z^5 coefficient of the normalisation
is visible, and the two normalisations differ by the scalar right factor 1 + m5 u^4 zhat^5 (m5 = [z^5] P^stored_11(0, z)/c0).
A case is (F, Gm, frame, normalisation); F, Gm as in un_gauge.system: f = t_c' F(u), g = -(varpi_0/(1-u)) Gm(u), b = (f, B/(f^2 g^2), f,
g, g), B = -1/(u(1-u)^5) (so tau, v(u) are those of the quintic).  F = Gm = 1 is the formal quintic in that conifold frame.
API: CASES; build(name) -> (res, info);  rows(res, prime, g, K) -> (Qh plain lists, nat);  row_checks(...).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from fq_fp import QRing
from fq_series import Ser
from un_gauge import pert_gauge, rows_from_gauge, check_G
from un_core import Arith, lift_rows, symplectic_check, right_scalar
from un_engine import get_nat
Q = QRing(); F = Fr

SYS = {
    "quintic": ((1,), (1,), ("0", "0")),
    "quintic_ab": ((1,), (1,), ("1/3", "-2/7")),
    "pertA": ((1, F(1, 3), F(-2, 7), F(5, 11)), (1, F(-3, 5), F(1, 2), F(2, 3)), ("0", "0")),          # un_e2 case 5
    "pertB": ((1, F(-7, 4), 3, F(-1, 9)), (1, F(11, 6), F(-5, 3), F(4, 13)), ("0", "0")),               # un_e2 case 6
    "lam32": ((F(3, 2), F(1, 3)), (1, F(-3, 5), F(1, 2)), ("0", "0")),                                  # un_e2 case 7, f(0) = 3/2
    "pert_ab": ((1, F(2, 5), F(-1, 3)), (1, F(1, 7), F(3, 4)), ("1/3", "-2/7")),                       # un_e2 case 8
}
CASES = {   # name: (system, normalisation, description)
    "a": ("quintic", "canon", "formal quintic, frame (0,0), CANONICAL normalisation P_11(0,z) = 1 (must reproduce the natural run)"),
    "b": ("quintic", "stored", "formal quintic, frame (0,0), STORED normalisation P_00(0,z) = 1 (P_11(0,z) = 1 - c5 z^5 - ..)"),
    "c1": ("pertA", "canon", "perturbed system f(1 + u/3 - 2u^2/7 + 5u^3/11), g(1 - 3u/5 + u^2/2 + 2u^3/3), CANONICAL"),
    "c2": ("lam32", "canon", "perturbed system with f(0) = 3/2: f(3/2 + u/3), g(1 - 3u/5 + u^2/2), CANONICAL (P_11(0,z) = 2/3)"),
    "d": ("pertA", "stored", "perturbed system of c1 in the STORED normalisation"),
    "c3": ("pertB", "canon", "perturbed system f(1 - 7u/4 + 3u^2 - u^3/9), g(1 + 11u/6 - 5u^2/3 + 4u^3/13), CANONICAL"),
    "c4": ("pert_ab", "canon", "perturbed system in frame (1/3,-2/7): f(1 + 2u/5 - u^2/3), g(1 + u/7 + 3u^2/4), CANONICAL"),
    "e": ("quintic_ab", "canon", "formal quintic, frame (1/3,-2/7), CANONICAL"),
}
_GAUGE = {}


def gauge_of(system, KG=12, NU=24):
    g = _GAUGE.get(system)
    if g is None:
        Fp, Gp, ab = SYS[system]; g = pert_gauge(KG, NU, Fp, Gp, ab); _GAUGE[system] = g
    return g


def p11_series(p):
    """P_11(0, z) = Sum_{k = 0 mod 5} p_{k,1}(0) z^k  as {k: Fraction}"""
    return {k: (p[k][1].coef(0) if p[k][1].prec > 0 else None) for k in range(0, len(p), 5)}


def normalise(res, mode):
    """returns (res', multiplier m(z) = {k: Fraction}) with p' = p m(z)"""
    p = res['p']; KG = len(p) - 1; P11 = p11_series(p)
    assert all(x is not None for x in P11.values()), "P_11(0, z) not available"
    if mode == "stored": return res, {0: Fr(1)}
    assert mode == "canon"
    inv = {0: 1 / P11[0]}
    for k in range(5, KG + 1, 5):
        inv[k] = -sum(P11[j] * inv[k - j] for j in range(5, k + 1, 5)) / P11[0]
    m = {k: P11[0] * x for k, x in inv.items()}
    pc = []
    for k in range(KG + 1):
        row = []
        for b in range(5):
            acc = None
            for j, mj in m.items():
                if j > k or mj == 0: continue
                t = p[k - j][b].scale(mj); acc = t if acc is None else acc.add(t)
            row.append(acc)
        pc.append(row)
    out = dict(res); out['p'] = pc
    return out, m


def build(name):
    system, mode, desc = CASES[name]
    res0 = gauge_of(system); res, m = normalise(res0, mode)
    P11 = p11_series(res['p']); lam = Fr(SYS[system][0][0])
    info = dict(name=name, system=system, mode=mode, desc=desc, frame=SYS[system][2], lam=lam,
                P11={k: str(x) for k, x in P11.items()}, mult={k: str(x) for k, x in m.items()},
                P11_stored={k: str(x) for k, x in p11_series(res0['p']).items()})
    return res, info


def rows(res, prime, g=4, K=None):
    K = 2 * g - 2 if K is None else K
    nat = get_nat(prime, g, ("0", "0"), K)
    return rows_from_gauge(res, nat['Rh'], nat['R'], nat['KZ'], K), nat


def gauge_checks(res, info, check, log):
    hol = all(x is None or x >= 0 for inf in res['info'] for x in inf['vals'])
    check("   Theorem U: every p_{k,b}, k <= %d, of the (stored) gauge is a power series in u" % (len(res['p']) - 1), hol)
    reg = all(s is None or (not s.c) or s.v >= 0 for row in res['p'] for s in row)
    check("   the gauge in the %s normalisation is regular at u = 0" % info['mode'], reg)
    ok, c, bad, nchk, rng = check_G(res, len(res['p']) - 1, 6)
    log("   P_11(0, z) in this normalisation: %s" % " + ".join("(%s) z^%d" % (info['P11'][k], k) for k in sorted(info['P11'])))
    check("   (G): P Y^C_model(v(u)) = P_11(0, z) Y^C_system on z^%d..z^%d mod u^6 (%d coefficients of z, exact over Q)" % (rng[0], rng[1], nchk), ok, bad[:3])
    const = all(Fr(info['P11'][k]) == 0 for k in info['P11'] if k > 0)
    log("   => positive powers of z in P Y^C_model: %s" % ("NONE (canonical: P_11(0, z) is constant)" if const else "present (P_11(0, z) is not constant)"))
    return const


def row_checks(Qh, nat, info, prime, check, log, g=4, K=None):
    """symplecticity in the window; for the quintic: regression against the project rows (canon) or the scalar-factor picture (stored)"""
    K = 2 * g - 2 if K is None else K; KZ = nat['KZ']; p = nat['p']; R = nat['R']; ar = Arith(p, K, 1); QL = lift_rows(ar, Qh)
    sb = symplectic_check(ar, QL, KZ, nat['eta_inv'])
    check("   genus %d: the rows Ghat Rhat(tau zhat) are symplectic mod (zhat^%d, u^%d)" % (g, KZ + 1, K), not sb, sb[:4])
    if info['system'] in ("quintic", "quintic_ab"):
        ref = get_nat(prime, g, info['frame'], K)['Qh']; RL = lift_rows(ar, ref)
        if info['mode'] == "canon":
            same = all(Qh[a][i][n] == ref[a][i][n] for a in range(5) for i in range(5) for n in range(KZ + 1))
            check("   genus %d REGRESSION: rows from the canonical gauge = the project rows in frame %s (25 x %d x %d residues)" % (g, info['frame'], KZ + 1, K), same)
        else:
            m5 = Fr(info['P11_stored'][5]) / Fr(info['P11_stored'][0])
            sc = right_scalar(ar, RL, {0: ar.one(), 5: ar.mono(R.of(m5), 4)}, KZ)
            same = all(QL[a][i][n] == sc[a][i][n] for a in range(5) for i in range(5) for n in range(KZ + 1))
            check("   genus %d: stored rows = project rows x (1 + m5 u^4 zhat^5), m5 = %s (scalar right factor; 25 x %d x %d residues)" % (g, m5, KZ + 1, K), same)
            diff = sum(1 for a in range(5) for i in range(5) for n in range(KZ + 1) if Qh[a][i][n] != ref[a][i][n])
            log("      residue slots (a, i, zhat^n) where stored and project rows differ: %d (all at zhat^n, n >= 5)" % diff)
    return QL
