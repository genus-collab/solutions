# -*- coding: utf-8 -*-
"""gi_hat_core.py -- the vanishing covector in the HATTED conventions of the graph sum (GAUGE_IDENTITY.md s.5).  Exact, stdlib only.

Flat hatted index b in Z/5 (class Hhat^b), colour i in Z/5, zhat = z/s^4, u = s^5, c = (0, -1, 3, 2, 1), eps = (1, -1, -1, 1, 1).
  y_b(t) := Y^C_b(v = 1, z = t)                      (the model's vanishing solution at v = 1; y_b in t^-(c_b+1) Q[[t^-5]])
          = (-1/t, -1, -1/(24 t^4), -1/(6 t^3), 1/(2 t^2)) + (relative corrections in t^-5)
  hatted model covector (flat lower-index components):   Yhat^mod_b(zhat; u) = u tau^(-1/2) y_b(tau zhat)
        (= s^3 sigma^(1-c_b) Y^C_b(v(u), z), sigma = s tau^(-1/4); homogeneity of the model: deg v = 5, deg z = 4)
  hatted quintic covector:   Yhat^Q_a(zhat; u) = s^(4-c_a) Y^C_{Q,a}(u, s^4 zhat) = Sum_n u^((4 - c_a - 4n)/5) zhat^-n y^(n)_a(u)
  identity (G) in hatted form:    Ghat(u, zhat) Yhat^mod = Yhat^Q,      Ghat = hatted gauge (Qhat = Ghat Rhat(tau zhat))
  colour form:  f_i(x) = Sum_b (Psi0^-1)_{ib} Yhat^mod_b(x),  Psi0[b][i] = eps_b zeta^(i(b+1)),  (Psi0^-1)[i][b] = (1/5) eps_b zeta^(-i(b+1));
        G(u, zhat) f has no positive power of zhat  (G = Psi(u)^-1 Ghat Psi0 the colour x colour factor of R~ = G M_0(tau zhat)).
  functional on the half-edge / leaf insertions of a vertex:  Pi(tau_m(e~_i)) = phi_{i,m} := (-1)^m [x^-m] f_i(x).
Fourier dictionary (verify_theory/v3_fourier.py):  Xhat(m, m') = (1/25) Sum_ij zeta^(m i) G_ij zeta^(m' j)
        = (1/5) eps_a eps_b [Ghat(u,0)^-1 Ghat(u,zhat)]_{ab},   a = m - 1,  b = -m' - 1  (mod 5);  the entry (a, b) is Hhat^(3-a) (x) Hhat^b.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import Q, Ser, ic_terms, MEXP, SGN
CEXP = (0, -1, 3, 2, 1); EPS = (1, -1, -1, 1, 1)


def y_coeffs(nmax):
    """y[b] = {n: coefficient of t^-n in y_b(t)}, n <= nmax"""
    out = []
    for b in range(5):
        d = {}
        for (n, j, c) in ic_terms(MEXP[b], nmax): d[n] = d.get(n, Fr(0)) + c * SGN[b]
        out.append(d)
    return out


def tau_series(NU):
    """tau(u) = 1 / Sum_n 4 u^n/(5n + 4) as a Ser over Q, O(u^NU)"""
    return Ser.from_list(Q, [Fr(4, 5 * n + 4) for n in range(NU)], 0, NU).inv()


def tau_neg_powers(tau, nmax, NU):
    """[tau^0, tau^-1, ..., tau^-nmax] truncated at O(u^NU)"""
    ti = tau.inv().truncate(NU); out = [Ser.const(Q, 1, NU)]
    for n in range(1, nmax + 1): out.append(out[-1].mul(ti).truncate(NU))
    return out


def flat_entry_of_fourier(m, mp):
    """(a, b, sign): Xhat(m, m') = (1/5) sign * [flat hatted entry (a, b)]"""
    a = (m - 1) % 5; b = (-mp - 1) % 5
    return a, b, EPS[a] * EPS[b]


def fourier_of_flat_entry(a, b):
    return (a + 1) % 5, (-b - 1) % 5


def apply_term_to_y(entries, y, nmax):
    """entries: {(k, a, b): coefficient} = Sum coefficient * zhat^k E_ab (flat hatted matrix units, u-independent).
    Returns {a: {power of zhat: coefficient}} of (X y)(zhat) with y the leading covector y(zhat) (tau = 1)."""
    out = {}
    for (k, a, b), c in entries.items():
        for n, yc in y[b].items():
            if n > nmax: continue
            d = out.setdefault(a, {}); d[k - n] = d.get(k - n, Fr(0)) + c * yc
    return {a: {e: v for e, v in d.items() if v != 0} for a, d in out.items()}


def positive_part(vec):
    return {(a, e): v for a, d in vec.items() for e, v in d.items() if e > 0}


def bivector_name(a, b):
    nm = lambda k: "1" if k == 0 else ("H" if k == 1 else "H^%d" % k)
    return "%s (x) %s" % (nm((3 - a) % 5), nm(b))


def edge_degree_bound(N):
    """zhat-degree bound of [u^N] of 1 - G(-psi)^T G(-psi') (two factors of G): max over n1 + n2 = N of d(n1) + d(n2), d(n) = floor((5n+2)/4)"""
    d = lambda n: (5 * n + 2) // 4
    return max(d(n1) + d(N - n1) for n1 in range(N + 1))
