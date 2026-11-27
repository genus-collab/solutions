# -*- coding: utf-8 -*-
"""gi_hat_table.py -- the gauge identity in the hatted Fourier conventions, over Q (GAUGE_IDENTITY.md s.5, s.6).

 (A) REGRESSION of the dictionary on the linear-response data of record (natural_form/extra/nat_kregular_linear.log, NATURAL_FORM
     s.1 item 6): a left dressing X(zhat) (a few Fourier components zhat^k Xhat(m, m')) is ANNIHILATING iff [X y]_(>0) = 0, y the
     leading vanishing covector.  Prediction: annihilating => invisible.  The Euler-dilaton block of G_2 must be annihilating.
 (B) the table natural_form/referee/G_fourier_NZ16_K12.json (the ACTUAL G of the formal quintic, from-scratch rows of a referee,
     exact rationals, zhat <= 16, u <= 11, frame (0,0)):
       (a) every coefficient [u^N zhat^e], e >= 1, of  Gt(u, zhat) y(tau zhat)  vanishes   (Gt = Psi0 G Psi0^-1);
       (b) the rest equals the hatted quintic vanishing solution Yhat^Q / (g_a(u) u tau^(-1/2)), g_a = tau^((1-c_a)/4) p_{0,a}(u);
       (c) leaf identity: Pi(T') = -t_c + u tau^(-3/2)/xhat  (in the form [zhat^-1] W_0 = -t_c xhat tau^(1/2)/u).
Run: py -I -B gi_hat_table.py      log: gi_hat_table.log
"""
import os, sys, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from gi_core import Q, Ser, PROJ, gauge_frame, quintic_vanishing
from gi_hat_core import (CEXP, EPS, y_coeffs, tau_series, tau_neg_powers, flat_entry_of_fourier, apply_term_to_y, positive_part,
                         bivector_name)
OUT = open(os.path.join(HERE, "gi_hat_table.log"), "w", encoding="ascii"); FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); OUT.write(s + "\n"); OUT.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail != "" else ""))
    if not ok: FAILED.append(name)


# (k, {(m, m'): coefficient}, status on record, label)
A = Fr(-127, 15750)
SINGLE = [
    (1, {(0, 1): 1, (1, 0): 1}, "INVISIBLE", "Pi-type zhat (H^3 (x) H^4 + H^4 (x) H^3)"),
    (1, {(2, 4): 1, (4, 2): 1}, "INVISIBLE", "Pi-type zhat (1 (x) H^2 + H^2 (x) 1)"),
    (1, {(3, 3): 1}, "VISIBLE", "zhat H (x) H"),
    (2, {(0, 2): 1, (2, 0): -1}, "INVISIBLE", "zhat^2 (H^4 (x) H^2 - H^2 (x) H^4)"),
    (2, {(3, 4): 1, (4, 3): -1}, "VISIBLE", "zhat^2 (H (x) 1 - 1 (x) H)"),
    (3, {(0, 3): 1, (3, 0): 1}, "VISIBLE", "zhat^3 (H^4 (x) H + H (x) H^4)"),
    (3, {(1, 2): 1, (2, 1): 1}, "INVISIBLE", "zhat^3 (H^3 (x) H^2 + H^2 (x) H^3)"),
    (3, {(4, 4): 1}, "VISIBLE", "zhat^3 1 (x) 1"),
    (4, {(0, 4): 1, (4, 0): -1}, "VISIBLE", "zhat^4 (H^4 (x) 1 - 1 (x) H^4)"),
    (4, {(1, 3): 1, (3, 1): -1}, "VISIBLE", "zhat^4 (H^3 (x) H - H (x) H^3)"),
    (5, {(0, 0): 1}, "VISIBLE", "zhat^5 H^4 (x) H^4"),
    (5, {(1, 4): 1, (4, 1): 1}, "VISIBLE", "zhat^5 (H^3 (x) 1 + 1 (x) H^3)"),
    (5, {(2, 3): 1, (3, 2): 1}, "VISIBLE", "zhat^5 (H^2 (x) H + H (x) H^2)"),
    (6, {(0, 1): 1, (1, 0): -1}, "INVISIBLE", "zhat^6 (H^4 (x) H^3 - H^3 (x) H^4)"),
    (6, {(2, 4): 1, (4, 2): -1}, "VISIBLE", "zhat^6 (H^2 (x) 1 - 1 (x) H^2)"),
]


def entries_of(k, comps):
    out = {}
    for (m, mp), c in comps.items():
        a, b, sg = flat_entry_of_fourier(m, mp)
        out[(k, a, b)] = out.get((k, a, b), Fr(0)) + 5 * sg * Fr(c)
    return out


def part_A():
    log("=" * 100); log("(A) regression of the dictionary on single-term left dressings (status on record: nat_kregular_linear.log, both primes)")
    y = y_coeffs(12)
    log("   leading covector y(zhat), flat components b = 0..4 (classes 1, H, H^2, H^3, H^4): %s" % [
        " + ".join("(%s) zhat^-%d" % (c, n) for n, c in sorted(d.items())[:2]) for d in y])
    agree = 0; exceptions = []
    for k, comps, status, label in SINGLE:
        ent = entries_of(k, comps)
        pos = positive_part(apply_term_to_y(ent, y, 12))
        ann = not pos
        names = ", ".join("%s E_%d%d [%s]" % (c, a, b, bivector_name(a, b)) for (kk, a, b), c in sorted(ent.items()))
        log("   zhat^%d %-22s flat: %-60s annihilates y: %-5s on record: %s" % (k, sorted(comps), names, ann, status))
        if ann == (status == "INVISIBLE"): agree += 1
        else: exceptions.append((label, "annihilating" if ann else "not annihilating", status))
    check("every ANNIHILATING single term is INVISIBLE on record (the implication predicted by the mechanism)",
          all(st == "INVISIBLE" for k, comps, st, lb in SINGLE if not positive_part(apply_term_to_y(entries_of(k, comps), y, 12))))
    check("the three terms listed VISIBLE in the brief (zhat H(x)H, zhat^2 (H(x)1 - 1(x)H), zhat^3 1(x)1) are NOT annihilating",
          all(positive_part(apply_term_to_y(entries_of(k, comps), y, 12)) for k, comps in ((1, {(3, 3): 1}), (2, {(3, 4): 1, (4, 3): -1}), (3, {(4, 4): 1}))))
    check("the Pi-type terms and zhat^2 (H^4(x)H^2 - H^2(x)H^4) ARE annihilating",
          all(not positive_part(apply_term_to_y(entries_of(k, comps), y, 12)) for k, comps in ((1, {(0, 1): 1, (1, 0): 1}), (1, {(2, 4): 1, (4, 2): 1}), (2, {(0, 2): 1, (2, 0): -1}))))
    log("   agreement annihilating <=> invisible: %d of %d single terms; exceptions: %s" % (agree, len(SINGLE), exceptions))
    # Euler-dilaton block of G_2 (values of the production table): Xhat(3,3) zhat, Xhat(3,4) zhat^2, Xhat(4,3) zhat^2, Xhat(4,4) zhat^3
    ED = {}
    for k, comps in ((1, {(3, 3): A}), (2, {(3, 4): A, (4, 3): -A}), (3, {(4, 4): -A})): ED.update(entries_of(k, comps))
    pos = positive_part(apply_term_to_y(ED, y, 12))
    log("   Euler-dilaton block of [u^2] G (alpha = %s): flat entries %s" % (A, {("zhat^%d" % k, "E_%d%d" % (a, b)): str(c) for (k, a, b), c in sorted(ED.items())}))
    check("the Euler-dilaton block alpha zhat (H - zhat)(x)(H + zhat) of G_2 ANNIHILATES y (all signs/transposes of the dictionary)", not pos, pos)
    EDc = {}
    for k, comps in ((1, {(3, 3): A}), (2, {(3, 4): -A, (4, 3): A}), (3, {(4, 4): -A})): EDc.update(entries_of(k, comps))
    posc = positive_part(apply_term_to_y(EDc, y, 12))
    check("control: the opposite relative sign zhat (H + zhat)(x)(H - zhat) (visible on record: -1/150, -7/7200) does NOT annihilate y", bool(posc),
          {key: str(v) for key, v in posc.items()})


def load_table():
    tab = json.load(open(os.path.join(PROJ, "natural_form", "referee", "G_fourier_NZ16_K12.json"), encoding="utf-8"))
    G = {}
    for key, val in tab["G_fourier"].items():
        un, zk, mm = key.split(); m, mp = mm.strip("()").split(",")
        G[(int(un[1:]), int(zk[1:]), int(m), int(mp))] = Fr(val)
    return G, tab["deg"]


def part_B():
    log("=" * 100); log("(B) the actual G (referee table over Q, frame (0,0), zhat <= 16, u <= 11)")
    G, deg = load_table(); NU = 12; NY = 34
    log("   %d Fourier components; zhat-degrees of [u^n] G, n = 0..11: %s" % (len(G), deg))
    Gt = {}                                                 # (k, a, b) -> list of 12 Fractions
    for (n, k, m, mp), val in G.items():
        a, b, sg = flat_entry_of_fourier(m, mp)
        Gt.setdefault((k, a, b), [Fr(0)] * NU)[n] += 5 * sg * val
    check("Gt(u, zhat = 0) = 1 and Z/5 grading a = b + k", all((a - b - k) % 5 == 0 for (k, a, b) in Gt) and all(
        Gt[(0, a, a)] == [Fr(1)] + [Fr(0)] * (NU - 1) for a in range(5)))
    GtS = {key: Ser.from_list(Q, c, 0, NU) for key, c in Gt.items()}
    tau = tau_series(NU); tn = tau_neg_powers(tau, NY, NU); y = y_coeffs(NY)
    W = {}
    for (k, a, b), g in GtS.items():
        if not g.c: continue
        for n, yc in y[b].items():
            t = g.mul(tn[n]).scale(yc).truncate(NU)
            d = W.setdefault(a, {}); d[k - n] = t if (k - n) not in d else d[k - n].add(t)
    kmax = max(k for (k, a, b) in Gt)
    # (a)
    tested = 0; bad = []; genuine = 0
    for a in range(5):
        for e in range(1, kmax + 1):
            if (e - (a - 1)) % 5: continue
            for N in range(NU):
                if e > (5 * N + 2) // 4: continue            # beyond the degree of [u^N] G: empty slot
                tested += 1
                val = W[a][e].coef(N) if e in W.get(a, {}) else Fr(0)
                if val != 0: bad.append((a, e, N, str(val)))
                cnt = sum(1 for (k, aa, b), g in Gt.items() if aa == a and (k - e) in y[b] for n1 in range(N + 1) if g[n1] != 0 and tn[k - e].coef(N - n1) != 0)
                if cnt >= 2: genuine += 1
    check("(a) ALL %d graded coefficients [u^N zhat^e] (e >= 1, e <= floor((5N+2)/4), N <= 11) of Gt y(tau zhat) vanish; %d of them are sums of >= 2 non-zero products"
          % (tested, genuine), not bad, bad[:4])
    # the slot that pins c_5: a = 1, e = 5, N = 4
    terms = [("zhat^%d E_1%d" % (k, b), str(g[n1]), "y_%d[zhat^-%d] tau^-%d [u^%d]" % (b, k - 5, k - 5, 4 - n1), str(y[b][k - 5] * tn[k - 5].coef(4 - n1)))
             for (k, aa, b), g in sorted(Gt.items()) if aa == 1 and (k - 5) in y[b] for n1 in range(5) if g[n1] != 0 and tn[k - 5].coef(4 - n1) != 0]
    log("   the slot (a, e, N) = (1, 5, 4) [where P -> P (1 + c z^5) would enter as -c]: terms (G entry, value, y factor, value): %s" % terms)
    # (b)
    res = gauge_frame(0, 40); p0 = res['p'][0]; C = res['C']; lt = tau.log()
    YQ = quintic_vanishing(res['bq'], 30, NY)
    xh = C.varpi0.truncate(NU)
    check("p_{0,0}(u) = tau^(-1/4)/varpi_0  (so g_0 = tau^(1/4) p_{0,0} = 1/xhat)", p0[0].truncate(NU).sub(lt.scale(Fr(-1, 4)).exp().mul(xh.inv())).truncate(NU).iszero())
    pref = [lt.scale(Fr(1 - CEXP[a], 4) - Fr(1, 2)).exp().mul(p0[a].truncate(NU)).truncate(NU) for a in range(5)]      # g_a tau^(-1/2)
    tested = 0; bad = []
    for a in range(5):
        for j in range(0, NY - 12):
            if (-j - (a - 1)) % 5: continue
            ex = 4 - CEXP[a] - 4 * j; assert ex % 5 == 0
            rhs = YQ[-j][a].shift(ex // 5) if (-j) in YQ else Ser.zero(Q, NU)      # Yhat^Q_a [zhat^-j]
            lhs = (W[a][-j] if -j in W.get(a, {}) else Ser.zero(Q, NU)).mul(pref[a]).shift(1)      # g_a u tau^(-1/2) W_a
            for N in range(NU):
                if N >= rhs.prec or N >= lhs.prec: continue
                tested += 1
                if lhs.coef(N) != rhs.coef(N): bad.append((a, -j, N, str(lhs.coef(N)), str(rhs.coef(N))))
    check("(b) g_a u tau^(-1/2) [zhat^-j] (Gt y(tau zhat))_a = [zhat^-j] Yhat^Q_a for all a, j = 0..%d, u^0..u^11 (%d coefficients, quintic side from the periods only)"
          % (NY - 13, tested), not bad, bad[:4])
    # (c)
    tc = C.tc.truncate(NU + 1)
    lhs = W[0][-1]
    rhs = tc.shift(-1).mul(xh).mul(lt.scale(Fr(1, 2)).exp()).neg().truncate(NU)
    check("(c) LEAF: [zhat^-1] (Gt y(tau zhat))_0 = -t_c xhat tau^(1/2)/u through u^11, i.e. Pi(T') = -t_c + u tau^(-3/2)/xhat", lhs.sub(rhs).truncate(NU - 1).iszero() and lhs.sub(rhs).iszero())
    vhat = lt.scale(Fr(-3, 2)).exp().mul(xh.inv()).shift(1).truncate(NU)
    d = tc.sub(vhat).truncate(8)
    log("   vhat = u tau^(-3/2)/xhat = %s" % vhat.fmt(4)); log("   t_c - vhat = -Pi(T') = %s" % d.fmt(5))
    check("t_c - vhat = u^2/30 + O(u^3) in frame (0,0): one leaf gives (2-2g) c_g u^(1-2g) u^2/30 = -(g-1) c_g/15 u^(3-2g)", d.coef(2) == Fr(1, 30) and d.coef(1) == 0)


if __name__ == "__main__":
    t0 = time.time(); part_A(); part_B()
    log("GI_HAT_TABLE %s   (%.0f s)" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED, time.time() - t0)); OUT.close()
