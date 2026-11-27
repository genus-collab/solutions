# -*- coding: utf-8 -*-
"""rf_gauge.py -- REFEREE: the normalisation c^act(z) = P^act_11(0,z) of the gauge of the QRR-normalised WKB function, by three
routes, with the QRR polynomials taken from the referee's large-radius computation (rf_qexp), not from (HPF):
 (F) the worker's formula  c_5q = p_(5q,0)(0) + p'_(5q-1,1)(0)  (POLY s.3), re-implemented;
 (R) rows 1..4 of P rebuilt from row 0 by the gauge equation z P' = A P - P At (row a: P_(a+1,b) = (z P_ab' + bt_(b-1) P_(a,b-1))/b_a),
     closing of the cycle (row 5 = row 0), (R1), (R2), and c read off directly as the entry (1,1) at u = 0;
 (G) Theorem G (GI s.3) in row 0: Sum_b P_0b(u,z) Y^C_b(v(u),z) = c(z) (Y^C_Q)_0(u,z), tested with c = 1 at every u^N, N <= NU, and
     every z-order of the complete window.
Row 0 from the hatted identity row (UR s.4, re-derived): (1/C0) Phi = Sum_k w^k g_k(u) F_(-k)(tau w), g_k = u^mu tau^((1-c_b)/4) p_(k,b),
b = -k mod 5, mu = (4k + c_b)/5.  Controls: stored normalisation (P_00(0,z) = 1), no tau-weights, wrong mu table (4q + r), C0 of
the other frame.  Run: py -I -B rf_gauge.py Q NU    (frames (0,0) and (1/3,-2/7), primes 2^61-1 and 2^89-1)
"""
import os, sys, time, re
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import P61, P89, CEXP, md, Ser, frame_data, model_rows
from rf_qexp import hatted_polys
from rf_vanish import model_vanishing, quintic_vanishing

FRAMES = ((Fr(0), Fr(0)), (Fr(1, 3), Fr(-2, 7)))
M5 = Fr(-43397, 807975000)


def extract(S, polys, C0, tau, F, NK, variant=''):
    """row 0 of P: {k: p_(k,-k)(u)} and the valid lengths"""
    p, K = S.p, S.K; iC0 = S.inv(C0); tp = [S.one()]
    for j in range(NK): tp.append(S.mul(tp[-1], tau))
    twi = {b: S.powr(tau, Fr(CEXP[b] - 1, 4)) for b in range(5)}
    g = {}; pk = {}; prec = {}
    for N in range(NK + 1):
        lhs = S.mul(iC0, (polys[N] + [0] * K)[:K])
        for k in range(N):
            c = F[(-k) % 5][N - k]
            if c: lhs = S.add(lhs, S.mul(g[k], tp[N - k]), (-c) % p)
        bN = (-N) % 5; g[N] = S.scale(lhs, F[bN][0])
        mu = (4 * N + CEXP[bN]) // 5
        if variant == 'mu' and N % 5 == 4: mu += 1
        x = S.shift(g[N], -mu)                                        # asserts divisibility (regularity of P_0.)
        pk[N] = x if variant == 'tw' else S.mul(x, twi[bN]); prec[N] = K - mu
    return pk, prec


def c_formula(pk, Q, p):
    return [(pk[5 * q][0] + (pk[5 * q - 1][1] if q else 0)) % p for q in range(Q + 1)]


def rebuild(S, pk, prec, b, bt, NK):
    """rows 0..5 of P as dicts k -> (series, valid length); row 5 must equal row 0"""
    ib = [S.inv(x) for x in b]                                        # ib[1] = 1/(u b_1)
    R = [{k: (pk[k], prec[k]) for k in pk}]
    for a in range(5):
        new = {}
        for k in range(NK + 1):
            col = (a + 1 - k) % 5; num = S.z(); L = 10 ** 9
            if k - 1 in R[a]:
                d, Ld = R[a][k - 1]; d = S.deriv(d); Ld -= 1
                if a == 1: d = S.shift(d, 1); Ld += 1
                num = d; L = Ld
            X, LX = R[a][k]
            if (col - 1) % 5 == 1:
                if a == 1: t = S.mul(bt[1], X)
                else: t = S.mul(bt[1], S.shift(X, -1)); LX -= 1
            else:
                t = S.mul(bt[(col - 1) % 5], X)
                if a == 1: t = S.shift(t, 1); LX += 1
            num = S.add(num, t); L = min(L, LX)
            new[k] = (S.mul(num, ib[a]), L)
        R.append(new)
    return R


def checks_R(R, NK, Q):
    close = all(R[5][k][0][:min(R[5][k][1], R[0][k][1])] == R[0][k][0][:min(R[5][k][1], R[0][k][1])] for k in range(NK + 1))
    r1 = True; r2 = True
    for a in range(5):
        for k in range(NK + 1):
            col = (a - k) % 5; x = R[a][k][0]
            if (a == 2 and col != 2) or (col == 1 and a != 1): r1 = r1 and x[0] == 0
            if a == 2 and col == 1: r2 = r2 and x[0] == 0 and x[1] == 0
    r1 = r1 and all(R[2][5 * q][0][0] == R[1][5 * q][0][0] for q in range(Q + 1))
    return close, r1, r2, [R[1][5 * q][0][0] for q in range(Q + 1)]


def theorem_G(S, pk, prec, YC, ys, NK, NU):
    """Sum_b P_0b Y^C_b == (Y^C_Q)_0 at u^0..u^NU, all complete z-orders; returns (#identities, failures)"""
    nmax = max(n for d in YC for n, y in d.items() if any(y[:NU + 1]))
    nid = 0; bad = []; ms = []; nz = 0
    for m in range(-nmax, NK + 1):
        if m % 5 != 4: continue                                     # other z-orders vanish on both sides by the Z/5 grading
        lhs = S.z(); complete = True
        for bb in range(5):
            for n, y in YC[bb].items():
                k = m + n
                if k < 0 or k % 5 != (-bb) % 5: continue
                if k > NK or prec[k] < NU + 1: complete = False; break
                lhs = S.add(lhs, S.mul(pk[k], y))
            if not complete: break
        if not complete: continue
        rhs = ys[-m][0] if (m <= 0 and -m < len(ys)) else S.z(); ms.append(m)
        for N in range(NU + 1):
            nid += 1; nz += rhs[N] != 0
            if lhs[N] != rhs[N]: bad.append((m, N))
    return nid, bad, (min(ms), max(ms), nz)


def worker_pact(path):
    vals = {}
    if not os.path.exists(path): return vals
    for line in open(path, encoding='utf-8'):
        mm = re.search(r"q=(\d+)\s+\[z\^\d+\] c\^act = \S+\s+p\^act_\(5q,0\)\(0\) = (-?\d+)(?:/(\d+))?", line)
        if mm: vals[int(mm.group(1))] = Fr(int(mm.group(2)), int(mm.group(3) or 1))
    return vals


def run(p, Q, NU, out):
    NK = 5 * Q; K = (4 * NK) // 5 + NU + 8
    S = Ser(p, K); polys, term = hatted_polys(p, NK); F = model_rows(p, NK)
    assert all(term)
    fds = [frame_data(S, al, be) for al, be in FRAMES]
    wp = worker_pact(os.path.join(HERE, "..", "logs", "pp_run_norm_Q14.log")); ok = True
    for fi, fd in enumerate(fds):
        tag = "frame (%s,%s)" % FRAMES[fi]
        pk, prec = extract(S, polys, fd['C0'], fd['tau'], F, NK)
        cF = c_formula(pk, Q, p)
        R = rebuild(S, pk, prec, fd['b'], fd['bt'], NK); close, r1, r2, cR = checks_R(R, NK, Q)
        YC = model_vanishing(S, fd['v'], NU); ys = quintic_vanishing(S, fd['b'], NU)
        nid, bad, mr = theorem_G(S, pk, prec, YC, ys, NK, NU)
        okF = cF == [1] + [0] * Q; okR = cR == cF
        line = ("p=%d %s: (F) c^act = 1 through z^%d: %s; (R) closing %s, (R1) %s, (R2) %s, c from entry (1,1) = (F): %s; "
                "(G) %d identities (%d non-zero) u^0..u^%d, z^%d..z^%d (z = 4 mod 5), failures %d"
                % (p, tag, 5 * Q, okF, close, r1, r2, okR, nid, mr[2], NU, mr[0], mr[1], len(bad)))
        out.append(line); ok = ok and okF and okR and close and r1 and r2 and not bad and nid > 0
        if fi == 0:
            nw = [q for q in range(Q + 1) if q in wp]
            same = all(pk[5 * q][0] == md(wp[q], p) for q in nw)
            out.append("   p^act_(5q,0)(0) = worker's logged values for q in %s..%s (%d values): %s; p_(5,0)(0) = 43397/807975000: %s"
                       % (nw[0], nw[-1], len(nw), same, pk[5][0] == md(Fr(43397, 807975000), p)))
            ok = ok and same
            # control 1: stored normalisation P_00(0,z) = 1, i.e. P^st = P^act k(z), k = 1/P^act_00(0,z)
            P00 = [pk[5 * q][0] for q in range(Q + 1)]; kk = [1] + [0] * Q
            for n in range(1, Q + 1): kk[n] = -sum(P00[j] * kk[n - j] for j in range(1, n + 1)) % p
            pst = {n: [0] * K for n in pk}
            for n in pk:
                for q in range(n // 5 + 1):
                    if kk[q]: pst[n] = S.add(pst[n], pk[n - 5 * q], kk[q])
            cst = c_formula(pst, Q, p); nid2, bad2, _ = theorem_G(S, pst, prec, YC, ys, NK, NU)
            out.append("   CONTROL stored normalisation: [z^5] c = m5 = -43397/807975000: %s; c != 1: %s; (G) failures %d of %d (first %s)"
                       % (cst[1] == md(M5, p), cst != [1] + [0] * Q, len(bad2), nid2, bad2[:1]))
            ok = ok and cst[1] == md(M5, p) and len(bad2) > 0
            for var in ('tw', 'mu'):
                pv, prv = extract(S, polys, fd['C0'], fd['tau'], F, NK, var)
                nv, bv, _ = theorem_G(S, pv, prv, YC, ys, NK, NU)
                try:
                    Rv = rebuild(S, pv, prv, fd['b'], fd['bt'], NK); cl = checks_R(Rv, NK, Q)[0]; msg = "closing %s" % cl
                except AssertionError as e:
                    cl = False; msg = "rebuild assertion (%s)" % e
                out.append("   CONTROL variant '%s': %s, (G) failures %d of %d -> %s"
                           % (var, msg, len(bv), nv, "fails as it should" if (not cl and bv) else "UNEXPECTED"))
                ok = ok and not cl and len(bv) > 0
    # control 3: C0 of frame (1/3,-2/7) with the system of frame (0,0)
    try:
        pv, prv = extract(S, polys, fds[1]['C0'], fds[0]['tau'], F, NK)
        Rv = rebuild(S, pv, prv, fds[0]['b'], fds[0]['bt'], NK); cl = checks_R(Rv, NK, Q)[0]
        msg = "closing %s" % cl; failed = not cl
    except AssertionError as e:
        msg = "assertion (%s)" % e; failed = True
    out.append("   CONTROL mismatched frames: %s -> %s" % (msg, "fails as it should" if failed else "UNEXPECTED"))
    return ok and failed


if __name__ == "__main__":
    Q = int(sys.argv[1]) if len(sys.argv) > 1 else 4; NU = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    t0 = time.time(); out = []
    print("rf_gauge: c^act through z^%d, Theorem-G window u^0..u^%d, QRR polynomials from rf_qexp (large radius)" % (5 * Q, NU))
    ok = True
    for p in (P61, P89):
        ok = run(p, Q, NU, out) and ok
    for l in out: print(l)
    print("RESULT rf_gauge: %s  [%.1f s]" % ("PASS" if ok else "FAIL", time.time() - t0))
