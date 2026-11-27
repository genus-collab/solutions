# -*- coding: utf-8 -*-
"""rf_frob_run.py -- REFEREE driver for rf_frob.py: the regular gauge from regularity alone, by linear algebra over Q.
 (1) homogeneous system (no normalisation): the solution space must be the c(z^5) directions only (Theorem U uniqueness, numerically);
 (2) canonical normalisation C_11(z) = 1: C = P^can(0, z) determined slot by slot; compared with the PROJECT gauge ur_gauge (read-only,
     via rf_ur) normalised by its own P_11(0, z);
 (3) the higher u-orders P^(n) obtained from C by the u-recursion, compared with ur_gauge at every complete slot.
Run: py -I -B rf_frob_run.py KC NMAX al be      log: rf_frob_run_<KC>_<frame>.log
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_core import systems
from rf_frob import recursion, solve
from rf_ur import stored_gauge, p11_zero, normalise


def main():
    KC = int(sys.argv[1]); NMAX = int(sys.argv[2]); al = Fr(sys.argv[3]); be = Fr(sys.argv[4])
    tag = "%d_ab_%s_%s" % (KC, str(al).replace("/", "o").replace("-", "m"), str(be).replace("/", "o").replace("-", "m"))
    OUT = open(os.path.join(HERE, "rf_frob_run_%s.log" % tag), "w", encoding="ascii"); FAIL = []
    def log(s): print(s, flush=True); OUT.write(s + "\n"); OUT.flush()
    def check(name, ok, det=""):
        log(("PASS " if ok else "FAIL ") + name + ((": " + str(det)) if det != "" else ""))
        if not ok: FAIL.append(name)
    t0 = time.time()
    log("rf_frob_run: frame (%s, %s), unknown C(z) = P(0, z) mod z^%d, u-recursion to n = %d (own periods, own recursion)" % (al, be, KC + 1, NMAX))
    S = systems(NMAX + 4, al, be)
    idx, lst, P, reach = recursion(S['b'], S['bt'], KC, NMAX)
    log("   %d unknowns; recursion done in %.0f s" % (len(lst), time.time() - t0))
    rows = []; graded = True; per_n = []
    for n in range(1, NMAX + 1):
        cnt = 0
        for (a, c, e), form in P[n].items():
            if (e - (a - c)) % 5: graded = False
            if e < 0 and e + reach[n][(a, c)] <= KC and form:
                rows.append((form, 0)); cnt += 1
        per_n.append(cnt)
    check("grading of every P^(n) slot (e = a - b mod 5)", graded)
    log("   complete regularity equations per u-order n = 1..%d: %s (total %d); reach_n (max over entries): %s"
        % (NMAX, per_n, len(rows), [max(reach[n].values()) for n in range(1, NMAX + 1)]))
    norm_slots = [idx[(1, 1, k)] for k in range(0, KC + 1, 5)]
    order = sorted([i for i in range(len(lst)) if i not in norm_slots], key=lambda i: (lst[i][2], i)) + norm_slots
    val, bad = solve(rows, len(lst), order)
    free = [i for i in range(len(lst)) if i not in val]
    log("   (1) homogeneous system: rank %d of %d unknowns, inconsistent rows %d; free unknowns (a, b, k): %s"
        % (len(val), len(lst), bad, [lst[i] for i in free]))
    kfree = min([lst[i][2] for i in free if i not in norm_slots] or [KC + 1])
    check("(1) below z^%d the only free unknowns are the normalisation slots C_11 at z^(5q) (solution space = C c(z^5))" % kfree,
          all(i in norm_slots for i in free if lst[i][2] < kfree) and all(i in free for i in norm_slots if lst[i][2] < kfree))
    nrows = rows + [({idx[(1, 1, k)]: Fr(1)}, 1 if k == 0 else 0) for k in range(0, KC + 1, 5)]
    val2, bad2 = solve(nrows, len(lst), order)
    det = {lst[i]: v for i, (v, vf) in val2.items() if not vf}
    kdet = min([k for (a, b, k) in lst if (a, b, k) not in det] or [KC + 1])
    log("   (2) canonical normalisation: %d of %d unknowns determined, all slots with k < %d determined; inconsistent rows %d"
        % (len(det), len(lst), kdet, bad2))
    check("(2a) the normalised system is consistent", bad2 == 0)
    p, prec = stored_gauge(max(KC, 10), max(3 * NMAX + 10, 40), al, be)
    pc, ci = normalise(p, p11_zero(p))
    mism = [(a, b, k) for (a, b, k), x in det.items() if k <= len(pc) - 1 and pc[k][b][0] != x]
    log("   C diagonal at z^5 (b = 0..4): %s" % [str(det.get(((b + 5) % 5, b, 5), "?")) for b in range(5)])
    log("   C_30 (pi) at z^3, z^8: %s, %s" % (det.get((3, 0, 3), "?"), det.get((3, 0, 8), "?")))
    check("(2b) C = P^can(0, z) of the project gauge (ur_gauge, normalised here) in all %d determined slots k < %d" % (len(det), kdet), not mism, mism[:5])
    # (3) higher u-orders from the determined C
    tested = 0; bad3 = []
    for n in range(1, NMAX + 1):
        if n >= prec: break
        for (a, c, e), form in P[n].items():
            if e < 0 or e + reach[n][(a, c)] > KC or e > len(pc) - 1: continue
            if any(lst[i] not in det for i in form): continue
            x = sum(det[lst[i]] * y for i, y in form.items())
            tested += 1
            if x != pc[e][c][n]: bad3.append((n, a, c, e))
    check("(3) P^(n)[z^e] from the u-recursion = ur_gauge canonical, %d complete slots (n <= %d, e >= 0)" % (tested, NMAX), not bad3, bad3[:5])
    log("RF_FROB_RUN %s: %s (%.0f s)" % (tag, "ALL PASS" if not FAIL else "FAIL %s" % FAIL, time.time() - t0)); OUT.close()


if __name__ == "__main__":
    main()
