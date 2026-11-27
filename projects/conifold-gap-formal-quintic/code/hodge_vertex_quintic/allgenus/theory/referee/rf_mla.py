# -*- coding: utf-8 -*-
"""rf_mla.py -- referee re-derivation of the ML_Lambda test (REDUCTION_DRAFT.md s.2.4) with own frames (general v, homogeneity checked),
own list enumeration and contraction, the referee's own M_0 from the LG recursion fed into the project's MVertex6, at a THIRD prime
(p3 = the least prime = 1 mod 10 above 2^80, not used by the worker), plus statistics the worker did not report:
  - genuine cancellations: per singular list, the number of non-zero flat terms of its multilinear expansion (0 = trivially zero);
  - control C1 (worker's): frames -> leading monomials;  C3: ONE correction coefficient of one frame doubled (element leaves T_v);
  - C4: frames built from a WRONG vanishing series (sign of the (-1/(4z))^i factor flipped);  C5: frames built with Y^C(z) in place
    of Y^C(-z) in the pairing.  Every control should produce non-zero singular correlators (the test can fail).
Run: py -I -B rf_mla.py regress | h NT      Log: rf_mla_regress.log, rf_mla_g<h>_n<NT>.log (appended).
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_mla_core import Corr, frames_upto, frame, singular_lists, w, lg_rows
from rf_flat_core import yc_monomials


def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


def p3():
    n = (1 << 80) + 1
    while not (n % 10 == 1 and is_prime(n)): n += 1
    return n


def regress(C, log):
    g = lg_rows(8)
    tab0 = [Fr(1), Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000)]
    tab4 = [Fr(1), Fr(-3, 20), Fr(21, 160), Fr(-5229, 80000), Fr(-1411641, 6400000), Fr(138080943, 128000000)]
    log("  LG rows g^0, g^4 (n <= 5) = LG_MODEL.md table: %s, %s" % (g[0][:6] == tab0, g[4][:6] == tab4))
    H = {(0, 1): (Fr(1), Fr(0))}
    one = {(0, 0): (Fr(1), Fr(0))}; t1 = {(1, 0): (Fr(1), Fr(0))}
    fr = lambda k, a: frame(k, a)
    checks = [("F_2 = <>_2", C.flat(2, []), Fr(1, 1200)), ("F_3 = <>_3", C.flat(3, []), Fr(1, 25200)),
              ("<HHH>_0", C.flat(0, [(0, 1)] * 3), Fr(-5)), ("<H>_1", C.flat(1, [(0, 1)]), Fr(-1, 12)),
              ("<tau_1(1)>_1", C.flat(1, [(1, 0)]), Fr(5, 24)), ("<HH>_1", C.flat(1, [(0, 1)] * 2), Fr(1, 12)),
              ("<H>_2 = F_2'", C.flat(2, [(0, 1)]), Fr(-2, 1200)), ("<frame(1,0)>_1", C.multi(1, [fr(1, 0)])[0], Fr(1, 8)),
              ("<frame(3,3)>_2", C.multi(2, [fr(3, 3)])[0], Fr(3, 3200)), ("<frame(4,2)>_2", C.multi(2, [fr(4, 2)])[0], Fr(-41, 9600)),
              ("<frame(1,0)>_2", C.multi(2, [fr(1, 0)])[0], Fr(0)), ("<tau_1(H) tau_0(1)>_2 = <H>_2", C.flat(2, [(1, 1), (0, 0)]), Fr(-2, 1200))]
    ok = True
    for lab, val, exp in checks:
        good = (val - C.R.of(exp)) % C.p == 0; ok = ok and good
        log("  %-34s = %-12s expected %-10s %s" % (lab, C.fmt(val), exp, "ok" if good else "MISMATCH"))
    log("  frames at general v: f(1,0) = %s; f(2,0) = %s; f(4,2) = %s; f(5,3) = %s" % (fr(1, 0), fr(2, 0), fr(4, 2), fr(5, 3)))
    return ok


def run(C, h, NT, log):
    t0 = time.time(); wmax = 2 * h - 2 + NT + 2
    FR = frames_upto(wmax)
    Ybad = [(b, j, n, cf * (-1) ** (n - j)) for (b, j, n, cf) in yc_monomials(60)]
    FR4 = frames_upto(wmax, Y=Ybad); FR5 = frames_upto(wmax, sign=False)
    st = dict(sing=0, zero=0, fail=[], skip=0, triv=0, canc=0, maxnz=0, dimtriv=0, flatzero=0, c1=0, c3=0, c3n=0, c4=0, c4n=0, c5=0, reg=0, regnz=0)
    tlast = time.time()
    for combo, m, e in singular_lists(h, NT, FR):
        if time.time() - tlast > 300:
            tlast = time.time()
            log("     [progress %.0f s] singular tested %d, vanish %d, fail %d, skipped %d, genuine %d, C1 %d, C3 %d/%d, C5 %d/%d"
                % (tlast - t0, st['sing'], st['zero'], len(st['fail']), st['skip'], st['canc'], st['c1'], st['c3'], st['c3n'], st['c5'], st.get('c5n', 0)))
        vecs = [FR[c] for c in combo] + [{(0, 1): (Fr(1), Fr(0))}] * m
        try: val, nzt, nallow, nnontriv = C.multi(h, vecs, stats=True)
        except ValueError: st['skip'] += 1; continue
        if e >= 0:
            st['reg'] += 1; st['regnz'] += bool(val); continue
        st['sing'] += 1
        if val: st['fail'].append((combo, m, e, C.fmt(val)))
        else: st['zero'] += 1
        if nzt == 0:
            st['triv'] += 1
            if nnontriv == 0: st['dimtriv'] += 1
            else: st['flatzero'] += 1
        else: st['canc'] += 1; st['maxnz'] = max(st['maxnz'], nzt)
        lead = [{c: (Fr(1), Fr(0))} for c in combo] + [{(0, 1): (Fr(1), Fr(0))}] * m
        st['c1'] += bool(C.multi(h, lead)[0])
        # C3: double the lowest-z correction coefficient of the first frame that has corrections
        for i, c in enumerate(combo):
            corr = sorted(kc for kc in FR[c] if kc[1] == 1)
            if corr:
                v3 = dict(FR[c]); x, ee = v3[corr[0]]; v3[corr[0]] = (2 * x, ee)
                st['c3n'] += 1; st['c3'] += bool(C.multi(h, vecs[:i] + [v3] + vecs[i + 1:])[0]); break
        if any(FR4[c] != FR[c] for c in combo):
            st['c4n'] += 1; st['c4'] += bool(C.multi(h, [FR4[c] for c in combo] + vecs[len(combo):])[0])
        if any(FR5[c] != FR[c] for c in combo):
            st['c5n'] = st.get('c5n', 0) + 1; st['c5'] += bool(C.multi(h, [FR5[c] for c in combo] + vecs[len(combo):])[0])
    log("  genus %d, <= %d insertions (frames of weight <= %d): %d singular lists tested, %d vanish, %d FAIL, %d skipped (deficit > 6)"
        % (h, NT, wmax, st['sing'], st['zero'], len(st['fail']), st['skip']))
    log("     genuine cancellations: %d lists with >= 1 non-zero flat term (max %d terms); %d with every flat term 0, of which %d "
        "trivially (each flat term killed by dimension or by the string equation alone) and %d with a dimensionally allowed flat "
        "correlator that vanishes" % (st['canc'], st['maxnz'], st['triv'], st['dimtriv'], st['flatzero']))
    log("     controls (number of singular lists NON-zero): C1 leading monomials %d; C3 one coefficient doubled %d of %d; "
        "C4 wrong series %d of %d; C5 Y^C(z) pairing %d of %d;  regular lists (e = 0, 1): %d of %d non-zero"
        % (st['c1'], st['c3'], st['c3n'], st['c4'], st['c4n'], st['c5'], st.get('c5n', 0), st['regnz'], st['reg']))
    if st['fail']: log("     FAILURES (first 5): %s" % st['fail'][:5])
    log("     %s (%.0f s)" % ("PASS" if not st['fail'] and st['sing'] else "FAIL", time.time() - t0))
    return st


def main():
    tag = "regress" if sys.argv[1] == "regress" else "g%s_n%s" % (sys.argv[1], sys.argv[2])
    LOG = open(os.path.join(HERE, "rf_mla_%s.log" % tag), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    p = p3(); C = Corr(p)
    log("=" * 110); log("rf_mla %s, p3 = %d (= 1 mod 10) [%s]" % (" ".join(sys.argv[1:]), p, time.strftime("%Y-%m-%dT%H:%M:%S")))
    if sys.argv[1] == "regress":
        log("RF_MLA REGRESS: %s" % ("ALL PASS" if regress(C, log) else "FAIL"))
    else:
        h = int(sys.argv[1]); NT = int(sys.argv[2])
        st = run(C, h, NT, log)
        log("RF_MLA genus %d: %s" % (h, "PASS" if not st['fail'] else "FAIL"))
    LOG.close()


if __name__ == "__main__":
    main()
