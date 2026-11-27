# -*- coding: utf-8 -*-
"""th_mla.py -- direct test of the Model Lemma ML_Lambda (REDUCTION_DRAFT.md s.2) on the frozen model (v = 1), exact mod p.

ML_Lambda(h), base-point form: for n >= 1 insertions e_1..e_n in T_v (the frame vectors frame(k, a), a != 1, of th_core) and m >= 0
insertions of tau_0(H), the ancestor correlator <e_1 .. e_n tau_0(H)^m>_h(v) = C v^e,  e = 2 - 2h - m + Sum_i w(e_i),
w(frame(k, a)) = (4k + c_a - 4)/5, has C = 0 whenever e < 0.  (Correlators with non-integer e vanish by the Z/5 selection rule and
are not counted.)  At v = 1 the test is C = <..>_h(1) = 0.
Controls: (C1) the same insertion lists with every frame replaced by its leading monomial z^k H^a (not in T_v): the number of
non-zero singular correlators is printed (the H-corrections of the frames are what makes ML_Lambda non-trivial); (C2) regular
lists (e = 0, 1) with frames: the number of non-zero values is printed (the frames are not trivially annihilated).
Correlators: th_core.Model (MVertex6 on the project symbol M_0, psi-deficit <= 6; lists outside that range are SKIPPED and counted).
Run: py -I -B th_mla.py 89|61 [genera, e.g. 0,1,2,3] [max total insertions per genus, e.g. 5,5,4,3]
Log: th_mla_<prime>.log (appended).
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from th_core import Model, frame, weight5, P1, P2

CEXP = (0, -1, 3, 2, 1)


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"
    genera = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0, 1, 2, 3]
    ntots = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [5, 5, 4, 3]
    p = P2 if prime == "89" else P1
    LOG = open(os.path.join(HERE, "th_mla_%s.log" % prime), "a", encoding="ascii")

    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()

    t0 = time.time(); M = Model(p)
    log("=" * 110)
    log("th_mla: ML_Lambda on the frozen model, p = 2^%s-1, genera %s, max insertions %s  [%s]" % (prime, genera, ntots, time.strftime("%Y-%m-%dT%H:%M:%S")))
    H = {(0, 1): Fr(1)}
    allok = True
    for h, NT in zip(genera, ntots):
        kmax = 3 * h + 3
        A = [(k, a) for k in range(kmax + 1) for a in (0, 2, 3, 4)]
        FR = {ka: frame(*ka) for ka in A}
        LEAD = {ka: {ka: Fr(1)} for ka in A}
        st = dict(sing=0, sing_zero=0, sing_fail=0, skip=0, ctrl_nonzero=0, reg=0, reg_nonzero=0)
        fails = []; ctrl_ex = []; reg_ex = []
        th = time.time()
        for ntot in range(1, NT + 1):
            if 2 * h - 2 + ntot <= 0: continue
            for m in range(0, ntot):
                n = ntot - m
                for combo in itertools.combinations_with_replacement(A, n):
                    e5 = 5 * (2 - 2 * h) - 5 * m + sum(weight5(k, a) for k, a in combo)
                    if e5 % 5: continue
                    e = e5 // 5
                    if e > 1: continue
                    try:
                        val = M.mcorr(h, [FR[c] for c in combo] + [H] * m)
                    except ValueError:
                        st['skip'] += 1; continue
                    lab = " ".join("f(%d,%d)" % c for c in combo) + (" H^%d" % m if m else "")
                    if e < 0:
                        st['sing'] += 1
                        if val % p == 0: st['sing_zero'] += 1
                        else:
                            st['sing_fail'] += 1; fails.append((lab, e, M.fmt(val)))
                        try:
                            cv = M.mcorr(h, [LEAD[c] for c in combo] + [H] * m)
                            if cv % p:
                                st['ctrl_nonzero'] += 1
                                if len(ctrl_ex) < 4: ctrl_ex.append((lab.replace("f(", "lead("), e, M.fmt(cv)))
                        except ValueError:
                            pass
                    else:
                        st['reg'] += 1
                        if val % p:
                            st['reg_nonzero'] += 1
                            if len(reg_ex) < 3: reg_ex.append((lab, e, M.fmt(val)))
        ok = st['sing_fail'] == 0 and st['sing'] > 0
        allok = allok and ok
        log("-" * 110)
        log("genus %d (insertions <= %d, frames k <= %d): %d singular lists (e < 0) tested, %d vanish, %d FAIL; %d skipped (psi-deficit > 6)   (%.0f s)"
            % (h, NT, kmax, st['sing'], st['sing_zero'], st['sing_fail'], st['skip'], time.time() - th))
        log("   control C1 (frames -> leading monomials): %d of the singular lists are NON-zero, e.g. %s" % (st['ctrl_nonzero'], ctrl_ex))
        log("   control C2 (regular lists, e in {0, 1}): %d of %d non-zero, e.g. %s" % (st['reg_nonzero'], st['reg'], reg_ex))
        if fails: log("   FAILURES (first 10):", fails[:10])
        log(("PASS " if ok else "FAIL ") + "ML_Lambda(%d) on all tested singular lists" % h)
    log("TH_MLA %s: %s   (%.0f s)" % (prime, "ALL PASS" if allok else "FAILURES", time.time() - t0))
    LOG.close()


if __name__ == "__main__":
    main()
