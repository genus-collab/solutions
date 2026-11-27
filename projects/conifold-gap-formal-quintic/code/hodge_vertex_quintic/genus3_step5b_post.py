# -*- coding: utf-8 -*-
"""genus3_step5b_post.py -- post-processing of step 5b (two-sided covariant-scalar dressing): the second direction of the 2-dimensional O(u)
solution family, (b_5, b_5') -> (b_5 - t_2, b_5' + t_2), moves the central element zhat^5 . 1 from the right to the left factor and leaves M
unchanged, so every [u^2] quadric depends on t_1 only (checked: the t_2-dependent coefficients vanish); the resultant-in-t_2 test of
genus3_step5b.py is then vacuous (all resultants vanish identically) and the correct common-zero test is the univariate one: gcd of the
quadratics in t_1 over F_p, discriminants, roots.  Run:  py -I genus3_step5b_post.py --prime 89|61  (appends to the log, writes state keys
obs5b_uni_<graph>, obs5b_uni_gcd_degree).
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_fp import FpRing, P1, P2
from fq3_rec import Rec
from fq3_steps import State
import genus3_step5b as g5


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); a = ap.parse_args()
    p = P2 if a.prime == "89" else P1; R = FpRing(p); rec = Rec(R)
    state = State(os.path.join(HERE, "g3_state_%s.json" % a.prime))
    LOG = open(os.path.join(HERE, "logs", "genus3_ou_%s.log" % a.prime), "a", encoding="utf-8")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    def check(name, ok, detail=""): log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    log("-" * 110); log("(5b post) the two-sided family is effectively one-dimensional: univariate analysis of the [u^2] quadrics in t_1")
    ns = state.get('sys5b_ns')
    check("second direction of the O(u) family = (0, 0, -1, 0, 0, 1): the central zhat^5 . 1 moved between the factors (M unchanged)", ns[1] == [0, 0, p - 1, 0, 0, 1])
    Q = {i: state.get('obs5b_u2_%d' % i) for i in (3, 6, 7, 8, 17)}
    check("all t_2-dependent coefficients of the five quadrics vanish", all(c[k] == 0 for c in Q.values() for k, m in enumerate(g5.MONOS) if m[1] > 0))
    uni = {}
    for i, c in Q.items():
        a0, a1, a2 = c[0], c[1], c[3]; uni[i] = [a0, a1, a2]; state.put('obs5b_uni_%d' % i, uni[i])
        disc = (a1 * a1 - 4 * a2 * a0) % p; leg = pow(disc, (p - 1) // 2, p)
        roots = []
        if leg == 1:
            r = pow(disc, (p + 1) // 4, p); roots = [((-a1 + s * r) * pow(2 * a2, p - 2, p)) % p for s in (1, -1)]
        log("      graph %2d: Q(t_1) = %s + %s t_1 + %s t_1^2; discriminant a square mod p: %s; roots %s" % (i, rec.fmt(a0), rec.fmt(a1), rec.fmt(a2), leg == 1, [rec.fmt(x) for x in roots]))
    gg = None
    for i in (6, 7, 8, 17):
        h = g5.pgcd(uni[3], uni[i], p); log("      gcd(Q_3, Q_%d): degree %d" % (i, len(h) - 1)); gg = h if gg is None else g5.pgcd(gg, h, p)
    state.put('obs5b_uni_gcd_degree', len(gg) - 1)
    check("KEY (two-sided covariant-scalar class, genus 2 + 3 O(u) solutions): the five [u^2] quadrics in t_1 have NO common zero over the algebraic closure (every gcd(Q_3, Q_j) is a constant)", len(gg) == 1 and gg[0] != 0)
    LOG.close()


if __name__ == "__main__":
    main()
