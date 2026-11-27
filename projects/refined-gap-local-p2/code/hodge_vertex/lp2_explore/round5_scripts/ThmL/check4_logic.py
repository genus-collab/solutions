# Check 4: mechanical run of the vanishing argument (s.6 of the note), with the uniform convention (U parts >= 1).
# For all G, D <= 2G-3, all T at level D and ALL parts k0 of T, R = T - k0, enumerates every (alpha, pi, assignment of R up to multiset, (g_B,k_B))
# with Sum e_B = 2G, and verifies:
#  (1) count Sum (delta_B - k_B) = D_R + k0 - 2G = D+1-2G
#  (2) some block has delta_B < k_B
#  (3) that block is killed by: falling factorial (0 <= delta < k), or g<G and D_U <= 2g-3 (H1'), or g=G and D_U<D (H2');
#      the only exception is the target (alpha=(k0), g=G, k=1, T_B=R), which occurs exactly once.
import sys
from tl_lib import *
from mvx import sub_multisets, multiset_minus

Gmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rmax = int(sys.argv[2]) if len(sys.argv) > 2 else 99

def assignments(Rt, blocks):
    """distribute multiset Rt among labelled blocks (multiset-level)."""
    if blocks == 1:
        yield [Rt]
        return
    for sub in sub_multisets(Rt):
        rest = multiset_minus(Rt, sub)
        rest = tuple(sorted(rest, reverse=True))
        for a in assignments(rest, blocks - 1):
            yield [sub] + a

total_terms = targets = survivors = viol = 0
cases = 0
from functools import lru_cache
for G in range(2, Gmax + 1):
    tG = 0
    for D in range(1, 2 * G - 2):
        for lam in partitions(D):
            T = tuple(p + 1 for p in lam)
            if len(T) > rmax: continue
            for k0 in sorted(set(T)):
                R = list(T); R.remove(k0); R = tuple(R)
                DR = sum(R) - len(R)
                assert DR == D - k0 + 1
                cases += 1
                ntgt = 0
                for alpha in partitions(k0):
                    al = list(alpha)
                    for pi in set_partitions(list(range(len(al)))):
                        m = len(pi)
                        Bs = [[al[i] for i in b] for b in pi]
                        for assign in assignments(R, m):
                            # enumerate (g,k) per block with total e = 2G
                            def rec(i, left, cur):
                                global total_terms, survivors, viol
                                nonlocal_ntgt = 0
                                if i == m:
                                    if left != 0: return 0
                                    total_terms_inc()
                                    s = 0
                                    killed = False
                                    tgt = False
                                    for (B, Tb, g, k) in cur:
                                        U = tuple(Tb) + tuple(B)
                                        DU = sum(U) - len(U)
                                        d = DU - 2 * g + 2
                                        s += d - k
                                    if s != D + 1 - 2 * G: viol_inc()
                                    has = False
                                    for (B, Tb, g, k) in cur:
                                        U = tuple(Tb) + tuple(B)
                                        DU = sum(U) - len(U)
                                        d = DU - 2 * g + 2
                                        if d < k:
                                            has = True
                                            if 0 <= d: killed = True
                                            elif g < G:
                                                assert g >= 2
                                                if DU <= 2 * g - 3: killed = True
                                            elif g == G:
                                                if DU < D: killed = True
                                                elif alpha == (k0,) and tuple(Tb) == R and k == 1 and len(B) == 1:
                                                    tgt = True
                                    if not has: viol_inc()
                                    if tgt and not killed:
                                        return 1
                                    if not killed:
                                        print("SURVIVOR", G, D, T, k0, alpha, cur)
                                        surv_inc()
                                    return 0
                                B = Bs[i]; Tb = assign[i]; sz = len(B)
                                cnt = 0
                                g = 0
                                while True:
                                    kmin = 1 if sz % 2 == 1 else 0
                                    e0 = 2 * g - 2 + sz + kmin
                                    if e0 > left: break
                                    for k in range(kmin, left - (2 * g - 2 + sz) + 1):
                                        e = 2 * g - 2 + sz + k
                                        assert e >= 0
                                        cnt += rec(i + 1, left - e, cur + [(B, Tb, g, k)])
                                    g += 1
                                return cnt
                            def total_terms_inc():
                                global total_terms
                                total_terms += 1
                            def viol_inc():
                                global viol
                                viol += 1
                            def surv_inc():
                                global survivors
                                survivors += 1
                            ntgt += rec(0, 2 * G, [])
                if ntgt != 1:
                    print("TARGET COUNT", G, D, T, k0, ntgt)
                    viol += 1
                targets += ntgt
    print("G =", G, "cumulative cases", cases, "terms", total_terms, "survivors", survivors, "violations", viol, flush=True)
print("cases", cases, "targets", targets, "terms", total_terms, "survivors", survivors, "violations", viol)
