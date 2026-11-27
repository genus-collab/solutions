# -*- coding: utf-8 -*-
"""fq_wsearch.py -- search for a weight multiset v realising prescribed power sums p_k = Sum_a v_a^{-k} (k = 1, 3, 5) over F_p
(genus2_symbol step 1 (b): the Hodge/Gamma hypothesis  b_k = -B_{2k}/(2k(2k-1)) p_{2k-1}(1/v)  for the symbol's odd log coefficients).

Candidates: (1) the named multisets: the tangent+twist weights wt_i = {zeta^i - zeta^j} U {-5 zeta^i}, their negatives, all
sub-multisets, and their multiples by 1, 5, 1/5 (and -1); (2) a pool search: Galois-orbit blocks (so that the power sums are
rational, as they must be when the data are rational) B = {c zeta^{mi}(1 - zeta^{mj})}_m (size 4), {c zeta^m}_m (size 4,
power sums -c^{-k}), {c (1 - zeta^m)}_m (size 4) and rational singletons {c}, for c in a list of small rationals; all multisets of
blocks of total size <= 6 are enumerated by a meet-in-the-middle over (p_1, p_3, p_5) triples (block sums are additive).
Elements with v_a = 0 are excluded.  Everything is mod p; a match is exact mod p (both primes in the driver).
API: named_candidates(R) -> [(name, [elements])]; power_sums(R, v, ks); search_blocks(R, target {1:,3:,5:}, max_size=6) -> matches.
"""
import itertools
from fractions import Fraction as Fr


def power_sums(R, v, ks=(1, 3, 5)):
    p = R.p
    return tuple(sum(R.pow(R.inv(x % p), k) for x in v) % p for k in ks)


def wt_weights(R, i):
    lams = [R.z(k) for k in range(5)]
    return [(lams[i] - lams[j]) % R.p for j in range(5) if j != i] + [(-5 * lams[i]) % R.p]


def named_candidates(R):
    p = R.p; out = []
    for i in range(5):
        w = wt_weights(R, i)
        for cname, c in (("1", 1), ("-1", p - 1), ("5", 5), ("-5", p - 5), ("1/5", R.inv(5)), ("-1/5", (-R.inv(5)) % p)):
            full = [x * c % p for x in w]; n = len(full)          # 5 weights: 4 tangent + 1 twist
            for r in range(1, n + 1):
                for sub in itertools.combinations(range(n), r):
                    out.append(("%s*wt_%d%s" % (cname, i, "" if r == n else " sub%s" % (sub,)), [full[k] for k in sub]))
    return out


def blocks(R):
    """Galois-orbit blocks: (name, elements) with rational power sums"""
    p = R.p; z = R.z; out = []
    rats = [Fr(1), Fr(-1), Fr(2), Fr(-2), Fr(3), Fr(-3), Fr(4), Fr(-4), Fr(5), Fr(-5), Fr(1, 2), Fr(-1, 2), Fr(1, 5), Fr(-1, 5),
            Fr(2, 5), Fr(-2, 5), Fr(5, 2), Fr(-5, 2), Fr(1, 3), Fr(-1, 3), Fr(3, 5), Fr(-3, 5), Fr(1, 4), Fr(-1, 4), Fr(1, 10), Fr(-1, 10),
            Fr(10), Fr(-10), Fr(6), Fr(-6), Fr(1, 6), Fr(-1, 6), Fr(5, 3), Fr(-5, 3), Fr(3, 2), Fr(-3, 2), Fr(2, 3), Fr(-2, 3)]
    for c in rats:
        cc = R.of(c)
        out.append(("{%s}" % c, [cc]))
        out.append(("{%s zeta^m}" % c, [cc * z(m) % p for m in range(1, 5)]))
        for j in range(1, 5):
            for i in range(5):
                orb = [cc * z(m * i) % p * (1 - z(m * j)) % p for m in range(1, 5)]
                out.append(("{%s zeta^(m*%d)(1 - zeta^(m*%d))}" % (c, i, j), orb))
    # dedupe by sorted element tuple
    seen = {}; res = []
    for name, els in out:
        key = tuple(sorted(els))
        if key in seen: continue
        seen[key] = name; res.append((name, els))
    return res


def search_blocks(R, target, max_size=6, ks=(1, 3, 5)):
    """multisets of blocks of total size <= max_size with power sums = target (tuple over ks).  Blocks have size 1 (rationals)
    or 4 (Galois orbits); with max_size <= 7 at most one orbit block fits.  The singleton part is found by meet in the middle
    (dictionary of all singleton multisets of size <= ceil(remaining/2))."""
    p = R.p; B = blocks(R)
    info = [(len(els), power_sums(R, els, ks), name) for name, els in B]
    singles = [b for b in range(len(info)) if info[b][0] == 1]; orbits = [b for b in range(len(info)) if info[b][0] > 1]
    def single_combos(limit):
        out = {}
        def rec(start, size, sums, chosen):
            out.setdefault(sums, []).append(tuple(chosen))
            if size == limit: return
            for b in range(start, len(singles)):
                ns = tuple((x + y) % p for x, y in zip(sums, info[singles[b]][1]))
                rec(b, size + 1, ns, chosen + [singles[b]])
        rec(0, 0, tuple(0 for _ in ks), [])
        return out
    tgt = tuple(x % p for x in target); matches = set()
    def find_singles(need, limit):
        """all singleton multisets of size <= limit with sums = need"""
        h1 = single_combos((limit + 1) // 2); h2 = single_combos(limit // 2) if limit // 2 else {tuple(0 for _ in ks): [()]}
        res = set()
        for sums1, lst1 in h1.items():
            rest = tuple((t - x) % p for t, x in zip(need, sums1))
            for c2 in h2.get(rest, []):
                for c1 in lst1:
                    if len(c1) + len(c2) <= limit: res.add(tuple(sorted(c1 + c2)))
        return res
    for s in find_singles(tgt, max_size): matches.add(s)
    for b in orbits:
        if info[b][0] > max_size: continue
        need = tuple((t - x) % p for t, x in zip(tgt, info[b][1]))
        for s in find_singles(need, max_size - info[b][0]): matches.add(tuple(sorted((b,) + s)))
    out = []
    for m in sorted(matches, key=lambda t: (sum(info[b][0] for b in t), t)):
        out.append((sum(info[b][0] for b in m), [info[b][2] for b in m]))
    return out, len(B)
