# Fresh tree localisation on Mbar_{0,n}(P2,d), twisted by e^{-1}(R pi_* f^*N). Evidence only.
from fractions import Fraction as Fr
from math import factorial as fac
from itertools import combinations


class Th:
    def __init__(s, al, eps, tw=True):
        s.a = [Fr(x) for x in al]
        s.e = [Fr(x) for x in eps]
        s.tw = tw
        a = s.a
        s.T = [(a[(i+1) % 3]-a[i])*(a[(i+2) % 3]-a[i]) for i in range(3)]
        s.Wp = [s.T[i]*(3*a[i] if tw else 1) for i in range(3)]
        s.glob = 1/(s.e[0]*s.e[1]) if tw else Fr(1)
        s.etai = [s.glob/s.Wp[i] for i in range(3)]  # eta(e_i,e_i)

    def H(s, p):
        return lambda i: (-s.a[i])**p

    def edge(s, p, c, d):
        a = s.a
        w = (a[c]-a[p])/d
        den = Fr(1)
        for k in range(3):
            for m in range(d+1):
                if (k == p and m == 0) or (k == c and m == d):
                    continue
                den *= a[k]-a[p]-m*w
        num = Fr(1)
        if s.tw:
            for l in range(1, 3*d):
                num *= 3*a[p]+l*w
        return num/den/d


def subsets(l):
    l = list(l)
    for r in range(len(l)+1):
        for c in combinations(l, r):
            yield frozenset(c)


def trees(col, marks, deg):
    """rooted trees: (col,S,children), children sorted tuple of (d,tree)"""
    out = set()
    for S in subsets(marks):
        rest = frozenset(marks)-S
        for ch in kids(col, rest, deg):
            out.add((col, tuple(sorted(S)), ch))
    return out


def kids(col, marks, deg):
    if deg == 0:
        return {()} if not marks else set()
    out = set()
    for d in range(1, deg+1):
        for c2 in range(3):
            if c2 == col:
                continue
            for D in range(d, deg+1):
                for M in subsets(marks):
                    for t in trees(c2, M, D-d):
                        for r in kids(col, frozenset(marks)-M, deg-D):
                            out.add(tuple(sorted(((d, t),)+r, key=repr)))
    return out


def autf(tree):
    _, _, ch = tree
    f = 1
    seen = {}
    for x in ch:
        seen[x] = seen.get(x, 0)+1
        f *= autf(x[1])
    for v in seen.values():
        f *= fac(v)
    return f


def prod_fac(l):
    r = 1
    for x in l:
        r *= fac(x)
    return r


def vertex_int(n, slots):
    tot = n-3

    def rec(i, left, acc, exps):
        if i == len(slots)-1:
            return acc*slots[i](left)*Fr(fac(tot))/prod_fac(exps+[left])
        r = Fr(0)
        for k in range(left+1):
            c = slots[i](k)
            if c != 0:
                r += rec(i+1, left-k, acc*c, exps+[k])
        return r
    return rec(0, tot, Fr(1), [])


def ev(th, tree, marks, parent=None):
    col, S, ch = tree
    a = th.a
    flags = []
    for d, t in ch:
        flags.append((a[t[0]]-a[col])/d)
    if parent is not None:
        d, pc = parent
        flags.append((a[pc]-a[col])/d)
    val = len(flags)
    n = val+len(S)
    r = th.Wp[col]**(val-1)
    for m in S:
        r *= marks[m][0](col)
    if n >= 3:
        slots = [marks[m][1] for m in S]+[(lambda w: (lambda k: 1/w**(k+1)))(w) for w in flags]
        r *= vertex_int(n, slots)
    elif val == 1 and len(S) == 1:
        r *= marks[S[0]][2](-flags[0])
    elif val == 1 and len(S) == 0:
        r *= flags[0]
    elif val == 2 and len(S) == 0:
        r /= flags[0]+flags[1]
    else:
        raise ValueError
    for d, t in ch:
        r *= th.edge(col, t[0], d)*ev(th, t, marks, (d, col))
    return r


def corr(th, marks, d):
    """<marks>_{0,n,d}; marks[i]=(cl,g,h); rooted at mark 0"""
    n = len(marks)
    tot = Fr(0)
    for col in range(3):
        for S in subsets(range(1, n)):
            S0 = tuple(sorted({0} | S))
            for ch in kids(col, frozenset(range(1, n))-S, d):
                t = (col, S0, ch)
                val = len(ch)
                if val+len(S0) < 3 and d == 0:
                    continue
                tot += ev(th, t, marks)/autf(t)
    return tot*th.glob


def pm(cl, k=0):
    return (cl, (lambda j: Fr(1) if j == k else Fr(0)), (lambda x: x**k))


def main():
    AL = (Fr(2, 7), Fr(5, 11), Fr(-3, 13))
    th0 = Th(AL, (3, 7), tw=False)
    print('untwisted <pt,pt>_1, <pt^5>_2:', corr(th0, [pm(th0.H(2))]*2, 1), corr(th0, [pm(th0.H(2))]*5, 2))
    th = Th(AL, (3, 7))
    ee = th.e[0]*th.e[1]
    print('eps1eps2<H,H,H>_d/d^3, d=1,2,3:', [str(corr(th, [pm(th.H(1))]*3, d)*ee/d**3) for d in (1, 2, 3)])
    es = [(lambda i: (lambda j: Fr(1) if i == j else Fr(0)))(i) for i in range(3)]
    etainv = [1/th.etai[i] for i in range(3)]
    e1 = sum(th.a)
    e3 = th.a[0]*th.a[1]*th.a[2]
    print('eta_ab*3eps1eps2 vs eq:etaflat:',
          [[str(sum(th.H(a)(i)*th.H(b)(i)*th.etai[i] for i in range(3))*3*ee) for b in range(3)] for a in range(3)],
          str(1/e3), str(e1))
    H = [th.H(p) for p in range(3)]
    bad = 0
    cnt = 0
    for d in (0, 1, 2, 3):
        for k in (0, 1, 2):
            for ic, ia, ib in [(0, 1, 2), (1, 1, 1), (1, 2, 2), (2, 2, 2), (1, 1, 0), (2, 1, 1)]:
                L = corr(th, [pm(H[ib], k+1), pm(H[ic]), pm(H[ia])], d)
                R = Fr(0)
                for d1 in range(1, d+1):
                    d2 = d-d1
                    for i in range(3):
                        R += corr(th, [pm(H[ib], k), pm(es[i])], d1)*etainv[i]*corr(th, [pm(es[i]), pm(H[ic]), pm(H[ia])], d2)
                cnt += 1
                if L != R:
                    bad += 1
                    print('TRR3 FAIL', d, k, ic, ia, ib, L, R)
    print('TRR, n=3, d<=3:', cnt, 'cases, failures', bad)
    bad = 0
    cnt = 0
    for d in (0, 1, 2):
        for k in (0, 1, 2):
            for ic, ia, ib, it in [(1, 1, 2, 1), (2, 1, 1, 2), (1, 2, 2, 2), (0, 1, 1, 1)]:
                L = corr(th, [pm(H[ib], k+1), pm(H[ic]), pm(H[ia]), pm(H[it])], d)
                R = Fr(0)
                for d1 in range(d+1):
                    d2 = d-d1
                    for i in range(3):
                        x = corr(th, [pm(H[ib], k), pm(es[i]), pm(H[it])], d1)
                        y = corr(th, [pm(es[i]), pm(H[ic]), pm(H[ia])], d2)
                        R += x*etainv[i]*y
                        if d1 > 0:
                            x = corr(th, [pm(H[ib], k), pm(es[i])], d1)
                            y = corr(th, [pm(es[i]), pm(H[ic]), pm(H[ia]), pm(H[it])], d2)
                            R += x*etainv[i]*y
                cnt += 1
                if L != R:
                    bad += 1
                    print('TRR4 FAIL', d, k, ic, ia, ib, it, L, R)
    print('TRR, n=4 (one t), d<=2:', cnt, 'cases, failures', bad)
    L = corr(th, [pm(H[1], 1), pm(H[1]), pm(H[1])], 2)
    R = Fr(0)
    for d1, d2 in [(1, 1), (2, 0)]:
        for i in range(3):
            R += corr(th, [pm(H[1]), pm(es[i])], d1)*th.T[i]*corr(th, [pm(es[i]), pm(H[1]), pm(H[1])], d2)
    print('control, untwisted diagonal differs (must be True):', L != R)
    # control: equivariant shift of psi by a weight must break the relation
    L2 = L+Fr(1, 3)*corr(th, [pm(H[1]), pm(H[1]), pm(H[1])], 2)
    R = Fr(0)
    for d1, d2 in [(1, 1), (2, 0)]:
        for i in range(3):
            R += corr(th, [pm(H[1]), pm(es[i])], d1)*etainv[i]*corr(th, [pm(es[i]), pm(H[1]), pm(H[1])], d2)
    print('control, psi shifted by a weight differs (must be True):', L2 != R, ' unshifted equal:', L == R)
    bad = 0
    for d in (1, 2):
        for k in (0, 1, 2, 3):
            for ia, ib in [(0, 1), (1, 2), (2, 2), (2, 0)]:
                L = corr(th, [pm(H[ib], k), pm(H[ia]), pm(H[1])], d)
                R = d*corr(th, [pm(H[ib], k), pm(H[ia])], d)
                if k > 0:
                    R += corr(th, [pm((lambda ib: (lambda i: H[ib](i)*H[1](i)))(ib), k-1), pm(H[ia])], d)
                if L != R:
                    bad += 1
                L = corr(th, [pm(H[ib], k), pm(H[ia]), pm(H[0])], d)
                R = corr(th, [pm(H[ib], k-1), pm(H[ia])], d) if k > 0 else Fr(0)
                if L != R:
                    bad += 1
    print('descendant divisor and string, d<=2, failures', bad)
    ok = True
    for m in (1, 2, 3):
        L = corr(th, [pm(H[1], m-1), pm(H[2])]+[pm(H[1])]*m, 0)
        R = sum(H[1](i)*H[2](i)*(-th.a[i])**m*th.etai[i] for i in range(3))
        ok = ok and L == R
    print('d=0 two-point with H-insertions equal eta(a H^m,b):', ok)

    def mul(x, y):
        r = [Fr(0)]*4
        for p in range(4):
            for q in range(4-p):
                r[p+q] += x[p]*y[q]
        return r

    def expser(x):
        r = [Fr(1), Fr(0), Fr(0), Fr(0)]
        term = [Fr(1), Fr(0), Fr(0), Fr(0)]
        for n in range(1, 4):
            term = [c/n for c in mul(term, x)]
            r = [r[j]+term[j] for j in range(4)]
        return r
    for z in (Fr(13, 3), Fr(-7, 2)):
        for i in range(3):
            a = th.a
            Id = [Fr(1)]
            for d in range(1, 4):
                num = Fr(1)
                for m in range(3*d):
                    num *= 3*a[i]-m*z
                den = Fr(1)
                for l in range(3):
                    for m in range(1, d+1):
                        den *= a[l]-a[i]+m*z
                Id.append(num/den)
            g = [Fr(0)]+[3*Fr(fac(3*d-1), fac(d)**3)*(-1)**d for d in (1, 2, 3)]
            Qy = expser(g)
            Q = mul([Fr(0), Fr(1), Fr(0), Fr(0)], Qy)
            Qp = [[Fr(1), Fr(0), Fr(0), Fr(0)]]
            for d in range(1, 4):
                Qp.append(mul(Qp[-1], Q))
            g1 = (lambda k: z**(-k-2))
            h1 = (lambda x: 1/(z*(z-x)))
            ei = es[i]
            J = [Fr(1), Fr(0), Fr(0), Fr(0)]
            for d in range(1, 4):
                c = corr(th, [(ei, g1, h1)], d)*etainv[i]
                J = [J[j]+c*Qp[d][j] for j in range(4)]
            pref = expser([-a[i]*x/z for x in g])
            Rhs = mul(pref, J)
            print('J=I, z=', z, 'P', i, ':', [Id[j] == Rhs[j] for j in range(4)])


main()
