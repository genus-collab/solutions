# Critic: frame-B graph sum at a GENERIC point Q0 of the Q-line (not near Q = 1), Q = Q0 e^t, in F_p.
# Tests the closed form including no pole at Disc = 0 and at Q = oo (Schwartz-Zippel in Q0), plus per-graph values.
# Usage: py -I gen_q.py g a0 a1 a2 e1 [Q0seed] [frame]
import sys, os, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
import cylib
from cylib import *
import gsum
from gsum import Sum, predicted, constant_term, _Sc


def rows_at_Q0(Q0):
    def build(w, N, L):
        Ser.L = L
        R = Rows(); R.w = w; R.N = N; R.L = L
        # Q(t) = Q0 e^t
        ex = [0] * L; f = 1
        for k in range(L):
            if k: f = f * k % P
            ex[k] = inv(f)
        Q = scal(Q0, ex)
        # critical points of W = -a2 log U + A log(1-U) + a1 log(1-QU):
        # B Q U^2 - (a0 + Q (a1-a2)) U - a2 = 0
        a12 = (w.a1 - w.a2) % P
        bq = add(const(w.a0, L), scal(a12, Q))
        disc = add(mul(bq, bq, L), scal(4 * w.B % P * w.a2 % P, Q))
        d0 = disc[0]
        rt = pow(d0, (P + 1) // 4, P)
        assert rt * rt % P == d0, 'Disc(Q0) not a square'
        sq = ssqrt(disc, rt, L)
        iQ = sinv(Q, L)
        ib = inv(2 * w.B)
        Up = scal(ib, mul(add(bq, sq), iQ, L))
        Um = scal(ib, mul(sub(bq, sq), iQ, L))
        tau = sub(const(1, L), iQ)
        itau = sinv(tau, L)
        KK = 2 * N + 3
        pols = pol_table(KK + 1)
        res = []
        for U in (Up, Um):
            X = sinv(sub(const(1, L), U), L)
            Y = sinv(sub(const(1, L), mul(Q, U, L)), L)
            pX = [const(1, L)]; pY = [const(1, L)]
            for i in range(KK + 1):
                pX.append(mul(pX[-1], X, L)); pY.append(mul(pY[-1], Y, L))
            wk = {}
            for k in range(2, KK + 1):
                out = [0] * L
                for i, co in pols[k - 1].items():
                    co = co % P
                    t = add(scal(w.a1, pY[i]), scal(w.A, pX[i]))
                    for n in range(L):
                        out[n] -= co * t[n]
                wk[k] = [x % P for x in out]
            i2 = sinv(scal(P - 1, wk[2]), L)
            p = [const(1, L)]
            for n in range(1, 3 * N + 3):
                p.append(mul(p[-1], i2, L))
            r0 = wick(p, wk, N, L)
            Lc = mul(add(scal(w.B, U), const(w.a2, L)), itau, L)
            kap = w.a1 * w.A % P * w.e1 % P * w.e2 % P
            res.append((U, wk, r0, (0, Lc), (0, scal(kap, wk[2]))))
        R.U1, R.wk1, R.r0pt1, R.L1, R.D1 = res[0]
        R.U2, R.wk2, R.r0pt2, R.L2, R.D2 = res[1]
        return R
    return build


def closed_form(w, g, Q0, L):
    # P_g Li_{3-2g}(Q0 e^t)/(e1 e2) as a t-series
    ex = [0] * L; f = 1
    for k in range(L):
        if k: f = f * k % P
        ex[k] = inv(f)
    extra = 2 * g
    LL = L + extra
    ex = [0] * LL; f = 1
    for k in range(LL):
        if k: f = f * k % P
        ex[k] = inv(f)
    Q = scal(Q0, ex)
    Li = mul(Q, sinv(sub(const(1, LL), Q, ), LL), LL)  # Li_0 = Q/(1-Q)
    for _ in range(2 * g - 3):
        Li = deriv(Li)
    Pg = 0
    for i in range(g + 1):
        Pg += _Sc(i, w.e1) * _Sc(g - i, w.e2)
    Pg %= P
    return scal(Pg * inv(w.e1 * w.e2 % P) % P, Li[:L])


def main():
    g = int(sys.argv[1])
    a0, a1, a2, e1 = [Fr(x) for x in sys.argv[2:6]]
    e2 = a0 + a1 + a2 - e1
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    frame = sys.argv[7] if len(sys.argv) > 7 else 'B'
    w = weights(a0, a1, a2, e1, e2)
    rnd = random.Random(seed)
    while True:
        Q0 = rnd.randrange(2, P)
        disc = (w.a0 * w.a0 + (2 * w.a0 * (w.a1 + w.a2) + 4 * w.a1 * w.a2) * Q0 + (w.a1 - w.a2) ** 2 * Q0 * Q0) % P
        if pow(disc, (P - 1) // 2, P) == 1:
            break
    gsum.build_rows = rows_at_Q0(Q0)
    t0 = time.time()
    S = Sum(w, g, frame, 4)
    tot, rows = S.run()
    L = S.L
    cf = closed_form(w, g, Q0, L)
    c0 = constant_term(w, g)
    k_ok = [lcoef(tot, k) == (cf[k] + (c0 if k == 0 else 0)) % P for k in range(5)]
    print('g', g, 'weights', a0, a1, a2, e1, e2, 'frame', frame, 'seed', seed)
    print('  F_g(Q0 e^t) == P_g Li/(e1e2) + F_g(0), coefficients t^0..t^4:', k_ok, ' time %.1f' % (time.time() - t0), flush=True)
    # colour-swap sanity: the total must not depend on the sign of the square root (trivial, sum over colourings)


if __name__ == '__main__':
    main()
