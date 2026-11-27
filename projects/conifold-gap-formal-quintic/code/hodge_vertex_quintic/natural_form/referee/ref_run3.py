# -*- coding: utf-8 -*-
"""ref_run3.py -- REFEREE item 3: attempts to BREAK the claim  Rt(zhat; u) = G(u, zhat) M_0(tau(u) zhat).

 (a) non-triviality: perturb one coefficient of tau, use tau = 1, use a 'neighbouring' natural function; count the nonzero
     Fourier components of [u^n zhat^k] G above the claimed degree (they must ALL vanish for the true tau).
 (b) M_0 on the LEFT:  G' = M_0(tau zhat)^(-1) Rt;  is there ANY rescaling c(u) with M_0(c zhat)^(-1) Rt polynomial at order u^1 ?
 (c) flat-colour form (left factor Psi(u) kept): Qh(zhat; u) M_0(tau zhat)^(-1) -- degrees (must equal those of G).
 (d) symplecticity of G: G(-zhat)^T G(zhat) = 1.
 (e) log G: degrees, number of nonzero Fourier components vs G, support rule W, the H (x) H = (3,3) component of the zhat-linear
     part at every u-order, all Fourier components through u^3.
Run: py -I ref_run3.py --prime 89 [--NZ 16] [--K 12]
"""
import os, sys, time, argparse
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ref_fp import P61, P89, P521, of, inv, fmt, mid, mzero, mmul, madd, mscale, miszero, bmul, bclean, bflip_T, blog, smul, sinv
from ref_hat import HatRows, tau_pred, M0_of, G_of, Mtau, degrees, fourier

PR = {'89': P89, '61': P61, '521': P521}
W = (0, 4, 3, 2, 1)


def ceil_div(a, b): return -((-a) // b)


def nfour(X, zeta, p, pred):
    return sum(len(fourier(A, zeta, p)) for key, A in X.items() if pred(key))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="89"); ap.add_argument("--NZ", type=int, default=16)
    ap.add_argument("--K", type=int, default=12); a = ap.parse_args(); p = PR[a.prime]; NZ = a.NZ; K = a.K
    H = HatRows(p, NZ, K); z = H.zeta; Rt = H.Rt; M0, M0i = M0_of(Rt, NZ, p); F = lambda x: fmt(x, p, z)
    tau = [of(x, p) for x in tau_pred(K)]; G = G_of(Rt, M0i, tau, p, NZ, K); deg = degrees(G, K)
    dclaim = [(5 * n + 2) // 4 for n in range(K)]
    print("=" * 100); print("ref_run3: p = 2^%s-1, NZ = %d, K = %d; checks %s" % (a.prime, NZ, K, H.checks))
    print("true tau: degrees %s; nonzero Fourier components above floor((5n+2)/4): %d" % (deg, nfour(G, z, p, lambda key: key[0] > dclaim[key[1]])))
    # ---------------------------------------------------------------- (a)
    print("(a) NON-TRIVIALITY")
    slots = sum(5 * (NZ - dclaim[n]) for n in range(1, K))
    print("    Z/5-allowed Fourier components above the claimed degree, u^1..u^%d, zhat <= %d: %d  (all vanish for the true tau)" % (K - 1, NZ, slots))
    for m0 in (1, 2, 5):
        t2 = list(tau); t2[m0] = (t2[m0] + 1) % p; G2 = G_of(Rt, M0i, t2, p, NZ, K)
        cnt = nfour(G2, z, p, lambda key: key[1] == m0 and key[0] > dclaim[m0])
        print("    tau_%d -> tau_%d + 1: degree of [u^%d] G = %s; nonzero Fourier components at u^%d above degree %d: %d of %d" %
              (m0, m0, m0, degrees(G2, K)[m0], m0, dclaim[m0], cnt, 5 * (NZ - dclaim[m0])))
    alts = {"tau = 1 (no rescaling)": [1] + [0] * (K - 1)}
    ser = lambda fr: [of(x, p) for x in fr]
    # (1-u)^(4/9): same tau_1 = -4/9
    c = [Fr(1)]
    for k in range(1, K): c.append(c[-1] * (Fr(4, 9) - (k - 1)) / k * (-1))
    alts["tau = (1-u)^(4/9) (same tau_1)"] = ser(c)
    # 1/2F1(1, 4/5; 9/5; u) with the wrong lower parameter: 1/Sum (4/5)_n/(14/5)_n u^n ... use 1/Sum 9 u^n/(5n+9)*... simple neighbour:
    alts["tau = 1/Sum_n 9 u^n/(5n+9)"] = sinv(ser([Fr(9, 5 * k + 9) for k in range(K)]), p)
    alts["tau = (u/t_c)^(4/5) (flat coordinate instead of the canonical one)"] = None
    tcu = [H.tc[k + 1] for k in range(K)]                              # t_c/u
    r = sinv(tcu, p); lg = [0] * K                                    # log(u/t_c)
    for n in range(1, K): lg[n] = (r[n] - sum(k * lg[k] * r[n - k] for k in range(1, n)) * inv(n, p)) % p
    lg = [x * of(Fr(4, 5), p) % p for x in lg]; ex = [1] + [0] * (K - 1)
    for n in range(1, K): ex[n] = sum(k * lg[k] * ex[n - k] for k in range(1, n + 1)) % p * inv(n, p) % p
    alts["tau = (u/t_c)^(4/5) (flat coordinate instead of the canonical one)"] = ex
    for name, t2 in alts.items():
        G2 = G_of(Rt, M0i, t2, p, NZ, K)
        print("    %-70s degrees of [u^n] G: %s" % (name, degrees(G2, K)))
    # ---------------------------------------------------------------- (b)
    print("(b) M_0 ON THE LEFT")
    Gl = bmul(Mtau(M0i, tau, p, NZ, K), Rt, p, NZ, K)
    print("    G' = M_0(tau zhat)^-1 Rt with the same tau: degrees %s" % degrees(Gl, K))
    # order u^1 with a free c_1: [u^1] G' = M0^-1 Rt_1 - c_1 M0^-1 zhat M0'  -> c_1 from each entry of the top two orders
    Rt1 = [Rt.get((n, 1), mzero()) for n in range(NZ + 1)]
    A = [mzero() for _ in range(NZ + 1)]; B = [mzero() for _ in range(NZ + 1)]
    for i in range(NZ + 1):
        for j in range(NZ + 1 - i):
            A[i + j] = madd(A[i + j], mmul(M0i[i], Rt1[j], p), p)
            if j: B[i + j] = madd(B[i + j], mmul(M0i[i], M0[j], p), p, j)
    for k in (NZ, NZ - 1):
        cands = set(A[k][i][j] * inv(B[k][i][j], p) % p for i in range(5) for j in range(5) if B[k][i][j])
        print("    order u^1, zhat^%d: distinct values of c_1 required by the 25 entries: %d  (1 would be needed)" % (k, len(cands)))
    # ---------------------------------------------------------------- (c)
    Qm = {}
    for n in range(NZ + 1):
        for k in range(K):
            Qm[(n, k)] = [[H.Q[aa][i][n][k] for i in range(5)] for aa in range(5)]
    Gf = bmul(bclean(Qm), Mtau(M0i, tau, p, NZ, K), p, NZ, K)
    print("(c) flat-colour rows Qh(zhat; u) M_0(tau zhat)^-1 (= Psi(u) G): degrees %s" % degrees(Gf, K))
    # ---------------------------------------------------------------- (d)
    S = bmul(bflip_T(G, p), G, p, NZ, K)
    print("(d) G(-zhat)^T G(zhat) = 1 mod (zhat^%d, u^%d): %s" % (NZ + 1, K, S == {(0, 0): mid()}))
    # ---------------------------------------------------------------- (e)
    t0 = time.time(); Lg = blog(G, p, NZ, K, K); dl = degrees(Lg, K)
    print("(e) log G (%.0f s): degrees %s" % (time.time() - t0, dl))
    print("    in sp (coefficient of zhat^k symmetric for k odd, antisymmetric for k even): %s" %
          all(A == mscale([list(r) for r in zip(*A)], 1 if k % 2 else p - 1, p) for (k, n), A in Lg.items()))
    for n in range(1, K):
        if dclaim[n] > NZ: continue
        print("    u^%-2d: nonzero Fourier components  G: %3d   log G: %3d   | (3,3) = H(x)H component of the zhat-linear part: G %s, log G %s" % (
            n, nfour(G, z, p, lambda key: key[1] == n), nfour(Lg, z, p, lambda key: key[1] == n),
            F(fourier(G.get((1, n), mzero()), z, p).get((3, 3), 0)), F(fourier(Lg.get((1, n), mzero()), z, p).get((3, 3), 0))))
    viol = 0; nonsharp = 0; tot = 0; nmin = {}
    for (k, n), Am in Lg.items():
        for (m, mp), v in fourier(Am, z, p).items():
            key = (((-m - 1) % 5, (-mp - 1) % 5), k); nmin[key] = min(n, nmin.get(key, K))
    for aa in range(5):
        for bb in range(5):
            for k in range(1, NZ + 1):
                if (aa + bb + 2 + k) % 5: continue
                b = max(1, ceil_div(4 * k - 2 + W[aa] + W[bb], 5))
                if b >= K: continue
                tot += 1; nm = nmin.get(((aa, bb), k))
                if nm is None or nm > b: nonsharp += 1
                elif nm < b: viol += 1
    print("    rule W for log G: violations %d, non-sharp %d, of %d (component, k) pairs" % (viol, nonsharp, tot))
    for n in range(1, 4):
        for k in range(1, dl[n] + 1):
            if (k, n) in Lg: print("    [u^%d zhat^%d] log G: %s" % (n, k, {key: F(v) for key, v in sorted(fourier(Lg[(k, n)], z, p).items())}))


if __name__ == "__main__":
    main()
