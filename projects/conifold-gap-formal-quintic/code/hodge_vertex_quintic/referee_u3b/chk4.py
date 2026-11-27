# -*- coding: utf-8 -*-
"""chk4.py -- referee NOTE: is the consistency of the rank-6 N_2 system structural?
At the solution points with N_2-rank 6 the whole row grad q_17(Delta) = 2 A_17 Delta vanishes (15 circulant coordinates), so N_2 exists iff
[u^3] G_17(N_1, N_2 = 0) = 0 (chain graph g1-g1-g1).  Test: graph 17 ALONE, true u-dependent vertices (fq3_par.job_hat, one process), at
random circulant Delta in ker A_17 that are NOT solutions of anything else; controls outside ker A_17; and the stored solution points.
Run: py -I chk4.py --prime 61|89      (nothing is written to the project folder)"""
import sys, os, json, random, time, argparse
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, PROJ)


def kernel(M, n, p):
    M = [list(r) for r in M]; piv = []; rk = 0
    for c in range(n):
        pr = None
        for r in range(rk, len(M)):
            if M[r][c] % p: pr = r; break
        if pr is None: continue
        M[rk], M[pr] = M[pr], M[rk]; inv = pow(M[rk][c], p - 2, p); M[rk] = [x * inv % p for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] % p:
                f = M[r][c]; M[r] = [(x - f * y) % p for x, y in zip(M[r], M[rk])]
        piv.append(c); rk += 1
    ker = []
    for c in range(n):
        if c in piv: continue
        v = [0] * n; v[c] = 1
        for r, pc in enumerate(piv): v[pc] = (-M[r][c]) % p
        ker.append(v)
    return rk, ker


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--prime", default="61"); ap.add_argument("--nn", type=int, default=12); a = ap.parse_args()
    from fq3_circsys import CircSys
    from fq3_steps import frame_of
    from fq3_par import job_hat
    from fq3_m import sp_mats
    S = CircSys(a.prime); p = S.p; rng = random.Random(4242); t0 = time.time()
    A17 = S.Q[17]; rk, ker = kernel(A17, 15, p)
    print("prime 2^%s-1: rank of the 15 x 15 matrix of q_17 on the circulant class: %d (kernel dimension %d); m_1 vanishes on the kernel: %s"
          % (a.prime, rk, len(ker), all(sum(S.m1[j] * v[j] for j in range(15)) % p == 0 for v in ker)), flush=True)
    rkK, kerK = kernel(S.L2 + S.J3, 15, p)
    AK = [[sum(kerK[x][i] * A17[i][j] * kerK[y][j] for i in range(15) for j in range(15)) % p for y in range(len(kerK))] for x in range(len(kerK))]
    print("   rank of q_17 restricted to K_4^circ (dimension %d): %d" % (len(kerK), kernel(AK, len(kerK), p)[0]))
    def rnd(B):
        co = [rng.randrange(p) for _ in B]
        return [sum(c * b[j] for c, b in zip(co, B)) % p for j in range(15)]
    cs = json.load(open(os.path.join(PROJ, "g3_cubic_%s.json" % a.prime)))
    AB2 = [(x + 2 * y) % p for x, y in zip(cs['Pi_A'], cs['Pi_B'])]
    pts = [("random Delta in ker A_17 (point 1)", rnd(ker)), ("random Delta in ker A_17 (point 2)", rnd(ker)), ("random Delta in ker A_17 (point 3)", rnd(ker)),
           ("control: random Delta with m_1 = 0 on K_4^circ", None), ("control: random circulant Delta", [rng.randrange(p) for _ in range(15)]),
           ("Delta* (running circstar pass)", cs['circ_star']), ("A + 2B (circpi)", AB2), ("genuine c1_0 = 0 point (circb)", cs['circ_b_solutions'][0]),
           ("extra[3] (b = -zeta^2, not Pi)", cs['circ_b_extra'][3]), ("extra[0] (b = +zeta^2)", cs['circ_b_extra'][0]), ("B of Pi", cs['Pi_B'])]
    # control on K_4^circ with m_1 = 0: impose m_1 on the kernel basis
    mk = [sum(S.m1[j] * v[j] for j in range(15)) % p for v in kerK]; j0 = max(i for i in range(len(kerK)) if mk[i])
    co = [rng.randrange(p) for _ in kerK]; co[j0] = (-sum(mk[i] * co[i] for i in range(len(kerK)) if i != j0)) * pow(mk[j0], p - 2, p) % p
    pts[3] = (pts[3][0], [sum(co[i] * kerK[i][j] for i in range(len(kerK))) % p for j in range(15)])
    frame = frame_of(p, (0, 0), 34, 4); out = {}
    for name, v in pts:
        N1 = sp_mats(S.to_sp(v), p); g = S.quad_grad(17, v)
        idx, coefs, info = job_hat((frame, ('M0exp', [N1], 'right'), 17, a.nn))
        res = dict(grad_zero=not any(g), u01_zero=(coefs[0] == 0 and coefs[1] == 0), u2_is_q17=(coefs[2] == S.quad(17, v)), u2_zero=(coefs[2] == 0), u3_zero=(coefs[3] == 0),
                   maxdeg=info['maxdeg'])
        out[name] = res
        print("   %-48s grad q_17 = 0: %-5s [u^0] = [u^1] = 0: %s; [u^2] = q_17(Delta): %s (zero: %-5s); [u^3] G_17(N_1, 0) = 0: %-5s (eps-degree %s, nodes %d, %.0f s)"
              % (name, res['grad_zero'], res['u01_zero'], res['u2_is_q17'], res['u2_zero'], res['u3_zero'], info['maxdeg'], a.nn, time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "chk4_%s.json" % a.prime), "w"), indent=1)
    print("CHK4 done (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
