# -*- coding: utf-8 -*-
"""fq_mreport.py -- validation blocks of genus2_symplectic.py (steps 1-3), with the reference values of
gsum_scratch/pergraph_polar_out.txt (frames F5, F4), FORMAL_QUINTIC_genus2.md s.4 (total) and SCALAR_SCAN_genus2.md s.3.
Each block takes (log, check) callables and returns data for the cross-prime comparison."""
import itertools, random, time, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_rows import edge_coeffs
from fq_hodge import wk
from fq_phi import PhiVertexTable, phi_frame
from gsum_q import to_sseries, genus2_sums, ser_fmt
from fq_sser import SympM, rows_prime, EdgeLeaf
from fq_mvertex import MVertex, tensor_get
from fq_msum import GraphSummer, sum_series

F5_REF = {1: (Fr(1807, 3072000), Fr(-7, 8000), Fr(67, 288000)), 2: (Fr(-811, 256000), Fr(763, 192000), Fr(-31, 480000)),
          3: (Fr(-57, 64000), Fr(89, 20000), Fr(-609, 100000)), 4: (Fr(843, 1024000), Fr(-53, 48000), Fr(1, 32000)),
          5: (Fr(3, 51200), Fr(-67, 192000), Fr(371, 480000)), 6: (Fr(87, 32000), Fr(-143, 40000), Fr(149, 100000)),
          7: (Fr(9, 12800), Fr(-59, 16000), Fr(199, 40000))}
TOTAL_REF = {-10: Fr(1, 1200), -5: Fr(-7, 6000), 0: Fr(97, 72000), 5: Fr(-6809, 18000000), 10: Fr(-2551, 15000000)}
F4_REF = {1: (Fr(1807, 3072000), Fr(-37, 122880), Fr(-7, 8000), Fr(-119, 221184), Fr(0), Fr(67, 288000)),
          2: (Fr(-811, 256000), Fr(343, 153600), Fr(763, 192000), Fr(89, 55296), Fr(-19, 3840), Fr(-31, 480000)),
          3: (Fr(-57, 64000), Fr(-29, 12800), Fr(89, 20000), Fr(-7, 4608), Fr(1, 150), Fr(-609, 100000)),
          4: (Fr(843, 1024000), Fr(-11, 122880), Fr(-53, 48000), Fr(455, 221184), Fr(1, 2880), Fr(1, 32000)),
          5: (Fr(3, 51200), Fr(-71, 30720), Fr(-67, 192000), Fr(-245, 55296), Fr(373, 57600), Fr(371, 480000)),
          6: (Fr(87, 32000), Fr(1, 3840), Fr(-143, 40000), Fr(1, 2304), Fr(-1, 800), Fr(149, 100000)),
          7: (Fr(9, 12800), Fr(19, 7680), Fr(-59, 16000), Fr(11, 4608), Fr(-7, 960), Fr(199, 40000))}
F4_ORD = (-10, -6, -5, -2, -1, 0)
TYPES = [(2, ()), (2, (2,)), (2, (3,)), (2, (4,)), (2, (2, 2)), (2, (2, 3)), (2, (2, 2, 2)), (1, (0,)), (1, (1,)), (1, (0, 0)), (1, (0, 1)),
         (1, (1, 1)), (1, (0, 2)), (1, (0, 0, 0)), (1, (0, 0, 2)), (1, (0, 0, 2, 2)), (0, (0, 0, 0)), (0, (0, 0, 0, 0)), (0, (0, 0, 0, 1)),
         (0, (0, 0, 0, 0, 0)), (0, (0, 0, 0, 2))]


def step1(S, rows, graphs, log, check):
    R = S.R; p = R.p; ring = rows.ring
    log("   x = Cf(C0)/L = %s" % S.x.fmt(20))
    check("Delta_i^{1/2} = zeta^{4i} x squares to Cf(Delta_i) for all colours", S.delta_ok)
    vt = S.val_table(4); cexp = (0, -1, 3, 2, 1)
    check("val Q^{(i)}_{a,n} = -4n + (0,-1,3,2,1)_a - 1, colour-independent, n <= 4", all(vt[(a, n)] == [-4 * n + cexp[a] - 1] for a in range(5) for n in range(5)), vt)
    check("symplecticity Sum_ab Q_a^{(i)}(z) eta^{ab} Q_b^{(j)}(-z) = delta_ij to z^4 in s-series", not S.symplectic_check(4))
    I = SympM.identity(R); EL = EdgeLeaf(S, rows_prime(S, I), I); t0 = time.time(); bad = 0
    for i in range(5):
        for j in range(5):
            for (m, n), el in edge_coeffs(rows, i, j, 2).items():
                ref = to_sseries(el, ring, S.imu).mul(S.xinv).mul(S.xinv).scale(R.inv(S.zh[i] * S.zh[j] % p)); mine = EL.E(m, n)[i][j]
                if not ref.eq(mine, min(ref.prec, mine.prec) - 1): bad += 1
    check("M = 1: E'_{ij}(m,n) = Cf(fq_rows.edge_coeffs)/(Delta_i Delta_j)^{1/2} for all colours, m, n <= 2 (%.1f s)" % (time.time() - t0), bad == 0, bad)
    bad = 0
    for i in range(5):
        for d in (2, 3, 4):
            ref = Ser(R, 0, [])
            for k, v in rows.Dil(i, d).items(): ref = ref.add(Ser.mono(R, -k, v, 60))
            ref = ref.mul(S.xinv).scale(R.inv(S.zh[i])); mine = EL.T(d)[i]
            if not ref.eq(mine, min(ref.prec, mine.prec) - 1): bad += 1
    check("M = 1: T'^{(i)}_d = Dil_i(d)/Delta_i^{1/2} (fq_rows.Dil), d = 2..4", bad == 0, bad)
    GS = GraphSummer(S, I); t0 = time.time(); res = GS.all_graphs(graphs); log("   M = 1 graph sums in %.1f s" % (time.time() - t0))
    for idx in sorted(res): log("   graph %d: %s" % (idx, ser_fmt(res[idx], R, kmax=0)))
    check("M = 1: colour sums of the 7 graphs = point frame F5 (pergraph_polar_out.txt) at s^-10, s^-5, s^0; no other terms through s^0",
          all(res[i].coef(k) == R.of(F5_REF[i][j]) for i in res for j, k in enumerate((-10, -5, 0))) and
          all(res[i].coef(k) == 0 for i in res for k in range(-10, 1) if k not in (-10, -5, 0)))
    T = sum_series(list(res.values())); log("   total: %s" % ser_fmt(T, R, kmax=10))
    check("M = 1: total = (1/1200) s^-10 - (7/6000) s^-5 + 97/72000 - (6809/18000000) s^5 - (2551/15000000) s^10, u-series through s^10 (prec O(s^%d))" % T.prec,
          T.prec >= 11 and all(T.coef(k) == R.of(TOTAL_REF.get(k, 0)) for k in range(-10, 11)))
    pc = [GS.graph(graphs[0], colmask=(c,)) for c in range(5)]
    check("M = 1: vertex graph per colour = (1807/15360000) s^-10 - (7/40000) s^-5 + 67/1440000 + ...",
          all(s.coef(-10) == R.of(Fr(1807, 15360000)) and s.coef(-5) == R.of(Fr(-7, 40000)) and s.coef(0) == R.of(Fr(67, 1440000)) for s in pc))
    s40 = GS.graph(graphs[3], colmask=(0, 0)); s41 = GS.graph(graphs[3], colmask=(0, 1)); w = (R.z(2) + R.z(3)) % p
    check("M = 1: graph 4 at (0,0) = -31/15360000 s^-10 + 1/240000 s^-5; at (0,1): 1/40000 - w/30000, -1/30000 + 11w/240000 (SCALAR_SCAN table)",
          s40.coef(-10) == R.of(Fr(-31, 15360000)) and s40.coef(-5) == R.of(Fr(1, 240000)) and
          s41.coef(-10) == (R.of(Fr(1, 40000)) - w * R.inv(30000)) % p and s41.coef(-5) == (R.of(Fr(-1, 30000)) + 11 * w * R.inv(240000)) % p)
    return res, T


def step2(R, log, check, rng):
    p = R.p; t0 = time.time()
    MV = MVertex(R, SympM.identity(R)); bad = 0
    for (g, d) in TYPES:
        T = MV.tensor(g, d); n = len(d)
        for idx in itertools.product(range(5), repeat=n):
            exp = 0
            if sum(d) == 3 * g - 3 + n and len(set(idx)) <= 1: exp = R.of(wk(g, d)) * (R.z(4 * idx[0] * (2 * g - 2 + n)) if n else 1) % p
            if n == 0 and sum(d) == 3 * g - 3: exp = 0            # Sum_j zeta^{8j} wk = 0 (phase sum)
            if tensor_get(T, idx) != exp: bad += 1
    check("Omega^1 = omega: diagonal zeta^{4i(2g-2+n)} WK(g; d) for %d (g, d) types" % len(TYPES), bad == 0, bad)
    pts = [(Fr(-1, 12), Fr(1, 360), "F4 point"), (Fr(-3, 20), Fr(-7, 2500), "true class"), (rng.randrange(p), rng.randrange(p), "random (a1, a2)")]
    for a1, a2, lab in pts:
        a1 = R.of(a1); a2 = R.of(a2)
        MV = MVertex(R, SympM.diagonal(R, [a1 * R.z(-i) % p for i in range(5)], [a2 * R.z(-3 * i) % p for i in range(5)])); bad = 0
        for (g, d) in TYPES:
            T = MV.tensor(g, d); n = len(d)
            for idx in itertools.product(range(5), repeat=n):
                exp = 0
                if n == 0: exp = sum(PhiVertexTable(R, a1, a2, i).val(g, d) * R.z(4 * i * (2 * g - 2)) for i in range(5)) % p
                elif len(set(idx)) <= 1: exp = PhiVertexTable(R, a1, a2, idx[0]).val(g, d) * R.z(4 * idx[0] * (2 * g - 2 + n)) % p
                if tensor_get(T, idx) != exp: bad += 1
        check("Omega^M for M = diag(exp(a1 zeta^-i z + a2 zeta^-3i z^3)) at %s = zeta^{4i(2g-2+n)} PhiVertexTable (fq_phi), monochrome" % lab, bad == 0, bad)
    a1s = [rng.randrange(p) for i in range(5)]; a2s = [rng.randrange(p) for i in range(5)]
    MV = MVertex(R, SympM.diagonal(R, a1s, a2s)); bad = 0
    for (g, d) in TYPES:
        n = len(d)
        if n == 0: continue
        T = MV.tensor(g, d)
        for i in range(5):
            if tensor_get(T, (i,) * n) != PhiVertexTable(R, a1s[i], a2s[i], 0).val(g, d) * R.z(4 * i * (2 * g - 2 + n)) % p: bad += 1
    check("Omega^M for a non-covariant diagonal M (independent a1_i, a2_i) = per-colour PhiVertexTable", bad == 0, bad)
    M = SympM.random(R, rng); MV = MVertex(R, M); bad = 0
    for (g, d) in [(2, (2, 3)), (1, (0, 1)), (1, (0, 0, 2)), (0, (0, 0, 0, 1)), (2, (2, 2, 3))]:
        n = len(d); T = MV.tensor(g, d)
        for perm in itertools.permutations(range(n)):
            T2 = MV.tensor(g, tuple(d[perm[l]] for l in range(n)))
            for idx in itertools.product(range(5), repeat=n):
                if tensor_get(T, idx) != tensor_get(T2, tuple(idx[perm[l]] for l in range(n))): bad += 1
    check("random full M: Omega^M symmetric under simultaneous permutation of (colour, psi-exponent) legs", bad == 0, bad)
    unit = [R.z(-4 * k) for k in range(5)]; bad = 0
    for (g, d) in [(1, (1,)), (2, (2,)), (2, (3,)), (1, (1, 1)), (1, (0, 2)), (0, (0, 0, 1)), (2, (2, 3)), (1, (0, 0, 2)), (0, (0, 0, 0, 0))]:
        n = len(d); T = MV.tensor(g, d); Ts = MV.tensor(g, d + (0,)); Td = MV.tensor(g, d + (1,))
        for idx in itertools.product(range(5), repeat=n):
            s = sum(unit[k] * tensor_get(Ts, idx + (k,)) for k in range(5)) % p; exp = 0
            for l in range(n):
                if d[l] >= 1:
                    dl = list(d); dl[l] -= 1; exp = (exp + tensor_get(MV.tensor(g, tuple(dl)), idx)) % p
            if s != exp: bad += 1
            if sum(unit[k] * tensor_get(Td, idx + (k,)) for k in range(5)) % p != (2 * g - 2 + n) * tensor_get(T, idx) % p: bad += 1
    check("random full M: string and dilaton equations of Omega^M with the unit 1 = Sum_k Delta_k^{-1/2} e~_k (flat unit; 9 types)", bad == 0, bad)
    log("   step 2 checks in %.1f s" % (time.time() - t0))


def nring_graph_series(frame, deco, ring, imu, R):
    """colour-summed per-graph Cf series of a gsum_q Frame (independent N-ring route)"""
    Ssum = genus2_sums(frame, deco); out = {}
    for (idx, col), el in Ssum.items():
        s = to_sseries(el, ring, imu); out[idx] = s if idx not in out else out[idx].add(s)
    return out


def step3(S, rows, deco, graphs, log, check, rng):
    R = S.R; p = R.p; ring = rows.ring; data = {}
    a1 = R.of(Fr(-1, 12)); a2 = R.of(Fr(1, 360))
    M4 = SympM.diagonal(R, [a1 * R.z(-i) % p for i in range(5)], [a2 * R.z(-3 * i) % p for i in range(5)])
    res4 = GraphSummer(S, M4).all_graphs(graphs)
    check("M = diag(phi_i) at (-1/12, 1/360) covariant: colour sums of the 7 graphs = frame F4 at s^-10, s^-6, s^-5, s^-2, s^-1, s^0 (pergraph_polar_out.txt)",
          all(res4[i].coef(k) == R.of(F4_REF[i][j]) for i in res4 for j, k in enumerate(F4_ORD)) and all(res4[i].coef(k) == 0 for i in res4 for k in range(-10, 1) if k not in F4_ORD))
    T4 = sum_series(list(res4.values()))
    # random covariant diagonal point vs fq_phi.phi_frame through gsum_q (N ring)
    ra = (rng.randrange(p), rng.randrange(p)); t0 = time.time()
    Mr = SympM.diagonal(R, [ra[0] * R.z(-i) % p for i in range(5)], [ra[1] * R.z(-3 * i) % p for i in range(5)])
    resr = GraphSummer(S, Mr).all_graphs(graphs)
    ref = nring_graph_series(phi_frame(rows, ra[0], ra[1], 'phirand'), deco, ring, S.imu, R)
    check("random covariant diagonal (a1, a2) = (%d, %d): colour sums per graph = fq_phi.phi_frame via gsum_q (N ring), through s^10 (%.1f s)" % (ra[0], ra[1], time.time() - t0),
          all(resr[i].eq(ref[i], 11) for i in resr))
    # random non-covariant diagonal vs fq_diag_lib.diag_frame
    from fq_diag_lib import diag_frame
    pt = tuple(rng.randrange(p) for _ in range(10)); t0 = time.time()
    Md = SympM.diagonal(R, list(pt[:5]), list(pt[5:]))
    resd = GraphSummer(S, Md).all_graphs(graphs)
    refd = nring_graph_series(diag_frame(rows, pt, 'diagrand'), deco, ring, S.imu, R)
    check("random non-covariant diagonal (10 independent a1_i, a2_i): colour sums per graph = fq_diag_lib.diag_frame via gsum_q, through s^10 (%.1f s)" % (time.time() - t0),
          all(resd[i].eq(refd[i], 11) for i in resd))
    data['T4'] = T4
    return data
