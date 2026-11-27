# -*- coding: utf-8 -*-
"""cyc_ref_frame.py N -- referee cyc_formal: C5, the Theorem and the Corollary of CYC_PROOF.md in conifold frames, exact over Q.

Frames: y_0 = mu (f_0 + alpha g + beta f_2) (the parametrisation of fq_conifold, pb_core, GI, UR), y_1 = nu g, t = y_1/y_0;
(alpha, beta, mu, nu) = (0,0,1,1), (1/3,-2/7,1,1), (7/11,-3/5,1,1), (-5/2,4,-3,7/2); and the critic's frame y_0 = f_0 + f_1/3 - 2 f_2/7.
  F1  C5: adapted symplectic basis (y_0, y_1, y_2, y_3') by linear algebra over Q; its Bmat is the signed permutation matrix;
      tau_0 = tau_1 = tau_2 = 0, u y_0^2 Dt^3 tau_3 = 1; Y := th^2 f_2 = tau_3; th^2 f_3 = -t Y.
  F2  Theorem: M = (1/y_0) D (1/Dt) D (u y_0^2 Dt^2) D (1/Dt) D (1/y_0), expanded as an operator, equals PF4 coefficientwise;
      u y_0^2 Dt^3 th^2(y/y_0) = c_2 - c_3 t on the adapted basis (= 1 on y_2, -t on y_3').
  F3  Corollary: O := (1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D, c = (Dt, 1/(u y_0^2 Dt^2), Dt, y_0, y_0):
      O o (1/y_0) = (1/y_0) o D o PF4; prod c_a = 1/u.  In the project normalisation (mu = nu = 1): UR/GI's
      b = (t', -1/(u y_0^2 (1-u)^3 t'^2), t', -y_0/(1-u), -y_0/(1-u)) equals -c/(1-u); prod b = -1/(u(1-u)^5).
  Controls: t = f_1/y_0 and t = (f_1 + f_2/2)/y_0 (non-isotropic); middle factor with Dt^3; c_1 with Dt^3; c_2 <-> c_3;
  C5's span statement for t = f_1/y_0 (while tau_1 = tau_2 = 0 still hold there: they do not use isotropy).
Usage: py -I -B cyc_ref_frame.py N    (stored as logs/cyc_ref_frame_N.log)
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from cyc_ref_core import Ser, ONE, U
from cyc_ref_ops import Op, REl
from cyc_ref_jet import det_inv, matT_B_mat
from cyc_ref_frob import PF4, basis, bmat, frame, combo

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NF = [0]
IU = Ser(-1, [Fr(1)])
DOP = Op.Dk(1)


def check(name, ok, control=False):
    if control:
        print(('PASS control (fails as it must): ' if not ok else 'FAIL control did not fail: ') + name)
        NF[0] += bool(ok)
    else:
        print(('PASS ' if ok else 'FAIL ') + name)
        NF[0] += (not ok)


def Mop(y0, t, mid_pow=2):
    e = t.D(); ie = e.inv(); iy = y0.inv()
    mid = U * y0 * y0
    for _ in range(mid_pow):
        mid = mid * e
    return Op.fn(iy) * (DOP * (Op.fn(ie) * (DOP * (Op.fn(mid) * (DOP * (Op.fn(ie) * (DOP * Op.fn(iy))))))))


def Ochain(ic):                              # ic = (1/c_0, .., 1/c_4): (1/c_4) D (1/c_3) D (1/c_2) D (1/c_1) D (1/c_0) D
    O = Op.fn(ic[0]) * DOP
    for a in (1, 2, 3, 4):
        O = Op.fn(ic[a]) * (DOP * O)
    return O


def prec_of(op):
    return min(s.N for s in op.A.values())


def solve2(r0, r1, b0, b1):                  # x in Q^4 with r0.x = b0, r1.x = b1 (first non-singular 2x2 minor)
    for i in range(4):
        for j in range(i + 1, 4):
            d = r0[i] * r1[j] - r0[j] * r1[i]
            if d:
                x = [Fr(0)] * 4
                x[i] = (b0 * r1[j] - b1 * r0[j]) / d; x[j] = (r0[i] * b1 - r1[i] * b0) / d
                return x
    raise ValueError('degenerate')


def in_span_1_t(x, t, M):                    # x = a + b t below u^M ?  (t = c u + ..., c != 0)
    a = x.coef(0); b = x.coef(1) / t.coef(1)
    return (x - Ser(0, [a]) - t.sc(b)).zero_below(M)


def taus(fs, Bi, ie, kmax=3):                # tau_k = Sum Bi_ij f_i th^k f_j, th = (1/Dt) D
    out = []
    for k in range(kmax + 1):
        acc = None
        for i in range(4):
            for j in range(4):
                if Bi[i][j]:
                    x = fs[j]
                    for _ in range(k):
                        x = x.D().mulS(ie)
                    s = (fs[i] * x).sc(Fr(Bi[i][j])); acc = s if acc is None else acc + s
        out.append(acc)
    return out


t0 = time.time()
bas = basis(N)
B = bmat(bas['Y'])[0]
gvec = [0, 1, bas['gamma'], 0]
J = [[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]]
FR = [(0, 0, 1, 1, 'project'), (Fr(1, 3), Fr(-2, 7), 1, 1, 'project'), (Fr(7, 11), Fr(-3, 5), 1, 1, 'project'),
      (Fr(-5, 2), 4, -3, Fr(7, 2), 'project'), (Fr(1, 3), Fr(-2, 7), 1, 1, 'critic')]
for (al, be, mu, nu, par) in FR:
    al, be, mu, nu = Fr(al), Fr(be), Fr(mu), Fr(nu)
    tag = '%s frame (%s, %s), mu = %s, nu = %s' % (par, al, be, mu, nu)
    y0, y1, t = frame(bas, al, be, mu, nu, par)
    mid = gvec if par == 'project' else [0, 1, 0, 0]
    v0 = [mu * (int(i == 0) + al * mid[i] + be * int(i == 2)) for i in range(4)]
    v1 = [nu * x for x in gvec]
    Bv = lambda a, b: sum(a[i] * B[i][j] * b[j] for i in range(4) for j in range(4))
    r0 = [sum(v0[i] * B[i][j] for i in range(4)) for j in range(4)]
    r1 = [sum(v1[i] * B[i][j] for i in range(4)) for j in range(4)]
    v3 = solve2(r0, r1, 1, 0); v2 = solve2(r0, r1, 0, 1)
    c = Bv(v2, v3); v2 = [x - c * y for x, y in zip(v2, v0)]
    Bn = matT_B_mat([[v[i] for v in (v0, v1, v2, v3)] for i in range(4)], B)
    check('F1 %s: B(y_0, y_1) = 0; adapted basis with Bmat = signed permutation (Bmat^-1 = -Bmat)' % tag,
          Bv(v0, v1) == 0 and Bn == J)
    ys = [combo(bas, v) for v in (v0, v1, v2, v3)]
    check('F1 %s: the adapted basis starts with y_0 and y_1' % tag,
          (ys[0] - REl.of(y0)).zero_below(N) and (ys[1] - REl.of(y1)).zero_below(N))
    iy0 = y0.inv(); e = t.D(); ie = e.inv()
    fs = [y.mulS(iy0) for y in ys]
    yuk = U * y0 * y0 * e * e * e
    tk = taus(fs, [[-x for x in r] for r in J], ie)
    P = min(s.N for x in tk for s in x.F.values())
    check('F1 %s: tau_0 = tau_1 = tau_2 = 0 and u y_0^2 Dt^3 tau_3 = 1 below u^%d' % (tag, P - 1),
          all(tk[k].zero_below(P) for k in range(3)) and (tk[3].mulS(yuk) - REl.of(ONE)).zero_below(P - 1))
    th2 = lambda x: x.D().mulS(ie).D().mulS(ie)
    Yf = th2(fs[2])
    check('F1 %s: Y = th^2 f_2 = tau_3 and th^2 f_3 = -t Y' % tag,
          (Yf - tk[3]).zero_below(P - 1) and (th2(fs[3]) + Yf.mulS(t)).zero_below(P - 1))
    M = Mop(y0, t)
    p = prec_of(M)
    check('F2 %s: M = PF4 coefficientwise (D^0..D^4, leading u) below u^%d' % (tag, p),
          M.diff(PF4, p) and set(M.A) <= {0, 1, 2, 3, 4})
    im = [th2(x) for x in fs[2:]]
    check("F2 %s: u y_0^2 Dt^3 th^2(y/y_0) = c_2 - c_3 t on y_2, y_3' (= 1 and -t)" % tag,
          all(s.zero_below(s.N) for x in im for k, s in x.F.items() if k) and (yuk * im[0].F[0] - ONE).zero_below(P - 2)
          and (yuk * im[1].F[0] + t).zero_below(P - 2))
    ic = [ie, U * y0 * y0 * e * e, ie, iy0, iy0]
    pc = ic[0] * ic[1] * ic[2] * ic[3] * ic[4]
    O = Ochain(ic) * Op.fn(iy0); R = Op.fn(iy0) * (DOP * PF4)
    po = min(prec_of(O), prec_of(R))
    check('F3 %s: prod c_a = 1/u; O o (1/y_0) = (1/y_0) o D o PF4 below u^%d (closing relation <=> PF5)' % (tag, po),
          (pc - U).zero_below(pc.N) and O.diff(R, po))
    if mu == 1 and nu == 1:
        tp = t.D().sc(-1) * Ser(0, [Fr(1), Fr(-1)]).inv(N)          # t' = dt/du = -Dt/(1-u)
        i1u = Ser(0, [Fr(1), Fr(-1)]).inv(N)
        b = [tp, (IU * (y0 * y0 * tp * tp).inv()).sc(-1) * i1u * i1u * i1u, tp, (y0 * i1u).sc(-1), (y0 * i1u).sc(-1)]
        cv = [e, IU * (y0 * y0 * e * e).inv(), e, y0, y0]
        ok_b = all((b[a] + cv[a] * i1u).zero_below(min(b[a].N, N - 6)) for a in range(5))
        pb = b[0] * b[1] * b[2] * b[3] * b[4]
        target = (IU * i1u * i1u * i1u * i1u * i1u).sc(-1)
        check('F3 %s: UR/GI b_a = -c_a/(1-u) for a = 0..4; prod b = -1/(u (1-u)^5)' % tag,
              ok_b and (pb - target).zero_below(min(pb.N, target.N)))
    if (al, be, par) == (0, 0, 'project'):
        check('F3 frame (0,0): varpi_0 = 1 + O(u^3), t_c = u + 7/10 u^2 + 41/75 u^3 + O(u^4) (UR s.1, GI)',
              [y0.coef(k) for k in range(3)] == [1, 0, 0] and [t.coef(k) for k in range(4)] == [0, 1, Fr(7, 10), Fr(41, 75)])
# controls, project frame (1/3, -2/7)
y0, y1, t = frame(bas, Fr(1, 3), Fr(-2, 7))
f = bas['f']; iy0 = y0.inv()
for name, tt in (('t = f_1/y_0', f[1] * iy0), ('t = (f_1 + f_2/2)/y_0', (f[1] + f[2].sc(Fr(1, 2))) * iy0)):
    M = Mop(y0, tt)
    check('F2 control, frame (1/3,-2/7), %s (non-isotropic pair): M = PF4' % name, M.diff(PF4, prec_of(M)), control=True)
M = Mop(y0, t, mid_pow=3)
check('F2 control: middle factor u y_0^2 Dt^3 instead of Dt^2: M = PF4', M.diff(PF4, prec_of(M)), control=True)
e = t.D(); ie = e.inv()
R = Op.fn(iy0) * (DOP * PF4)
for name, ic in (('c_1 = 1/(u y_0^2 Dt^3)', [ie, U * y0 * y0 * e * e * e, ie, iy0, iy0]),
                 ('c_2 <-> c_3', [ie, U * y0 * y0 * e * e, iy0, ie, iy0])):
    O = Ochain(ic) * Op.fn(iy0)
    check('F3 control: %s: O o (1/y_0) = (1/y_0) D PF4' % name, O.diff(R, min(prec_of(O), prec_of(R))), control=True)
tt = f[1] * iy0; e = tt.D(); ie = e.inv()
det, Bi = det_inv(B)
fs = [y.mulS(iy0) for y in bas['Y']]
tg = taus(fs, Bi, ie, kmax=2)
P = min(s.N for x in tg for s in x.F.values())
check('F1 observation: for t = f_1/y_0 (non-isotropic) tau_0 = tau_1 = tau_2 = 0 still hold (Frobenius basis and Binv)',
      all(x.zero_below(P) for x in tg))
x2 = (U * y0 * y0 * e * e * e) * fs[2].D().mulS(ie).D().mulS(ie).F[0]
check('F1 control, t = f_1/y_0: u y_0^2 Dt^3 th^2(f_2/y_0) lies in span(1, t)', in_span_1_t(x2, tt, x2.N), control=True)
print('CYC_REF_FRAME %d: %s   (%.1f s)' % (N, 'ALL CHECKS PASSED' if NF[0] == 0 else '%d FAILURES' % NF[0], time.time() - t0))
