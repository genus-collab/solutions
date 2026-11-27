# -*- coding: utf-8 -*-
"""fq3_par.py -- multiprocessing workers for the genus-3 production runs (genus3_ou.py): hatted sums per outer graph and frozen
directional derivatives per sp direction, with the per-process setup (rows, conifold frame, symbol, M_0) cached in module globals.

Every job is specified by picklable plain data:
  frame = (p, alpha_str, beta_str, prec_abs, K)         (Fractions as strings; K = u-precision of the hatted theory)
  mspec = ('const', mu)  constant M = exp(mu) (mu = 7 matrices as lists)  |  ('M0',)  the tautological M_0  |  ('id',)
        | ('M0exp', Ns, placement)  M = M_0 exp(Sum_k u^k N_k) (right) / exp(..) M_0 (left), Ns = list of sp series (7 matrices each)
        | ('M0exp2', Ls, Rs)  two-sided M = exp(Sum u^k L_k) M_0 exp(Sum u^k R_k)
        | ('M0exp_eps', N, eps)  constant M_0 exp(eps N) (right) / ('eps_M0', N, eps) exp(eps N) M_0 (left)   [frozen theory]
Jobs:
  job_hat((frame, mspec, idx, NN))      -> (idx, [u^0..u^{K-1}] coefficients of G_Gamma for the outer graph idx, info)
       constant M: MVertex6 directly; u-dependent: UVertex interpolation over NN eps-nodes (NN = 21 checks the degree bound 18).
  job_frozen((frame, mspec, idxs))      -> {idx: F_p} frozen sums (u = 0) of the listed outer graphs for a CONSTANT mspec.
  job_jac((frame, N, idxs, nodes, placement, tag)) -> (tag, {idx: [values at the eps-nodes]}) of the frozen sums with M_0 exp(eps N) (right) or
       exp(eps N) M_0 (left): the caller interpolates (fq_hat.vander_inv / interp_fp) to get the eps^1 coefficient D_N L_Gamma.
run_pool(func, jobs, nproc, log) runs the jobs with multiprocessing.Pool (imap_unordered) and returns the list of results.
"""
import os, sys, time, pickle
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)

_SETUP = {}


def setup(frame):
    key = tuple(frame)
    st = _SETUP.get(key)
    if st is not None: return st
    from fq_fp import FpRing
    from fq_rows import Rows
    from fq_conifold import Conifold
    from fq_sser import SRows
    from fq_symbol import Symbol
    from fq3_msum import HatRowsU
    from fq3_m import M0_from_symbol
    from fq_leadsolve import FrozenRows
    from fq3_graphs import outer_graphs3
    p, a, b, prec_abs, K = frame
    R = FpRing(p); rows = Rows(R, Nz=6); C = Conifold(R, Nu=10, alpha=Fr(a), beta=Fr(b)); S = SRows(rows, C, prec_abs=prec_abs)
    sym = Symbol(S); H = HatRowsU(S, K); M0, mu, chk = M0_from_symbol(sym, R); fr = FrozenRows(sym)
    graphs = {G['idx']: G for G in outer_graphs3()}
    st = dict(R=R, rows=rows, C=C, S=S, sym=sym, H=H, M0=M0, mu0=mu, M0checks=chk, frozen=fr, graphs=graphs, K=K)
    _SETUP[key] = st
    return st


def make_M(st, mspec):
    """returns (M object, kind) with kind 'const' (MConst) or 'u' (MU)"""
    from fq3_mvertex import MConst
    from fq3_m import MU, sp_scale
    from fq_symbol import mat_series_mul, mat_series_exp
    R = st['R']; p = R.p; M0 = st['M0']
    t = mspec[0]
    if t == 'id': return MConst.identity(R), 'const'
    if t == 'M0': return M0, 'const'
    if t == 'const': return MConst.from_log(R, mspec[1]), 'const'
    if t == 'M0exp': return MU(R, M0, mspec[1], st['K'], mspec[2]), 'u'
    if t == 'M0exp2':
        from fq3_m import MUprod
        return MUprod(R, M0, mspec[1], mspec[2], st['K']), 'u'
    if t in ('M0exp_eps', 'eps_M0'):
        N = mspec[1]; eps = mspec[2] % p
        E = mat_series_exp(sp_scale(N, eps, p), R, 6); Ei = mat_series_exp(sp_scale(N, (-eps) % p, p), R, 6)
        if t == 'M0exp_eps': M = mat_series_mul(M0.M, E, p, 6); Mi = mat_series_mul(Ei, M0.Minv, p, 6)
        else: M = mat_series_mul(E, M0.M, p, 6); Mi = mat_series_mul(M0.Minv, Ei, p, 6)
        return MConst(R, M, Mi, check=False), 'const'
    raise ValueError(mspec)


def job_hat(args):
    frame, mspec, idx, NN = args
    from fq3_mvertex import MVertex6
    from fq3_msum import GraphSummer6, UVertex
    t0 = time.time(); st = setup(frame); R = st['R']; G = st['graphs'][idx]
    M, kind = make_M(st, mspec)
    if kind == 'const': vt = MVertex6(R, M)
    else: vt = UVertex(R, M, st['K'], NN)
    s = GraphSummer6(st['H'], M, vt).graph(G)
    coefs = [s.coef(k) if k < s.prec else None for k in range(st['K'])]
    info = dict(time=time.time() - t0, prec=s.prec, maxdeg=getattr(vt, 'maxdeg', None), E=G['E'])
    return idx, coefs, info


def job_frozen(args):
    frame, mspec, idxs = args
    from fq3_msum import frozen_sums6
    st = setup(frame); M, kind = make_M(st, mspec); assert kind == 'const'
    return frozen_sums6(st['frozen'], M, [st['graphs'][i] for i in idxs])


def job_jac(args):
    frame, N, idxs, nodes, placement, tag = args
    from fq3_msum import frozen_sums6
    st = setup(frame); out = {i: [] for i in idxs}
    for eps in nodes:
        M, _ = make_M(st, ('M0exp_eps' if placement == 'right' else 'eps_M0', N, eps))
        res = frozen_sums6(st['frozen'], M, [st['graphs'][i] for i in idxs])
        for i in idxs: out[i].append(res[i])
    return tag, out


def run_pool(func, jobs, nproc, log=None, label=""):
    import multiprocessing as mp
    t0 = time.time(); out = []
    if nproc <= 1:
        for j in jobs: out.append(func(j))
        return out
    with mp.Pool(processes=min(nproc, len(jobs))) as pool:
        for k, r in enumerate(pool.imap_unordered(func, jobs)):
            out.append(r)
            if log and ((k + 1) % max(1, len(jobs) // 8) == 0 or k + 1 == len(jobs)):
                log("      %s %d/%d jobs done (%.0f s)" % (label, k + 1, len(jobs), time.time() - t0))
    return out
