# -*- coding: utf-8 -*-
"""referee evaluation of [u^0..u^3] G_Gamma for M = M_0 exp(u N_1 + u^2 N_2 + u^3 N_3), frozen or u-dependent vertex tensors"""
import sys, os, time
sys.dont_write_bytecode = True
PROJ = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); sys.path.insert(0, PROJ)
from fq_fp import P1, P2
from fq3_steps import frame_of
from fq3_par import setup
from fq3_m import MU, sp_mats, NSP
from fq3_mvertex import MVertex6
from fq3_msum import GraphSummer6, UVertex


def evaluate(tag, Ns, idxs, mode='frozen', NN=19, K=4, verbose=False):
    """Ns: list of 75-vectors (N_1, N_2, ..); returns ({idx: [u^0..u^(K-1)]}, maxdeg)"""
    p = P2 if tag == "89" else P1
    st = setup(frame_of(p, (0, 0), 34, K)); R = st['R']
    M = MU(R, st['M0'], [sp_mats(n, p) for n in Ns], K, 'right')
    if mode == 'frozen':
        vt = st.get('mv0')
        if vt is None: vt = MVertex6(R, st['M0']); st['mv0'] = vt
    else:
        vt = UVertex(R, M, K, NN, degree_check=False)
    GS = GraphSummer6(st['H'], M, vt); out = {}
    for i in idxs:
        t = time.time(); s = GS.graph(st['graphs'][i]); out[i] = [s.coef(k) if k < s.prec else None for k in range(K)]
        if verbose: print("      graph %d: %.0f s" % (i, time.time() - t), flush=True)
    return out, getattr(vt, 'maxdeg', None)
