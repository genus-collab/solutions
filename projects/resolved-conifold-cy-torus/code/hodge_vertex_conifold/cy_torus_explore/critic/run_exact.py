# Critic: run the other code (exact_test/cylib.py) at a fresh sample for a cross-check of row_audit.
import sys,os,time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'exact_test'))
from cylib import *
from fractions import Fraction as Fr
W = cyw(Fr(-3,8),Fr(13,6),Fr(11,9),Fr(709,360))
print({k:(v) for k,v in W.items()} if 0 else 'ok', flush=True)
for g,T,zmax,N in ((2,3,5,30),(3,2,7,34),(4,2,9,42)):
    t0=time.time()
    DH=Decorated(g,True)
    rows=Rows(W,zmax,N=N)
    E=(W['a0']+W['a1']+W['a2'])%P
    e1,e2=W['e1'],W['e2']
    wB=[E,-e1%P,-e2%P]
    exp=expected_total(g,W,T); cm=const_map(g,W); pol=expected_total(g,W,0)
    fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
    gs=GraphSum(DH,fr,T).run(percol=True)
    tot=gs.total()
    diffs=[]
    for k in sorted(set(tot)|set(exp)):
        x=(tot.get(k,0)-exp.get(k,0))%P
        if k==0: x=(x-cm)%P if cm is not None else 0
        if x: diffs.append(k)
    bg=gs.by_graph()
    sg=[i for i,d in bg.items() if any(k<0 for k in d['c'])]
    vg=[i for i,d in bg.items() if not d['ed']][0]
    print('g',g,'total ok',not diffs,'singular stable graphs',sg,'of',len(bg),'vertex graph index',vg,'%.0fs'%(time.time()-t0),flush=True)
