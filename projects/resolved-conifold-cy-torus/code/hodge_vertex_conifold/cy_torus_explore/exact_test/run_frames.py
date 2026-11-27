# run_frames.py g [T] : all frames at all samples; writes results_g<g>.txt (per decorated graph) and summary_g<g>.txt
import sys,os,time,json,random
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
import cylib
g=int(sys.argv[1]); T=int(sys.argv[2]) if len(sys.argv)>2 else 2
zmax=3*g-3; N=26+6*(g-2)
OUT=os.path.dirname(os.path.abspath(__file__))
fd=open(os.path.join(OUT,"decorated_g%d.txt"%g),"w")
fs=open(os.path.join(OUT,"summary_g%d.txt"%g),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
t00=time.time()
DP=Decorated(g,False); DH=Decorated(g,True)
pr("genus",g,"zmax",zmax,"rel prec N",N,"T",T,"point",DP.counts(),"hodge",DH.counts(),"enumeration",round(time.time()-t00,1),"s","p=",P)
random.seed(20261102)
samples={}
samples['F(x=3,5) fibrewise']=(3,5,3,5)
samples['F(5,11) fibrewise']=(5,11,5,11)
samples['S1']=(3,4,5,11)
samples['S1b same e']=(2,7,5,11)
samples['S2']=(7,13,2,29)
samples['S3']=(17,-5,19,-8)
samples['S4 rational']=(Fr(3,7),Fr(-11,5),Fr(13,9),Fr(17,4))
samples['S5 a0=0 swapped a1=e2,a2=e1']=(11,5,5,11)
samples['S6 a0=0, a1 != e1,e2']=(3,13,5,11)
samples['S7 u=0 (a1=E/2)']=(8,3,5,11)
samples['S8 random F_p']=None
def frames(W):
    E=(W['a0']+W['a1']+W['a2'])%P
    a1,a2,e1,e2=W['a1'],W['a2'],W['e1'],W['e2']
    wB=[E,-e1%P,-e2%P]; wa1=[E,-a1%P,(a1-E)%P]; wa2=[E,-a2%P,(a2-E)%P]
    fr={}
    fr['i point']=(DP,{1:None,2:None},-1)
    fr['ii B=(E,-e1,-e2)']=(DH,{1:wB,2:wB},-1)
    fr['iii (E,-a1,a1-E)']=(DH,{1:wa1,2:wa1},-1)
    fr['iii (E,-a2,a2-E)']=(DH,{1:wa2,2:wa2},-1)
    fr['iv c1=B c2=point']=(DH,{1:wB,2:None},-1)
    fr['iv c1=B c2=(E,-a1,a1-E)']=(DH,{1:wB,2:wa1},-1)
    fr['iv c1=(E,-a1,a1-E) c2=B']=(DH,{1:wa1,2:wB},-1)
    fr['iv c1=point c2=B']=(DH,{1:None,2:wB},-1)
    # wrong classes (controls)
    fr['ctl modeA (-E,e1,e2)']=(DH,{1:[-x%P for x in wB],2:[-x%P for x in wB]},-1)
    fr['ctl wrong (E+1,-e1,-e2)']=(DH,{1:[(E+1)%P,-e1%P,-e2%P],2:[(E+1)%P,-e1%P,-e2%P]},-1)
    fr['ctl sign s=+1 (total must fail)']=(DH,{1:wB,2:wB},1)
    return fr
def runW(name,tup):
    if tup is None:
        vals=[random.randrange(2,P) for _ in range(4)]
        W=dict(a1=vals[0],a2=vals[1],e1=vals[2],e2=vals[3]); W['a0']=(W['e1']+W['e2']-W['a1']-W['a2'])%P
        tup=tuple(vals)
    else: W=cyw(*tup)
    t0=time.time()
    rows=Rows(W,zmax,N=N)
    E=(W['a0']+W['a1']+W['a2'])%P
    pr("\n=== sample",name,tup if name!='S8 random F_p' else '(random residues)',"a0=",fmt(W['a0']),"E=",fmt(E),"; Delta_1 t^2 -> ",fmt(rows.Delta[1][1][0]),"(= -E^3 e1 e2:",(rows.Delta[1][1][0]+E*E*E*W['e1']*W['e2'])%P==0,"); rows built",round(time.time()-t0,2),"s")
    exp=expected_total(g,W,T)
    cm=const_map(g,W)
    pol=expected_total(g,W,0)
    for nm,(D,wc,s) in frames(W).items():
        t0=time.time()
        fr=Frame(rows,wc,s=s,gmax=g+1)
        gs=GraphSum(D,fr,T).run(percol=True)
        tot=gs.total()
        ok=True; diffs={}
        for k in sorted(set(tot)|set(exp)):
            x=(tot.get(k,0)-exp.get(k,0))%P
            if k==0: x=(x-cm)%P if cm is not None else 0
            if x: ok=False; diffs[k]=fmt(x)
        bg=gs.by_graph()
        tally={}
        edge_sing=0; edge_sum_sing=0; maxviol=0; oddlead=0; col_sing_edge=0; dec_cols_sing_in_regular=0
        lines=[]
        for r in gs.results:
            neg=[k for k,y in r['c'].items() if k<0]
            po=-min(neg) if neg else 0
            b=r['bound']
            tally[(po,b)]=tally.get((po,b),0)+1
            if po%2: oddlead+=1
            if po>b: maxviol+=1
            if r['ed'] and po>0: edge_sing+=1
            colsing=sum(1 for c,d in r['pc'].items() if any(k<0 for k in d))
            if r['ed']: col_sing_edge+=colsing
            if po==0 and colsing: dec_cols_sing_in_regular+=1
            fd.write("%s | %s | graph %d gv %s key %s | pole %d bound %d | polar %s | colourings singular %d\n"%(name,nm,r['graph'],r['gv'],r['key'],po,b,{k:fmt(y) for k,y in sorted(r['c'].items()) if k<0},colsing))
        sg=[]
        for i,d in sorted(bg.items()):
            neg={k:y for k,y in d['c'].items() if k<0}
            if neg: sg.append(i)
        edge_graph_sing=[i for i in sg if bg[i]['ed']]
        # vertex graph polar
        vg=[i for i,d in bg.items() if not d['ed']][0]
        vp={k:y for k,y in bg[vg]['c'].items() if k<0}
        gapok = vp=={k:y for k,y in pol.items() if k<0}
        pr(" %-34s total ok: %-5s | stable graphs singular (edge graphs) %d/%d | edge decorated graphs singular %d | pole>bound %d | odd lead %d | singular colourings in edge graphs %d | regular dec with singular colourings %d | vertex graph polar = gap: %s | tally(pole,bound) %s | %.1fs"%(nm,ok,len(edge_graph_sing),len(bg)-1,edge_sing,maxviol,oddlead,col_sing_edge,dec_cols_sing_in_regular,gapok,sorted(tally.items()),time.time()-t0))
        if not ok: pr("     total differs at:",dict(list(diffs.items())[:4]))
for n,t in samples.items():
    runW(n,t)
pr("done",round(time.time()-t00,1),"s")
