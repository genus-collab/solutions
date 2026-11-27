# run_g4.py: genus 4, frame B (and the colour-1-only variant) at generic CY weights; also one wrong frame on a subset is skipped for time
import sys,os,time,json,random
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
g=4; T=2; zmax=9; N=42
OUT=os.path.dirname(os.path.abspath(__file__))
fs=open(os.path.join(OUT,"summary_g4.txt"),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
t00=time.time()
DH=Decorated(g,True)
pr("genus 4 hodge",DH.counts(),"enum",round(time.time()-t00,1),"s")
random.seed(44)
vals=[random.randrange(2,P) for _ in range(4)]
W8=dict(a1=vals[0],a2=vals[1],e1=vals[2],e2=vals[3]); W8['a0']=(W8['e1']+W8['e2']-W8['a1']-W8['a2'])%P
samples={'F(3,5)':cyw(3,5,3,5),'S1':cyw(3,4,5,11),'S3':cyw(17,-5,19,-8),'S8 random F_p':W8}
for name,W in samples.items():
    t0=time.time()
    rows=Rows(W,zmax,N=N)
    pr("=== sample",name,"rows",round(time.time()-t0,1),"s")
    E=(W['a0']+W['a1']+W['a2'])%P
    a1,a2,e1,e2=W['a1'],W['a2'],W['e1'],W['e2']
    wB=[E,-e1%P,-e2%P]; 
    exp=expected_total(g,W,T); cm=const_map(g,W)
    pol=expected_total(g,W,0)
    for nm,wc in (('ii B',{1:wB,2:wB}),):
        t0=time.time()
        fr=Frame(rows,wc,gmax=g+1)
        gs=GraphSum(DH,fr,T).run(percol=True)
        tot=gs.total()
        diffs={}
        for k in sorted(set(tot)|set(exp)):
            x=(tot.get(k,0)-exp.get(k,0))%P
            if k==0: x=(x-cm)%P if cm is not None else 0
            if x: diffs[k]=fmt(x)
        bg=gs.by_graph()
        tally={}; viol=0; odd=0; edge_sing=0; sg=[]
        for r in gs.results:
            neg=[k for k,y in r['c'].items() if k<0]; po=-min(neg) if neg else 0; b=r['bound']
            tally[(po,b)]=tally.get((po,b),0)+1
            if po>b: viol+=1
            if po%2: odd+=1
            if r['ed'] and po: edge_sing+=1
        for i,d in bg.items():
            if any(k<0 for k in d['c']): sg.append(i)
        vg=[i for i,d in bg.items() if not d['ed']][0]
        vp={k:y for k,y in bg[vg]['c'].items() if k<0}
        pr(nm,"total ok:",not diffs,diffs if diffs else "","| singular stable graphs:",sg,"of",len(bg),"| edge dec singular",edge_sing,"| pole>bound",viol,"| odd lead",odd,"| vertex graph polar=gap:",vp=={k:y for k,y in pol.items() if k<0},"| tally",sorted(tally.items()),"| %.0fs"%(time.time()-t0))
pr("done",round(time.time()-t00,1))
