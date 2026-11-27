# analysis2.py: coloured stable graphs, alpha dependence dec by dec, actual five-weight frames, non-CY probes, wrong-constant controls
import sys,os,time,random,itertools
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
OUT=os.path.dirname(os.path.abspath(__file__))
fs=open(os.path.join(OUT,"analysis2_out.txt"),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
Dec={}
def DD(g,h):
    if (g,h) not in Dec: Dec[(g,h)]=Decorated(g,h)
    return Dec[(g,h)]
NN={2:26,3:32,4:40}
random.seed(99)
vals=[random.randrange(2,P) for _ in range(4)]
Wr=dict(a1=vals[0],a2=vals[1],e1=vals[2],e2=vals[3]); Wr['a0']=(Wr['e1']+Wr['e2']-Wr['a1']-Wr['a2'])%P
# ---------- A. coloured stable graphs in frame B
pr("A. coloured stable graphs in frame B: sum over labelled decorations at fixed colouring")
for name,W in (('S1',cyw(3,4,5,11)),('S3',cyw(17,-5,19,-8)),('random',Wr)):
    for g in (2,3,4):
        t0=time.time()
        rows=Rows(W,3*g-3,N=NN[g])
        E=(W['a0']+W['a1']+W['a2'])%P; wB=[E,-W['e1']%P,-W['e2']%P]
        fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
        D=DD(g,True)
        nc=0;sing=0;singedge=0;ncedge=0;dec_col_viol=0;dec_col_sing=0;dec_col=0;sing_genus_le1=0
        for G in D.graphs:
            ai=inv(G['aut'])
            bound=sum(max(2*h-2,0) for h in G['gv'])
            for col in itertools.product((1,2),repeat=G['V']):
                tot={}
                for key,decs in G['orbs'].items():
                    sub={}
                    for dec in decs:
                        x=eval_dec(G,dec,col,fr,-1)
                        if x is None: continue
                        for i,y in enumerate(x[1]):
                            k=x[0]+i
                            if k<0 and y: sub[k]=(sub.get(k,0)+y)%P
                    for k,y in sub.items(): tot[k]=(tot.get(k,0)+y*ai)%P
                    # per coloured decorated orbit (orbit sum at this labelled colouring)
                    neg=[k for k,y in sub.items() if y]
                    dec_col+=1
                    if neg:
                        dec_col_sing+=1
                        if -min(neg)>bound: dec_col_viol+=1
                        if bound==0: sing_genus_le1+=1
                tot={k:y for k,y in tot.items() if y}
                nc+=1
                if G['ed']: ncedge+=1
                if tot:
                    sing+=1
                    if G['ed']: singedge+=1
        pr("  ",name,"g",g,": coloured stable graphs",nc,"(with edge",ncedge,"); singular:",sing,", singular with an edge:",singedge,"; coloured dec orbits",dec_col,"singular",dec_col_sing,"; above Sum(2g_v-2)+:",dec_col_viol,"; singular with all g_v<=1:",sing_genus_le1,"; %.0fs"%(time.time()-t0))
# ---------- B. actual five-weight frames
pr("\nB. frames whose vertex class is the actual fixed-point class, and related")
for name,tup in (('S1',(3,4,5,11)),('S3',(17,-5,19,-8))):
    W=cyw(*tup)
    a0,a1,a2,e1,e2=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']; E=(a0+a1+a2)%P
    wp=[a0,a1,a2,-e1%P,-e2%P]; wm=[-a0%P,(a0+a1)%P,(a0+a2)%P,-e1%P,-e2%P]; wB=[E,-e1%P,-e2%P]
    for g in (2,3):
        rows=Rows(W,3*g-3,N=NN[g]); ex=expected_total(g,W,1); cm=const_map(g,W)
        for nm,wc in (('actual (wt_U1 at c1, wt_U2 at c2) as (p+,p-)',{1:wp,2:wm}),('actual swapped',{1:wm,2:wp}),('c1=B c2=wt-',{1:wB,2:wm}),('c1=B c2=wt+',{1:wB,2:wp}),('both = nw (a1,a0+a2,-e1,-e2)',{1:[a1,(a0+a2)%P,-e1%P,-e2%P],2:[a1,(a0+a2)%P,-e1%P,-e2%P]})):
            fr=Frame(rows,wc,gmax=g+1)
            gs=GraphSum(DD(g,True),fr,1).run()
            tot=gs.total()
            ok=all(((tot.get(k,0)-ex.get(k,0)-(cm if k==0 else 0))%P)==0 for k in set(tot)|set(ex))
            bg=gs.by_graph()
            se=sum(1 for d in bg.values() if d['ed'] and any(k<0 for k in d['c']))
            pr("  ",name,"|",nm,"| g",g,"total ok",ok,"| singular edge graphs",se,"of",len(bg)-1)
# ---------- C. controls that must fail
pr("\nC. controls")
W=cyw(3,4,5,11); a0,a1,a2,e1,e2=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']; E=(a0+a1+a2)%P; wB=[E,-e1%P,-e2%P]
for g in (2,3):
    ex=expected_total(g,W,2); cm=const_map(g,W)
    def tot_ok(gs):
        tot=gs.total()
        return all(((tot.get(k,0)-ex.get(k,0)-(cm if k==0 else 0))%P)==0 for k in set(tot)|set(ex))
    rows0=Rows(W,3*g-3,N=NN[g],ffc=[1]+[0]*(3*g-3))
    fr=Frame(rows0,{1:wB,2:wB},gmax=g+1); gs=GraphSum(DD(g,True),fr,2).run()
    bg=gs.by_graph();
    pr("   g",g,"no QRR constant (ffc=1): total ok",tot_ok(gs),"| singular edge graphs",sum(1 for d in bg.values() if d['ed'] and any(k<0 for k in d['c'])))
    # constant with e1 -> e1+1
    W2=dict(W); W2['e1']=(W['e1']+1)%P
    rows1=Rows(W,3*g-3,N=NN[g],ffc=ffc_closed(W2,3*g-3))
    fr=Frame(rows1,{1:wB,2:wB},gmax=g+1); gs=GraphSum(DD(g,True),fr,2).run()
    bg=gs.by_graph();
    pr("   g",g,"wrong constant (e1->e1+1 in ffc): total ok",tot_ok(gs),"| singular edge graphs",sum(1 for d in bg.values() if d['ed'] and any(k<0 for k in d['c'])))
    rows=Rows(W,3*g-3,N=NN[g])
    fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
    wX=[(E+1)%P,-e1%P,-e2%P]
    fr.vt={1:VertexTable(wX),2:VertexTable(wX)}
    gs=GraphSum(DD(g,True),fr,2).run()
    pr("   g",g,"rows twisted by h_w, vertex class w'=(E+1,-e1,-e2): total ok",tot_ok(gs))
    fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
    gs=GraphSum(DD(g,True),fr,2).run()
    pr("   g",g,"correct frame B (positive control): total ok",tot_ok(gs))
# ---------- D. alpha dependence dec by dec at fixed (e1,e2)
pr("\nD. alpha dependence of the frame-B polar coefficients, decorated graph by decorated graph, at fixed (e1,e2)=(5,11)")
fixed=(5,11)
alist=[(3,4),(2,7),(5,11),(8,3),(11,5),(3,13),(1,6),(9,-2),(4,-7),(13,1),(6,6),(7,2),(-3,5),(10,9),(12,-1),(2,1),(Fr(1,2),3),(5,Fr(1,3)),(Fr(7,2),Fr(9,5)),(-1,-2)]
def solve(A,b):
    m=len(A); n=len(A[0]) if A else 0
    R=[row[:]+[b[i]] for i,row in enumerate(A)]
    piv=[];r=0
    for c in range(n):
        p=None
        for i in range(r,m):
            if R[i][c]%P: p=i;break
        if p is None: continue
        R[r],R[p]=R[p],R[r]
        iv=inv(R[r][c])
        R[r]=[x*iv%P for x in R[r]]
        for i in range(m):
            if i!=r and R[i][c]%P:
                f=R[i][c]
                R[i]=[(x-f*y)%P for x,y in zip(R[i],R[r])]
        piv.append(c);r+=1
    for i in range(r,m):
        if R[i][n]%P: return None
    sol=[0]*n
    for i,c in enumerate(piv): sol[c]=R[i][n]
    return sol,len(piv)
for g in (2,3):
    res=[]
    for (a1_,a2_) in alist:
        W=cyw(a1_,a2_,*fixed)
        E=(W['a0']+W['a1']+W['a2'])%P; wB=[E,-W['e1']%P,-W['e2']%P]
        rows=Rows(W,3*g-3,N=NN[g])
        fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
        gs=GraphSum(DD(g,True),fr,-1).run()
        u=(E-2*W['a1'])*(E-2*W['a2'])%P
        res.append((u,W['a0'],{(r['graph'],r['key']):r['c'] for r in gs.results}))
    keys=list(res[0][2].keys())
    for pw in range(-(2*g-2),0):
        # monomials in u, a0^2: test increasing total weight (u, a0^2 have degree 2)
        basis_by_deg=lambda d:[(i,j) for i in range(d+1) for j in range(d+1-i)]
        tally={}
        for key in keys:
            ys=[r[2][key].get(pw,0) for r in res]
            if all(y==0 for y in ys):
                tally['zero']=tally.get('zero',0)+1; continue
            found=None
            for d in range(0,3*g):
                bas=basis_by_deg(d)
                if len(bas)+3>len(res): break
                A=[[pow(r[0],i,P)*pow(r[1]*r[1]%P,j,P)%P for (i,j) in bas] for r in res]
                sol=solve(A,ys)
                if sol is not None:
                    # record which monomials occur and if u divides
                    found=(d,all(s==0 for s,(i,j) in zip(sol[0],bas) if i==0))
                    break
            k=('deg %d in (u,a0^2), divisible by u: %s'%found) if found else 'not polynomial in (u,a0^2)'
            tally[k]=tally.get(k,0)+1
        pr("   g",g,"t^%d:"%pw,sorted(tally.items()))
# ---------- E. off the CY locus
pr("\nE. off the CY locus: a0+a1+a2 = Ea != e1+e2 = Ee (the theorem's hypothesis fails; no closed form for the total used)")
for tup in ((2,3,4,5,11),(3,3,4,5,11),(1,7,2,5,3)):
    W=weights(*tup); a0,a1,a2,e1,e2=W['a0'],W['a1'],W['a2'],W['e1'],W['e2']
    Ea=(a0+a1+a2)%P; Ee=(e1+e2)%P
    for g in (2,3):
        rows=Rows(W,3*g-3,N=NN[g])
        pointtot=GraphSum(DD(g,False),Frame(rows,{1:None,2:None},gmax=g+1),1).run().total()
        pr("   weights",tup,"g",g,"t^2 Delta_1 =",fmt(rows.Delta[1][1][0]),"-Ea^3 e1 e2 =",fmt(-Ea**3*e1*e2%P))
        for nm,w in (('(Ea,-e1,-e2)',[Ea,-e1%P,-e2%P]),('(Ee,-e1,-e2)',[Ee,-e1%P,-e2%P])):
            fr=Frame(rows,{1:w,2:w},gmax=g+1)
            gs=GraphSum(DD(g,True),fr,1).run()
            tot=gs.total()
            bg=gs.by_graph()
            sg=sorted(i for i,d in bg.items() if any(k<0 for k in d['c']))
            r01ok=all(s_coef(fr.r[1][0][m],0)==pow(-1,m,P)*pow(inv(Ea),m,P)%P for m in range(3*g-2))
            vg=[i for i,d in bg.items() if not d['ed']][0]
            pr("      ",nm,"total = point-frame total:",tot==pointtot,"| r'_{0,1}(z;0)=Ea/(Ea+z):",r01ok,"| singular stable graphs:",len(sg),"of",len(bg),sg[:8],"| vertex polar powers:",sorted(k for k in bg[vg]['c'] if k<0),"| singular edge graphs:",sum(1 for i in sg if bg[i]['ed']))
pr("done")
