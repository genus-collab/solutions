import pickle
from fractions import Fraction as F
exec(open('opsym.py').read().split("import pickle")[0])
Q,Rec=pickle.load(open('ops.pkl','rb'))
def divlin(p,c0):  # divide p(m) by (m+c0); return quotient, remainder (polys in z,E)
    deg=max([k[0] for k in p]+[0])
    coef=[{} for _ in range(deg+1)]
    for (a,b,c),x in p.items(): coef[a][(b,c)]=coef[a].get((b,c),0)+x
    q=[{} for _ in range(deg)]
    cur={}
    for a in range(deg,0,-1):
        cur={k:coef[a].get(k,0)+cur.get(k,0) for k in set(coef[a])|set(cur)}
        cur={k:v for k,v in cur.items() if v}
        q[a-1]=cur
        cur={k:-c0*v for k,v in cur.items()}
    rem={k:coef[0].get(k,0)+cur.get(k,0) for k in set(coef[0])|set(cur)}
    rem={k:v for k,v in rem.items() if v}
    Qp={}
    for a,d in enumerate(q):
        for (b,c),x in d.items(): Qp[(a,b,c)]=x
    return Qp,rem
D={}
R=dict(Q)
for j in [1,0,-1,-2]:
    s=j+1
    p=R.get(s,{})
    # p = d_j(m)*27*(m+j)
    if j==0:
        # rec_1(m)=27 m: divide by m
        qq,rem=divlin(p,0)
    else:
        qq,rem=divlin(p,j)
    assert not rem, ("nonzero remainder",j,rem)
    d=pscal(F(1,27),qq)
    D[j]=d
    R=oadd(R,{k:pscal(F(-1),v) for k,v in omul({j:d},Rec).items()})
print("remaining after subtraction:", R)
def show(p):
    return " + ".join("%s*m^%d z^%d E^%d"%(v,a,b,c) for (a,b,c),v in sorted(p.items()))
for j in sorted(D): print("d_%d ="%j, show(D[j]))
chk=oadd(Q,{k:pscal(F(-1),v) for k,v in omul(D,Rec).items()})
print("Q - D o Rec == 0:", chk=={})
pickle.dump(D,open('D.pkl','wb'))
