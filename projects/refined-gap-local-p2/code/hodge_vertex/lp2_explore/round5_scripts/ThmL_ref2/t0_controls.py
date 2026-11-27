from hodge import *
def chk(name,a,b):
    print(("ok  " if a==b else "FAIL"),name,a,b)
chk("<tau_4>_2",wk(2,(4,)),Fr(1,1152)); chk("<tau_2 tau_3>_2",wk(2,(2,3)),Fr(29,5760))
chk("<tau_7>_3",wk(3,(7,)),Fr(1,82944)); chk("<tau_10>_4",wk(4,(10,)),Fr(1,7962624))
chk("lam1^3 g2",corr(2,(),(1,1,1)),Fr(1,2880))
chk("<lam1 tau_3>_2",corr(2,(3,),(1,)),Fr(1,480))
chk("<lam2 tau_2>_2",corr(2,(2,),(1,1))/2,Fr(7,5760))
# lambda_3 = ch1^3/6 + 2 ch3
chk("<lam3 tau_4>_3",corr(3,(4,),(1,1,1))/6+2*corr(3,(4,),(3,)),Fr(31,967680))
chk("lam3^2 g3",corr(3,(),(1,)*6)/36+Fr(2,3)*corr(3,(),(1,1,1,3))+4*corr(3,(),(3,3)),Fr(0))
# lambda_4 psi^6 in genus 4 : (2^7-1)/2^7 |B_8|/8! = 127/128 * (1/30)/40320
# lambda_4 = e4 of chern roots: from exp(ch1 s + 0 + 2 ch3 s^3): s^4: ch1^4/24 + 2 ch1 ch3
chk("<lam4 tau_6>_4",corr(4,(6,),(1,1,1,1))/24+2*corr(4,(6,),(1,3)),Fr(127,128)*Fr(1,30)/40320)
tau=Fr(2,7); s=tau*(1+tau)
chk("X_1(1)",X(1,(1,),tau),-(1+1/s)/24)
# Omega expansion control against direct Lambda product in genus 2
def direct_g2(ds,tau):
    ts=[Fr(1),tau,-1-tau]
    P={0:Fr(1)}
    for t in ts:
        Q={}
        for k,v in P.items():
            Q[k]=Q.get(k,0)+v*t*t
            Q[k+1]=Q.get(k+1,0)-v*t
            Q[k+2]=Q.get(k+2,0)+v*Fr(1,2)
        P=Q
    sv=(-tau*(1+tau))**2
    return sum(v*corr(2,ds,(1,)*k) for k,v in P.items())/sv
chk("M_2(5)",M(2,(5,),tau),sum(Fr(5)**d*direct_g2((d,),tau) for d in range(5)))
# Omega_g genus-3 direct: Lambda^v(t) = t^3 - l1 t^2 + l2 t - l3 ; l1=ch1,l2=ch1^2/2,l3=ch1^3/6+2ch3
def direct_g3(ds,tau):
    ts=[Fr(1),tau,-1-tau]
    # poly dict: (n1,n3)-> coef
    P={(0,0):Fr(1)}
    for t in ts:
        Q={}
        def add(k,v): Q[k]=Q.get(k,0)+v
        for (a,b),v in P.items():
            add((a,b),v*t**3)
            add((a+1,b),-v*t**2)
            add((a+2,b),v*t*Fr(1,2))
            add((a+3,b),-v*Fr(1,6))
            add((a,b+1),-2*v)
        P=Q
    sv=(-tau*(1+tau))**3
    return sum(v*corr(3,ds,(1,)*a+(3,)*b) for (a,b),v in P.items())/sv
chk("M_3(4)",M(3,(4,),tau),sum(Fr(4)**d*direct_g3((d,),tau) for d in range(8)))
for t in [Fr(2,7),Fr(-5,3),Fr(3,5),Fr(-1,2),Fr(11,4)]:
    print("tau",t,"X2(2)",X(2,(2,),t),"X2(1)",X(2,(1,),t),"X2(1,1)",X(2,(1,1),t),"X2(1,2)",X(2,(1,2),t),"X2(3) (out)",X(2,(3,),t),"X2(2,2) (out)",X(2,(2,2),t),"X2()",X(2,(),t))
