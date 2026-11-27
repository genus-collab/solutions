# appended: Lemma 1 closed forms check
def Sm(mm):
    out={}
    for l in range(0,CAP-1):
        zl=Z0(); ip=iprod(max(l+1-mm,0))
        for k in range(NZ+1-l): zl[k+l]=ip[k]*r[l]
        out[l]=zl
    return S(out)
S1=Sm(1); S2_=Sm(2); iv=zf([1],-1)
Ph1=mul(iv,P1); Ph2=mul(iv,P2)
cf1=sc(E,mul(mul(C1i,mul(thv,iv)),S1))
A=mul(C1i,thv)
tl=mul(th(A),inv(A))  # theta log(C1i theta v)
inner=add(sc(E,mul(thv,mul(mul(iv,iv),S2_))), mul(add(zz(tl),sc(-1,d2)),mul(S1,iv)))
cf2=sc(E,mul(mul(Del,mul(c1,thv)),inner))
D1=add(Ph1,sc(-1,cf1)); D2=add(Ph2,sc(-1,cf2))
print("Lemma1 row1 residual v-orders:",[k for k in sorted(D1) if k<CAP-3],"row2:",[k for k in sorted(D2) if k<CAP-4])
