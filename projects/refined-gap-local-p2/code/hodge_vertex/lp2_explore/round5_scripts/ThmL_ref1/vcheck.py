import math
def parts(n,m=None):
    m=n if m is None else m
    if n==0:
        yield ()
        return
    for k in range(min(n,m),0,-1):
        for r in parts(n-k,k):
            yield (k,)+r
lam=0.37
s=lambda t: math.sin(t*lam/2)
maxerr=0
for n in range(1,11):
    for nu in parts(n):
        l=len(nu)
        V=1.0
        for a in range(l):
            for b in range(a+1,l):
                V*=s(nu[a]-nu[b]+b-a)/s(b-a)
        for i in range(l):
            for v in range(1,nu[i]+1):
                V/=2*s(v-(i+1)+l)
        conj=[sum(1 for x in nu if x>j) for j in range(nu[0])]
        H=1.0
        for i in range(l):
            for j in range(nu[i]):
                H/=2*s(nu[i]-j+conj[j]-i-1)
        maxerr=max(maxerr,abs(V/H-1))
print("max rel err (3) vs hook form with 2^|nu|:",maxerr)
