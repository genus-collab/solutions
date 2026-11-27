from rows import *
def solve_mod(A,b):
    # least squares free: A m x n, exact solve with extra equations check; returns x or None
    m=len(A); n=len(A[0]) if m else 0
    M=[row[:]+[b[i]] for i,row in enumerate(A)]
    piv=[];r=0
    for c in range(n):
        p=None
        for i in range(r,m):
            if M[i][c]%P: p=i;break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        iv=inv(M[r][c])
        M[r]=[x*iv%P for x in M[r]]
        for i in range(m):
            if i!=r and M[i][c]%P:
                f=M[i][c]
                M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c);r+=1
    for i in range(r,m):
        if M[i][n]%P: return None
    x=[0]*n
    for i,c in enumerate(piv): x[c]=M[i][n]
    return x

