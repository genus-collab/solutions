def multisets(D,minp=2,maxr=9):
    """multisets of parts>=2 with sum(part-1)=D, at most maxr parts"""
    def rec(D,maxpart,left):
        if D==0:
            yield ()
            return
        if left==0: return
        for k in range(min(D+1,maxpart),1,-1):
            for r in rec(D-(k-1),k,left-1):
                yield (k,)+r
    return list(rec(D,10**6,maxr))
