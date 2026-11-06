import math
def crit(k):
    # unique critical point in each I_j=(-j-1,-j), j=0..k-2: root of g(x)=sum 1/(x+i), g strictly decreasing on I_j
    out=[]
    for j in range(k-1):
        lo,hi=-j-1+1e-15,-j-1e-15
        for _ in range(200):
            mid=(lo+hi)/2
            g=sum(1.0/(mid+i) for i in range(k))
            if g>0: lo=mid
            else: hi=mid
        c=(lo+hi)/2
        logM=sum(math.log(abs(c+i)) for i in range(k))
        out.append((j,c,logM))
    return out
bad=[]
for k in range(3,201):
    L=[t[2] for t in crit(k)]            # log M_j, j=0..k-2
    n=len(L)
    sym=all(abs(L[j]-L[n-1-j])<1e-6*max(1,abs(L[j])) for j in range(n))
    mid=(k-2)/2
    dec=all(L[j]>L[j+1]-1e-9 for j in range(0,int(math.floor(mid))))   # decreasing from j=0 to the centre
    inc=all(L[j]<L[j+1]+1e-9 for j in range(int(math.ceil(mid)),n-1))  # increasing from the centre to j=k-2
    strict=all(L[j]>L[j+1] for j in range(0,int(math.floor(mid)))) and all(L[j]<L[j+1] for j in range(int(math.ceil(mid)),n-1))
    if not(sym and strict): bad.append(k)
print("k in 3..200 where symmetric + strictly V-shaped (|f(c_j)| increasing outward) FAILS:",bad)
