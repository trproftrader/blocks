import math, sys, collections
def iroot_ceil(v,n):
    r=int(round(v**(1.0/n)))
    while r**n<v: r+=1
    while r>1 and (r-1)**n>=v: r-=1
    return r
def maxexp(p):
    best=(p,1)
    for n in range(2,64):
        if 2**n>p: break
        r=int(round(p**(1.0/n)))
        for t in (r-1,r,r+1):
            if t>=2 and t**n==p: best=(t,n)
    return best
C=60
sol=collections.defaultdict(list)
for k,B in [(2,10**13),(3,10**15),(4,10**15),(5,10**15),(6,10**15),(7,10**15),(8,10**15),(9,10**15)]:
    x=1
    while True:
        b=math.prod(range(x,x+k))
        if b>B: break
        seen=set()
        for n in range(2,int(math.log2(b+C))+1):
            lo,hi=b+1,b+C
            r=iroot_ceil(lo,n)
            while r**n<=hi:
                seen.add(r**n); r+=1
        for p in seen:
            z,n=maxexp(p); c=p-b
            if k==4 and c==1: continue
            if k==2 and n==2 and False: pass
            sol[c].append((k,x,z,n))
        x+=1
tot=0
for c in sorted(sol):
    L=sorted(set(sol[c]))
    tot+=len(L)
print("total solutions with 1<=c<=%d, b<=1e13/1e15:"%C,tot)
for c in sorted(sol):
    print("c=%d:"%c," ".join("b%d(%d)+%d=%d^%d"%(k,x,c,z,n) for k,x,z,n in sorted(set(sol[c]))))
