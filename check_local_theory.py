import math, random
from math import gcd
def factor(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1 if d==2 else 2
    if n>1: f[n]=f.get(n,0)+1
    return f
def bperp(b,l):  # part of b from primes != l, not 1 mod l
    r=1
    for p,e in factor(b).items():
        if p!=l and p%l!=1: r*=p**e
    return r
def blk(k,x): return math.prod(range(x,x+k))
# 1) known g=3 solutions
for k,x,V,U in [(4,25818,56551,763231),(4,196441,3879095,11566799)]:
    b=blk(k,x); assert U**3-V**3==b
    d=U-V; q=(U**3-V**3)//d; bp=bperp(b,3)
    print("k",k,"x",x,"d|b",b%d==0,"bperp|d",d%bp==0,"log bperp/log b=%.3f"%(math.log(bp)/math.log(b)),"<1/3 ok",bp<b**(1/3),"gcd(d,q)=",gcd(d,q),"primes of q mod3:",sorted({p%3 for p in factor(q)}))
# 2) g=4 solutions: primes 3 mod 4 part <= D < sqrt(b)
for k,x,V,U in [(3,40,11,17),(3,88,7,29)]:
    b=blk(k,x); assert U**4-V**4==b
    D=U*U-V*V; q=U*U+V*V
    p3=1
    for p,e in factor(b).items():
        if p%4==3: p3*=p**e
    print("g=4 x",x,"p3 part",p3,"D",D,"sqrt b %.1f"%b**.5,"D%p3==0",D%p3==0,"gcd(D,q)",gcd(D,q),"q primes mod4",sorted({p%4 for p in factor(q)}))
# 3) unitary count of w^2-v^2=b, coprime
def cnt(b):
    return sum(1 for v in range(1,int(b**.5)+2) if (lambda w2: math.isqrt(w2)**2==w2 and gcd(v,math.isqrt(w2))==1)(b+v*v)) if False else None
def brute(b):
    c=0
    for d in range(1,math.isqrt(b)+1):
        if b%d==0:
            e=b//d
            if d<e and (d+e)%2==0:
                w=(d+e)//2; v=(e-d)//2
                if v>=1 and gcd(v,w)==1: c+=1
    return c
bad=0
for b in list(range(1,3000))+[blk(k,x) for k in (2,3,4,5) for x in range(1,60)]:
    om=len(factor(b//4)) if b%4==0 else None
    pred=0 if b%8 else (2**(om-1) if b//4>1 else 0)
    if b%8==0 and b//4==1: pred=0
    if brute(b)!=pred: bad+=1; print("mismatch",b,brute(b),pred)
print("unitary-count formula mismatches:",bad)
# 4) density of x with bperp<b^(1/3), k=4
import sys
X=60000; k=4
spf=list(range(X+k+2))
for i in range(2,int((X+k+1)**.5)+1):
    if spf[i]==i:
        for j in range(i*i,X+k+2,i):
            if spf[j]==j: spf[j]=i
ok=0; fr=[]
for x in range(1,X+1):
    lb=0.0; lp=0.0
    for t in range(k):
        n=x+t
        while n>1:
            p=spf[n]; n//=p; lb+=math.log(p)
            if p!=3 and p%3!=1: lp+=math.log(p)
    fr.append(lp/lb)
    if lp<lb/3: ok+=1
fr.sort()
print("k=4,x<=%d: fraction with bperp<b^(1/3): %.5f (%d of %d); median log bperp/log b = %.3f"%(X,ok/X,ok,X,fr[len(fr)//2]))
