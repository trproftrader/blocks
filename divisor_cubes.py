import sys
from math import isqrt, gcd
X=int(sys.argv[1]); K=int(sys.argv[2])
spf=list(range(X+K+2))
for i in range(2,int((X+K+1)**.5)+1):
    if spf[i]==i:
        for j in range(i*i,X+K+2,i):
            if spf[j]==j: spf[j]=i
def fac(n,f):
    while n>1:
        p=spf[n]; f[p]=f.get(p,0)+1; n//=p
def divisors(f):
    ds=[1]
    for p,e in f.items():
        ds=[d*p**i for d in ds for i in range(e+1)]
    return ds
def ispp(n):
    if n<2: return False
    for e in range(2,64):
        r=round(n**(1.0/e))
        for c in (r-1,r,r+1):
            if c>=2 and c**e==n: return True
    return False
out=[]
for x in range(1,X+1):
    f={}
    b=1
    for t in range(K):
        fac(x+t,f); b*=x+t
    for d in divisors(f):
        if d**3>=b: continue
        q=b//d
        disc=12*q-3*d*d
        if disc<=0: continue
        s=isqrt(disc)
        if s*s!=disc or (s-3*d)%6: continue
        y=(s-3*d)//6
        if y<1: continue
        z=y+d
        assert z**3-y**3==b
        out.append((x,y,z,gcd(gcd(b,y),z)))
cop=[o for o in out if o[3]==1]
print(f"k={K} x<={X}: solutions of z^3-y^3=b_k(x): {len(out)}; coprime: {len(cop)}")
for x,y,z,g in cop:
    print(f"b_{K}({x}) + {y}^3 = {z}^3   (z^3={z**3:.3e}) y_pp={ispp(y)} z_pp={ispp(z)}")
