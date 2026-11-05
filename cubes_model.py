import sys, math, numpy as np
def run(k,X):
    spf=list(range(X+k+2))
    for i in range(2,int((X+k+1)**.5)+1):
        if spf[i]==i:
            for j in range(i*i,X+k+2,i):
                if spf[j]==j: spf[j]=i
    E=0.0
    for x in range(1,X+1):
        f={}; b=1
        for t in range(k):
            n=x+t; b*=n
            while n>1:
                p=spf[n]; f[p]=f.get(p,0)+1; n//=p
        ds=[1]
        for p,e in f.items():
            pw=[p**i for i in range(e+1)]
            ds=[d*q for d in ds for q in pw]
        d=np.array(ds,dtype=float)
        bf=float(b)
        m=d**3<bf
        d=d[m]
        disc=12*bf/d-3*d*d
        E+=float(np.sum(1.0/(2*np.sqrt(disc))))
    return E
for k,X,obs in [(4,30000,10),(4,200000,11),(5,100000,3),(6,8000,1)]:
    print(f"k={k} x<={X}: expected (naive, all reduced solutions) {run(k,X):.2f}   observed {obs}")
