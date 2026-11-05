import sys, math, csv, re, collections
import numpy as np
# usage: python3 heur_fact2.py [solutions.csv] k [k ...]   (default csv: solutions_N1e16.csv)
args=sys.argv[1:]
CSV=args.pop(0) if args and args[0].endswith('.csv') else 'solutions_N1e16.csv'
N=10**16; T=14
def spf_sieve(M):
    spf=list(range(M+1))
    for i in range(2,int(M**.5)+1):
        if spf[i]==i:
            for j in range(i*i,M+1,i):
                if spf[j]==j: spf[j]=i
    return spf
def run(k):
    Xmax=1
    while math.prod(range(Xmax,Xmax+k))<=N: Xmax+=1
    Xmax-=1
    spf=spf_sieve(Xmax+k+2)
    E=np.zeros((T+1,T+1))
    ts=np.arange(1,T+1,dtype=float)
    for x in range(1,Xmax+1):
        f={}
        b=1
        for t in range(k):
            n=x+t; b*=n
            while n>1:
                p=spf[n]; f[p]=f.get(p,0)+1; n//=p
        ds=[1]
        for p,e in f.items():
            pw=[p**i for i in range(e+1)]
            ds=[d*q for d in ds for q in pw]
        d=np.array(ds,dtype=np.int64); e_=b//d
        m=(d<e_)&(((d+e_)&1)==0)
        d=d[m]; e_=e_[m]
        w=(d+e_)//2; v=(e_-d)//2
        m=(v>=2)&(v<=10**8)&(w<=141421356)
        v=v[m]; w=w[m]
        if len(v)==0: continue
        m=np.gcd(v,w)==1
        v=v[m].astype(float); w=w[m].astype(float)
        if len(v)==0: continue
        PV=(v[None,:]**(1.0/ts[:,None]-1))/ts[:,None]*(v[None,:]>=2.0**ts[:,None])
        PW=(w[None,:]**(1.0/ts[:,None]-1))/ts[:,None]*(w[None,:]>=2.0**ts[:,None])
        E[1:,1:]+=PV@PW.T
    return Xmax,E
obs=collections.Counter()
for r in csv.DictReader(open(CSV)):
    if r['type']!='2.17' or '1^M' in r['equation']: continue
    m=re.match(r'b_(\d+)\((\d+)\) \+ (\d+)\^(\d+) = (\d+)\^(\d+)',r['equation'])
    k,x,y,mm,z,n=map(int,m.groups()); obs[(k,mm,n)]+=1
for k in [int(a) for a in args]:
    Xmax,E=run(k)
    print(f"\n=== k={k}  (x<= {Xmax}) ===  signature (m,n)=(2a,2c): observed / expected(factorisation model)")
    tot_o=tot_e=0; rows=[]
    for a in range(1,T+1):
        for c in range(1,T+1):
            if not (1/k+1/(2*a)+1/(2*c)<1-1e-12): continue
            o=obs.get((k,2*a,2*c),0); e=E[a,c]
            tot_o+=o; tot_e+=e
            if o>0 or e>=0.15: rows.append((2*a,2*c,o,e))
    for m_,n_,o,e in sorted(rows,key=lambda r:-(r[2]+r[3])): print(f"  (k,m,n)=({k},{m_},{n_})  obs {o:3d}  exp {e:7.2f}")
    odd=sum(v for (kk,mm,n),v in obs.items() if kk==k and not(mm%2==0 and n%2==0))
    print(f"  TOTAL even-even hyperbolic: observed {tot_o}, expected {tot_e:.1f}   (observed with an odd exponent: {odd})")
