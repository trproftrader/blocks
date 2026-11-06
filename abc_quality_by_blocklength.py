import csv,re,math,collections,sys
from math import gcd
sys.setrecursionlimit(10000)
def factor(n):
    f=set();d=2
    while d*d<=n:
        while n%d==0: f.add(d); n//=d
        d+=1 if d==2 else 2
    if n>1: f.add(n)
    return f
cache={}
def pf(n):
    if n not in cache: cache[n]=factor(n)
    return cache[n]
def parse(t):
    m=re.match(r'b_(\d+)\((\d+)\)',t)
    if m:
        e,x=int(m.group(1)),int(m.group(2)); s=set(); v=1
        for i in range(e): s|=pf(x+i); v*=x+i
        return v,s,e
    m=re.match(r'(\d+)\^(\d+|M)',t); b=int(m.group(1))
    return (b**(1 if m.group(2)=='M' else int(m.group(2)))), (pf(b) if b>1 else set()), 0
res=collections.defaultdict(list)
for r in csv.DictReader(open('solutions_N1e16.csv')):
    parts=re.split(r' \+ | = ',r['equation'])
    T=[parse(t) for t in parts]
    (A,sa,ka),(B,sb,kb),(C,sc,kc)=T
    g=gcd(A,B); a,b,c=A//g,B//g,C//g
    if c<=2: continue            # trivial 1+1=2
    P=sa|sb|sc
    rad=1
    for p in P:
        if any(v%p==0 for v in (a,b,c)): rad*=p
    q=math.log(c)/math.log(rad) if rad>1 else float('nan')
    kmin=min([k for k in (ka,kb,kc) if k] or [0])
    res[(r['type'],'k>=3' if kmin>=3 else ('k=2' if kmin==2 else 'pure'))].append((q,r['equation'],g,rad,c))
for key in sorted(res):
    v=sorted(res[key],reverse=True)
    print(key,"n=%d max q=%.3f  #q>=1: %d"%(len(v),v[0][0],sum(1 for x in v if x[0]>=1)),"|",v[0][1],"g=",v[0][2])
