import csv,re,math,collections
from math import gcd
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
        return v,s
    m=re.match(r'(\d+)\^(\d+|M)',t)
    b=int(m.group(1)); e=1 if m.group(2)=='M' else int(m.group(2))
    return b**e,(pf(b) if b>1 else set())
def radred(vals,primes):
    r=1
    for p in primes:
        if any(v%p==0 for v in vals): r*=p
    return r
rows=[]
for r in csv.DictReader(open('/home/claude/solutions_N1e16.csv')):
    parts=re.split(r' \+ | = ',r['equation'])
    (A,sa),(B,sb),(C,sc)=[parse(t) for t in parts]
    g=gcd(A,B); A2,B2,C2=A//g,B//g,C//g
    P=sa|sb|sc
    rad=radred([A2,B2,C2],P)
    q=math.log(C2)/math.log(rad) if C2>1 and rad>1 else float('nan')
    # unreduced
    rad0=math.prod(P); q0=math.log(C)/math.log(rad0)
    rows.append((q,q0,r['type'],r['equation'],g,A,B,C))
rows.sort(key=lambda r:-r[0])
print("reduced abc quality q = log(C/g)/log rad(ABC/g^3-part):")
for q,q0,t,eq,g,A,B,C in rows[:14]: print(f"  q={q:.3f} (unreduced {q0:.3f}) g={g}  {t}  {eq}")
by=collections.defaultdict(list)
for q,q0,t,eq,g,A,B,C in rows: by[t].append(q)
print()
for t,v in sorted(by.items()): print(t,"n=%d max q=%.3f #q>1: %d"%(len(v),max(v),sum(1 for x in v if x>1)))
# theta statistic for 2.17: log b / log C vs abc lower bound 1-1/m-1/n
print()
import statistics
th=[];gap=[]
for r in csv.DictReader(open('/home/claude/solutions_N1e16.csv')):
    if r['type']!='2.17' or '1^M' in r['equation']: continue
    mm=re.match(r'b_(\d+)\((\d+)\) \+ (\d+)\^(\d+) = (\d+)\^(\d+)',r['equation'])
    k,x,y,m,z,n=map(int,mm.groups())
    b=math.prod(range(x,x+k)); C=int(r['sum'])
    t=math.log(b)/math.log(C); lb=1-1/m-1/n
    th.append((t,lb,r['equation'])); gap.append(t-lb)
print("2.17: log b_k(x)/log C : min %.3f median %.3f; abc lower bound 1-1/m-1/n: min %.3f median %.3f"%(min(t for t,_,_ in th),statistics.median(t for t,_,_ in th),min(l for _,l,_ in th),statistics.median(l for _,l,_ in th)))
print("smallest theta - bound:",min(gap)," median gap:",statistics.median(gap))
# radical ratio rad(b)/b
rr=[]
for r in csv.DictReader(open('/home/claude/solutions_N1e16.csv')):
    if r['type']!='2.17' or '1^M' in r['equation']: continue
    mm=re.match(r'b_(\d+)\((\d+)\)',r['equation']); k,x=int(mm.group(1)),int(mm.group(2))
    s=set()
    for i in range(k): s|=pf(x+i)
    b=math.prod(range(x,x+k)); rr.append(math.log(math.prod(s))/math.log(b))
print("2.17: log rad(b)/log b : min %.3f median %.3f max %.3f"%(min(rr),statistics.median(rr),max(rr)))
