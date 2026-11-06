# Checks Theorem 1(d): for coprime U>V>=1 and an odd prime l with l | U^l - V^l we have l | U-V and
# v_l(U^l-V^l) = v_l(U-V)+1.  Also checks it on the output of gdiv (lines "k=.. g=3 x=.. V=.. U=..").
# usage: python3 check_theorem1d.py [gdiv_output.txt ...]
import sys,re,math,random
def v(p,n):
    c=0
    while n%p==0: n//=p; c+=1
    return c
bad=0
for l in (3,5,7,11):
    for _ in range(20000):
        U=random.randint(2,10**6); V=random.randint(1,U-1)
        if math.gcd(U,V)!=1: continue
        b=U**l-V**l
        if b%l==0 and not ((U-V)%l==0 and v(l,b)==v(l,U-V)+1): bad+=1
print("random violations:",bad)
for fn in sys.argv[1:]:
    for line in open(fn):
        m=re.match(r'k=(\d+) g=3 x=(\d+) V=(\d+) U=(\d+)',line)
        if m:
            k,x,V,U=map(int,m.groups()); b=math.prod(range(x,x+k))
            print(fn,k,x,"ok" if v(3,b)==v(3,U-V)+1>=2 else "VIOLATION")
