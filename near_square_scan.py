from math import isqrt
def run(k,X):
    best=1e99; out=[]
    for x in range(1,X+1):
        b=1
        for i in range(k): b*=x+i
        s=isqrt(b); d=min(b-s*s,(s+1)**2-b)
        if d==0: continue
        r=d/b**(0.5-1.0/k)
        if r<best: best=r; out.append((x,d,round(r,3)))
    return out
for k,X in [(3,400000),(5,200000),(7,100000)]:
    o=run(k,X); print("k=%d x<=%d running minima of |b-y^2|/b^(1/2-1/k):"%(k,X),o[-4:],flush=True)
