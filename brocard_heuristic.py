import math
def E(Bmax,ncap=12):
    tot=0.0; by={}
    for k in range(2,40):
        if k==4: continue
        x=1; kk=0.0
        while True:
            b=math.prod(range(x,x+k))
            if b>Bmax: break
            for n in range(2,ncap):
                if k==2 and n==2: continue
                kk+=(1/n)*b**(1/n-1)
            x+=1
            if x>2*10**6: break
        by[k]=kk; tot+=kk
    return tot,by
for B in (1e12,1e16,1e24):
    t,by=E(B); print("Bmax=%g expected ~ %.2f"%(B,t), {k:round(v,2) for k,v in by.items() if v>0.03},flush=True)
