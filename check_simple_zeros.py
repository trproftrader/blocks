import sympy as sp, numpy as np, math
x=sp.symbols('x')
res=[]
for k in range(2,41):
    f=sp.Poly(sp.prod([x+i for i in range(k)]),x)
    g=f+1
    # number of distinct roots and multiplicities via squarefree factorisation
    sqf=sp.sqf_list(g)[1]
    mult=sorted(m for _,m in sqf for _ in [0])
    # simple zeros = sum of degrees of factors with multiplicity 1
    simple=sum(sp.degree(p,x) for p,m in sqf if m==1)
    res.append((k,simple,[ (sp.degree(p,x),m) for p,m in sqf if m>1]))
for r in res: 
    if r[1]<r[0]: print("k=%d simple zeros=%d  repeated factors(deg,mult)=%s"%r)
print("all k in 2..40 with fewer than 3 simple zeros:",[r[0] for r in res if r[1]<3])
print("k with fewer than k simple zeros (i.e. repeated roots):",[r[0] for r in res if r[1]<r[0]])
# monotonicity of critical-value magnitudes M_j (max |f| on (-j-1,-j)), k up to 60
bad=[]
for k in range(3,61):
    c=np.poly1d(np.poly(-np.arange(k)))
    cr=sorted(c.deriv().roots.real)
    M=[abs(c(t)) for t in cr]   # ordered left to right = j descending
    n=len(M); half=n//2
    inc_right=all(M[i]<M[i+1] for i in range(n//2, n-1)) if n>1 else True
    sym=all(abs(M[i]-M[n-1-i])<=1e-6*M[i] for i in range(n))
    if not (sym and inc_right): bad.append(k)
print("k in 3..60 failing symmetric/outward-increasing check:",bad)
print("k=6 critical values:",[round(float(np.poly1d(np.poly(-np.arange(6)))(t)),4) for t in sorted(np.poly1d(np.poly(-np.arange(6))).deriv().roots.real)])
