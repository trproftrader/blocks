import sympy as sp
x,t=sp.symbols('x t')
for k in range(2,11):
    f=sp.prod([x+i for i in range(k)])
    R=sp.resultant(sp.diff(f,x), sp.expand(f-t), x)      # polynomial in t whose roots are the critical values
    R=sp.Poly(R,t)
    facs=sp.factor_list(R)[1]
    lin=[(p.as_expr(),m) for p,m in facs if p.degree()==1]
    print("k=%d: rational critical values (with multiplicity):"%k, [(sp.solve(p,t)[0],m) for p,m in lin], " irreducible degs:",[p.degree() for p,m in facs])
