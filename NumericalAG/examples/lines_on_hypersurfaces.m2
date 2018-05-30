needsPackage "MonodromySolver"
        
sparseFamily = method(Options=>{ParameterSymbol=>W})
sparseFamily PolySystem := o -> PS -> (
    polys := flatten entries PS.PolyMap;
    ind := flatten apply(#polys,i-> -- indices for parameters                                              
    	apply(exponents polys#i, t->(i,t))                                                                 
    	);                                                                   
    R := PS.PolyMap.ring;                   
    AR := CC[apply(ind,i->(o.ParameterSymbol)_i)][gens R];                                                                   
    polysP := for i to #polys-1 list -- system with parameteric coefficients and same support              
    sum(exponents polys#i, t->((o.ParameterSymbol)_(i,t))_AR*AR_(t));	
    polySystem transpose matrix {polysP}
    )  

sliceProjectiveFamily = method(Options=>{ParameterSymbol=>H})
sliceProjectiveFamily PolySystem := o -> P -> (
    PS:=specializeSystem(point {toList((numgens coefficientRing ring P):1_CC)},P); 
    c:=P.NumberOfVariables-P.NumberOfPolys;
    R:=CC[gens ring P];
    lForms:=sparseFamily(polySystem apply(apply(c-1,i->random(1,R)) | {random(1,R)-1},p->sub(p,ring P)),ParameterSymbol=>o.ParameterSymbol);
    S:=CC[gens coefficientRing ring lForms | gens coefficientRing ring P][gens ring P];
    polySystem(apply(equations lForms,e->sub(e,S))|apply(equations P,e->sub(e,S)))
    )




end

restart
needs "twentyseven.m2"

n=3
d=2*n-3
R=CC[x_0..x_n]
f=random(d,R)
P=sparseFamily polySystem {f}
peek P
gens ring P
gens coefficientRing ring P

S=coefficientRing ring P
R'=S[toList(p_0..p_n) | toList(q_0..q_n) | {s,t}]
F = map(R',ring P,
            s*matrix{{p_0..p_n}}+
            t*matrix{{q_0..q_n}}
            )
R=S[take(gens R',2*(n+1))]
gens R
cFx=sub(last coefficients(F (first equations P),Variables=>{s,t}),R)
P=polySystem flatten entries sub(cFx,R)
peek P

--slice to get a 0-dimensional family
Q=sliceProjectiveFamily P
peek Q

--solve a random system in this family
setRandomSeed 0
elapsedTime (sys,sols)=solveFamily Q;
-- why doesn't semicolon suppress "failure" message?
#sols


--check that all solutions are ok?
resids=apply(sols,s->sum apply(sys,p->norm evaluate(matrix{{p}},s)));
max resids


pluck = exteriorPower(2,matrix{{p_0..p_n},{q_0..q_n}})
pluckerCoords = apply(sols,p->evaluate(pluck,p));
min(apply(pluckerCoords,p->abs p_(0,0)))
affinePluckerCoords = apply(pluckerCoords,p->(1/p_(0,0))*p);
netList sort apply(pluckerCoords,p->p_(0,1))