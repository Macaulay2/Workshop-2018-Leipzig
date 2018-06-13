restart
needsPackage "Graphs"
needsPackage "MonodromySolver"
G=completeGraph(7)

gaussCC = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

normalPoint = m -> point matrix{apply(m,i->realPart gaussCC() + 0*ii)}

--todo: handle case where v0 not a source vertex
powerEquations = G -> (
    v0:=first vertices G;
    V:=drop(vertices G,1);
    n=#V;
    b=symbol b;
    c=symbol c;
    S=CC[apply(edges G,edge->b_(sort keys edge))|apply(V,i->c_i)];
    x=symbol x;
    y=symbol y;
    R=S[apply(V,i->x_i)|apply(V,i->y_i)];
    X={1}|take(gens R,{0,n-1});
    Y={0}|take(gens R,{n,2*n-1});
    polySystem(apply(V,i->sum apply(toList neighbors(G,i),j->b_(sort {i,j})_R * (X#i_R*Y#j_R-X#j_R*Y#i_R)))
	         | apply(V,i->x_i_R^2+y_i_R^2-c_i_R))
	     )

-- root count for complete (Hamiltonian?) graphs
cptRC = n -> binomial(2*(n-1),n-1)

end
restart
needs "powerNetwork.m2"


--complete graphs: max rc
setRandomSeed 0
for n from 3 to 7 do (
    G=completeGraph n;
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(n-1));--if equal, monodromy worked
    )
