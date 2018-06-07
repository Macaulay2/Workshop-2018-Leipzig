restart
needsPackage "Graphs"
needsPackage "MonodromySolver"
G=completeGraph(7)

powerEquations = G -> (
    V=drop(vertices G,1);
    n=#V;
    b=symbol b;
    c=symbol c;
    S=CC[apply(edges G,edge->b_(sort keys edge))|apply(V,i->c_i)];
    x=symbol x;
    y=symbol y;
    R=S[apply(V,i->x_i)|apply(V,i->y_i)];
    polySystem(apply(V,i->b_{0,i}_R*y_i_R+sum apply(remove(toList neighbors(G,i),0),j->b_(sort {i,j})_R * (x_i_R*y_j_R-x_j_R*y_i_R)))
	         | apply(V,i->x_i_R^2+y_i_R^2-c_i_R))
	     )
end
restart
needs "powerNetwork.m2"

setRandomSeed 0
for n from 4 to 7 do (
    G=completeGraph n;
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(n-1));
    )

