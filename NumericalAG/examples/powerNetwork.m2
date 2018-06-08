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
    polySystem(apply(V,i->b_{v0,i}_R*y_i_R+sum apply(remove(toList neighbors(G,i),v0),j->b_(sort {i,j})_R * (x_i_R*y_j_R-x_j_R*y_i_R)))
	         | apply(V,i->x_i_R^2+y_i_R^2-c_i_R))
	     )
end
restart
needs "powerNetwork.m2"

G=graph hashTable {0=>{1,2,3,4,5,6,7,8},1=>{0,2},2=>{0,1},
                 3=>{0,4,5},4=>{0,3,5},5=>{0,3,4},
		 6=>{0,7,8},7=>{0,6,8},8=>{0,6,7}}

powerEquations G

--complete graphs
setRandomSeed 0
for n from 5 to 5 do (
    G=completeGraph n;
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    --sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    --print(#sols- length V.PartialSols,2^(n-1));
    )
length V.PartialSols
T=track(V.SpecializedSystem,specializeSystem(normalPoint numgens coefficientRing ring P,P),points V.PartialSols)
isReal first T
apply(T,t->imaginaryPart t)
