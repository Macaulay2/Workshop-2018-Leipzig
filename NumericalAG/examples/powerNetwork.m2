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
end
restart
needs "powerNetwork.m2"

cptRC = n -> binomial(2*(n-1),n-1)


--cut vertex implies all solutions are trivial?
G=graph hashTable {0=>{1,2,3,6},1=>{0,2,3,6},2=>{0,1},
                 3=>{0,1,4,5},4=>{3,5},5=>{3,4},
		 6=>{0,1,7,8},7=>{6,8},8=>{6,7}}
P=powerEquations G
peek P
(V,npaths)=monodromySolve(P,NumberOfNodes=>3)
length V.PartialSols
cptRC #(toList vertices G)

--complete graphs: max rc
setRandomSeed 0
for n from 3 to 7 do (
    G=completeGraph n;
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    --sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    --print(#sols- length V.PartialSols,2^(n-1));
    )

--cycle graphs: max rc
setRandomSeed 0
for n from 3 to 7 do (
    G=cycleGraph n;
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(n-1));
    )

--barbell graphs: cut vertex, 0 nontrivial
setRandomSeed 0
for n from 3 to 3 do (
    G=barbellGraph n;
    v:=#(toList vertices G);
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    print(length V.PartialSols | "\n");
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(v-1));
    )

--star graphs: cut vertex
setRandomSeed 0
for n from 3 to 4 do (
    G=starGraph n;
    v:=#(toList vertices G);
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    print(length V.PartialSols | "\n");
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(v-1));
    )

--path graphs: cut vertex
setRandomSeed 0
for n from 3 to 4 do (
    G=pathGraph n;
    v:=#(toList vertices G);
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    print(length V.PartialSols | "\n");
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(v-1));
    )

--path graphs: cut vertex
setRandomSeed 0
for n from 3 to 4 do (
    G=pathGraph n;
    v:=#(toList vertices G);
    P=powerEquations G;
    assert(P.NumberOfPolys==P.NumberOfVariables);
    (V,npaths)=monodromySolve(P);
    print(length V.PartialSols | "\n");
    sols=solveSystem specializeSystem(point random(CC^1,CC^(numgens coefficientRing ring P)),P);
    print(#sols- length V.PartialSols,2^(v-1));
    )

