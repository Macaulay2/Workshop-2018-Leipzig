needsPackage "EdgeIdeals"


--The function gives a boolean for the definition of splitting vertex:
--INPUTS : The ideal I, a vertex
--OUTPUT : boolean

isSplittingVertex = ( G , v ) -> (
    return ( degreeVertex (G , v ) == 0 and  
	edges ( inducedGraph ( G , 
		sort toList set vertices ( G ) - set { v } ) ) == {} ) 
    )

--examplse: 
--R=QQ[x_1..x_5]
--genI= x_1*{x_2,x_3,x_4,x_5}
--I=ideal genI
--v=x_1
--isSplittingVertex(I,v)
--false

-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------

--The function gives a boolean for the definition of splitting edges:
--INPUTS : The ideal I, a edge in the form of list of its vertices
--OUTPUT : boolean


isSplittingEdge= (G,e)->(
    N := set neighbors(G,e_0);
    N1 := N + set{e_0};
    M := set neighbors(G,e_1);
    M1 := M + set{e_1};
    return (isSubset(N,M1) or isSubset(M,N1));
	    )

-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------
-------------------------------------------------------------------------

--The function gives a list of spliting edges
--INPUTS : an edge ideal I
--OUTPUT : list

splittingEdges = (G) -> (
    select (edges G, e -> isSplittingEdge(G,e) )
    )

--The function that checks Betti Splitting:

--INPUTS : The ideals J and K in R
--OUTPUT : True if two ideals J and K in the ring R 
--are Betti Splitting fot the ideal I=J+K and False otherwise

isBettiSplitting=(J,K)->(
    R := ring J;
    I := J+K;
    SetI := set flatten entries mingens I;
    SetJ := set flatten entries mingens J;
    SetK := set flatten entries mingens K;
    if  SetJ*SetK===set{} and SetJ+SetK===SetI then (
	a := betti res (I*R^1);
    	b := betti res (J*R^1);
    	c := betti res (K*R^1);
    	d := betti res (intersect(J,K)*R^1);    	    	    	    
    	h := apply(keys d, (i,l,j)->(i+1,l,j));
    	L := unique  (keys a|keys b|keys c|h);
    	A := apply(L,(i,j,l) -> a_(i,j,l)==b_(i,j,l)+c_(i,j,l)+d_(i-1,j,l));
	return unique A == {true};)
    else return false
    )


--The function that gives vertex splitting and edge splitting:

--INPUTS : The ideal I, a vertex or edge u comming from I
--OUTPUT : Ideals J and k where I=J+k and J and K are splliting for I

makeSplitting = (G,v)->( 
    I := edgeIdeal  G
    genI := first entries mingens I;
    genJ := delete (v,first entries mingens ( I + ideal (v) ));
    genK := toList set genI - set genJ;
    J := ideal genJ;
    K := ideal genK;
    return (J,K)    
    )

