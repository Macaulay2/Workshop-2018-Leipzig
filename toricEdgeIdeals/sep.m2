restart
loadPackage "EdgeIdeals"


-- Goal: Writing a fucntion to check if a given pair of ideals (J,K) form a betti splitting for I=J+K

-------------------------------------------------------------------
-- isBettiSplitting function
-------------------------------------------------------------------

-----------------------------------------------------------	
--The function returns true or false depending on whether I=J+K is a betti splitting
--INPUTS : monomial ideals J and K
--OUTPUT : Boolean
-----------------------------------------------------------


isBettiSplitting = method()

isBettiSplitting = (J,K)->(
    I := J+K;
    SetI := set flatten entries mingens I;
    SetJ := set flatten entries mingens J;
    SetK := set flatten entries mingens K;
    if (SetJ*SetK===set{} and SetJ+SetK===SetI)  then (
       	a := betti res (I*R^1);
    	b := betti res (J*R^1);
    	c := betti res (K*R^1);
    	d := betti res (intersect(J,K)*R^1);    	    	    	    
	h := apply(keys d, (i,l,j)->(i+1,l,j));
    	L := unique  (keys a|keys b|keys c|h);
	A := apply(L,(i,j,l) -> a_(i,j,l) == b_(i,j,l)+c_(i,j,l)+d_(i-1,j,l));
	for i from 0 to #A-1 do ( if all(A, l ->  l == true) then return true );
	false);
    false
    )


---------------------------------------------------------
--Tests of IsBettiSplitting function: 
---------------------------------------------------------
---------------------------------------------------------
--Test1:Removing a splitting edge and result is True

R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_4*x_1,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6)
J=ideal(x_1*x_2)
K=ideal(x_2*x_3,x_3*x_4,x_4*x_1,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6)

isBettiSplitting(J,K)


---------------------------------------------------------
---------------------------------------------------------
--Test2:Removing a NON splitting edge and result is False


R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)
J=ideal(x_1*x_6)
K=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_4*x_6,x_5*x_6,x_3*x_5)

isBettiSplitting(J,K)


-- checking the betti numbers to make sure our function is correct

a := peek betti res (I*R^1)
b := peek betti res (J*R^1)
c := peek betti res (K*R^1)
d := peek betti res (intersect(J,K)*R^1)

--problem occurs at betti (3, {5}, 5) of I

---------------------------------------------------------
---------------------------------------------------------
--Test3:Removing a NON splitting edge and result is True


R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)
J=ideal(x_2*x_4)
K=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)

isBettiSplitting(J,K)

---------------------------------------------------------
---------------------------------------------------------
--Test4: Removing a spllitting vertex and the result is True

R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)
J = ideal (x_2*x_3 , x_1*x_3 , x_3*x_4 , x_3*x_5)
K = ideal (x_1*x_2 , x_2*x_4 , x_1*x_5 , x_1*x_6 , x_4*x_6 , x_5*x_6)

isBettiSplitting(J,K)


-----------------------------------------------------------
-----------------------------------------------------------
-- isSplittingVertex function 
-----------------------------------------------------------
-----------------------------------------------------------

-----------------------------------------------------------	
--The function that checks if a given vertex is a splitting vertex in the graph associated to I
--INPUTS : edge ideal I and a vertex a in the form of a list
--OUTPUT : Boolean
-----------------------------------------------------------

isSplittingVertex = method()

isSplittingVertex = ( I , a ) -> (
    G := graph I;
    V := vertices G;
    d := degreeVertex(G,a);
    H := inducedGraph ( G , delete(a,V));
    if (d == 0 or set edges H === set {}) then return false;
    true
    )

---------------------------------------------------------
--Tests of IsSplittingVertex function: 
---------------------------------------------------------
---------------------------------------------------------

-- Test 1: we take the center of a star which is not a splitting vertex

R=QQ[x_1..x_10]

I=ideal(x_1*x_2,x_1*x_3,x_1*x_4,x_1*x_5)

isSplittingVertex(I,x_1) --center->false
isSplittingVertex(I,x_2) --true
isSplittingVertex(I,x_6) --isolated vertex->false
-----------------------------------------------------------
-----------------------------------------------------------
-- Test 2: we generate a random graph to check examples

R=QQ[x_1..x_10]
G=randomGraph(R,5)
I=edgeIdeal G
vertices G
degreeVertex(G,x_3)
isSplittingVertex(I,x_3)

--example random edge ideal I=monomialIdeal (x_1*x_6 , x_2*x_4,x_4*x_6,x_5*x_8,x_5*x_9)

                     

-----------------------------------------------------------
-----------------------------------------------------------
-- isSplittingEdge function 
-----------------------------------------------------------
-----------------------------------------------------------

-----------------------------------------------------------	
--The function that checks if a given edge is a splitting edge in a graph
--INPUTS : edge ideal I and an edge in the form of a list
--OUTPUT : Boolean
-----------------------------------------------------------

isSplittingEdge = method()

isSplittingEdge = (I,e)->(
    G := graph I;
    E := set edges G;
    N := set neighbors(G,e_0);
    N1 := N+set{e_0};
    M := set neighbors(G,e_1);
    M1 := M+set{e_1};
    if not isSubset(set{e},E) then error"second argument must be an edge";
    return (isSubset(N,M1) or isSubset(M,N1))
	    )


---------------------------------------------------------
--Tests of IsSplittingEdge function: 
---------------------------------------------------------
---------------------------------------------------------
--Test1: splitting edge and non-edge

R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_4*x_1,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6)



isSplittingEdge(I,{x_1,x_2}) -- true
isSplittingEdge(I,{x_3,x_5}) -- gives an error when the seocnd argument is not an edge

---------------------------------------------------------
---------------------------------------------------------
--Test2: not a splitting edge and non-edge


R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)


isSplittingEdge(I,{x_1,x_6}) -- false
isSplittingEdge(I,{x_1,x_4}) -- not an edge

-----------------------------------------------------------
-----------------------------------------------------------
-- numSplittingEdges function 
-----------------------------------------------------------
-----------------------------------------------------------


-----------------------------------------------------------	
--The function gives a list of splitting edges
--INPUTS : a graph G
--OUTPUT : list
-----------------------------------------------------------

numSplittingEdges = method()

numSplittingEdges = (G) -> (
    if class G === Graph then (
	return # select(edges G, e -> isSplittingEdge(I,e)))
    else error "the entry is not a graph"
    )

---------------------------------------------------------
--Tests of splittingEdges function: 
---------------------------------------------------------
---------------------------------------------------------
--Test1: all edges are splitting edges

I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_4*x_1,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6)
G= graph I

select(edges G, e -> isSplittingEdge(I,e))

numSplittingEdges(G)


---------------------------------------------------------
---------------------------------------------------------
--Test2: not all the edges are splitting edges


R=QQ[x_1..x_6]
I=ideal(x_1*x_2,x_2*x_3,x_3*x_4,x_1*x_3,x_2*x_4,x_1*x_5,x_1*x_6,x_4*x_6,x_5*x_6,x_3*x_5)
G= graph I

select(edges G, e -> isSplittingEdge(I,e))

numSplittingEdges(G)

