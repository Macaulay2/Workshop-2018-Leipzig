--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--    This is draft to build a new package to work on tensors.
--    It is based on the previous file "Tensors.m2" written during
--    	  the Macaulay2 Workshop in Boise (2015).
--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-- PREAMBLE -----------------------------------------------------
-- -*- coding: utf-8 -*-
newPackage(
    "Tensors",
    Version => "0.1",
    Date => "5 June 2018",
    Authors => {}, -- TODO
    Headline => "some tensor constructions",
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

-- EXPORT LIST --------------------------------------------------
export {
    -- Types
    "TensorSpace",
    "Tensor",
    -- methods
    "tensorSpace",
    "makeTensor",
    "kronecker",
    "factorsTensor",
    "mergeTensor",
    "glAction",
    "flattening",
    "symmetrize",
    "slice",
    "contraction",
    -- symbols
    "dims", "coeff", "baseRing", "tensorBasis"
    }

protect dims
protect coeff
protect baseRing
protect tensorBasis

-- DOCUMENTATION ------------------------------------------------
beginDocumentation()
doc ///
  Key
    Tensors
  Headline
     some tensor constructions
  Description
   Text
    {\em Tensors} is work in progress.
///

-- CODE ---------------------------------------------------------

------------------------------------------------------------------------
-- CONSTRUCTIONS OF TENSOR SPACES AND TENSORS
------------------------------------------------------------------------

-- Definition of new types for tensor spaces (TensorSpace) and tensors (Tensor)
TensorSpace = new Type of HashTable
Tensor = new Type of HashTable

-- construction of a TENSOR SPACE. The attributes in input are:
--    R = base ring
--    L = list of lists of variables, one for each vector space forming 
--    	  the tensor space V_1 * ... * V_s
--    D = optional list of integers {d_1,...,d_s}, if the V_i's are 
--    	  the d_i-th homogeneous part of a symmetric (or exterior) algebra;
--    	  Default value: D = {1,...,1}.
--    A = optional list of booleans indicating which one of the V_i's has to be
--    	  an exterior algebra instead of a symmetric algebra;
--    	  Default value: Anti = {false,...,false}
tensorSpace = method();
tensorSpace (Ring, VisibleList, VisibleList) := (R,Xs,N) -> (
    d := #N;
    if d == 0 then (
	Tmod := R;
    	) else (
	Tmod = R[(Xs#0)_(0,0)..(Xs#0)_(0,N_0-1)];
	for i from 1 to d-1 do (
	    Tmod = Tmod ** R[(Xs#i)_(i,0)..(Xs#i)_(i,N_i-1)];
	    );
	);
    new TensorSpace from hashTable{
	baseRing => R,
	dims => N,
	tensorBasis => first entries basis(toList(#N:1),Tmod)
	}
    )
tensorSpace (Ring,Symbol,VisibleList) := (R,X,N) -> tensorSpace(R, #N:X, N)

-- Definition of the way a tensor space 'looks like' when printed in output
expression (TensorSpace) := V -> (
    N := V#dims;
    if #N == 0 then (
        return toString(V#baseRing);
        );
    expr := toString(V#baseRing)|"[";
    for i to #N-2 do (
	expr = expr | toString(N_i) | "x";
	);
    expr = expr | toString(last(N)) | "]";
    return expression expr
)

net (TensorSpace) := V -> net expression V

-- function to construct a TENSOR
makeTensor = method();
makeTensor (VisibleList,TensorSpace) := (L,V) -> (
    if (#L != #(V#tensorBasis)) then (
	return "error: coefficients do not match the dimension"
	);
    new Tensor from hashTable{
	coeff => toList(L) / (i -> sub(i, V#baseRing)),
	tensorSpace => V
	}
    )

expression (Tensor) := T -> (
    Tspace := T#tensorSpace;
    Tcoeff := T#coeff;
    i0 := 0;
    while ((T#coeff)_i0 == 0_(Tspace#baseRing) and i0 < product(Tspace#dims)-1) do (
	    i0 = i0+1;
	    );
    expr := expression toString(Tcoeff_i0 * (Tspace#tensorBasis)_i0);
    for i from i0+1 to #(Tspace.tensorBasis)-1 do (
	if Tcoeff_i != 0_(Tspace.baseRing) then (
	    expr = expression expr + expression (toString(Tcoeff_i * (Tspace#tensorBasis)_i));
	);
    );
    return expression (expr)
)

--orderTensor = method()
--orderTensor (Tensor) := T -> (
--    return #((T#tensorSpace)#dims)
--    )

--orderTensorSpace = method()
--orderTensorSpace (TensorSpace) := V -> (
--    return #(V#dims)
--    )

net (Tensor) := T -> net expression T
Tensor#{Standard,AfterPrint} = T -> (
    << endl;
    << toString(class(T)) | " in " | net(T#tensorSpace)
    << endl;
    )
------------------------------------------------------------------------
-- ALGEBRA OF TENSORS
------------------------------------------------------------------------

-- access to tensor basis
TensorSpace _ Sequence := (V,s) -> (
    N := V#dims;
    if #s != #(V#dims) or  any(toList(0..#s-1), i -> s_i > ((V#dims)_i)-1 or s_i < 0) then (
	return "error: the sequence does not match the format"
	);
    d := #s;
    ind := s#0;
    for i in 1..<#s do ind = ind*(N#i) + s#i;
    ind = ind + 1;
    I := ((ind-1):0) | (1:1) | ((product(N)-ind):0);
    return makeTensor(I,V)
    )
TensorSpace _ ZZ := (V,z) -> (
    s := append((),z);
    return V_s
    )

-- indexed tensor
Tensor _ Sequence := (T,s) -> (
    V := T#tensorSpace;
    N := V#dims;
    if #s != #(V#dims) or  any(toList(0..#s-1), i -> s_i > (V#dims)_i-1 or s_i < 0) then (
	return "error: the sequence does not match the format"
	);
    ind := s#0;
    for i in 1..<#s do ind = ind*(N#i) + s#i;
    return (T#coeff)_ind
    )
Tensor _ ZZ := (T,z) -> (
    s := append((),z);
    return T_s
    )

-- tensor product

TensorSpace ** TensorSpace := (V,W) -> (
    if V#baseRing =!= W#baseRing then (
	return "error: base rings are different"
    );
    N := V#dims | W#dims;
    R := if #(V#dims) == 0 then W#baseRing else if #(W#dims) == 0 then V#baseRing else ring (first V#tensorBasis) ** ring (first W#tensorBasis);
    new TensorSpace from hashTable{
	baseRing => V#baseRing,
	dims => N,
	tensorBasis => first entries basis(toList(#N:1),R)
	}
    )

pickSymbol = method();
pickSymbol (TensorSpace) := V -> (
    d := #(V#dims);
    M := decompose ideal first V#tensorBasis;
    return for i to d-1 list (
	(baseName((M#(d-1-i))_0))#0
	)
    )

factorsTensor = method()
factorsTensor(TensorSpace) := V -> (
    L := pickSymbol V;
    return for i in 0..#(V#dims)-1 list (
	tensorSpace(V#baseRing,L_i,{(V#dims)#i})
	)
    )

--factorsTensor = method()
--factorsTensor(TensorSpace) := V -> (
--	return for i in 0..#(V#dims)-1 list (
--	    new TensorSpace from hashTable{
--		baseRing => V#baseRing,
--		dims => (V#dims)#i,
--		tensorSpaceOrder => 1,
--		tensorSpaceDim => (V#dims)#i,
--		tensorBasis => first entries basis(apply(toList((#(V#dims)):0),j->if j==i then 1 else 0),R)
--		}
--	    )
--	)
--
--	    tensorSpace(V#baseRing,first(V#tensorBasis),{(V#dims)#i}))
--	)

mergeTensor = method()
mergeTensor(TensorSpace) := V -> (
   return tensorSpace(V#baseRing,first pickSymbol V,{product for i in V#dims list i})
   )

kronecker = method()
kronecker(TensorSpace,TensorSpace) := (V,W) -> (
   if V#baseRing =!= W#baseRing then (
	return "error: base rings are different"
	);
  if #(V#dims) =!= #(W#dims) then (
	return "error: the number of factors of the given tensor spaces are not equal"
	);
   L := factorsTensor(V);
   M := factorsTensor(W);
   Z := mergeTensor((L#0)**(M#0));
   for i in 1..#L-1 do Z=Z**mergeTensor((L#i)**(M#i));
   return Z
   )

-- flattening = method()
-- flattening (TensorSpace,List) := (V,L) -> (
--     F := factorsTensor(V);
--     L1 := {};
--     L2 := {};
--     for i in 0..<#F do (
-- 	if member(i,L) then (
-- 	    L1 = append(L1,F#i)
-- 	 ) else (
-- 	    L2 = append(L2,F#i)
-- 	    )
-- 	);
--     A := L1#0;
--     for j from 1 to #L1-1 do A = A**(L1#j);
--     B := L2#0;
--     for j from 1 to #L2-1 do B = B**(L2#j);
--     return mergeTensor(A) ** mergeTensor(B)
--     )

-- flattening (Tensor,List) := (T,L) -> (
--     compL := for i in 1..#(T#tensorSpace.dims) list (if not member(i,L) then i)
--     )

--kronecker = method()
--kronecker(TensorSpace,TensorSpace) := (V,W) -> (
--    if V#baseRing =!= W#baseRing then (
--	return "error: base rings are different"
--    );
--    if #(V#dims) =!= #(W#dims) then (
--	return "error: the number of factors of the given tensor spaces are not equal"
--    );
--	return tensorSpace(V#baseRing,symbol x, for i in 0..#(V#dims)-1 list (((V#dims)#i)*((W#dims)#i)))
--	)

Tensor ** Tensor := (T,U) -> (
    M := flatten for i in T#coeff list for j in U#coeff list i*j;
    R := T#tensorSpace ** U#tensorSpace;
    return makeTensor(M,R)
	)

--kronProduct(Tensor,Tensor) := (T,U) -> (
--    M = flatten for i in T#coeff list for j in U#coeff list i*j;
--    R = T#tensorSpace ** U#tensorSpace;
--    return makeTensor(M,R)
--	)

Tensor ^** ZZ := (T,n) -> (
    if n == 0 then return 1_((T#tensorSpace).baseRing);
    U := T;
    for i from 1 to n-1 do U = U**T;
    U
    )

--- Tensor operations

TensorSpace == TensorSpace := (W,V) -> (
   if  W.baseRing === V.baseRing and W.dims == V.dims then true
   else false
    )
    
Tensor == Tensor := (T,T') -> (
    if T'#tensorSpace == T#tensorSpace and T'#coeff == T#coeff then true
    else false 
    )

Tensor + Tensor := (T,T') -> (
     if not  T'#tensorSpace ==  T#tensorSpace then error "Tensor+Tensor not from the same TensorSpace";
     makeTensor(T#coeff + T'#coeff, T'#tensorSpace)
     )
 
Thing * Tensor := (r,T) -> (
    return makeTensor(sub(r,(T#tensorSpace).baseRing)*(T#coeff), T#tensorSpace)
    )

Tensor * Thing := (T,r) -> (
    return makeTensor(sub(r,(T#tensorSpace).baseRing)*(T#coeff), T#tensorSpace)
    )

- Tensor := T -> (-1)*T
 
 
Tensor - Tensor := (T,T') -> (
     return T + (-T')
     )


--- Group actions 

glAction = method()
glAction (List,Tensor) := (G,T) -> (
    V := T#tensorSpace;
    N := apply(V.dims,i->i-1);
    d := #N;
    coeffT := for J in (d:0)..toSequence(N) list (
	 sum for I in (d:0)..toSequence(N) list (
	     T_I * product for k to d-1 list (
		(G_k)_(J_k,I_k)
		)
	    )
	);
    return makeTensor(coeffT,V)
    )
glAction (Matrix,Tensor) := (G,T) -> (
    d := #(T#tensorSpace.dims);
    GG := toList(d:G);
    return glAction(GG,T)
    )


--- Symmetrize 

symmetrize = method()
symmetrize (Tensor) := (T) -> (
    V := T#tensorSpace;
    N := apply(V#dims,i->i-1);
    d :=  #N; 
    L := for J in (d:0)..toSequence(N) list (
      if toList(J) == sort(toList(J)) then J else continue);
    P := for J in L list permutations(toList(J));
    S := for J in P list set J;
    return sum flatten for O in P list (
	repr = O_0;
	(1/(d!))*(sum for I in O list T_(toSequence I))*(sum for I in toList(set O) list V_(toSequence I))
		    )
		)
	    
IsSymmetric = method()
IsSymmetric (Tensor) := (T) -> (
    if symmetrize(T) == T then true else false
    )


---- Slices and contractions 


slice = method();
slice (Tensor, List) := (T,L) -> (
    left := apply(L, a -> if a === null then 0 else a);
    right := apply(L, 0..<#L, (a,i) -> if a === null then (T#tensorSpace#dims#i)-1 else a);
    dims' := for i to #(T#tensorSpace#dims)-1 list if L#i === null then T#tensorSpace#dims#i else continue;
    Xs := pickSymbol(T#tensorSpace);
    Xs' := for i to #(T#tensorSpace#dims)-1 list if L#i === null then Xs#i else continue;
    V := tensorSpace(T#tensorSpace#baseRing, Xs', dims');
    coeff' := (left..right) / (J -> T_(toSequence(J)));
    return makeTensor(coeff', V);
    )

contraction = method();
contraction (Tensor, List, List) := (T,K,L) -> (
    D := T#tensorSpace#dims;
    KD := apply(K, k->D#k);
    LD := apply(L, k->D#k);
    if KD != LD then error "dimension mismatch";
    f := i-> (
            sliceList := new MutableList from (#D:null);
            scan(#K, j->(sliceList#(K#j) = i#j; sliceList#(L#j) = i#j));
            return slice(T, toList sliceList);
            );
    Tslices := apply((#K:0)..<(toSequence KD), f);
     sum toList Tslices
     );
contraction (Tensor, ZZ, ZZ) := (T,k,l) -> contraction(T,{k},{l});
contraction (Tensor,Tensor,List,List) := (T,U,K,L) -> (
    Td := T#tensorSpace#dims;
    Ud := U#tensorSpace#dims;
    KD := apply(K, k->Td#k);
    LD := apply(L, k->Ud#k);
    if KD != LD then error "dimension mismatch";
    slices := apply((#K:0)..<(toSequence KD), i-> (
            TsliceList := new MutableList from (#Td:null);
            UsliceList := new MutableList from (#Ud:null);
            scan(#K, j->(TsliceList#(K#j) = i#j; UsliceList#(L#j) = i#j));
            Tslice := slice(T, toList TsliceList);
            Uslice := slice(U, toList UsliceList);
            -- if #K == #Td or #K == #Ud then Tslice*Uslice else Tslice**Uslice
            Tslice**Uslice
            ));
     sum toList slices
     );
contraction (Tensor,Tensor,ZZ,ZZ) := (T,U,k,l) -> contract(T,U,{k},{l})


-- TESTS --------------------------------------------------------

TEST ///
    V = tensorSpace(QQ,symbol X,{2,2,2})
    W = tensorSpace(QQ,Y,{3,3,3})
    T1 = makeTensor(1..8,V)
    T2 = makeTensor(1..8,V)
    assert(2*T1 == T1+T2)
    assert(T1 == T2)
    assert(class V_(1,1,1) === Tensor)
    assert(class T1_(0,0,1) === T1#tensorSpace#baseRing)
    V**W

    T3 = slice(T1, {,0,})
    assert(T3#tensorSpace#dims == {2,2})
    assert(T3#coeff == {1,2,5,6})
    T4 = slice(T1, {1,,0})
    assert(T4#coeff == {5,7})

    T5 = contraction(T1, {0}, {1})
    assert(T5#coeff == {8,10})

    W = tensorSpace(QQ,symbol Y,{2,2,2})

    assert(contraction(V_(0,0,0), W_(0,0,0), {0,1,2}, {0,1,2}) == makeTensor({1}, tensorSpace(QQ, {}, {})))
    assert(contraction(V_(0,0,0), W_(0,0,0), {0,1}, {1,2}) == makeTensor({1,0,0,0}, tensorSpace(QQ, {symbol X, symbol Y}, {2,2})))
///

end--------------------------------------------------------------

uninstallPackage "Tensors"
restart
installPackage "Tensors"
check "Tensors"
viewHelp "Tensors"


-- equations of rank 1 (symmetric) tensor P2xP2xP2
ringT = QQ[a_0..a_8,b_0..b_8,c_0..c_8,Z_(0,0,0)..Z_(2,2,2)]
G1 = sub(genericMatrix(QQ[a_0..a_8],3,3),ringT)
G2 = sub(genericMatrix(QQ[b_0..b_8],3,3),ringT)
G3 = sub(genericMatrix(QQ[c_0..c_8],3,3),ringT)
use ringT
P222 = tensorSpace(ringT,symbol X,{3,3,3})
T = P222_(0,0,0)
generic222 = makeTensor(Z_(0,0,0)..Z_(2,2,2), P222)
orbitT = glAction({G1,G2,G3},T)	    	    -- about 80 secs :(
orbitTsym = glAction(G1,T)	      	      	    -- about 80 secs :(
Isym = ideal (generic222 - orbitTsym)#coeff
I = ideal (generic222 - orbitT)#coeff
eliminate(Isym,toList(a_0..a_8))
eliminate(I,toList(a_0..a_8 | b_0..b_8 | c_0..c_8))

