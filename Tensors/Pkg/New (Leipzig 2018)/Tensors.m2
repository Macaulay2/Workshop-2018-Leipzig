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
    Date => "June 2018",
    Authors => {
	{Name => "Alessandro Oneto",
	    Email => "alessandro.oneto@upc.edu",
	    HomePage => "https://sites.google.com/view/alessandrooneto/home"
	    },
	{Name => "Emanuele Ventura",
	    Email => "emanueleventura.sw@gmail.com, eventura@math.tamu.edu",
	    HomePage => "https://sites.google.com/site/emanueleventurasw/"
	    },
	{Name => "Iman Bahmani Jafarloo",
	    Email => "ibahmani@unito.it",
	    HomePage => "http://calvino.polito.it/~imanbj/#home"
	    },
	{Name => "Nasrin Altafi",
	    Email => "nasrinar@kth.se"
	    },
        {Name => "Francesco Galuppi",
	    Email => "Francesco.Galuppi@mis.mpg.de",
	    HomePage => "https://www.mis.mpg.de/nlalg/members/francesco-galuppi.html"
	    },
        {Name => "Luca Sodomaco",
	    Email => "luca.sodomaco@unifi.it"
	    },
	{Name => "Markus Wageringel",
	    Email => "markus.wageringel@uni-osnabrueck.de"
	    }
	}, -- TODO
    Headline => "A package on Tensors",
    AuxiliaryFiles => false,
    DebuggingMode => true,
    Reload => true
    )

-- EXPORT LIST --------------------------------------------------
export {
    -- types
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
    "pickSymbol",
    "entriesTensor",
    -- symbols
    "dims", "coeff", "baseRing", "tensorBasis"
    }

protect dims
protect coeff
protect baseRing
protect tensorBasis

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

-- construction of a TENSOR. The attributes in input are:
--    L = a list of coefficients
--    V = a TensorSpace
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
    while ((T#coeff)_i0 == 0_(Tspace#baseRing) and i0 < product(Tspace#dims, identity)-1) do (
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

net (Tensor) := T -> net expression T
Tensor#{Standard,AfterPrint} = T -> (
    << endl;
    << toString(class(T)) | " in " | net(T#tensorSpace)
    << endl;
    )

-- Basic features of Tensors

--orderTensor = method()
--orderTensor (Tensor) := T -> (
--    return #((T#tensorSpace)#dims)
--    )

--orderTensorSpace = method()
--orderTensorSpace (TensorSpace) := V -> (
--    return #(V#dims)
--    )

-- 'entries' returns the list of coefficients in a nested list
entriesTensor = method()
entriesTensor (Tensor) := T -> (
    f := (coeffs, ds) -> (
        if #ds <= 1 then coeffs
        else (
            ds' := drop(ds, 1);
            apply(pack(#coeffs // ds#0, coeffs), cs -> f(cs, ds'))
            )
        );
    f(T#coeff, T#tensorSpace#dims)
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
    I := ((ind-1):0) | (1:1) | ((product(N, identity)-ind):0);
    return makeTensor(I,V)
    )
TensorSpace _ ZZ := (V,z) -> (
    s := append((),z);
    return V_s
    )

-- access to tensor indices
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
    N := toSequence(V#dims) | toSequence(W#dims);
    R := if #(V#dims) == 0 then W#baseRing else if #(W#dims) == 0 then V#baseRing else ring (first V#tensorBasis) ** ring (first W#tensorBasis);
    new TensorSpace from hashTable{
	baseRing => V#baseRing,
	dims => N,
	tensorBasis => first entries basis(toList(#N:1),R)
	}
    )

-- given a tensor space, it returns the list of symbols used to name all the basis of each single space
pickSymbol = method();
pickSymbol (TensorSpace) := V -> (
    d := #(V#dims);
    M := decompose ideal first V#tensorBasis;
    return for i to d-1 list (
	(baseName((M#(d-1-i))_0))#0
	)
    )

-- given a tensor space, it returns the factors of the tensor space
factorsTensor = method()
factorsTensor(TensorSpace) := V -> (
    L := pickSymbol V;
    return for i in 0..#(V#dims)-1 list (
	tensorSpace(V#baseRing,L_i,{(V#dims)#i})
	)
    )

-- identifies a tensor space to a vector space of same dimension
mergeTensor = method()
mergeTensor(TensorSpace) := V -> (
   return tensorSpace(V#baseRing,first pickSymbol V,{product for i in V#dims list i})
   )

-- kronecker product of tensor spaces
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

-- flattening

-- to BE FIXED INDICES
-- flattening = method()
-- flattening (TensorSpace,List) := (V,L) -> (
--        F := factorsTensor(V); 
--        L1 := {};
--        L2 := {};
--        for i in 0..<#F do (
-- 	  if member(i,L) then (
--           L1 = append(L1,F#i)
--       ) else (
--       	    L2 = append(L2,F#i)
-- 	    	    )
--       );
--        A := L1#0;
--        for j from 1 to #L1-1 do A = A**(L1#j);
--        B := L2#0;
--        for j from 1 to #L2-1 do B = B**(L2#j);
--        if member(0,L) then return mergeTensor(A) ** mergeTensor(B)
--        else return mergeTensor(B) ** mergeTensor(A)
--        )

-- flattening (Tensor,List) := (T,L) -> (
-- R := flattening(T#tensorSpace,L);
-- return makeTensor(T#coeff, R)
-- )

flattening = method();
flattening (TensorSpace,List) := (V,L) -> (
    d := V#dims;
    X := pickSymbol(V);
    Lc := toList(0..#d-1);
    for i in L do Lc = delete(i,Lc);
    N := {product for i in L list d_i, product for i in Lc list d_i};
    I := for i in L list (
	    (i,0)..(i,d_i)
	);
    Ic := for i in Lc list (
	    (i,0)..(i,d_i)
	    );
    indL := toList(toSequence(for i in L list (i,0))..toSequence(for i in L list (i,d_i)));
    indLc := toList(toSequence(for i in Lc list (i,0))..toSequence(for i in Lc list (i,d_i)));
    Tmod := V#baseRing[for i in indL list (X_0)_i, for j in indLc list (X_0)_j];
    new TensorSpace from hashTable{
	baseRing => V#baseRing,
	dims => N,
	tensorBasis => flatten for i in indL list for j in indLc list value((X_0)_i)*value((X_0)_j)
	}
    )

-- tensor product of tensors
Tensor ** Tensor := (T,U) -> (
    M := flatten for i in T#coeff list for j in U#coeff list i*j;
    R := T#tensorSpace ** U#tensorSpace;
    return makeTensor(M,R)
	)


kroneckerProduct = method()
-- to BE FIXED INDICES
kroneckerProduct(Tensor,Tensor) := (T,U) -> (
    M := flatten for i in T#coeff list for j in U#coeff list i*j;   
    R := kronecker(T#tensorSpace, U#tensorSpace);
    return makeTensor(M,R)
        )



-- tensor power of a tensor
Tensor ^** ZZ := (T,n) -> (
    if n == 0 then return 1_((T#tensorSpace).baseRing);
    U := T;
    for i from 1 to n-1 do U = U**T;
    U
    )

-- equality between on TensorSpace
TensorSpace == TensorSpace := (W,V) -> (
    -- Why does equality for VisibleList not work?
    if  W.baseRing === V.baseRing and toSequence(W.dims) == toSequence(V.dims) then true
    else false
    )
    
-- equality between on Tensor    
Tensor == Tensor := (T,T') -> (
    if T'#tensorSpace == T#tensorSpace and T'#coeff == T#coeff then true
    else false 
    )

-- sum of tensors
Tensor + Tensor := (T,T') -> (
     if not  T'#tensorSpace ==  T#tensorSpace then error "Tensor+Tensor not from the same TensorSpace";
     makeTensor(T#coeff + T'#coeff, T'#tensorSpace)
     )
 
-- multiplication by a scalar
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
 

-- natural group actions of the general linear groups over tensor spaces
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


-- Symmetrize 

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
	(1/(d!))*(sum for I in O list T_(toSequence I))*(sum for I in toList(set O) list V_(toSequence I))
		    )
		)
	    
	     
IsSymmetric = method()
IsSymmetric (Tensor) := (T) -> (
    if symmetrize(T) == T then true else false
    )


-- Symmetric space 

-- to BE CONTINUED
--symm (TensorSpace) = method()
--symm (TensorSpace) := (V) -> (
   --F:=factorsTensor V;
   --if #(set F) != 1 then return "error: vector spaces need to have the same dimensions ";
    --V := T#tensorSpace;
    --N := apply(V#dims,i->i-1);
    --d :=  #N; 
    --L := for J in (d:0)..toSequence(N) list (
    --if toList(J) == sort(toList(J)) then J else continue);
    --P := for J in L list permutations(toList(J));
    --P':= for J in P list set J  
--)







-- slices and contractions 
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

-- DOCUMENTATION ------------------------------------------------
beginDocumentation()
document{
  Key => Tensors,
  Headline => "a package on tensors",
  
    EM "Tensors", " is work in progress started during the ", EM "Macaulay2 Workshop", 
    " at Max Plank Institute of Leipzig (Germany), June 4th -- 8th, 2018.",
    
    BR{}, BR{},
    BOLD "Overview:",
    
    BR{}, BR{},
    "The goal of this project is to provide a list of functions that might be useful for 
    the every-day life of a researcher on tensors.",
    
    BR{}, BR{},
    "In this package, we defined two new types:",
    UL{
    	LI {TO TensorSpace, " which is an ", TO HashTable, " containing:",
	    UL{
		LI {TO Ring, ": the base ring of coefficients"},
		LI {TO List, ": the list of dimensions of the vector spaces forming the tensor space"},
		LI {TO List, ": the basis of the tensor space"},
		},
	    },
	LI {TO Tensor, "which is an ", TO HashTable, " containing:",
	    UL{
		LI {TO List, ": a list of coefficients"},
		LI {TO TensorSpace, ": the ambient tensor space}"},
		},
	    }
	},

    BOLD "User main functions:",
    
    UL{
	LI {TO tensorSpace, " -- to construct a ", TO TensorSpace},
        LI {TO makeTensor, " -- to construct a ", TO Tensor},
	}
}

doc ///
Key
    TensorSpace
Headline
    the type of tensor spaces
Description
    Text
    	Given vector spaces $V_1,\ldots,V_s$ of dimensions $n_1,\ldots,n_s$, respectively,
	we construct the vector space $V_1\otimes\ldots\otimes V_s$ of dimension 
	$n_1\cdots n_s$.
	
	A TensorSpace is a @TO HashTable@ with attributes:
	
	$\bullet$ @TO baseRing@: @TO Ring@, the ring of coefficients;
	
	$\bullet$ @TO dims@: @TO List@, the ordered list of dimensions of the vector spaces;	
	
	$\bullet$ @TO tensorBasis@: @TO List@, the list of basis elements ordered lexicographically.
    Example
    	tensorSpace(QQ,X,{2,2,2})
	peek oo
SeeAlso
    tensorSpace
///


document{
Key => {makeTensor}, 
Headline => "Element of a tensor space",
Usage => "makeTensor(V, L)",
Inputs => {
    "V" => { "a tensor space"},
    "L" => { "a list of elements in the base ring of V" }
    },
Outputs => {
    " The tensor in V whose coefficients are elements of L with respet to the basis of V"},
EXAMPLE {
    "V = tensorSpace(QQ, symbol X, {2,2,2})",
    "T = makeTensor(toList(1..8), V)"
    }
}

document{
Key => {entriesTensor}, 
Headline => "Entries of a tensor",
Usage => "entriesTensor(T)",
Inputs => {"T" => { "a tensor"} },
Outputs => {
    "A nested list representing the entries of the tensor T."},
EXAMPLE {
    "V = tensorSpace(QQ, symbol X, {2,2,2})",
    "T = makeTensor(toList(1..8),V)",
    "entriesTensor T"
    }
}

--document{
--Key=> {TensorSpace == TensorSpace},
--Headline => "Identity between tensor spaces",
--Usage => "V == W",
--Inputs => {
--    "V" => {"a tensor space"},
--    "W" => {"a tensor space"}
--    },
--Outputs => {
--    "true if V and W have the same base ring and the same list of dimensions"},
--EXAMPLE {
--    "X = tensorSpace(QQ, symbol x, {2,2,2})",
--    "Y = tensorSpace(ZZ, symbol y, {2,2,2})",
--    "Z = tensorSpace(QQ, symbol z, {2,2,3})",
--    "W = tensorSpace(QQ, symbol w, {2,2,2})",
--    "X == Y, X == Z, X == W"
--    }
--}

--document{
--Key=> {Tensor == Tensor},
--Headline => "Identity between tensors",
--Usage => "T1 == T2",
--Inputs => {
--    "T1" => {"a tensor"},
--    "T2" => {"a tensor"}
--    },
--Outputs => {
--    "true if T1 and T2 are in the same tensor space and the lists of coefficients are the same"},
--EXAMPLE {
--    "X = tensorSpace(QQ, symbol x, {2,2})",
--    "T1 = makeTensor(X,{1,2,3,4})",
--    "T2 = makeTensor(X,{1,2,3,5})",
--    "T1 == T1, T1 == T2"
--   }
--}

--document{
--Key=> {TensorSpace ** TensorSpace},
--Headline => "Product of two tensor spaces",
--Usage => "V ** W",
--Inputs => {
--    "V" => {"a tensor space"},
--    "W" => {"a tensor space"}
--    },
--Outputs => {
--    "a new tensor space V ** W, which is the tensor product between V and W, whose base ring is the same as the one of V and W"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "W = tensorSpace(QQ, symbol Y, {3,4,5,6})",
--    "Z = V ** W"    
--    }
--}

document{
Key=> {pickSymbol},
Headline => "symbol for the basis elements of a tensor space",
Usage => "pickSymbol(V)",
Inputs => {
    "V" => {"a tensor space"}
    },
Outputs => {
    "The symbol used to write the basis elements of V"},
EXAMPLE {
    "V = tensorSpace(QQ, symbol X, {2,2})",
    "W = tensorSpace(QQ, symbol Y, {2,2})",
    "pickSymbol(V**W)"   
    }
}
    
document{
Key=> {factorsTensor},
Headline => "factors of a given tensor space",
Usage => "factorsTensor(V)",
Inputs => {
    "V" => {"a tensor space"}
    },
Outputs => {
    "The list of vector spaces that are the factors of V"},
EXAMPLE {
    "V = tensorSpace(QQ, symbol X, {2,2})",
    "L = factorsTensor(V)",
--    "(L#0) ** (L#1) == V"   -- This is not an example, should go in TESTS 
    }
}

document{
Key=> {mergeTensor},
Headline => "",
Usage => "mergeTensor(V)",
Inputs => {
    "V" => {"a tensor space"}
    },
Outputs => {
    ""},
EXAMPLE {
    "V = tensorSpace(QQ, symbol X, {2,2})",
    "W = mergeTensor(V)",
    }
}

--document{
--Key=> {Tensor + Tensor},
--Headline => "sum of tensors",
--Usage => "T1 + T2",
--Inputs => {
--    "T1" => {"a tensor"},
--    "T2" => {"a tensor"}
--    },
--Outputs => {
--    "the sum of T1 and T2 in the tensor space"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T1 = V_(0,0,0)",
--    "T2 = V_(1,1,0)",
--    "T = T1 + T2"
--    }
--}

--document{
--Key=> {Thing * Tensor},
--Headline => "left product of a tensor by a scalar",
--Usage => "a * T",
--Inputs => {
--    "a" => {"an element of the base ring of T"},
--    "T" => {"a tensor"}
--    },
--Outputs => {
--    "left product of T by an element of the base ring of T"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "T' = 2 * T"
--    }
--}

--document{
--Key=> {Tensor * Thing},
--Headline => "right product of a tensor by a scalar",
--Usage => "T * a",
--Inputs => {
--    "T" => {"a tensor"},
--    "a" => {"an element of the base ring of T"}
--    },
--Outputs => {
--    "right product of T by an element of the base ring of T"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "T' = T * 2"
--    }
--}

--document{
--Key=> {-Tensor},
--Headline => "the opposite of a tensor",
--Usage => "-T",
--Inputs => {
--    "T" => {"a tensor"}
--    },
--Outputs => {
--    "the opposite of T"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "T' = -T"
--    }
--}

--document{
--Key=> {Tensor-Tensor},
--Headline => "the difference of two tensors",
--Usage => "T1-T2",
--Inputs => {
--    "T1" => {"a tensor"},
--    "T2" => {"a tensor"}
--    },
--Outputs => {
--    "the difference of T1 and T2"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T1 = V_(0,0,0)",
--    "T2 = V_(0,0,1)",
--    "T = T1 - T2"
--    }
--}

--document{
--Key => {TensorSpace _ Sequence}, 
--Headline => "Extract a basis element from a tensor space",
--Usage => "V _ S ",
--Inputs => {
--    "V" => { "a tensor space"},
--    "S" => { "a sequence of nonnegative integers"}},
--Outputs => {
--    {"The basis element of V corresponding to S." } },
--"This function provides an easier to define a tensor.",
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "T' = makeTensor(V, {1,0,0,0,0,0,0,0})",
--    "T == T' "
--   },
--SeeAlso => makeTensor
--}

--document{
--Key => {Tensor _ Sequence}, 
--Headline => "Extract a coefficient of a tensor",
--Usage => "T _ S ",
--Inputs => {
--    "T" => { "a tensor"},
--    "S" => { "a sequence of nonnegative integers"}},
--Outputs => {
--    {"The coefficient of T with respect to the basis vector corresponding to S." } },
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = makeTensor(V, toList(1..8))",
--    "coeff = T_(1,1,1)",
--    "coeff == 8"
--    },
--SeeAlso => TensorSpace _ Sequence
--}

--document{
--Key => {TensorSpace _ ZZ}, 
--Headline => "Extract a coefficient of a tensor of order one",
--Usage => "T _ S ",
--Inputs => {
--    "T" => { "a tensor"},
--    "S" => { "a nonnegative integers"}},
--Outputs => {
--    {"The coefficient of T with respect to the basis vector corresponding to S." } },
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = makeTensor(V, toList(1..8))",
--    "coeff = T_(1,1,1)",
--    "coeff == 8"
--    },
--SeeAlso => TensorSpace _ Sequence
--}

--document{
--Key=> {glAction (List, Tensor)},
--Headline => "linear action on a tensor",
--Usage => "glAction(L,T)",
--Inputs => {
--    "L" => {"a list of matrices of the appropriate size over the same base ring of T"},
--    "T" => {"a tensor"}
--    },
--Outputs => {
--    "image of a tensor under the action of a list of matrices"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "m1 = matrix{{1,0},{0,1}}",
--    "m2 = matrix{{0,-1},{1,0}}",
--    "m3 = matrix{{1,2},{-3,5}}",
--    "L = {m1, m2, m3}",
--   "T' = glAction(L,T)"
--    },
--SeeAlso => glAction (Matrix, Tensor)
--}

--document{
--Key=> {glAction (Matrix, Tensor)},
--Headline => "linear action on a tensor",
--Usage => "glAction(M,T)",
--Inputs => {
--    "M" => {"a matrix of the appropriate size over the same base ring of T"},
--    "T" => {"a tensor"}
--    },
--Outputs => {
--   "image of a tensor under the action of a matrix"},
--EXAMPLE {
--    "V = tensorSpace(QQ, symbol X, {2,2,2})",
--    "T = V_(0,0,0)",
--    "M = matrix{{1,2},{-3,5}}",
--    "T' = glAction(M,T)"
--    },
--SeeAlso => glAction (List, Tensor)
--}

-- TESTS --------------------------------------------------------

TEST ///
    V = tensorSpace(QQ,symbol X,{2,2,2})
    W = tensorSpace(QQ,symbol Y,3:3)
    T1 = makeTensor(1..8,V)
    T2 = makeTensor(1..8,V)
    assert(2*T1 == T1+T2)
    assert(T1 == T2)
    assert(class V_(1,1,1) === Tensor)
    assert(class W_(1,1,1) === Tensor)
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

    -- contraction agrees with matrix product
    MS1 = tensorSpace(QQ, symbol Z, {2,3})
    MS2 = tensorSpace(QQ, symbol Z, {3,4})
    M1 = makeTensor(1..6, MS1)
    M2 = makeTensor(1..12, MS2)
    assert(matrix entriesTensor contraction(M1, M2, {1}, {0}) == (matrix entriesTensor M1) * (matrix entriesTensor M2))
///


end--------------------------------------------------------------

uninstallPackage "Tensors"
restart
installPackage "Tensors"
check "Tensors"
viewHelp "Tensors"

