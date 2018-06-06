--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--    This is draft to build a new package to work on tensors.
--    It is based on the previous file "Tensors.m2" written during
--    	  the Macaulay2 Workshop in Boise (2015).
--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
tensorSpace (Ring,Symbol,VisibleList) := (R,X,N) -> (
    d := #N;
    if d == 0 then (
	Tmod := R;
    	) else (
	Tmod = R[X_(0,0)..X_(0,N_0-1)];
	for i from 1 to d-1 do (
	    Tmod = Tmod ** R[X_(i,0)..X_(i,N_i-1)];
	    );
	);
    new TensorSpace from hashTable{
	baseRing => R,
	dims => N,
	tensorBasis => first entries basis(toList(#N:1),Tmod)
	}
    )

-- Definition of the way a tensor space 'looks like' when printed in output
expression (TensorSpace) := V -> (
    N := V#dims;
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
	coeff => toList(L),
	tensorSpace => V
	}
    )

expression (Tensor) := T -> (
    Tspace := T#tensorSpace;
    Tcoeff := T#coeff;
    expr := expression toString(Tcoeff_0*(Tspace#tensorBasis)_0);
    for i from 1 to #(Tspace.tensorBasis)-1 do (
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
------------------------------------------------------------------------
-- ALGEBRA OF TENSORS
------------------------------------------------------------------------

-- access to tensor basis
TensorSpace _ Sequence := (V,s) -> (
    if #s != #(V#dims) or  any(toList(0..#s-1), i -> s_i > ((V#dims)_i)-1 or s_i < 0) then (
	return "error: the sequence does not match the format"
	);
    N := V#dims;
    i := sum for j to #s-2 list s_j * product (for h from j+1 to #s-1 list N_h);
    i = i + last(s);
    I = apply(toList(0..(product N - 1)), j -> if j == i then 1 else 0);    
    return makeTensor(I,V)
    )
TensorSpace _ ZZ := (V,z) -> (
    s := append((),z);
    return V_s
    )

-- indexed tensor
Tensor _ Sequence := (T,s) -> (
    V := T#tensorSpace;
    if #s != #(V#dims) or  any(toList(0..#s-1), i -> s_i > (V#dims)_i-1 or s_i < 0) then (
	return "error: the sequence does not match the format"
	);
    N := V#dims;
    i := sum for j to #s-2 list s_j * product (for h from j+1 to #s-1 list N_h);
    i = i + last(s);
    return (T#coeff)_i
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
    N = V#dims | W#dims;
    R = ring (first V#tensorBasis) ** ring (first W#tensorBasis);
    
 new TensorSpace from hashTable{
	baseRing => V#baseRing,
	dims => N,
	tensorBasis => first entries basis(toList(#N:1),R)
	}
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
     if not  T'#tensorSpace ===  T#tensorSpace then error "Tensor+Tensor not from the same TensorSpace";
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
 
 
 
 
 
 
 
 
 
 
 --  TEST: