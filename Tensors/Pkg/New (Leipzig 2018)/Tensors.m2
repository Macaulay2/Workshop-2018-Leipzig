--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--    This is draft to build a new package to work on tensors.
--    It is based on the previous file "Tensors.m2" written during
--    	  the Macaulay2 Workshop in Boise (2015).
--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
tensorSpace (Ring,List,List,List) := (R,L,D,A) -> (
    d := #L;
    Tmod := R[];
    for i to d-1 do (
	Tmod = Tmod ** R[toList(L#i), SkewCommutative => A#i];
	);
    new TensorSpace from hashTable{baseField => R,
	tensorBasis => first entries basis({0}|D,Tmod),
	dims => for i in L list #i,
	degs => D,
	antiSym => A}
    )
tensorSpace (Ring,List,List) := (R,L,D) -> (
    A := toList(#L:false);
    return tensorSpace(R,L,D,A)
)
tensorSpace (Ring,List) := (R,L) -> (
    D := toList(#L:1);
    return tensorSpace(R,L,D)
)

-- Definition of the way a tensor space 'looks like' when printed in output
expression (TensorSpace) := V -> (
    expr := "";
    for i from 0 to #(V#dims)-2 do (
	if (V#antiSym)_i == true then (
	    expr = expr | "ext^"|toString((V#degs)_i)|"("|toString((V#dims)_i)|") * ";
	    ) else (
	    expr = expr | "sym^"|toString((V#degs)_i)|"("|toString((V#dims)_i)|") * ";
	    );
	);
    if last(V#antiSym) == true then (
	    expr = expr | "ext^"|toString(last(V#degs))|"("|toString(last(V#dims))|")";
	    ) else (
	    expr = expr | "sym^"|toString(last(V#degs))|"("|toString(last(V#dims))|")";
	    );
    return expression expr
)

net (TensorSpace) := V -> net expression V

-- function to construct a TENSOR
makeTensor = method();
makeTensor (List,TensorSpace) := (L,V) -> (
    if (#L != #(V#tensorBasis)) then (
	return "error: coefficients do not match the dimension"
	);
    new Tensor from hashTable{
	coeff => L,
	tensorSpace => V
	}
    )

expression (Tensor) := T -> (
    Tspace := T#tensorSpace;
    Tcoeff := T#coeff;
    expr := expression toString(Tcoeff_0*(Tspace#tensorBasis)_0);
    for i from 1 to #(Tspace.tensorBasis)-1 do (
	if Tcoeff_i != 0_(Tspace.baseField) then (
	    expr = expression expr + expression (toString(Tcoeff_i * (Tspace#tensorBasis)_i));
	);
    );
    return expression (expr)
)

--  EXAMPLE:
--    construction of a general tensor in sym^1(CC^2) ** sym^2(CC^2)
S = QQ[a_0..a_6];	       	       -- ring of coefficients
V = tensorSpace(S,{{x,y},{z,t}},{1,2}) -- tensor space
T = makeTensor(toList(a_1..a_6),V)     -- make the tensor
peek V	      	      	      	       -- to look at attributes of V
peek T	      	      	      	       -- to look at attributes of T
T#tensorSpace
V#dims
