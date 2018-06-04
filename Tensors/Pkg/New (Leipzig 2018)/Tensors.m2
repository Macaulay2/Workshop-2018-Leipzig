--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--    This is draft to build a new package to work on tensors.
--    It is based on the previous file "Tensors.m2" written during
--    	  the Macaulay2 Workshop in Boze (2015).
--~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TensorModule = new Type of HashTable
Tensor = new Type of HashTable

-- construction of tensors
tensorModule = method();
tensorModule (Ring,List,List,List) := (R,L,D,Anti) -> (
    d := #L;
    Tmod := R;
    for i to d-1 do (
	Tmod = Tmod ** R[L#i, SkewCommutative => Anti#i];
	);
    new TensorModule from hashTable{baseField => R,
	tensorBasis => first entries basis(D,Tmod),
	dims => for i in L list #i,
	degs => D,
	antiSym => Anti}
    )
tensorModule (Ring,List,List) := (R,L,D) -> (
    Anti := toList(#L:false);
    return tensorModule(R,L,D,Anti)
)
tensorModule (Ring,List) := (R,L) -> (
    D := toList(#L:1);
    return tensorModule(R,L,D)
)

expression (TensorModule) := V -> (
    expr := "";
    for i from 0 to #(V#dims)-2 do (
	if (V#antiSym)_i == true then (
	    expr = expr | "ext^"|toString((V#degs)_i)|"("|toString(V#baseField)|"^"|toString((V#dims)_i)|") * ";
	    ) else (
	    expr = expr | "sym^"|toString((V#degs)_i)|"("|toString(V#baseField)|"^"|toString((V#dims)_i)|") * ";
	    );
	);
    if last(V#antiSym) == true then (
	    expr = expr | "ext^"|toString(last(V#degs))|"("|toString(V#baseField)|"^"|toString(last(V#dims))|")";
	    ) else (
	    expr = expr | "sym^"|toString(last(V#degs))|"("|toString(V#baseField)|"^"|toString(last(V#dims))|")";
	    );
    return expression expr
)

net (TensorModule) := V -> net expression V

V = tensorModule(QQ,{{x,y},{z,t},{u,v}})

makeTensor = method();
makeTensor (List,TensorModule) := (L,V) -> (
    if (any(L, i -> ring(i) =!= V#baseField)) then (
    	return "error: coefficients are not in the base field"
	);
    if (#L != #(V#tensorBasis)) then (
	return "error: coefficients do not match the dimension"
	);
    new Tensor from hashTable{
	coeff => L,
	tensorModule => V
	}
    )

