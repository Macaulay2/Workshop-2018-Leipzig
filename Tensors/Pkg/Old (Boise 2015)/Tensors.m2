newPackage(
	"Tensors",
    	Version => "0.2.0", 
    	Date => "May 30, 2015",
    	Authors => {
	     {Name => "Andrew Critch", Email => "critch@math.berkeley.edu", HomePage => "http://www.acritch.com/"},
	     {Name => "Claudiu Raicu", Email => "claudiu@math.berkeley.edu", HomePage => "http://math.berkeley.edu/~claudiu/"},
	     {Name => "Hirotachi Abo"},
	     {Name => "Roberto Barrera"},
	     {Name => "Robert Krone"},
	     {Name => "Benjamin Reames"},
	     {Name => "Zach Teitler"}
	     },
	PackageImports => {
	     "PHCpack"
	     },
    	Headline => "tensors",
	AuxiliaryFiles => true,
	DebuggingMode => true 
    	)
 --Macaulay2-1.4/share/Macaulay2/Core/matrix1.m2 
 --needs to replaced for this package to work 

----------------------------------------
--Searchable comment legend:
--a.c. : andrew critch
----------------------------------------
export{
     "Tensor",
     "TensorModule",
     "makeTensor",
     "tensorModule",
     "tensorModuleProduct",
     "tensorDims",
     "genericTensor",
     "genericSymmetricTensor",
     "randomTensor",
     "indexedTensorProduct",
     "einsteinSum",
     "symmetrize",
     "isSymmetric",
     "Antisymmetrize",
     "tensorToPolynomial",
     "tensorToMultilinearForm",
     "polynomialToTensor",
     "multiplicationTensor",
     "matrixMultiplicationTensor",
     "tensorEigenvectors",
     "eigenDiscriminant",
     "tensorEigenvectorsCoordinates",
     "tensorFromSlices",
     "flattenTensor",
     "factorMap",
     "permutationSign",
     "associativeCartesianProduct"
     }

-------------------------
--Symbol methods
-------------------------

gs = getSymbol

cs=coreSymbol = method()     
cs String := nam -> (
     getGlobalSymbol(Core#"private dictionary",nam))

isSymbolic = x -> instance(x,Symbol) or instance(x,IndexedVariable)

--these are currently used for einstein summation,
--which needs to be rewritten

-----------------------
--Error methods
-----------------------
assertInstances=method()
assertInstances (List,Type) := (L,T) -> if not all(L,i->instance(i,T)) then (
     error ("expected list entries to be instances of "|toString(T)|"s"))
assertInstances (List,Type,String) := (L,T,context) -> if not all(L,i->instance(i,T)) then (
     error (context|" expected list entries to be instances of "|toString(T)|"s"))

allInstances = method()
allInstances (VisibleList,List) := (things,types) -> (
     all(things,i->not all(types,j->not instance(i,j)))
     )
allInstances (VisibleList,HashTable) := (things,type) -> (
     allInstances(things,{type})
     )

--------------------------------------------
--Load part 1 (minimize dependence on this)
--------------------------------------------
load "./Tensors/cartesian-list-methods.m2"

inserts=method()
inserts(VisibleList,VisibleList,VisibleList):=(locs,things,host)->(
     if not #locs===#things then error "#locations =!= #things to insert";
     for i in 0..<#locs do host=insert(locs#i,things#i,host);
     host
     )

find = method()
find (Thing,VisibleList) := (x,l) -> (
     position(l,i->i===x))

----------------------------------
--Part 2 of 3:
--Tensors and Tensor Modules
----------------------------------

----------------
--Tensor Modules
----------------
Tensor=new Type of Vector
Tensor.synonym="tensor"
Tensor.cache = new CacheTable
vector Tensor := t -> new Vector from t

TensorModule = new Type of Module
TensorModule.cache = new CacheTable
module TensorModule := M -> M.module
module Module := identity
------
--Using dimensions method previously defined for
--RNLs now for...
--if not class tensorDims === MethodFunction then (

tensorDims =
tensorDimensions = method()
tensorDims Module := M -> {rank ambient M}
tensorDims TensorModule := M -> M#(gs"dimensions")

tensorKeys = method(Dispatch=>Thing)
tensorKeys VisibleList := l ->  toList acp(apply(l,i->0..<i)) 
tensorKeys Tensor := t -> tensorKeys tensorDims t

--Printing TensorModules:
moduleSummary=M->(
     n := rank ambient M;
     if M.?generators then
     if M.?relations then << ", subquotient of " << ambient M
     else << ", submodule of " << ambient M
     else if M.?relations then << ", quotient of " << ambient M
     else if n > 0 then (
	  if not all(degrees M, d -> all(d, zero)) 
	  then << ", degrees " << if degreeLength M === 1 then flatten degrees M else degrees M;
	  );
     )

iftm=
isFreeTensorModule = method()
iftm TensorModule := M -> (
     try M.cache#(gs"isFree") else
     M.cache#(gs"isFree")=all(M#(gs"factors"),isFreeModule)
     )

TensorModule.synonym="tensor module"

net TensorModule := M -> (
     if isFreeTensorModule M then (
	  (net module M)|
	  "{"|(fold(M#(gs"dimensions")/toString,(i,j)->i|"x"|j))|"}"
	  ) else (
	  fold(apply(M#(gs"factors"),net@@module),(i,j)->i|" ** "|j)
	  )
     )

TensorModule#{Standard,AfterPrint} = M -> (
     << endl;				  -- double space
     n := rank ambient M;
     << concatenate(interpreterDepth:"o") << lineNumber << " : "
     << (if isFreeTensorModule M then "Free " else "")
     << ring M
     << "-TensorModule of order "|toString(#M#(gs"dimensions"))|
     ", dimensions "|toString(M#(gs"dimensions"));
     moduleSummary M;
     << endl;
     )

-------------------------------------
--Building tensor modules:
-------------------------------------
tm=
tensorModule = method()

--make a free module into a tensor module:
tensorModule (Ring,List) := (R,dims) -> (
     d:=product dims;
     new TensorModule of Tensor from (
	  new HashTable from (pairs R^d)|{
      	       gs"factors" =>  apply(dims,i->R^i),
     	       gs"dimensions" =>  dims,
	       symbol module => R^d}
     	  )
     )

--make a possibly non-free module into an order 1 tensor module, 
--for tensoring with other such modules to build higher-order
--non-free tensor modules:
tensorModule Module := M -> (
     if not isQuotientModule M then error "tensorModule(Module) expected a free module or quotient module";
      new TensorModule of Tensor from (
       	   new HashTable from (pairs M)|{
		gs"factors" =>  {M},
       	   	gs"dimensions" =>  {rank ambient M},
	        symbol module => M}
	   )
     )
tensorModule TensorModule := identity

--this is conceptually weird if M is not free and #L>1
tensorModule (Module,List) := (M,dims) -> (
     d:=product dims;
     if not rank ambient M == d then error "dimensions do not multiply to the number of entries";
     if not isQuotientModule M then error "tensorModule (Module,List) expected a quotient module";
     new TensorModule of Tensor from (
	   new HashTable from (pairs M)|{
	   	gs"factors" =>  {M},
       	   	gs"dimensions" =>  dims,
	        symbol module => M})
     )

tensorModule Tensor := T -> class T

--perhaps this should instead be
--t-> (classes := ancestors class t;
--     return classes#(position(classes,i->class i===TensorModule))
--     )

fm=--[INTERNAL]
factorModules=method()
factorModules TensorModule := T -> T#(gs"factors")
factorModules Module := M -> {M}

tensorDims Tensor := t -> tensorDims class t
dim (ZZ,Tensor) := (n,T) -> (tensorDims T)#n

--Tensor module from a list of modules to tensor product,
--which themselves may be tensor modules
tmp=
tensorModuleProduct=method(Dispatch=>Thing)
tensorModuleProduct Sequence := fctrs -> tensorModuleProduct toList fctrs
tensorModuleProduct List := fctrs -> (
     assertInstances(fctrs,Module,"tensorModuleProduct(List)");
     dims:=flatten(fctrs/tensorDims);
     f:=flatten(fctrs/factorModules);
     M:=fold(fctrs/module,(i,j)->i**j);
     T:=if all(fctrs,isFreeModule) then TensorModule else TensorModule;
      new T of Tensor from (
	   new HashTable from (pairs M)|{
	   	gs"factors" => f,
       	   	gs"dimensions" => dims,
	        symbol module => M})
     )

----------------------------
--Comparing tensor modules
----------------------------
TensorModule == TensorModule := (M,N) -> (M#(gs"factors") / module)==(N#(gs"factors") / module)

----------------------------
--TensorModule operations
----------------------------
TensorModule^ZZ := (M,n) -> tensorModuleProduct (n:M)
TensorModule**TensorModule := (M,N) -> tensorModuleProduct(M,N)

--permute the factors of a tensor module:
TensorModule @ List := (M,l) -> tensorModuleProduct M#(gs"factors")_l

-----------------------------
--Basic tensor methods
-----------------------------
--Get the ambient module of a tensor
module Tensor := t -> module class t;

--Convert a tensor back into a vector
vector Tensor := t -> new (module t) from t

--Extract an entry of a tensor
--by a multi-index

--fast access without error checking
tensorAccess = method()
tensorAccess (Tensor,Sequence) := (t,s) -> (
     dims := tensorDims t;
     if not #s == #dims then error "dimension mismatch";
     if not all(0..<#s,i->s#i<dims#i) then error "index out of range";
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     t_ind
     )

Tensor _ Sequence := tensorAccess

fta=
fastTensorAccess = method()
fta (Tensor,Sequence) := (t,s) -> (
     dims := tensorDims t;
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     t_ind
     )

------------------------------------------
--Making tensors without RNLs (previously TensorArrays)
------------------------------------------
tensor (TensorModule,Vector) := opts -> (M,v) -> new M from v
tensor (TensorModule,VisibleList) := opts -> (M,l) -> (
     new M from map(M,(ring M)^1,for i in l list {i}))

tensor Matrix := o -> M -> makeTensor entries M

tensor Vector := o -> V -> makeTensor entries V
--
makeTensor=method()
--a.c. fix "new M from" here...
makeTensor (VisibleList,VisibleList):=(dims,ents)->(
     R:=commonRing toList ents;
     M:=tensorModule(R,dims);
     tensor(M,ents)
     )
makeTensor (VisibleList,Function):=(dims,f)->(
     ents:=apply(tensorKeys dims,f);
     makeTensor(dims,ents))

Ring**Tensor := (r,t) -> error "not implemented yet"

Tensor/Function := (t,f) -> tensor(class t,apply(entries t,f))

substitute (Tensor,Thing) := (T,S) -> (
    E := entries T;
    E = apply(E, e->sub(e,S));
    M := tensorModule(ring first E,tensorDims T);
    tensor(M,E)
    )

----------------------------
--Access to basis elements
--by multi-index
----------------------------
TensorModule _ Sequence := (M,s) -> (
     dims := M#(gs"dimensions");
     if not #s == #dims then error "dimension mismatch";
     if not all(0..<#s,i->s#i<dims#i) then error "index out of range";
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     M_ind
     )

------------------------------
--Conversions between Tensors
--and RectangularRestedLists
------------------------------

makeTensor List := L -> (
     if not isrect(L) then error "makeTensor List expected a rectangular nested list";
     dims:=initialDimensions L;
     ents:=ultimate(flatten,L);
     makeTensor(dims,ents)
     )


rnl=
rectangularNestedList=method()
rnl(List,List):=(dims,L) -> (
     if not product dims == #L then error "dimensions mismatch";
     while #dims>1 do (
	  d:=last dims;
	  L = for i in 0..<round(#L/d) list take(L,{i*d,(i+1)*d-1});
	  dims = take(dims,{0,-2+#dims}));
     L)

tensorNet = method()
tensorNet Tensor := T -> (
     dims := tensorDims T;
     if #dims < 3 then return netList rnl(dims,entries T);
     colKeys := tensorKeys(remove(dims,0));
     rowKeys := 0..<dims#0;
     colWidth := j -> j => max apply(rowKeys,i->width net T_((1:i)|j));
     colWidths := hashTable apply(colKeys,colWidth);
     padding := I -> concatenate(colWidths#(remove(I,0)) - (width net T_I):" ");
     padEntry := I -> (net T_I)|(padding I);
     netList rnl(dims,apply(tensorKeys dims,padEntry))
     )

net Tensor := memoize tensorNet;

---------------------------
--Tensor operations
---------------------------
Tensor == Tensor := (v,w) -> (
    class v == class w and entries v == entries w
    )
Tensor + Tensor := (v,w) -> (
     if not class v === class w then error "Tensor+Tensor not from the same TensorModule";
     tensor(class v,(vector v)+(vector w))
     )
Tensor - Tensor := (v,w) -> (
     if not class v === class w then error "Tensor-Tensor not from the same TensorModule";
     tensor(class v,(vector v)-(vector w))
     )
RingElement * Tensor := (r,w) -> (
     if not ring r === ring w then error "RingElement*Tensor not over the same ring";
     tensor(class w,r*(vector w))
     )
Tensor * RingElement := (w,r) -> r*w
- Tensor := w -> (-1)*w
Tensor ** Tensor := (v,w) -> (
     M:=(class v)**(class w);
     tensor(M,(vector v)**(vector w))
     )
Tensor ^ ZZ := (t,n) -> fold(for i in 0..<n list t,(i,j)->i**j)

Tensor ^** ZZ := (T,n) -> (
    if n == 0 then return 1_(ring T);
    U := T;
    for i from 1 to n-1 do U = U**T;
    U
    )

--------------------------------
--Permuting the axes of a tensor
--------------------------------
Tensor @ List := (T,l) -> (
     assertInstances(l,ZZ,"tensor(Tensor,List)");
     dims:=tensorDims T;
     if not set l === set(0..<#dims) then error "
     Tensor @ List expected a permutation of 0..<#d, where
     d is the number of dimensions of the Tensor";
     l':=inversePermutation l;
     M:=(class T)@l;
     inds:=tensorKeys dims_l;
     ents:=apply(inds,i->T_(toSequence i_l'));
     tensor(M,ents)
     )

--------------------------------------
--Turn a free tensor into a function that
--accesses its entries
--------------------------------------
assertFreeTensor=method()
assertFreeTensor Tensor := t -> (
     if not isFreeModule class t then 
     error "expected a tensor in a free tensor module"
     )

tensorFunction = method()
tensorFunction Tensor := t -> (
     f:=method(Dispatch=>Thing);
     f Sequence := s -> t_s;
     f
     )

---------------------
--Tensor slices
---------------------
--use inserts function here!
Tensor_List := (t,l) -> (
     l':=toSequence select(l,i->not i===null);
     if #l'==#l then return t_(toSequence l);
     assertFreeTensor t;
     dims:=tensorDims t;
     blanks:=positions(l,i->i===null);
     odims:=dims_blanks;
     M:=class t;
     M':=tensorModuleProduct((M#(gs"factors"))_blanks);
     keylists:=toList \ tensorKeys odims;
     ents:=toList apply(keylists,i->t_(inserts(blanks,i,l')));
     tensor(M',ents)
     )
 
---------------------
--Contracting, symmetrizing, flattening
--------------------- 
symmetrize = method(Options=>{Antisymmetrize=>false})
symmetrize (Tensor) := o -> T -> symmetrize(T,toList (0..#(tensorDims T)-1),Antisymmetrize=>o.Antisymmetrize)
symmetrize (Tensor,List) := o -> (T,L) -> (
    d := #(tensorDims T);
    R := ring T;
    S := apply(permutations L, p->(
	    ind := new MutableList from (0..d-1);
	    scan(#L, j->(ind#(L#j) = p#j));
	    c := if o.Antisymmetrize then permutationSign p else 1;
	    c*T@(toList ind)
	    ));
    (1_R/((#L)!))*(sum S)
    )

isSymmetric = method(Options=>{Antisymmetrize=>false})
isSymmetric Tensor := o -> T -> (
    D := tensorDims T;
    if not all(#D, i->D#i == D#0) then return false;
    S := permutations(#D);
    all(S, s->(
	    c := if o.Antisymmetrize then permutationSign s else 1;
	    c*T@s == T
	    ))
    )

permutationSign = method()
permutationSign(List) := L -> (
    if (sum apply(#L-1, i-> number(drop(L,i+1),j-> j<=L#i)) % 2) == 0 then 1 else -1
    )

--utility function copies from Tensors/gentensors.m2
getIndex := (R,x) -> (
     M := try monoid R else error "expected a polynomial ring or quotient of one";
     if class x =!= R then error "expected an element of the ring";
     x = try baseName x else error "expected a variable of the ring";
     M.index#x)

genericSymmetricTensor = method(Options=>{Antisymmetrize=>false})
genericSymmetricTensor (Ring,ZZ,List) := o -> (R,i,dims) ->(
    n := dims#0;
    d := #dims;
    k := if o.Antisymmetrize then binomial(n,d) else binomial(n+d-1,d);
    varList := flatten entries genericMatrix(R,R_i,1,k);
    H := new MutableHashTable;
    s := 0;
    ents := for ind in (#dims:0)..<(toSequence dims) list (
	sind := sort toList ind;
	c := if o.Antisymmetrize then permutationSign toList ind else 1;
    	ent := if H#?sind then c*H#sind 
	else if o.Antisymmetrize and #(unique ind) < #ind then 0_R
	else (
	    s = s+1; 
	    c*varList#(s-1)
	    );
	H#(toList ind) = ent;
	ent
	);
    M := tensorModule(R,dims);
    new M from vector ents
    )
genericSymmetricTensor(Ring,List) := o -> (R,dims) -> genericSymmetricTensor(R,0,dims,Antisymmetrize=>o.Antisymmetrize)
genericSymmetricTensor(Ring,RingElement,List) := o -> (R,x,dims) -> 
     genericSymmetricTensor(R,getIndex(R,x),dims,Antisymmetrize=>o.Antisymmetrize)


contract (Tensor,List,List) := (T,K,L) -> (
    D := tensorDims T;
    KD := apply(K, k->D#k);
    LD := apply(L, k->D#k);
    if KD != LD then error "dimension mismatch";
    Tslices := apply((#K:0)..<(toSequence KD), i-> (
	    sliceList := new MutableList from (#D:null);
	    scan(#K, j->(sliceList#(K#j) = i#j; sliceList#(L#j) = i#j));
	    T_(toList sliceList)
	    ));
     sum toList Tslices
     );
contract (Tensor,ZZ,ZZ) := (T,k,l) -> contract(T,{k},{l})
contract (Tensor,Tensor,List,List) := (T,U,K,L) -> (
    Td := tensorDims T;
    Ud := tensorDims U;
    KD := apply(K, k->Td#k);
    LD := apply(L, k->Ud#k);
    if KD != LD then error "dimension mismatch";
    slices := apply((#K:0)..<(toSequence KD), i-> (
	    TsliceList := new MutableList from (#Td:null);
	    UsliceList := new MutableList from (#Ud:null);
	    scan(#K, j->(TsliceList#(K#j) = i#j; UsliceList#(L#j) = i#j));
	    Tslice := T_(toList TsliceList);
	    Uslice := U_(toList UsliceList);
	    if #K == #Td or #K == #Ud then Tslice*Uslice else Tslice**Uslice
	    ));
     sum toList slices
     );
contract (Tensor,Tensor,ZZ,ZZ) := (T,U,k,l) -> contract(T,U,{k},{l})

tensorFromSlices = method()
tensorFromSlices List := S -> (
    Sf := ultimate(flatten,S);
    D1 := initialDimensions S;
    fSf := first Sf;
    D2 := if instance(fSf,Tensor) then tensorDims fSf else {};
    R := ring fSf;
    M := tensorModule(R,D1|D2);
    tensor(M, if instance(fSf,Tensor) then flatten apply(Sf,entries) else Sf)
    )

flattenTensor = method()
flattenTensor (Tensor,List) := (T,L) -> (
   D := tensorDims T;
   LD := apply(L, k->D#k);
   Tslices := apply((#L:0)..<(toSequence LD), i-> (
	   sliceInd := new MutableList from (#D:null);
	   scan(#L, j->(sliceInd#(L#j) = i#j));
	   T_(toList sliceInd)
	   ));
   U := tensorFromSlices toList Tslices;
   d := #(tensorDims U);
   k := min L;
   w := toList apply(d,i->if i < k then i+1 else if i == k then 0 else i);
   U@w
   )

factorMap = method()
factorMap (Tensor,Matrix,ZZ) := (T,M,k) -> (
    D := tensorDims T;
    U := contract(T,tensor M,k,1);
    w := toList apply(#D,i->if i < k then i+1 else if i == k then 0 else i);
    U@w
    )


---------------------
--Tensor to and from Polynomials
--------------------- 
tensorToPolynomial = method()
tensorToPolynomial (Tensor,Symbol) := (T,x) -> (
    R := ring T;
    D := tensorDims T;
    if not all(#D, i->(D#i == D#0)) then error "tensor is not square";
    n := D#0;
    S := R[apply(n,i->x_i)];
    tensorToPolynomial(T,S,S_0)
    );
tensorToPolynomial (Tensor,Ring) := (T,S) -> tensorToPolynomial(T,S,S_0)
tensorToPolynomial (Tensor,Ring,RingElement) := (T,S,x) -> tensorToPolynomial(T,S,getIndex(S,x))
tensorToPolynomial (Tensor,Ring,ZZ) := (T,S,x) -> (
    R := ring T;
    D := tensorDims T;
    n := D#0;
    xTen := makeTensor(take(gens S,{x,x+n-1}));
    U := xTen ^** #D;
    L := toList (0..#D-1);
    contract(sub(T,S),U,L,L)
    );

tensorToMultilinearForm = method()
tensorToMultilinearForm (Tensor,Ring) := (T,S) -> tensorToMultilinearForm(T,S,S_0)
tensorToMultilinearForm (Tensor,Ring,RingElement) := (T,S,x) -> tensorToMultilinearForm(T,S,getIndex(S,x))
tensorToMultilinearForm (Tensor,Ring,ZZ) := (T,S,x) -> (
    D := tensorDims T;
    U := null;
    for n in D do (
	u := makeTensor(take(gens S, {x, x + n -1}));
	if U === null then U = u else U = U**u;
	x = x + n;
	);
    L := toList (0..#D-1);
    contract(sub(T,S),U,L,L)
    );
tensorToMultilinearForm (Tensor,Symbol) := (T,x) -> (
    R := ring T;
    D := tensorDims T;
    varList := flatten apply(#D, i->apply(D#i, j->x_(i,j)));
    S := R[varList];
    tensorToMultilinearForm(T,S)
    )

polynomialToTensor = method()
polynomialToTensor (RingElement ) := (f) -> (
    d := degree f;
    n := numgens ring f;
    M := tensorModule(QQ, toList(d_0:n));
    T := apply(terms f, m -> (
	    c := leadCoefficient m;
	    j := toSequence(flatten apply(#(exponents m)_0, 
		i -> toList(((exponents m)_0)_i:i)
		));
	    c*M_j
	    ));
    symmetrize(sum T)
    )

---------------------
--Eigenvector and signular value tuples
--------------------- 
tensorEigenvectors = method()
tensorEigenvectors (Tensor,ZZ,Ring,RingElement) := (T,k,S,x) -> tensorEigenvectors(T,k,S,getIndex(S,x))
tensorEigenvectors (Tensor,ZZ,Ring,ZZ) := (T,k,S,x) -> (
    R := ring T;
    D := tensorDims T;
    n := D#0;
    xTen := makeTensor(take(gens S,{x,x+n-1}));
    U := xTen^**(#D-1);
    L := toList (0..#D-1);
    V := contract(sub(T,S),U,drop(L,{k,k}),drop(L,-1));
    minors(2, matrix{entries V, entries xTen})
    );
tensorEigenvectors (Tensor,ZZ,Ring) := (T,k,S) -> tensorEigenvectors (T,k,S,S_0)
tensorEigenvectors (Tensor,ZZ,Symbol) := (T,k,x) -> (
    R := ring T;
    n := (tensorDims T)#0;
    S := R[apply(n,i->x_i)];
    tensorEigenvectors(T,k,S,S_0)
    );

SingularVectorTuples := (T,x)-> (
 K := ring T;
 d:=tensorDims(T);
 myVars:=for i from 0 to #d-1 list toList(x_(i,0)..x_(i,d#i-1));
 S:=K[flatten flatten myVars];
 myVars=for i from 0 to #d-1 list toList(x_(i,0)..x_(i,d#i-1));
 f:=tensorToMultilinearForm(T,S);
 L:=for k in 0..#d-1 list for ind in (k,0)..(k,d#k-1) list ind;
 L=sum for i in 0..(#L-1) list minors(2,contract(matrix({myVars#i}),f)||matrix({myVars#i}));
 myVarsId:=for i in myVars list ideal i;
 L=fold(saturate, join({L},myVarsId))
)


eigenDiscriminant = method()
eigenDiscriminant (ZZ,ZZ,Ring) := (n,d,Sa) -> (
    K := coefficientRing Sa;
    x := symbol x;
    vx := toList apply(n,i->x_i);
    vs := gens Sa;
    S := K[vs,vx];
    vx = take(gens S, -n);
    vs = take(gens S, n^d);
    T := genericTensor(S,toList (d:n));
    I := tensorEigenvectors(T,0,S,first vx);
    jj := diff(transpose matrix{vx},gens I);
    singI := minors(n-1,jj)+I;
    J := saturate(singI,ideal vx);
    sub(eliminate(vx,J),Sa)
    )

tensorEigenvectorsCoordinates = method()
tensorEigenvectorsCoordinates (Tensor,ZZ,Symbol) := (T,k,x) -> (
    n := (tensorDims T)#0;
    I := tensorEigenvectors(T,k,x);
    R := ring I;
    S := CC[toSequence entries vars R];
    J := sub(I,S);
    rr := (vars S | matrix{{1_S}})*transpose random(CC^1,CC^(n+1));
    L := J + ideal rr;
    F := first entries gens L;
    solveSystem(F)
    )

---------------------
--Multiplication tensors
--------------------- 
multiplicationTensor = method()
multiplicationTensor Ring := R -> (
    Bmatrix := basis R;
    B := flatten entries Bmatrix;
    K := coefficientRing R;
    V := tensorModule(K, {#B});
    L := for i from 0 to #B-1 list (
	for j from 0 to #B-1 list (
	    pVect := sub(last coefficients(B#i * B#j, Monomials=>Bmatrix), K);
	    pTens := makeTensor flatten entries pVect;
	    V_(1:i) ** V_(1:j) ** pTens
	    )
	);
    sum flatten L
    )


matrixMultiplicationTensor = method()
matrixMultiplicationTensor (Ring,ZZ,ZZ,ZZ) := (R,l,m,n) -> (
    M := tensorModule(R,{l*m,m*n,l*n});
    T := 0_M;
    for ijk in (0,0,0)..<(l,m,n) do (
	(i,j,k) := ijk;
	T = T + M_(i*m + j, j*n + k, i*n + k);
	);
    T
    )

--------
TEST///


///

-------------------
--Tensor marginals
--------------------
marg=
marginalize=method()
marg(Tensor,List) := (T,tosum) -> (
     assertFreeTensor T;
     dims:=tensorDims T;
     n:=#dims;
     if not all(tosum,i->instance(i,ZZ) and i<n) then 
      error "marginalize(Tensor,List) expected a list of integers less than the dimensions of the tensor";
     if #tosum===n then return sum entries T;
     tokeep := toList(0..<n)-set(tosum);
     keepkeys:=tensorKeys dims_tokeep;
     sumkeys:=tensorKeys dims_tosum;
     f := l -> sum apply(sumkeys,i->T_(inserts(tosum,i,l)));
     ents:=toList apply(keepkeys,f);
     M:=tensorModuleProduct((class T)#(gs"factors")_tokeep);
     tensor(M,ents)
     )


TEST ///
R=QQ[x 1]
M=tm(R,{2,2})
N=tm(R,{4})
assert(M==R^4)--they are equal as modules
assert(not M===R^4)
assert(not M==N)
assert(not M===N)
h=new MutableHashTable
h#M=1
h#N==1--unfortunately
///


--a.c. SORT THIS UPWARD:
diff(Tensor,RingElement) := (t,r) -> t/(i->diff(i,r))

---------------------
--Load part 3
---------------------
load "./Tensors/gentensors.m2"
load "./Tensors/indexedtensors.m2"

--

TEST  ///


///


load "./Tensors/tensors-documentation.m2"
end

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

