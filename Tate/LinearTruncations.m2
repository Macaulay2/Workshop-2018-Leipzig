///
restart
uninstallPackage"LinearTruncations"
restart
installPackage("LinearTruncations")
loadPackage("LinearTruncations",Reload=>true)
viewHelp "LinearTruncations"
peek loadedFiles
check "LinearTruncations"
///


newPackage(
    "LinearTruncations",
    Version => "0.1",
    Date => "June 6, 2018",
    Authors => {
	{ Name => "Navid Nemati", Email => "Navid.Nemati@imj-prg.fr"},
	{Name => "David Eisenbud", Email => "de@msri.org"}},
    Headline => "find the multigraded truncations that give linear resolutions",
    DebuggingMode => true
    )

export {
    "LL",
    "multigradedPolynomialRing",
    "coarseMultigradedRegularity",
    "isLinearComplex",
  --  "partialRegularities",
--    "findOtherLinearTruncations",
    "findAllLinearTruncations",        
    -- Options
    "CoefficientField"
    }

-- Code here

multigradedPolynomialRing = method(Options=>
    {CoefficientField=>ZZ/32003,
    Variables=>getSymbol "x"})

multigradedPolynomialRing(List) := opts -> n -> (
     kk := opts.CoefficientField;
     x:= opts.Variables; -- symbol x;
     t:= #n;
     xx:=flatten apply(t,i->apply(n_i+1,j->x_(i,j)));
     degs:=flatten apply(t,i->apply(n_i+1,k->apply(t,j->if i==j then 1 else 0)));
     kk[xx,Degrees=>degs])
    
multigradedPolynomialRing ZZ := opt -> n -> (multigradedPolynomialRing(toList(n:1)))

--prep for linearTruncations
LL = method()
LL (ZZ,ZZ) := (d,t) -> (
    if t==1 then {{d}} else
    flatten apply(d+1, i->
	apply(LL(d-i,t-1),ell ->flatten prepend(i,ell)))
    )

LL (ZZ,List) := (d,n) -> (
    L1 := LL(d,#n);
    select(L1, ell -> all(#n, i-> ell_i<= (1+n_i)));
    )

findMins = L->(
    if L == {} then return {};
    t := #(L_0);
    u := symbol u;
    P := ZZ/101[u_0..u_(t-1)];
    I := ideal apply(L, ell-> product(t, j-> P_j^(ell_j)));
    apply(flatten entries mingens I, m-> flatten exponents m)
    )

---------------

coarseMultigradedRegularity = method()
coarseMultigradedRegularity Module := M-> 
    coarseMultigradedRegularity res prune M

coarseMultigradedRegularity ChainComplex := F-> (
    t := degreeLength ring F;
   if t<=2 then 
   (
       if sum partialRegularities F <= regularity F 
       then 
       partialRegularities F+{ceiling((regularity F - sum partialRegularities F)/2),ceiling((regularity F - sum partialRegularities F)/2)}
       else
       partialRegularities F
       )
  else
  (
    range := toList(min F..max F-1);
    degsF := apply(range,i -> degrees (F_i));
    lowerbounds := flatten flatten apply(range, i->(
	    apply(degsF_i, d -> apply(LL(i,t), ell -> d-ell))
	    ));
    apply(t, i-> max apply(lowerbounds, ell->ell_i))
    )
    )

partialRegularities = method()
partialRegularities Module := M-> 
    partialRegularities res prune M

partialRegularities ChainComplex := F-> (
    t := degreeLength ring F;
    range := toList(min F..max F-1);
    apply(t, i-> max flatten apply(range, j-> (apply( (degrees F_j), ell-> ell_(i)-j))))
    )
isLinearComplex = method()
isLinearComplex ChainComplex := F->(
    if F == 0 then return true;
    t := degreeLength ring F;
    range := toList(min F..max F-1);
    dF := apply(range, i-> (degrees(F_i))/sum);
    mindF := min(dF_(min F));
    all(range, i->max(dF_i) === mindF+i-min F)
    )

linearTruncations = method(Options =>{Verbose =>false})
linearTruncations(Module):= o->M-> (
    t := degreeLength M;
    F := res prune M;
    r := coarseMultigradedRegularity F;
    d := regularity F;
    L0 := LL(d,t);
    candidates := set {};
    scan(L0, ell ->candidates =  candidates+set toList(ell..r));
    candidates = toList candidates;
    if o.Verbose == true then <<"candidates are"<< candidates<<endl;
    L := select(candidates, ell -> 
	isLinearComplex res truncate(ell,M));
    if o.Verbose == true then << "the candidates with linear truncations are" << L<<endl;
    findMins L
    )
linearTruncations(Module,List) := o-> (M, candidates) ->(
    L := select(candidates, c->(
        isLinearComplex res prune truncate(c, M)));
    findMins L)

isSameWeight = L ->(
        d := sum L_0;
        all(L, ell->sum ell == d))

isInterval = L-> (
    --L is a list of lists, of length 2.
    --the script reports whether they all have the same 
    --total weight, and, if so, whether they form an "interval"
    a := first L_0;
    all(#L, i -> first L_i == i+a)
    )
    
findOtherLinearTruncations = M->(
    --assume M is is gen in non-neg multi-degrees; and t = 2
    r := sum coarseMultigradedRegularity M;
    known := linearTruncations M;
    minfirst := min(known/first);
    minlast := min(known/last);
    candidates1 := apply(minfirst, i-> {i,r-i});
    candidates2 := apply(minlast, i-> {r-i,i});  
    linearTruncations(M,candidates1|candidates2)
    )

isInteresting = method()
isInteresting (Ideal,List) := (I,m) ->(
    --tests whether x^m is in range and not in the ideal
    T := ring I; -- ring of multi-degrees
    if numgens T != #m then 
        error "m should correspond to a monomial in ring I";
    mm =product(#m,i->(T_i)^(m_i));
    if mm % I ==0 then false else true)
T = ZZ/2[u,v]
I = ideal"u2,uv,v3"
m = {0,2}
isInteresting(I,m)

findAllLinearTruncations=method()
findAllLinearTruncations (List,Module) := (range,M) ->(
     u := symbol u;
     t := degreeLength ring M;
     T := ZZ/2[u_0..u_(t-1)];
     I := ideal 0_T;
     L := {};
     linTruncs := {};
     d0 := range_0;
     d1 := range_1;     
     for i from d0 to d1 do(
     L := LL(i,t);
     scan(L, ell -> (
	     if isInteresting(I,ell) and
	      isLinearComplex(res prune truncate(ell,M))
	      then (linTruncs = linTruncs|{ell};
--		  error"";
	      I = I + ideal(product(t,i->(T_i)^(ell_i)))
	      )
		  ));
	      );
	  linTruncs)

///
restart
uninstallPackage "LinearTruncations"
restart
installPackage "LinearTruncations"
debug needsPackage "LinearTruncations"

    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = ideal(x_1^3*x_3, x_2*x_3*x_4, x_3^4*x_4, x_4*x_2^2, x_1*x_4*x_3^3)
    M = S^1/I
    d = regularity M      
    c = sum coarseMultigradedRegularity M

linearTruncations M
findOtherLinearTruncations M
findAllLinearTruncations({d,c},M)
---
S = QQ[x_1..x_6,Degrees=>{{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}}]
I = ideal(x_1*x_4*x_6, x_1*x_3^2, x_3^2*x_4*x_5, x_2^2*x_5^2, x_1*x_4^2*x_5, x_1*x_2^2*x_4)
M = S^1/I

d = regularity M
c = sum coarseMultigradedRegularity M
linearTruncations M
findOtherLinearTruncations M
findAllLinearTruncations({d,c},M)
findAllLinearTruncations({d,c+3},M)

///     

-------------------------

-*
TEST ///
-- test code and assertions here
-- may have as many TEST sections as needed
///
*-

beginDocumentation()
doc ///
Key
  LinearTruncations
Headline
 Find the multigraded truncations that give linear resolutions.
Description
  Text
     Let k be a field and S=k[x_1,\dots ,x_n] be a \mathbb{Z}^r-graded polynomial ring over k.
      Let M be a finitely generated S-module,
      define M_{\geq d} := \oplus _{d'\geq d} M_d' be the truncation of M at d. We define the linear truncations of M to be 
    \lbrace d \in \mathbb{Z}^r \mid M_{\geq d} has a linear resolution \rbrace
   \subset \mathbb{Z}^r. The main purpose of this package is to find linear truncations of a module M. 
        It computes the multigraded truncations that give linear resolutions. 
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = ideal(x_1^3*x_3, x_2*x_3*x_4, x_3^4*x_4, x_4*x_2^2, x_1^3*x_2^3,x_3^3)
    M = S^1/I
    d = regularity M      
    c =  coarseMultigradedRegularity M
    L = findAllLinearTruncations({d,sum c},M)
    for i in L list isLinearComplex (res truncate(i,M))
  Text
    In the above example, the output is the minimal generators of linear truncations.  Although, in some cases, the linear truncations may have a total degree less than a regularity.
  Example
    S= QQ[x,y, Degrees=>{{1,0},{0,1}}];
    d = 2;
    I = ideal(x^d, y^(d+1), x*y);
    M = S^1/I;
    regularity M
    isLinearComplex res truncate({d-1,0},M)	
  Text
   For all d greater than 2, the regularity of M is  equal to d, and the truncation at (d-1,0) has a linear resolution. 
Caveat
  Note that, if the ring is standard $\mathbb{Z}$-graded then the answer is just the regularity. 
///

doc ///
Key
   LL
  (LL,ZZ,ZZ)
  (LL,ZZ,List)
Headline
  t-tuple numbers with sum equal to d
Usage
  LL(d,t)
  LL(d,L)
Inputs
  d: ZZ
    a number.
  t: ZZ
    a number.
  L: List
    a List. 
Outputs
  L: List
    a List of t-tuple pairs.
Description
  Text
    This function takes two numbers (d,t) and gives all the t-tuples pairs where the sum of each pair is equal to d.
  Example
    S = QQ[x,y,z,Degrees=>{{1,0}, {0,1},{0,1}}]
    t = degreeLength S
    d = 5
    LL(d,t)
  Text
    Here it gives all pairs with total degree d in a ring S.
///
///Caveat
  
SeeAlso
  
///
doc ///
Key
   coarseMultigradedRegularity
   (coarseMultigradedRegularity,Module)
   (coarseMultigradedRegularity,ChainComplex)
Headline
  this function gives a multidegree (but probably not minimal) where we get linear resolution after truncating the module.
Usage
  coarseMultigradedRegularity(M)
Inputs
  M: Module
    the module that we want to truncate.
  F: ChainComplex
    the minimal free resolution of the module that we want to truncate.   
Outputs
  L: List
    a list containing only one multidegree in the linear truncations.
Description
  Text
    This code is implemented by using Proposition 1.7 in 
    [David Eisenbud, Daniel Erman, and Frank-Olaf Schreyer. Tate resolutions for products of projective spaces. Acta Math. Vietnam., 40(1):5–36, 2015.]
    except in the bigraded case. In this case, the output is the partial regularities of M. 
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    I = ideal(x_1^3*x_3, x_2*x_3*x_4, x_3^4*x_4, x_4*x_2^2, x_1*x_4*x_3^3);
    M = S^1/I;
    regularity M 
    coarseMultigradedRegularity(S^1/I)
    betti res truncate({2,3}, S^1/I)
  Text
    We see that coarseMultigradedRegularity is not a minimal multidegree that we have linear resolution after truncating.
Caveat
  
SeeAlso
  regularity
///
doc ///
Key
   isLinearComplex
  (isLinearComplex, ChainComplex)
Headline
   isLinearComplex
Usage
  isLinearComplex (F)
Inputs
  F: ChainComplex
    a ChainComplex
Outputs
  L: Boolean
Description
  Text
    This function checks whether the complex is linear or not.
  Example
    S=QQ[x_1..x_4]
    I = ideal(x_1*x_2, x_1*x_3,x_1*x_4, x_2*x_3, x_3*x_4);
    F = res (S^1/I)
    betti F
    isLinearComplex F
    F' = res truncate(2,S^1/I)
    betti F'
    isLinearComplex F'
  Text
   In this example, the ideal I has linear resolution but F is not a linear complex.
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = power(ideal(vars S),4);
    F = res I
    betti F
    isLinearComplex F
    t = degreeLength S
    F' = res truncate({0,3},S^1/I)
    betti F'
    isLinearComplex F'    
  Text
    This function also works in the multigraded case.
Caveat
  
SeeAlso
  betti
///
doc ///
Key
    findAllLinearTruncations
   (findAllLinearTruncations,List,Module)
Headline
 Find all minimal multidegree where we get linear resolution after truncating at a specific range.
Usage
  findAllLinearTruncations(L,M)
Inputs
  M: Module
    the module that we want to truncate.
  L: List
    indicates the range that we are looking for multidegrees.
Outputs
  L: List
    a List of multigraded degree.
Description
  Text
    This function computes all the minimal multidegrees of linear truncations in a specific range.
  Example
    S = QQ[x_1..x_6,Degrees=>{{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}}]
    I = ideal(x_1*x_4*x_6, x_1*x_3^2, x_3^2*x_4*x_5, x_2^2*x_5^2, x_1*x_4^2*x_5, x_1*x_2^2*x_4)
    M = S^1/I;
    d = regularity M
    c = coarseMultigradedRegularity M
    findAllLinearTruncations({regularity M, sum coarseMultigradedRegularity M}, M)
  Text
    In this example, the minimal multidegrees have different total degrees.
  Example
    S = QQ[x_0,x_1,x_2,y_0..y_2,Degrees=>{{1,0},{1,0},{1,0},{0,1},{0,1},{0,1}}]
    I = ideal(x_0^2*y_0^2,x_1^2*y_1^2,x_2^2*y_2^2,(x_0+x_1+x_2)^2*(y_0+y_1+y_2)^2);
    B = intersect(ideal(x_0,x_1,x_2), ideal(y_0, y_1,y_2));
    J = saturate I;
    M = S^1/J;
    d = regularity M
    findAllLinearTruncations({9,10},M)
  Text
    M corresponds to 96 generic points in $\mathbb{P}^2\times \mathbb{P}^2$. The regularity of M is 11 but the truncation at the multidegree (6,4) has linear resolution. 
    Since there is no linear truncations with total degree 9, \{\{3, 7\}, \{4, 6\}, \{6, 4\}, \{7, 3\}\}
    is the set of minimal multidegrees.

Caveat
  By choosing different range the output will be different. 
SeeAlso
  coarseMultigradedRegularity
///
doc ///
Key
    multigradedPolynomialRing
   (multigradedPolynomialRing,List)
   (multigradedPolynomialRing,ZZ)
Headline
 Multigraded PolynomialRing
Usage
  multigradedPolynomialRing(L)
  multigradedPolynomialRing(n)
Inputs
  L: List
    indicates the dimension of the projective spaces that we want to have.
  n: ZZ
    indicates the number of projective spaces.    
Outputs
  S: PolynomialRing
    a multigraded polynomial ring.
Description
  Text
    it gives the coordinated ring of  product of projective spaces.
  Example
    multigradedPolynomialRing({2,3,5})
    multigradedPolynomialRing({2,3,5},CoefficientField => ZZ/5)
    multigradedPolynomialRing 4-- The result is the coordinate ring of product of 4 projective lines.
  Text
    There is an option to choose the coefficient field.
Caveat
    By default, the coefficient field is ZZ/32003 
SeeAlso
    CoefficientField
///
doc ///
Key
    CoefficientField
Headline
    Option for multigradedPolynomialRing
Caveat
    
SeeAlso
    multigradedPolynomialRing
///


--------------------------------------------------------
end--
--------------------------------------------------------
uninstallPackage "LinearTruncations"
restart
installPackage"LinearTruncations"
debug needsPackage "LinearTrunctions"
needsPackage "TateOnProducts"
needsPackage "RandomIdeals"

----------------------- Examples of the interval conjecture--------------
tally apply(10, i->(
I = ran apply(5,i-> 2+random (2+random 7));
M = S^1/I;
L = linearTruncations M;
coarseMultigradedRegularity M;
betti res prune truncate({5,3},M);
{isSameWeight L,isInterval L}
))
--if not isLinearComplex res truncate(r,M) then print toString M
--lin = linearTruncations M
S = multigradedPolynomialRing({2,2}, CoefficientField => ZZ/5)
S' = coefficientRing S[gens S]
M' = coker map(S'^1,,sub(presentation M, S'))
isHomogeneous M'
minimalBetti M'
elapsedTime tally apply(100, i->(
a = sort apply(2+random 8,i->-{1+random 5, 1+random 5});
<<a<<endl;
flush;
I = ideal random(S^1, S^a);
M = S^1/I;
L = linearTruncations M;
{isSameWeight L,isInterval L}
))

-----------------Example of cohomology table----------------
use S
M1 = cokernel matrix {{x_(1,0)^2*x_(1,1), x_(0,0)*x_(0,1)^3, x_(0,0)^2*x_(0,1)*x_(1,1)^2, x_(0,0)*x_(0,1)^2*x_(1,1)^2, x_(0,1)^3*x_(1,0)^2, x_(0,0)^3*x_(1,0)^3}}
M2 = cokernel matrix {{x_(0,0)^2*x_(1,1), x_(0,0)^2*x_(0,1)^2, x_(0,0)^2*x_(1,0)^3, x_(0,0)^3*x_(1,0)^2, x_(0,1)^3*x_(1,0)^2, x_(0,1)^2*x_(1,0)*x_(1,1)^3}}
linearTruncations M1
cohomologyMatrix(M1,{-3,-3},2*{3,3})
linearTruncations M2


--M = S^1/ran{2,3,4,5}
--linearTruncations M

--coarseMultigradedRegularity M

--interesting example, where the degree for linear res is {1,5}
I = ideal(x_(0,1)*x_(0,2),x_(0,2)^2*x_(1,1),
    x_(0,1)*x_(1,0)*x_(1,1)^2,x_(1,0)^4*x_(1,2))
M = S^1/ideal random(S^1, S^{2:{0,-2},2:-{2,3},2:-{2,0}})
d = regularity M
c = sum coarseMultigradedRegularity M
findAllLinearTruncations({d,c},M)
--linearTruncations M


cokernel | x_(0,1)^2x_(1,0) x_(0,0)^2x_(0,1)x_(1,0) x_(1,0)^3x_(1,1)^2 x_(0,1)x_(1,0)^3x_(1,1) x_(0,0)^4x_(0,1) x_(0,1)x_(1,0)^5

-----




----------------------Documentation for LinearTruncation------------ 
doc ///
Key
   linearTruncations
   (linearTruncations,Module)
   (linearTruncations,Module,List)
Headline
  compute the minimal pairs when we get linear resolution after truncation
Usage
  LinearTruncation(M)
  LinearTruncation(M,L)
Inputs
  M: Module
    the module that we want to truncate
  L: List
    a candidate list to chech if we get linear resolution after truncation or not
Outputs
  L: List
    a List of multigraded degree
Description
  Text
    this function computes the minimal multigraded degrees when we get linear resolution after truncating at that multigraded degree. 
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = ideal(x_1^3*x_3, x_2*x_3*x_4, x_3^4*x_4, x_4*x_2^2, x_1*x_4*x_3^3)
    t = degreeLength S
    L = LL(regularity (S^1/I), t )
    --linearTruncations(S^1/I)
    --linearTruncations(S^1/I,L)	  
    --findOtherLinearTruncations (S^1/I)
  Text
    Note that linearTruncation does not give all the minimal pairs.
  Example
    S = QQ[x_1..x_6,Degrees=>{{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}}]
    I = ideal(x_1*x_4*x_6, x_1*x_3^2, x_3^2*x_4*x_5, x_2^2*x_5^2, x_1*x_4^2*x_5, x_1*x_2^2*x_4)
    betti res I
  --  linearTruncations(S^1/I)
  Text
    In this example the total degree of the pairs are different
Caveat
 " ?????"
SeeAlso
  truncate
///




---------------Presentation-------------
restart
needsPackage "LinearTruncations"
---------------LinearTruncations---------------------

--export {
--    "multigradedPolynomialRing",
--    "coarseMultigradedRegularity",
--    "isLinearComplex",
--    "findAllLinearTruncations",        
--    }
--------------------------------------------------
S = multigradedPolynomialRing({1,1,1})-----P^1xP^1xP^1
I = ideal(S_0*S_3*S_5, S_0*S_2^2, S_2^2*S_3*S_4, S_1^2*S_4^2, S_0*S_3^2*S_4, S_0*S_1^2*S_3)
M = S^1/I
d = regularity M
coarseMultigradedRegularity(M)
isLinearComplex(res truncate(coarseMultigradedRegularity(M),M))
c = sum coarseMultigradedRegularity(M)
L = findAllLinearTruncations({d,c},M)
tally for i in L list (isLinearComplex(res truncate(i,M)))
-----------------------------------------------------

S = multigradedPolynomialRing({1,1})-----P^1xP^1
I = ideal(S_2^2*S_3, S_0^2*S_3^3, S_0^2*S_1^2*S_2, S_0^3*S_2*S_3, S_0*S_1*S_3^4, S_0^2*S_1^3*S_3^2)
M = S^1/I
d = regularity M
findAllLinearTruncations({6, 6},M)
coarseMultigradedRegularity(M)
findAllLinearTruncations({6, 11},M)
------------------------------------------------------
S = multigradedPolynomialRing({2,2})--------P^2xP^2
I = ideal (x_(0,0)*x_(0,1)^3, x_(0,0)*x_(0,1)*x_(1,0)*x_(1,1),  x_(1,0)^2*x_(0,1)*x_(1,1)^2, x_(0,0)*x_(1,0)*x_(1,1)^3, x_(0,0)*x_(0,2)^3*x_(1,1), x_(0,0)*x_(0,1)*x_(0,2)^2*x_(1,0)^2)
M = S^1/I
d = regularity M
findAllLinearTruncations({d,d},M)
findAllLinearTruncations({d,d+1},M)
findAllLinearTruncations({d,d+2},M)
c = coarseMultigradedRegularity M
findAllLinearTruncations({d,sum c},M)
viewHelp "LinearTruncations"
-------------------------------------------------An example of  linear Truncation  less than regularity!!!!--------
S = QQ[x,y,Degrees=>{{1,0},{0,1}}]
I = ideal(x^3, x*y,y^7 )
M = S^1/I
regularity M
isLinearComplex (res truncate({2,0},M))




if t<=2 then (
       	      if sum PR F> regularity F then
	      	  PR F
		      else
		  PR F+ {ceiling((regularity F - sum PR F)/2),ceiling((regularity F - sum PR F)/2)}  ) 
