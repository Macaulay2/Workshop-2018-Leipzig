newPackage(
    "LinearTruncations",
    Version => "0.1",
    Date => "June 6, 2018",
    Authors => {
	{ Name => "", Email => "", HomePage => "" }
	},
    Headline => "A new package",
    DebuggingMode => true
    )

export {
    "LL",
    "linearTruncations",
    "multigradedPolynomialRing",
    "coarseMultigradedRegularity",
    "isLinearComplex",
    "findOtherLinearTruncations",
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
    P := ZZ/101[vars(0..t-1)];
    I := ideal apply(L, ell-> product(t, j-> P_j^(ell_j)));
    apply(flatten entries mingens I, m-> flatten exponents m)
    )

---------------

coarseMultigradedRegularity = method()
coarseMultigradedRegularity Module := M-> 
    coarseMultigradedRegularity res prune M

coarseMultigradedRegularity ChainComplex := F-> (
    t := degreeLength ring F;
    range := toList(min F..max F-1);
    degsF := apply(range,i -> degrees (F_i));
    lowerbounds := flatten flatten apply(range, i->(
	    apply(degsF_i, d -> apply(LL(i,t), ell -> d-ell))
	    ));
    apply(t, i-> max apply(lowerbounds, ell->ell_i))
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
    a number 
  t: ZZ
    a number
  L: List
    a List 
Outputs
  L: List
    a List of t-tuple pairs.
Description
  Text
    This function takes a two numbers (d,t) and gives all the t-tuples paris where the sum of each of them is equal to d.
  Example
    S = QQ[x,y,z,Degrees=>{{1,0}, {0,1},{0,1}}]
    t = degreeLength S
    d = 5
    LL(d,t)
  Text
    it gives all pairs with total degree d in a ring S.
///
///Caveat
  ?????
SeeAlso
  ???????
///

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
    linearTruncations(S^1/I)
    linearTruncations(S^1/I,L)	  
  Text
    Note that linearTruncation does not give all the minimal pairs.
  Example
    S = QQ[x_1..x_6,Degrees=>{{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}}]
    I = ideal(x_1*x_4*x_6, x_1*x_3^2, x_3^2*x_4*x_5, x_2^2*x_5^2, x_1*x_4^2*x_5, x_1*x_2^2*x_4)
    betti res I
    linearTruncations(S^1/I)
  Text
    In this example the total degree of the pairs are different
Caveat
  ?????
SeeAlso
  truncate
///

doc ///
Key
   coarseMultigradedRegularity
  (coarseMultigradedRegularity,Module)
Headline
  this function gives a multidegree where we are sure that truncation at this degree has a linear resolution
Usage
  coarseMultigradedRegularity(M)
Inputs
  M: Module
    the module that we want to truncate
Outputs
  L: List
    a list containing just only one multidegree
Description
  Text
    this function gives a pair but probably non minimal when we truncate we get a linear resolution.
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = ideal(x_1^3*x_3, x_2*x_3*x_4, x_3^4*x_4, x_4*x_2^2, x_1*x_4*x_3^3)
    t = degreeLength S
    linearTruncations(S^1/I)	 
    coarseMultigradedRegularity(S^1/I)
  Text
    we see that coarseMultigradedRegularity is not a minimal multidegree.
Caveat
  ?????
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
    this function check whether the complex is linear or not.
  Example
    S=QQ[x_1..x_4]
    I = ideal(x_1*x_2, x_1*x_3,x_1*x_4, x_2*x_3, x_3*x_4)
    F = res (S^1/I)
    betti F
    isLinearComplex F
    F' = res truncate(2,S^1/I)
    betti F'
    isLinearComplex F'
  Text
   in this example the ideal I has linear resolution but the minimal free resolution of S^1/I is not a linear complex
  Example
    S=QQ[x_1..x_4,Degrees=>{{1,0},{1,0},{0,1},{0,1}}]
    I = power(ideal(vars S),4)
    F = res I
    betti F
    isLinearComplex F
    t = degreeLength S
    F' = res truncate({0,3},S^1/I)
    betti F'
    isLinearComplex F'    
  Text
    the function also works in the multigraded case
Caveat
  ???
SeeAlso
  betti
///

--------------------------------------------------------
end--
--------------------------------------------------------
uninstallPackage "LinearTruncations"
restart

installPackage"LinearTruncations"

needsPackage "TateOnProducts"

needsPackage "LinearTruncations"
needsPackage "RandomIdeals"

(S,E) = productOfProjectiveSpaces{2,2}
S' = coefficientRing S[gens S]
ran = L -> substitute(randomMonomialIdeal(L,S'), S)
I = ran{3,4,5,5,5,6};
M = S^1/I
linearTruncations M
low = -{3,3}
cohomologyMatrix(M,low, -2*low)
regularity M
coarseMultigradedRegularity M


netList apply(10, i->(
I = ran{3,3,7,9};
M = S^1/I;
{coarseMultigradedRegularity M, linearTruncations M}
))

tally apply(10, i->(
I = ran apply(5,i-> 2+random (2+random 7));
M = S^1/I;
L = linearTruncations M
coarseMultigradedRegularity M
betti res prune truncate({5,3},M)
{isSameWeight L,isInterval L}
))


if not isLinearComplex res truncate(r,M) then print toString M
--lin = linearTruncations M

use S
M1 = cokernel matrix {{x_(1,0)^2*x_(1,1), x_(0,0)*x_(0,1)^3, x_(0,0)^2*x_(0,1)*x_(1,1)^2, x_(0,0)*x_(0,1)^2*x_(1,1)^2, x_(0,1)^3*x_(1,0)^2, x_(0,0)^3*x_(1,0)^3}}
M2 = cokernel matrix {{x_(0,0)^2*x_(1,1), x_(0,0)^2*x_(0,1)^2, x_(0,0)^2*x_(1,0)^3, x_(0,0)^3*x_(1,0)^2, x_(0,1)^3*x_(1,0)^2, x_(0,1)^2*x_(1,0)*x_(1,1)^3}}
linearTruncations M1
cohomologyMatrix(M1,{-3,-3},2*{3,3})
linearTruncations M2


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


M = S^1/ran{2,3,4,5}
linearTruncations M

coarseMultigradedRegularity M

--interesting example, where the degree for linear res is {1,5}
I = ideal(x_(0,1)*x_(0,2),x_(0,2)^2*x_(1,1),
    x_(0,1)*x_(1,0)*x_(1,1)^2,x_(1,0)^4*x_(1,2))
M = S^1/ideal random(S^1, S^{2:{0,-2},2:-{2,3},2:-{2,0}})
linearTruncations M
coarseMultigradedRegularity M

cokernel | x_(0,1)^2x_(1,0) x_(0,0)^2x_(0,1)x_(1,0) x_(1,0)^3x_(1,1)^2 x_(0,1)x_(1,0)^3x_(1,1) x_(0,0)^4x_(0,1) x_(0,1)x_(1,0)^5

-----
--possible counterexample to the interval conjecture
S = multigradedPolynomialRing{2,2}	
m = 2
m = random(S^{3:{m,m}}, S^{{1,-1},{-1,1}});
I = minors(2,m);
M = S^1/I;
betti res M
coarseMultigradedRegularity M --{7,7}
linearTruncations M


----- An example where linearTruncation gives pairs with different total degree-----
S = QQ[x_1..x_6,Degrees=>{{1,0,0},{1,0,0},{0,1,0},{0,1,0},{0,0,1},{0,0,1}}]
S'=QQ[gens S]
I = ideal(x_1*x_4*x_6, x_1*x_3^2, x_3^2*x_4*x_5, x_2^2*x_5^2, x_1*x_4^2*x_5, x_1*x_2^2*x_4)
betti res I
linearTruncations(S^1/I)
coarseMultigradedRegularity(S^1/I)

viewHelp linearTruncations

///
restart
uninstallPackage"LinearTruncations"
restart
installPackage"LinearTruncations"
needsPackage "TateOnProducts"
needsPackage "RandomIdeals"
(S,E) = productOfProjectiveSpaces{2,2}
S' = coefficientRing S[gens S]
ran = L -> substitute(randomMonomialIdeal(L,S'), S)
I = ran{3,4,5,5,5,6};
M = S^1/I
linearTruncations M
findOtherLinearTruncations M
known = linearTruncations (M, Verbose => true)
minfirst = min(known/first)
minlast = min(known/last)
linearTruncations
code (linearTruncations,Module,List)
candidates1|candidates2

