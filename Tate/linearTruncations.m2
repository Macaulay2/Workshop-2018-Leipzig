newPackage(
              "linearTruncations",
              Version => "0.1", 
              Date => "",
              Authors => {{Name => "", 
                        Email => "", 
                        HomePage => ""}},
              Headline => "",
              DebuggingMode => false
              )

      export {}

      -- Code here

      beginDocumentation()

      doc ///
      Key
        linearTruncations
      Headline
      Description
        Text
        Example
      Caveat
      SeeAlso
      ///

      doc ///
      Key
      Headline
      Usage
      Inputs
      Outputs
      Consequences
      Description
        Text
        Example
        Code
        Pre
      Caveat
      SeeAlso
      ///

      TEST ///
      -- test code and assertions here
      -- may have as many TEST sections as needed
      ///

end--

resMax = method(Options => options coarseMultigradedRegularity)
resMax ChainComplex := o-> F -> (
    --we assume F starts in homol degree 0.
    el := length F;
    r := degreeLength ring F;
    D := apply((min F..max F-1), i-> unique degrees F_i);
    --replace with D = hashTable
    L := flatten apply(
	length D, i-> apply(
	    D_i, s -> s-toList(r:i)));
    regs := apply(r, p-> max(apply(L, q-> q_p)));
--    d := max(regularity F, sum regs);
    d := regularity F;
    e := d-sum regs;
    e' := floor(e/r);
    f := e-r*e';
    (regs, d, regs + toList(#regs:e') + (toList(f:1)|toList((#regs-f):0)))
    )

resMax1 = method(Options => options coarseMultigradedRegularity)
resMax1 ChainComplex := o-> F -> (
    --we assume F starts in homol degree 0.
    el := length F;
    r := degreeLength ring F;
    D := apply((min F..max F-1), i-> findMins unique degrees F_i);
    --replace with D = hashTable
    L := flatten apply(
	length D, i-> apply(
	    D_i, s -> s-toList(r:i)));
    regs := apply(r, p-> max(apply(L, q-> q_p)));
--    d := max(regularity F, sum regs);
    d := regularity F;
    e := d-sum regs;
    e' := floor(e/r);
    f := e-r*e';
    (regs, d, regs + toList(#regs:e') + (toList(f:1)|toList((#regs-f):0)))
    )

resMax Module := o->M->resMax(res prune M)
resMax1 Module := o->M->resMax(res prune M)

findMins = L->(
    t = #L_0;
    P = ZZ/101[vars(0..t-1)];
    I := ideal apply(L, ell-> product(t, j-> P_j^(ell_j)));
    apply(flatten entries mingens I, m-> flatten exponents m)
    )

coarseSet = M ->(
    (twistreg,d,reg) = resMax M;
    d' = d-sum twistreg;
    d'' = max(d',0);
    apply(d''+1, i-> {twistreg_0+i,twistreg_1+d'-i}))
coarseSet1 = M ->(
    (twistreg,d,reg) = resMax1 M;
    d' = d-sum twistreg;
    d'' = max(d',0);
    apply(d''+1, i-> {twistreg_0+i,twistreg_1+d'-i}))

LL = (d,t) -> (
    if t==1 then {{d}} else
    flatten apply(d+1, i->
	apply(LL(d-i,t-1),ell ->flatten prepend(i,ell)))
    )

findLins= M->(
    t := degreeLength ring M;
    r :=regularity M;
    L := LL(r,t);
    select(L,c->(
	F = res prune truncate(c,S^1/I);
    	all(toList(min F..max F-1),i-> 
	    max apply(degrees F_i, d->sum d) == i+sum c
	    )
	))
)
loadPackage ("TateOnProducts",Reload=>true)
(S,E) = productOfProjectiveSpaces{2,2}	

S' = coefficientRing S[gens S]
loadPackage"randomIdeals"
ran = L -> substitute(randomMonomialIdeal(L,S'), S)


I = ran{2,3,4,5}
--interesting example, where the degree for linear res is {1,5}
I = ideal(x_(0,1)*x_(0,2),x_(0,2)^2*x_(1,1),
    x_(0,1)*x_(1,0)*x_(1,1)^2,x_(1,0)^4*x_(1,2))
I = ideal random(S^1, S^{2:{0,-2},2:-{2,3},2:-{2,0}})
M = S^1/I
findLins M

resMax1 M
coarseSet M
coarseSet1 M

coarseSet(M)
resMax M
regularity M
--LL = toList({0,0}..{6,6})
LL(6,1)    
LL(6,2)    
LL(6,3)
	apply(
L = findMins select(LL,c->(
	F = res prune truncate(c,S^1/I);
    	all(toList(min F..max F-1),i-> 
	    max apply(degrees F_i, d->sum d) == i+sum c
	    )
	))

    
    
G = res truncate({2,4},M)
G.dd_3
F = res M
netList apply(toList(min F..max F-1),i->  unique degrees F_i)


high = {9,9};low = {0,0}
cohomologyMatrix(M,low,high)
cohomologyMatrix(T= cornerComplex(M,low,high), low, high)
cohomologyMatrix (cornerComplex(T,{1,5}), low, high)
cohomologyMatrix (cornerComplex(T,{1,4}), low, high)
cohomologyMatrix (lastQuadrantComplex(T,{1,4}), low, high)
lq = lastQuadrantComplex(T,{1,4})
betti oo
uq = firstQuadrantComplex(T,{1,4})
betti oo

