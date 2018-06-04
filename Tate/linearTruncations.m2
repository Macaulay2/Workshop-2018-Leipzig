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
restart
loadPackage("TateOnProducts",Reload=>true)
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

LL = method()
LL (ZZ,ZZ) := (d,t) -> (
    if t==1 then {{d}} else
    flatten apply(d+1, i->
	apply(LL(d-i,t-1),ell ->flatten prepend(i,ell)))
    )

LL (ZZ,List) := (d,n) -> (
    L1 = LL(d,#n);
    select(L1, ell -> all(#n, i-> ell_i\leq (1+n_i)));
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

multigradedPolynomialRing = method(Options=>
    {CoefficientField=>ZZ/32003,
    Variables=>{getSymbol "x"}})
multigradedPolynomialRing(List) := opts -> n -> (
     kk := opts.CoefficientField;
     x:= opts.Variables; -- symbol x;
     t:= #n;
     xx:=flatten apply(t,i->apply(n_i+1,j->x_(i,j)));
     degs:=flatten apply(t,i->apply(n_i+1,k->apply(t,j->if i==j then 1 else 0)));
     kk[xx,Degrees=>degs])
    
multigradedPolynomialRing ZZ := opt -> n -> (productOfProjectiveSpaces(toList(n:1)))

linearTruncations = method()
linearTruncations Module := M-> (
    t := degreeLength M;
    F := res prune M;
    range := toList(min F..max F-1);
    degsF := apply(range,i -> degrees (F_i));
    lowerbounds := flatten flatten apply(range, i->(
	    apply(degsF_i, d -> apply(LL(i,t), ell -> d-ell))
	    ));
    m := apply(t, i-> max apply(lowerbounds, ell->ell_i));
    r := regularity F;
    L1 := LL(r,t);
    L2 :=  apply(L1, ell-> toList(ell..m));
    L3 := toList sum(L2, ell->set ell);
    L4 := select(L3,c->(
	    F = res prune truncate(c,M);
    	all(range,i-> 
	    max apply(degrees F_i, d->sum d) == i+sum c
	    )));
    findMins L4
    )

M = S^1/ran{2,3,4,5}
linearTruncations M

linearTruncations(Module,List) := (M, candidates) ->(
    F := res prune M;
    L := select(candidates, c->
    	all(toList(min F..max F-1),i-> 
	    max apply(degrees F_i, d->sum d) == i+sum c
	    )
	);
    findMins L)

(S,E) = productOfProjectiveSpaces{2,2}	

S' = coefficientRing S[gens S]
loadPackage"randomIdeals"
ran = L -> substitute(randomMonomialIdeal(L,S'), S)

I = ideal random(S^1, S^{2:{0,-2},2:-{3,0}})

M = S^1/ran{2,3,4,5}
linearTruncations M
coarseMultigradedRegularity M
--interesting example, where the degree for linear res is {1,5}
I = ideal(x_(0,1)*x_(0,2),x_(0,2)^2*x_(1,1),
    x_(0,1)*x_(1,0)*x_(1,1)^2,x_(1,0)^4*x_(1,2))
M = S^1/ideal random(S^1, S^{2:{0,-2},2:-{2,3},2:-{2,0}})
linearTruncations M
