

needsPackage "FourierMotzkin";
needsPackage "Polyhedra";
needsPackage "NormalToricVarieties";


newPackage(
     "Chow",
     Version => "0.5",
     Date => "22 November  2009",
     Authors => {{
	       Name => "Diane Maclagan",
       	       Email => "D.Maclagan@warwick.ac.uk",
	       HomePage => "http://www.warwick.ac.uk/staff/D.Maclagan"}},
     Headline => "Chow computations for toric varieties",
     DebuggingMode => true
     )

export { 
     "AA",
     "SR",
     "nefCone",
     "nefCone2",
     "effCone",
     "isContainedCones",
     "nefEqualsIntersection",
     "chowRingBasis",
     "Mob",
     "Lcone",
     "Lconeb",
     "LconevsNef",
     "IntersectionRing"
     }


protect ChowRingBas 
protect ChowRingIdeal

---------------------------------------------------------------------------
-- CODE
---------------------------------------------------------------------------

 needsPackage "FourierMotzkin";
needsPackage "Polyhedra";
needsPackage "NormalToricVarieties";

--Is A contained in B?
--A, B lists
isContained = (A,B)->(
     all(A,i->member(i,B))
     )

--Index of lattice generated by rays of tau in the lattice generated by the 
--rays of the Normal ToricVariety X 
latticeIndex = (tau, X) ->(
      A:=matrix rays X;
      B:=A^tau;
      brank:=rank B;
      if brank < rank A then (
	   K:=gens kernel B;
	   C:=gens kernel transpose(A*K);
	   A=transpose(C)*A;
      );
      if not((rank A)==brank) then error("Something's wrong here");
      a:=(flatten entries mingens minors(brank,A))_0;
      b:=(flatten entries mingens minors(brank,B))_0;
      return(lift(b/a,ZZ));
)     
     
     
--Greg says this has been superseded - see orbits
--cones (ZZ, NormalToricVariety) :=  List => (i,X) ->(
-- cones (ZZ, NormalToricVariety) :=  (i,X) ->(
--      if not X.cache.?cones then X.cache.cones = new MutableHashTable;
--      if not X.cache.cones#?i then (
-- 	  if isSimplicial(X) then (
-- 	      X.cache.cones#i = unique flatten for sigma in max X  list subsets(sigma,i);
-- 	  )
--      	  else (     
-- 	       F:=fan X;
-- 	       Conesi:=cones(i,F);
-- 	       X.cache.cones#i=apply(Conesi,C->(
-- 		    RaysC := entries transpose rays C;
-- --This next bit to be changed when Rene fixes Polyhedra
-- RaysC2:=apply(RaysC,i->(apply(i,j->(lift(j,ZZ)))));
-- RaysC=RaysC2;		    
--      	       	    P:=positions(rays X, i->member(i,RaysC));
-- 		    P
-- 	       ));
--      	  );
--      );
--      X.cache.cones#i=sort X.cache.cones#i;
--      return X.cache.cones#i;
-- ) 




--Is cone with generators given by columns of M contained in cone generated
--by columns of N?
--Temporarily assuming full dimensional
--(ie this is a hack)
isContainedCones= (M,N) ->(
          Nfacets:=(fourierMotzkin(N))#0;
	  nonneg := flatten entries ((transpose M)* Nfacets);
	  if max(nonneg)>0 then return(false) else return(true);
);     

--Intersect the cones given by the columns of M and the columns of N
--Temporarily assuming full dimensional
--(ie this is a hack)
intersectCones=(M,N)-> (
        Nfacets:=(fourierMotzkin(N))#0;
        Mfacets:=(fourierMotzkin(M))#0;
	return( (fourierMotzkin(Nfacets | Mfacets))#0);
);	

--Chow
-- i is codim
AA=(i,X) -> ( 
     if i>dim(X) then error("i > dim(X)");
     if not(X.cache.?Chow) then X.cache.Chow = new MutableHashTable;
     if not(X.cache.Chow#?i) then (
     n:=dim X;
     -- Get the faces of dim i, i-1
     sigmaCodimi := orbits(X,n-i);
     tauCodimiminus1 := orbits(X,n-i+1);
     if #tauCodimiminus1 > 0 then (
         --Create the relations (Fulton-Sturmfels eqtn 1, p337)
     Relns:=apply(tauCodimiminus1, tau -> (
     	  Mtau:= entries transpose gens kernel (matrix rays X)^tau;	  
	  TauRelns:=apply(Mtau, u->(
     	       reln:= apply(sigmaCodimi,sigma->(
     	       	    relnsigma:=0;
		    if isContained(tau,sigma) then (
     	       	    	  j:=position(sigma, k->(not(member(k,tau))));
			  nvect:=(rays X)#(sigma#j);
      	       	    	  udotn:=0;
		      	  for k from 0 to #u-1 do
			       udotn=udotn+(u#k)*(nvect#k);
 		          nsigmamult:=latticeIndex(append(tau,sigma#j) ,X) // latticeIndex(tau,X);
			  relnsigma=udotn // nsigmamult;
		    )
	       	    else (	 
		      relnsigma=0;
		    );
	       	    relnsigma
	       ));
	       reln
     	  ));
     	  TauRelns
     ));
     Relns=flatten Relns;
     X.cache.Chow#i = prune coker transpose matrix Relns;
     )
     else X.cache.Chow#i = ZZ^(#sigmaCodimi);
     );
     X.cache.Chow#i
);	 
	 

--isCartier = D -> (
     
--)     


-- intersect (List, ZZ, List, NormalToricVariety) := List => (D, k, tau, X)  -> (
--      if isCartier(D) then (
--      n:=dim X;
--      Zorder := cones(n-k,X);
--      outputCones := cones(n-k+1, X);
--      --First rewrite D so that it is not supported on V(tau)
--      --????
--      --Then do the intersection 
--      DdotTau:=apply(outputCones, sigma -> (
--      	  if not(isContained(tau,sigma)) then 0
-- 	  else (
     	       	            	       	       
--      	  )
--      ));
--      return(DdotTau);
--      )
--      else (
-- 	  <<"D is not a Cartier divisor"<<endl;
-- --???should error trap properly
--      );
-- );     

--Intersect V(sigma) and V(tau) using the SR formulation.
--Not sure about the use of this.
-- intersect ( List, List, NormalToricVariety, ) := MutableHashTable =>( sigma, tau, X) -> (
--       if not isSimplicial X then error("Not implemented yet");
--       --We'll turn sigma into a product of torus-invariant divisors
--       --and do the intersection one-by-one
--       I:=SR(X);
--       R:=ring I;
--       m:=1_R;
--       for i in sigma do (
-- 	   m=m*R_i;
--       );
--       for i in tau do (
-- 	   m=m*R_i;
--       );
--       rem:= m % I;
--       rem   		   
-- );      





--Create SR ideal
SR=X->(
     if not X.cache.?ChowRingIdeal then (
	  z:=symbol z;
     	  R:=QQ[z_1..z_#(rays X)];
       	  I:= ideal apply(max X, sigma->(
	       	    mono:=1_R;
	       	    for j from 0 to #(rays X)-1 do 
		        if not(member(j,sigma)) then mono=mono*R_j;
	       	    mono
		    ));
     	  squaresIdeal:=ideal apply(gens R, xx->xx^2);       
     	  I=ideal flatten entries ((gens (squaresIdeal : I)) % squaresIdeal);
     	  I=I+ ideal apply(transpose rays X, a->(
	       genJ:=0_R;
	       for j from 0 to #a-1 do (
		    genJ=genJ+a#j*R_j;
	       );
	       genJ    
     	  ));
     X.cache.ChowRingIdeal=ideal mingens I;
     );
     X.cache.ChowRingIdeal
);          	       	    


--Create SR ideal
intersectionRing=X->(
     if not X.cache.?IntersectionRing then (
	  z:=symbol z;
     	  R:=QQ[z_1..z_#(rays X)];
       	  I:= ideal apply(max X, sigma->(
	       	    mono:=1_R;
	       	    for j from 0 to #(rays X)-1 do 
		        if not(member(j,sigma)) then mono=mono*R_j;
	       	    mono
		    ));
     	  squaresIdeal:=ideal apply(gens R, xx->xx^2);       
     	  I=ideal flatten entries ((gens (squaresIdeal : I)) % squaresIdeal);
     	  I=I+ ideal apply(transpose rays X, a->(
	       genJ:=0_R;
	       for j from 0 to #a-1 do (
		    genJ=genJ+a#j*R_j;
	       );
	       genJ    
     	  ));
     X.cache.IntersectionRing=R/(ideal mingens I);
     );
     X.cache.IntersectionRing
);          	       	    


--Compute a basis for the Chow ring
chowRingBasis=(X,i)->(
     if not X.cache.?ChowRingBas then
     	  X.cache.ChowRingBas = new MutableHashTable;
     I:=SR(X);
     R:=ring I;
     if not X.cache.ChowRingBas#?i then 	  
          X.cache.ChowRingBas#i=flatten entries lift(basis(dim X -i,R/I),R);
     return(X.cache.ChowRingBas#i);
);


--Code to compute the cone of nef cycles

--Currently returns a rather arbitrary  basis for the ith Chow group 
-- and then a matrix whose columns represent elements there
--generating the cone of nef i-cycles
--(Caveat: this is dimension i, so codimension n-i)

nefCone=(i,X)->(
     if not isSmooth(X) then error("Not implemented yet");
     n:=dim X;
     Conesi:=orbits(X,n-i);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     --Now create the multiplication map
     --First get a basis for AA_i
     chowBas:=chowRingBasis(X,i);
     mono:=1_R;
     for i in (max X)_0 do mono=mono*R_i;
     topBas1:=mono % I;
     Mat:=matrix unique apply(chowBas,m->(
	       apply(Conesi,sigma->(
			 mono:=1_R;
			 for j in sigma do mono=mono*R_j;
     	       	    	 --Assumes R has coefficients in QQ
     	       	    	 lift(((m*mono) % I)/topBas1, QQ)
     	       ))
     ));
--Temporarily assuming that cone is full-dimensional - is it always???
--<<"Got this far with nefCone"<<endl;
--<<rank source Mat <<"    "<<rank target Mat <<endl;
    matDual:=-1*(fourierMotzkin Mat)#0;
    return(matDual);
);


--Code to compute the cone nef^k_i, which is the cone of all
-- codimension k-cycles that intersect every effective i-dimensional
-- cycle in effective cycles.
--Currently returns a rather arbitrary  basis for the codim-k Chow group 
-- and then a matrix whose columns represent generators for this 
--cone

nefCone2=(k,i,X)->(
     if k>i then error("i must be at least k");
     if not isSmooth(X) then error("Not implemented yet");
     n:=#((rays X)#0);
     Conesk:=orbits(X,k);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     --Now create the multiplication map
     --First get a basis for AA_k
     if not X.cache.?ChowRingBas then
     	  X.cache.ChowRingBas = new MutableHashTable;
     if not X.cache.ChowRingBas#?(n-k) then 	  
          X.cache.ChowRingBas#(n-k)=flatten entries lift(basis(k,R/I),R);
     --We'll create a bas times Conesi matrix with (j,k) entry bas#j * Conesi #k
     --???need to edit from here.	  
     mono:=1_R;
     for i in (max X)_0 do mono=mono*R_i;
     topBas1:=mono % I;
     Mat:=matrix unique apply(X.cache.ChowRingBas#(n-k),m->(
	       apply(Conesk,sigma->(
			 mono:=1_R;
			 for j in sigma do mono=mono*R_j;
     	       	    	 --Assumes R has coefficients in QQ
     	       	    	 lift(((m*mono) % I)/topBas1, QQ)
     	       ))
     ));
--Temporarily assuming that cone is full-dimensional - is it always???
     matDual:=-1*(fourierMotzkin Mat)#0;
     return(matDual);
);


--Compute the Mob^k_i cone.  This is the cone
-- cap_{sigma in Sigma(n-i)} pos( V(tau) : tau in Sigma(k), tau cap sigma = emptyset)
-- k is the codimension
--Note:  In the ToricCycles.pdf notes it was suggested that this might equal 
--nef.  This seems not to be the case for ~20% of the smooth Fano 4-folds
Mob=(k,i,X)->(
     local Tauset; local Taulist;
     local mono; local tauCone;
     if k>i then error("i must be at least k");
     if not isSmooth(X) then error("Not implemented yet");
     n:=#((rays X)#0);
     Conesk:=orbits(X,n-k);
     Conesi:=orbits(X,i);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     K:=coefficientRing R;
     coeffMap:=map(K, R, apply(numgens R,m->1_K));
     --Now create the multiplication map
     --First get a basis for AA_k
     chowBask:=chowRingBasis(X,n-k);
     --We'll create a bas times Conesi matrix with (j,k) entry bas#j * Conesi #k
     --???need to edit from here.	  
     TotalFacets:=flatten apply(Conesi, sigma->(
	       Tauset=select(Conesk,tau->(all(tau,i->(not(member(i,sigma))))));
--To be changed:  Here we're assuming that the basis for the Chow ring is the 
-- one given by the Grobner basis for SR(I). 
	       Taulist=apply(Tauset,tau->(
			 mono=1_R;
		 	 for j in tau do mono=mono*R_j;
			 mono%I
	       ));
     	       contracted:=contract(matrix {chowBask}, transpose matrix {Taulist});
 	       tauCone = coeffMap transpose contracted;
	       transpose entries(-1*(fourierMotzkin(tauCone))#0)
     ));	       
     MobCone:=-1*(fourierMotzkin transpose matrix TotalFacets)#0;
     return(MobCone);
);


--Compute the L^k_a cone.  
--This is the cone 
--cap_{|sigma| leq n-k} pos(V(tau) : tau in Sigma(k), 
--tau cap sigma = emptyset, tau cup sigma in Sigma) + 
--span( V(tau) : tau in Sigma(k), tau cap sigma = emptyset, 
--tau cup sigma in Sigma}
--k is the codimension

Lcone=(k,i,X)->(
     local Tauset; local Taulist;
     local mono; local tauCone;
     local tauCone2; local tausigmainSigma;
     local tausigmaNot; local jFacets;
     if k>i then error("i must be at least k");
     if not isSmooth(X) then error("Not implemented yet");
     n:=#((rays X)#0);
     Conesk:=orbits(X,n-k);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     K:=coefficientRing R;
     coeffMap:=map(K, R, apply(numgens R,m->1_K));
     --Now create the multiplication map
     --First get a basis for AA_k
     chowBask:=chowRingBasis(X,n-k);
     --We'll create a bas times Conesi matrix with (j,k) entry bas#j * Conesi #k
     --???need to edit from here.	  
     TotalFacets:=flatten apply(n-i+1,j-> (
    	       Conesj:=orbits(X,n-j);
	       jFacets=flatten apply(Conesj, sigma->(
	       		 Tauset=select(Conesk,tau->(all(tau,i->(not(member(i,sigma))))));
--To be changed:  Here we're assuming that the basis for the Chow ring is the 
-- one given by the Grobner basis for SR(I). 
                         tausigmainSigma={};
	       		 tausigmaNot={};
     	       		 scan(Tauset,tau->(
 			 	   mono=1_R;
		 	 	   for j in tau do mono=mono*R_j;
			 	   if member(sort(sigma |tau), orbits(X,n-k-j)) then
			      	   tausigmainSigma=append(tausigmainSigma,mono%I)
			 	   else tausigmaNot=append(tausigmaNot,mono%I);
	       		 ));
     	       		 contracted:=contract(matrix {chowBask}, transpose matrix {tausigmainSigma});
 	       		 tauCone = coeffMap transpose contracted;
               		 contracted2:=contract(matrix {chowBask}, transpose matrix {tausigmaNot});
 	       		 tauCone2 = coeffMap transpose contracted2;
--<<tauCone<<endl<<tauCone2<<endl;			 
	       		 transpose entries(-1*(fourierMotzkin(tauCone,tauCone2))#0)
     	       ));
     	       jFacets     	       	       
     ));
     MobCone:=-1*(fourierMotzkin transpose matrix TotalFacets)#0;
     return(MobCone);
);




--Compute the L^k_b cone.  
--This is the cone 
--cap_{|sigma| leq n-k} pos(V(tau) : tau in Sigma(k), V(tau).V(sigma) is effective) +
-- span( V(tau) : V(tau).V(sigma)=0 in Sigma(k))
--k is the codimension

Lconeb=(k,i,X)->(
     local Tauset; local Taulist;
     local mono; local tauCone;
     local tauCone2; local tausigmainSigma;
     local tausigmaNot; local jFacets;
     local effkjFacets; local mono2;
     local remainderSigmaTau;
     if k>i then error("i must be at least k");
     if not isSmooth(X) then error("Not implemented yet");
     n:=#((rays X)#0);
     Conesk:=orbits(X,n-k);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     K:=coefficientRing R;
     coeffMap:=map(K, R, apply(numgens R,m->1_K));
     --Now create the multiplication map
     --First get a basis for AA_k
     chowBask:=chowRingBasis(X,n-k);
     --We loop through j from 0 to i, and compute for each j 
     --the facets of the cone pos(V(sigma) : V(sigma).V(tau) is effective for all
     -- tau in Sigma(j)) + span(V(sigma) : V(sigma).V(tau) =0 for all tau 
     -- in Sigma(j)
     TotalFacets:=flatten apply(n-i+1,j-> (
--<<"j   "<<j<<endl;	       
    	       Conesj:=orbits(X,n-j); 
               chowBaskj:=chowRingBasis(X,n-k-j);
	       --facets of effCone(n-k-j) are the rows of effkjFacets
	       effkjFacets=-1*(fourierMotzkin(effCone(n-k-j,X)))#0;
	       effkjFacets=transpose effkjFacets;
	       jFacets=flatten apply(Conesj, sigma->(
                         tausigmainSigma={};
	       		 tausigmaNot={};
--To be changed:  Here we're assuming that the basis for the Chow ring is the 
-- one given by the Grobner basis for SR(I). 
                         mono=1_R;
	       		 scan(sigma,j->(mono=mono*R_j));
   	       		 scan(Conesk,tau->(
 			 	   mono2=1_R;
		 	 	   for jj in tau do mono2=mono2*R_(jj);
				   remainderSigmaTau=((mono*mono2) % I);
				   if remainderSigmaTau==0_R then tausigmaNot=append(tausigmaNot,mono2%I)
				   else (
				   	remainderVector:=transpose coeffMap contract(matrix {chowBaskj}, remainderSigmaTau);
			 	   	if min(flatten entries(effkjFacets*remainderVector)) >=0 then 
			      	   	    tausigmainSigma=append(tausigmainSigma,mono2%I);
				   );	    
	       		 ));
--<<tausigmainSigma<<endl;		     
     	       		 contracted:=contract(matrix {chowBask}, transpose matrix {tausigmainSigma});
 	       		 tauCone = coeffMap transpose contracted;
               		 contracted2:=contract(matrix {chowBask}, transpose matrix {tausigmaNot});
 	       		 tauCone2 = coeffMap transpose contracted2;
--<<tauCone<<endl<<tauCone2<<endl;			 
	       		 transpose entries(-1*(fourierMotzkin(tauCone,tauCone2))#0)
     	       ));
--<<endl<<jFacets<<endl;
     	       jFacets     	       	       
     ));
     MobCone:=-1*(fourierMotzkin transpose matrix TotalFacets)#0;
     return(MobCone);
);
--???bug somewhere in what's above???



--Compute the L^k_b cone.  
--This is the cone 
--cap_{|sigma| leq n-k} pos(V(tau) : tau in Sigma(k), V(tau).V(sigma) is effective) +
-- span( V(tau) : V(tau).V(sigma)=0 in Sigma(k))
--k is the codimension

--This is only just started.
LconevsNef=(k,X)->(
     local Tauset; local Taulist;
     local mono; local tauCone;
     local tauCone2; local tausigmainSigma;
     local tausigmaNot; local jFacets;
     local effkjFacets; local mono2;
     local remainderSigmaTau;
     --I don't know what i was here - this is a cludge to make it run
--     if k>i then error("i must be at least k");
i:=k;
     if not isSmooth(X) then error("Not implemented yet");
     n:=dim X;
     Conesk:=orbits(X,n-k);
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     K:=coefficientRing R;
     coeffMap:=map(K, R, apply(numgens R,m->1_K));
     --Now create the multiplication map
     --First get a basis for AA_k
     chowBask:=chowRingBasis(X,n-k);
     --We loop through j from 0 to i, and compute for each j 
     --the facets of the cone pos(V(sigma) : V(sigma).V(tau) is effective for all
     -- tau in Sigma(j)) + span(V(sigma) : V(sigma).V(tau) =0 for all tau 
     -- in Sigma(j)
     TotalFacets:=flatten apply(n-i+1,j-> (
--<<"j   "<<j<<endl;	       
    	       Conesj:=orbits(X,n-j); 
               chowBaskj:=chowRingBasis(X,n-k-j);
	       --facets of effCone(n-k-j) are the rows of effkjFacets
	       effkjFacets=-1*(fourierMotzkin(effCone(n-k-j,X)))#0;
	       effkjFacets=transpose effkjFacets;
	       jFacets=flatten apply(Conesj, sigma->(
                         tausigmainSigma={};
	       		 tausigmaNot={};
--To be changed:  Here we're assuming that the basis for the Chow ring is the 
-- one given by the Grobner basis for SR(I). 
                         mono=1_R;
	       		 scan(sigma,j->(mono=mono*R_j));
   	       		 scan(Conesk,tau->(
 			 	   mono2=1_R;
		 	 	   for jj in tau do mono2=mono2*R_(jj);
				   remainderSigmaTau=((mono*mono2) % I);
				   if remainderSigmaTau==0_R then tausigmaNot=append(tausigmaNot,mono2%I)
				   else (
				   	remainderVector:=transpose coeffMap contract(matrix {chowBaskj}, remainderSigmaTau);
			 	   	if min(flatten entries(effkjFacets*remainderVector)) >=0 then 
			      	   	    tausigmainSigma=append(tausigmainSigma,mono2%I);
				   );	    
	       		 ));
--<<tausigmainSigma<<endl;		     
     	       		 contracted:=contract(matrix {chowBask}, transpose matrix {tausigmainSigma});
 	       		 tauCone = coeffMap transpose contracted;
               		 contracted2:=contract(matrix {chowBask}, transpose matrix {tausigmaNot});
 	       		 tauCone2 = coeffMap transpose contracted2;
--<<tauCone<<endl<<tauCone2<<endl;			 
	       		 transpose entries(-1*(fourierMotzkin(tauCone,tauCone2))#0)
     	       ));
--<<endl<<jFacets<<endl;
     	       jFacets     	       	       
     ));
     nefC:=transpose nefCone(n-k,X);
<<     nefC*(transpose matrix TotalFacets);
  --   MobCone:=-1*(fourierMotzkin transpose matrix TotalFacets)#0;
    -- return(MobCone);
);




--Compute the effective cone of i cycles in X
-- i is the dimension?
effCone=(i,X)->(
     if not isSmooth(X) then error("Not implemented yet");     
     n:=dim X;
     --Get SR ring
     I:=SR(X);
     R:=ring I;
     if not X.cache.?ChowRingBas then
     	  X.cache.ChowRingBas = new MutableHashTable;
     if not X.cache.ChowRingBas#?i then 	  
          X.cache.ChowRingBas#i=flatten entries lift(basis(n-i,R/I),R);
     --j=n-i
     Conesj:=orbits(X,i);
     EffMat:=transpose matrix apply(Conesj,sigma->(
    	  mono:=1_R;
	  for j in sigma do mono=mono*R_j;
     	  mono=mono % I;
	  apply(X.cache.ChowRingBas#i,m->(coefficient(m,mono)))
     ));	  
     return(EffMat);
);



--Compute the intersection of all simplices in eff containing nef
--???temporary working code???

nefEqualsIntersection=(i,X)->(
     d:=rank(AA(i,X));
     nefc:=nefCone(i,X);
     eff:=effCone(i,X);
     n:=rank source eff;
     Subs:=subsets(n,d);
     ConestoCheck:={};
     Indices:={};
     scan(Subs,s->(
	       As:=eff_s;
	       if abs(det(As))>0 then (
		      if isContainedCones(nefc,As) then (
		           ConestoCheck=append(ConestoCheck,As);
			   Indices=append(Indices,s);
		      );
     	       );
     ));
     A:=(fourierMotzkin(ConestoCheck#0))#0;
--<<A<<endl;     
     scan(#ConestoCheck-1, i->(
	       A = A | (fourierMotzkin(ConestoCheck#(i+1)))#0;
     ));  
     intersectAllCones:=(fourierMotzkin(A))#0;
--     <<nefc<<endl<<intersectAllCones<<endl;
     if (isContainedCones(nefc,intersectAllCones) and isContainedCones(intersectAllCones,nefc))
     	  then return(true) else return(false);
);     
  
     
     
   


---------------------------------------------------------------------------
-- DOCUMENTATION
---------------------------------------------------------------------------
beginDocumentation()

document { 
     Key => Chow,
     Headline => "intersection theory for normal toric varieties",
     "This is a subpackage for eventual inclusion into Greg Smith's
NormalToricVarieties package",
     PARA{},     
     "It contains routines to do compute the Chow ring and groups of a
normal toric variety, plus compute the nef and effective cones of
cycles." 
     }  

-- document { 
--      Key => {(cones, ZZ, NormalToricVariety)},
--      Headline => "the i-dimensional cones of the fan of X",
--      Usage => "cones(i,X)",
--      Inputs => {"i" => ZZ, "X" => NormalToricVariety},
--      Outputs => {List => " of lists; each entry is the index set of
-- rays in an i-dimensional cone of the fan of X."},
--      "This procedure produces a list of the i-dimensional cone of the fan of a toric variety X",
--      EXAMPLE lines ///
-- 	  X = projectiveSpace 2;
-- 	  cones(2,X)
-- 	  max X
-- 	  cones(1,X)
-- 	  rays X
-- 	  ///,
-- 	  "Note that an i-dimensional cone may have more than i extreme rays",
--      EXAMPLE lines ///
-- 	  X=normalToricVariety({{1,0,0},{1,2,0},{1,2,2},{1,0,2},{-1,-1,-1}}, {{0,1,2,3},{0,1,4}, {1,2,4},{2,3,4},{0,3,4}});
-- 	  cones(3,X)
-- 	  ///,
--      }


document { 
     Key => AA,
     Headline => "Chow rings for toric varieties",
     Usage => "AA(i,X)",
     Inputs => {"i" => ZZ, "X" => NormalToricVariety},
     Outputs => {"the codim-i Chow group $A^i(X)$, an abelian group (a  ZZ-module)"},
     "This procedure computes the ith Chow group of the
NormalToricVariety X.  It produces it as the cokernel of a matrix,
following the description given in Proposition 2.1 of Fulton-Sturmfels
\" Intersection Theory on toric varieties \" (Topology, 1996).",
PARA{}, "It is cached in X.cache.Chow#i. ", 
" ???say something about pruning map ", 
PARA{},
" These groups are all one-dimensional for projective space. ", 
EXAMPLE lines /// 
X = projectiveSpace 4 
rank AA(1,X) 
rank AA(2,X) 
rank AA(3,X)
///, 
"We next consider the blow-up of P^3 at two points." ,
EXAMPLE lines ///
X=normalToricVariety({{1,0,0},{0,1,0},{0,0,1},{-1,-1,-1},{1,1,1}, {-1,0,0}}, {{0,2,4},{0,1,4},{1,2,4},{1,2,5},{2,3,5},{1,3,5},{0,1,3},{0,2,3}})
AA(1,X) 
AA(2,X) 
/// 
}


document {
     Key => effCone,
     Headline => "the cone of effective T-invariant i-cycles",   
     Usage => "effCone(i,X)",
     Inputs => {"i" => ZZ, "X" => NormalToricVariety},
     Outputs => {Matrix => "whose columns are the generators for
the cone of effective i-cycles "},
     "This is currently only implemented for smooth toric varieties. ",
      "The columns should be given in a basis for the i-th Chow group
recorded in X.cache.ChowRingBas#i. ",
     EXAMPLE lines ///
     X = projectiveSpace 4
     effCone(2,X)
     ///,
     EXAMPLE lines ///
     X = hirzebruchSurface 1;
     effCone(1,X)
     ///,
}    


document {
     Key => nefCone,
     Headline => "the cone of nef T-invariant i-cycles",   
     Usage => "nefCone(i,X)",
     Inputs => {"i" => ZZ, "X" => NormalToricVariety},
     Outputs => {Matrix => "whose columns are the generators for
the cone of nef i-cycles "},
     "A cycle is nef if it intersects every effective cycle of
complementary dimension nonnegatively. ",
      "This is currently only implemented for smooth toric varieties. ",
      "The columns are given in a basis for the i-th Chow group
recorded in X.cache.ChowRingBas#i. ",
    EXAMPLE lines ///
    X=projectiveSpace 4
    nefCone(2,X)
    ///,
    EXAMPLE lines ///
    X=hirzebruchSurface 1;
    nefCone(1,X)
    ///,
}     

-- document { 
--      Key => isCartier,
--      Headline => " if a Weil divisor is a Cartier divisor ",
--      Usage => "isCartier(D)",
--      Inputs => { "D" => List, " of length the number of rays of X,
-- representing a Weil divisor"},
--      Outputs => {Boolean},
--      "???"
-- }     
     
     
-- document {
--      Key => (intersect, List, ZZ, List, NormalToricVariety),
--      Headline => "intersect a Cartier divisor with a k-cycle",
--      Usage => "intersect(D, k, Z, X)",
--      Inputs => {"D" => List, " of length the number of rays of X,
-- being the Cartier divisor representation of a Weil divisor", k => ZZ,
-- "Z" => List, " the ray indices of a cone of codimension k in the fan of X"},
--      Outputs => {List, " of length the number of cones of codimension
-- (k-1) in the fan of X, representing the (k-1)-cycle D.V(tau)"},
--      "   "
-- }     
     
     
     
document {
     Key => SR,
     Headline => "compute the Chow ring of a smooth toric variety",
     Usage => "SR(X)",
     Inputs => {"X" => NormalToricVariety},
     Outputs => {Ideal => "which defines the relations on the Chow
ring of X"},
     "The ring of the ideal has one generator for each ray of X, and
the ideal is the ideal given in the Stanley-Reisner presentation of
the cohomology ring of X. ", PARA{},
     "This assumes that X is smooth.  Eventually it will be
implemented for simplicial toric varieties. ",
    EXAMPLE lines ///
    X = projectiveSpace 2
    I = SR X
    R = ring I
    for i from 0 to 2 do <<hilbertFunction(i,R/I)<<endl
    ///,
   	PARA{},
 	  "Next we consider the blow-up of P^3 at 2 points. ",
    EXAMPLE lines ///
    X=normalToricVariety({{1,0,0},{0,1,0},{0,0,1},{-1,-1,-1},{1,1,1}, {-1,0,0}}, {{0,2,4},{0,1,4},{1,2,4},{1,2,5},{2,3,5},{1,3,5},{0,1,3},{0,2,3}})
    I = SR X
    R = ring I
    hilbertFunction(1,R/I)
    ///,
    "Note that the degree-one part of the ring has dimension the
Picard-rank, as expected. ", 
}

document {
     Key => isContainedCones,
     Headline => "decide if one cone is contained inside another",
     Usage => "isContainedCones(M,N)",
     Inputs => {"M"=>Matrix, "N"=>Matrix},
     Outputs => {Boolean},
     "Returns true if the cone generated by the columns of the matrix
M is contained in the cone generated by the columns of the matrix N.
",
"This currently assumes that both cones are full-dimensional, and is implemented in 
a somewhat hackish manner."
}
---------------------------------------------------------------------------
-- TEST
---------------------------------------------------------------------------

--Replace X by something more interesting
TEST ///
X=projectiveSpace 4
assert(rank AA(3,X) == rank AA(1,X))
assert(rank AA(3,X) == rank picardGroup X)
/// 

TEST ///
A=sort apply(3,i->random(5))
X=kleinschmidt(6,A)
R=QQ[x,y]
I=ideal(x^4,y^4)
for i from 0 to 6 do
     assert(rank AA(i,X) == hilbertFunction(i,R/I))
///

end

---------------------------------------------------------------------------
-- SCRATCH SPACE
---------------------------------------------------------------------------
 
restart
loadPackage "Chow"
--X is blow-up of P^3 at two points
raysX={{1,0,0},{0,1,0},{0,0,1},{-1,-1,-1},{1,1,1}, {-1,0,0}};
Sigma={{0,2,4},{0,1,4},{1,2,4},{1,2,5},{2,3,5},{1,3,5},{0,1,3},{0,2,3}};
X=normalToricVariety(raysX,Sigma);

--X is the blow-up of P^4 at 2 points

raysX={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{-1,-1,-1,-1},{1,1,1,1},{-1,0,0,0}};
Sigma={{0,1,2,4},{0,1,3,4}, {0,2,3,4}, {0,1,2,5}, {0,1,3,5},
{0,2,3,5}, {1,2,3,5}, {1,2,3,6}, {1,2,4,6}, {1,3,4,6}, {2,3,4,6}};
X=normalToricVariety(raysX,Sigma);

--Y is the blow-up of P^4 at 1 point
raysY={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{-1,-1,-1,-1},{1,1,1,1}};
Sigma2 = {{0,1,2,4}, {0,1,3,4}, {0,2,3,4}, {1,2,3,4}, {0,1,2,5}, {0,1,3,5},
     {0,2,3,5}, {1,2,3,5}};
Y=normalToricVariety(raysY,Sigma2);

-- document { 
--      Key => {(cones, (ZZ, NormalToricVariety)},
--      Headline => "the i dimension cones of the fan",
--      Usage => "cones(i,X)"
--      Inputs => {
-- 	  "i" => "a nonnegative integer",
-- 	  "X" => NormalToricVariety
-- 	  },
--      Outputs => {},

--      EXAMPLE lines ///
-- 	  PP1 = projectiveSpace 1;
-- 	  ///,
--      SeeAlso => {normalToricVariety, weightedProjectiveSpace,
-- 	  (ring,NormalToricVariety), (ideal,NormalToricVariety)}
--      }     


--stellarSubdivision bug/features
----doesn't add new ray at end
----when adding ray that is already there it doesn't realize it
----Also creates X.cache.cones (name clash)

--For example - try blow-up of P4 at 2 points



W=sort apply(5,i->random(7)+1);
while not all(subsets(W,4), s -> gcd s === 1) do 
      W=sort apply(5,i->random(7)+1);
X=resolveSingularities weightedProjectiveSpace(W);
summ=0;
I=SR(X);
R=ring I;
scan(5,i->(summ=summ+(hilbertFunction(i,R/I)-rank(AA(i,X)))^2;));
<<summ<<endl;


uninstallPackage "Chow"
restart
installPackage "Chow"
