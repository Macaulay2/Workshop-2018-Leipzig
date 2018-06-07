restart
------this function computes the Apolar ideal of a given F in S=QQ[x_1..x_n]
perp = (F) ->( 
    S := ring F;
    d := first degree F;
    I = trim sum for i from 1 to d+1 list ideal((sort basis(i,S))*syz diff(transpose (sort basis(d-i,S)),diff((sort basis(i,S)),F))))
----this function checks if given ideal is Apolar to a given form F
isApolar = (I,F) -> (perp(F)==trim I)
----this function computes the apolar ideal of a homogeneous ideal I
apolarIdeal = (I) -> (
    S := ring I;
    intersect(for i in flatten entries mingens I list perp(i))
    )
---
---this function gives the waring rank of a binary form F
rk = (F) ->(
    S := ring F;
    x := (gens S)_0;
    y := (gens S)_1;
    I := perp(F);
    g1 := first flatten entries gens I;
    g2 := last flatten entries gens I;
    m := (min {degree g1,degree g2})_0;
    g := if (degree g1)_0 == m then g1 else g2;
    if g!=0 and resultant(diff(x,g),g,x)!=0 or resultant(diff(y,g),g,y)!=0 then m else (max {degree g1,degree g2})_0   
    ) 
----this function computes the number of essential variables of a form F
numEssVars = (F) -> (
    S := ring F;
    d := first degree F;
    I := perp(F);
    m = hilbertFunction(1,I)
    )
----this function gives a set of essential variables of a form F
essVars = (F) -> (
    S := ring F;
    d := first degree F;
    flatten entries basis (1,S/perp(F))
)
----this function computes the Catalectican matrix of a form F
cat = (i,F) -> (
    R := ring(F);      d := first degree F;
    return diff(transpose basis(d-i,R), diff(basis(i,R),F));
    )
----this function computes the rank of a form F where 
----the number of essential variables is at most 2
rkGeneral = (F) -> (
    if numEssVars F == 2 then (
	R = QQ[essVars F];
	F = sub(F,R);
	rk F
	)   
    else (
    	if  numEssVars F == 1 then (1)
    	else error "there are more than two essential variables "
         )
    )


-----example
S = QQ[x,y,z]
F = x^2+y^2+z^2
F = (x+y)^2+z^2
F = x^2
essVar(F)
essVars(F)
numEssVars F
perp F
rk(F)
rkGeneral F

------

