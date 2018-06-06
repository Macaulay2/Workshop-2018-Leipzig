restart
------this function computes the Apolar ideal of a given F in S=QQ[x_1..x_n]
perp = (F) ->( 
    S := ring F;
    d := first degree F;
    I = trim sum for i from 1 to d+1 list ideal((sort basis(i,S))*syz diff(transpose (sort basis(d-i,S)),diff((sort basis(i,S)),F))))
----this function checks if given ideal is Apolar to a given form
isApolar = (I,F) -> (perp(F)==trim I)

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
----
numEssVars = (F) -> (
    S := ring F;
    d := first degree F;
    I := perp(F);
    m = hilbertFunction(1,I)
    )
essVars = (F) -> (
    S := ring F;
    d := first degree F;
    m :=numEssVars F;
    L = first entries ( basis(1,S) * syz diff(transpose (sort basis(d-1,S)),diff((sort basis(1,S)),F)));
    flatten entries basis (1,S/perp(F))
    )

rkGeneral = (F) -> (
    if numEssVars F == 2 then (
	R = QQ[essVars F];
	F = sub(F,R);
	rk F
	)
    else error "there are more than two essential variables "
    )
-----example
S = QQ[x,y]
F = x*y^2+y*x^2
F = x^3
F = x^2*y+x*y^2
F = random(5,S)
perp F
rk(F)
numEssVars(F)
S = QQ[x,y,z]
F = x^2+y^2+z^2
F = (x+y)^2+z^2
essVars F
numEssVars F
rkGeneral F

