close JuliaProcess -- if not already open
restart
needs "../julia.m2"
importJulia({"using HomotopyContinuation","import DynamicPolynomials: PolyVar"},JuliaProcess)

R'=CC[x,y]
f= x^2+y
g=y^2-pi*ii
P=polySystem {f,g};
sols=solveJulia P;
netList sols


R = QQ[x, y, a_1, a_2, b_1, b_2, r, s]

f = (x-a_1)^2 + (y-a_2)^2 - r^2
g = (x-b_1)^2 + (y-b_2)^2 - s^2

Jac = diff(matrix {{x, y}}, matrix {{f}, {g}})


I = ideal(f, g, det Jac)

J = eliminate({x,y}, I)

setRandomSeed 12
F = for i from 1 to 3 list 
	sub(J_0, {
	b_1 => random QQ, 
	b_2 => random QQ,
	s => random QQ});

--netList F

sols=solveJulia polySystem F;
realPoints sols

--netList sols
-- real solutions?
zer0=10.0^(-10);
r=select(sols,s->all(s.Coordinates,x->abs imaginaryPart x < zer0));
netList r
#r

