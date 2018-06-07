restart
R = QQ[x, y, a_1..a_6, b_1..b_6]

f = a_1*x^2 + a_2 *x*y + a_3*y^2 + a_4*x + a_5*y + a_6
g = b_1*x^2 + b_2 *x*y + b_3*y^2 + b_4*x + b_5*y + b_6 

Jac = diff(matrix {{x, y}}, matrix {{f}, {g}})

I = ideal(f, g, det Jac)

J = eliminate({x,y}, I)


p=(gens J)_(0,0)
m=last coefficients p
sum apply(flatten entries m,n->abs sub(n,CC))/(numrows m)
apply(flatten entries gens J,g->norm sub(first coefficients g,CC))

S'= CC[a_1..a_6,b_1..b_6]
S = CC[take(gens S',6)]

f = sub(J_0, S');


bs=apply(0..5,i->random CC)
F = for i from 1 to 5 list 
	sub(sub((1/10.0)*f, {
	b_1 => b#0,
	b_2 => b#1,
	b_3 => b#2,
	b_4 => b#3,	
	b_5 => b#4,
	b_6 => b#5}),S);

needs "../julia.m2"
importJulia({"using HomotopyContinuation","import DynamicPolynomials: PolyVar"},JuliaProcess)
P=polySystem F;
peek P

sols=solveJulia P;

nonSing=0
elapsedTime for s in sols do(
    print("nonsingular vs total conics so far: " | nonSing | " vs " | tot);
    p=s.Coordinates;
    R'=CC[x,y,z];
    f'=polySystem {p#0*x^2 + p#1 *x*y + p#2*y^2 + p#3*x*z + p#4*y*z + p#5*z^2};
    if #(solveSystem polySystem(ideal(f')+ ideal jacobian f'))==1 then nonSing=nonSing+1;
    tot=tot+1;
    )
    


L=apply(sols,x->x-> first x.Coordinates);
netList sort L
x=first L
