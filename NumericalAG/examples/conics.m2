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

F = for i from 1 to 5 list 
	sub(sub((1/100.0)*f, {
	b_1 => random CC,
	b_2 => random CC,
	b_3 => random CC,
	b_4 => random CC,	
	b_5 => random CC,
	b_6 => random CC}),S);
