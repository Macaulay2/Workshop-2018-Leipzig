formalLog = f -> (
    assert ((reverse flatten entries last coefficients(f))_0 == 1);
    d = #(coefficients f)_0;
    sum for i from 1 to d list (-1)^(i-1)/i * (f-1)^i
    )


-- TEST --

R = QQ[x,y]
S = QQ[x,y][t]/t^10

f = x*t + 1/2 * y^2 * t^2

g = exp(f)

formalLog(g)