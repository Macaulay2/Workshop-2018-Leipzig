
cumulantsToMoments = cumulants -> (
    assert(cumulants#0 == 0);
    S := QQ[t]/t^#cumulants;
    use S;
    p := exp(sum for i from 1 to #cumulants-1 list t^i * cumulants#i/i!);
    return for i from 0 to #cumulants-1 list i!*coefficient(t^i, p);
    )

momentsToCumulants = moments -> (
    assert(moments#0 == 1);
    S := QQ[t]/t^#moments;
    use S;
    p := formalLog(sum for i from 1 to #moments-1 list t^i * moments#i/i!);
    return for i from 0 to #moments-1 list i!*coefficient(t^i, p);
    )

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