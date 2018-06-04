restart
cumulantsToMoments = cumulants -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp)
    )

momentsToCumulants = moments -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, formalLog)
    )

powerSeriesTransform = (l, f) -> (
    S := QQ[t]/t^#l;
    use S;
    p := f(sum for i from 1 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p);
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

cumulants = {0, 1, 1}
cumulantsToMoments(cumulants)

moments = {1, 1, 2}
momentsToCumulants(moments)
