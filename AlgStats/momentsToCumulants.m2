restart
cumulantsToMoments = (cumulants,R) -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp, R)
    )

momentsToCumulants = (moments,R) -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, formalLog, R)
    )

powerSeriesTransform = (l, f, R) -> (
    S := R[t]/t^#l;
    use S;
    p := f(sum for i from 0 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p);
    )

formalLog = f -> (
    assert ((reverse flatten entries last coefficients(f))_0 == 1);
    d = #(coefficients f)_0;
    sum for i from 1 to d list (-1)^(i-1)/i * (f-1)^i
    )


-- TEST --

R = QQ[x,y]

cumulants = {0, 1, 1}
cumulantsToMoments(cumulants,R)

moments = {1, 1, 2}
momentsToCumulants(moments,R)

cumulants = {0, x, y^2, 0, 0, 0, 0}
l = cumulantsToMoments(cumulants,R)
l_2
