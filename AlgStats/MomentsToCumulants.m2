cumulantsToMoments = (cumulants,R) -> (
    assert(cumulants#0 == 0);
    return powerSeriesTransform(cumulants, exp, R)
    )

momentsToCumulants = (moments,R) -> (
    assert(moments#0 == 1);
    return powerSeriesTransform(moments, x -> formalLog(x,#moments), R)
    )

powerSeriesTransform = (l, f, R) -> (
    S := R[t]/t^#l;
    use S;
    p := f(sum for i from 0 to #l-1 list t^i * l#i/i!);
    return for i from 0 to #l-1 list i!*coefficient(t^i, p)
    )

formalExp = (f, d) -> (
    sum for k from 0 to d list 1/(k!) * f^k
)

formalLog = (f, d) -> (
    sum for k from 1 to d list (-1)^(k-1)/k * (f-1)^k
)

momentIdealToCumulants = (I,truncAbove) -> (
    idealTransform(I, truncAbove, formalExp, 0)  
    )

cumulantToMoments = (I, truncAbove) -> (
    idealTransform(I, truncAbove, formalLog, 1)
    )

-- Applies a change of variables to the ideal I
idealTransform = (I, truncAbove, f, substVal) -> (
    s := symbol s;
    t := symbol t;
    R2 := QQ[s_0..s_truncAbove];
    use R2;
    R2' := R2[t]/t^(truncAbove+1);
    use R2';
    p := sum for i from 0 to truncAbove list 1/i! * s_i * t^i;
    q := f(p,truncAbove+1);
    -- NOTE: here, we dehomogenize. May think about giving homog/nonhomog as an option.
    li := for i from 0 to truncAbove list i! * coefficient(t^i, q);
    li = for i from 0 to truncAbove list sub(li_i, s_0 => substVal);
    phi := map(R2, ring I,li);
    phi I
    )

-- momentsTensorFormat is a tensor (M_I)_I where M_I = m_I, I 
momentIdealToCumulantsMultivariate = (I, truncAbove, momentsTensorFormat) -> (
    return 0
)

end


-- TEST --

R = QQ[x,y] 

cumulants = {0,1,1}
cumulantsToMoments(cumulants,R)

moments = {1,1,2}
momentsToCumulants(moments, R)

cumulants = {0,x,y^2,0,0,0,0,0,0}
li = cumulantsToMoments(cumulants,R)
momentsToCumulants(li,R)

truncAbove = 4
I = momentIdealGaussian(1,truncAbove)
transpose gens gb I
C = momentIdealToCumulants(I,truncAbove)
M = cumulantIdealToMoments(C, truncAbove)
transpose gens gb M
transpose gens gb sub(I, {first gens ring I => 1})

truncAbove = 4
momentsTensorFormat = {truncAbove, truncAbove}
