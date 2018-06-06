restart
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
    sum for k from 0 to d-1 list 1/(k!) * f^k
)

formalLog = (f, d) -> (
    sum for k from 1 to d-1 list (-1)^(k-1)/k * (f-1)^k
)

-- TEST --

R = QQ[x,y] 

cumulants = {0,1,1}
cumulantsToMoments(cumulants,R)

moments = {1,1,2}
momentsToCumulants(moments, R)

cumulants = {0,x,y^2,0,0,0,0,0,0}
li = cumulantsToMoments(cumulants,R)
momentsToCumulants(li,R)

truncAbove = 3
I = momentIdealExponential(1,truncAbove)
R2 = QQ[k_0..k_truncAbove]
R2' = R2[t]/t^(truncAbove+1)
p = sum for i from 0 to truncAbove list 1/i! * k_i * t^i
q = formalExp(p,truncAbove+1)
li = for i from 0 to truncAbove list i! * coefficient(t^i, q)
li = for i from 0 to truncAbove list sub(li_i,{k_0 => 0})
phi = map(R2, QQ[m_0..m_truncAbove],li)
I = sub(I, source phi)
J = phi I
gens gb J
